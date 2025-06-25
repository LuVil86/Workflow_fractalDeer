

#set source R scripts wd
setwd("/media/loreto/Grande/ie-ofev-24-25/workflow_fractalDeer/PIPELINE_FRACTALDEER")

source("./getNSDValues.R")
source("./refineClusteringByNSD.R")
require(gridExtra)
require(multidplyr)
require(tidyverse)
require(lubridate)
require(adehabitatLT)
require(arulesViz)
require(gdata)
require(ggplot2)
require(ggtext)
require(mclust)
require(amt)
require(grid)
require(FactoMineR)
require(factoextra)

require(ordinal)
require(pals)
require(Hmisc)
require(adehabitatHR)
require(reshape2)


############ *** SHOWING DIFFERENCES BETWEEN DEERS WITH ONE OR TWO BEHAVIOURS FOR DIFFERENT trajectory MEASURES ***** ####################################################

tac2<-function(x, ...){
  x<-steps(x)
  ta<-x$ta_[-1]
  1/nrow(x) * sum(diff(cos(ta))^2 + diff(sin(ta))^2, na.rm=T)
}

##### generate path measures for residents ######
setwd("/media/loreto/Grande/ie-ofev-24-25/cerf_movement_patterns/plateau_selected/unimodal_cerf_plateau")
fullTrajStatList<-list()
k<-1
for(behaviourFile in list.files(pattern=".*REFINED.*.csv")){
  
  selectedStepSize<-strsplit(behaviourFile, split=" ")[[1]][5]
  selectedDeer<-strsplit(behaviourFile, split=" ")[[1]][1]
  
  selectedYear<-paste(strsplit(behaviourFile, split=" ")[[1]][2], strsplit(behaviourFile, split=" ")[[1]][3], sep="_")
  
  NSD<-c()
  SL<-c()
  
  
  selectedTraj<-readr::read_csv(behaviourFile)%>%    
    mutate(path_no=factor(new_path_no))%>%
    mutate(behaviour=ifelse(behaviour=="mc_1","in-patch", "in-matrix"))
  
  
  NSD<-c(NSD,selectedTraj%>%amt::make_track(.,x,y,t, crs=2056)%>%nsd()/1e6)
  SL<-c(SL,selectedTraj%>%amt::make_track(.,x,y,t, crs=2056)%>%step_lengths()/1000 )
  coordinates(selectedTraj)<-c("x","y")
  
  fullTrajStatList[[k]]<-data.frame(
    prenom=selectedDeer,
    deerYear=strsplit(selectedYear, split="_")[[1]][2],
    nbBehaviour=1,
    sumNSD=sum(NSD),
    meanNSD=mean(NSD),
    maxNSD=max(NSD, na.rm = T),
    sdNSD=sd(NSD, na.rm=T),
    sumSL=sum(SL, na.rm=T),
    meanSL=mean(SL, na.rm=T),
    sdSL=sd(SL, na.rm=T),
    STRA=as.numeric(selectedTraj%>%amt::make_track(.,x,y,t, crs=2056)%>%straightness()),
    UI=as.numeric(selectedTraj%>%amt::make_track(.,x,y,t, crs=2056)%>%intensity_use() ),
    MCP95=suppressWarnings(mcp(selectedTraj, percent=95, unin = "m", unout="ha")$area),
    TAC=as.numeric(selectedTraj%>%amt::make_track(.,x,y,t, crs=2056)%>%tac2() )
  )
  k<-k+1
} 

fullTrajStat_residents<-do.call(rbind, fullTrajStatList)


##### generate trajectory measures for migrants ######
setwd("/media/loreto/Grande/ie-ofev-24-25/cerf_movement_patterns/plateau_selected/bimodal_cerf_plateau")

fullTrajStatList<-list()
k<-1
for(behaviourFile in list.files(pattern=".*REFINED.*.csv")){
  
  selectedStepSize<-strsplit(behaviourFile, split=" ")[[1]][5]
  selectedDeer<-strsplit(behaviourFile, split=" ")[[1]][1]
  
  selectedYear<-paste(strsplit(behaviourFile, split=" ")[[1]][2], strsplit(behaviourFile, split=" ")[[1]][3], sep="_")
  
  NSD<-c()
  SL<-c()
  
  
  selectedTraj<-readr::read_csv(behaviourFile)%>%    
    mutate(path_no=factor(new_path_no))%>%
    mutate(behaviour=ifelse(behaviour=="mc_1","in-patch", "in-matrix"))
  
  
  NSD<-c(NSD,selectedTraj%>%amt::make_track(.,x,y,t, crs=2056)%>%nsd()/1e6)
  SL<-c(SL,selectedTraj%>%amt::make_track(.,x,y,t, crs=2056)%>%step_lengths()/1000 )
  coordinates(selectedTraj)<-c("x","y")
  fullTrajStatList[[k]]<-data.frame(
    prenom=selectedDeer,
    deerYear=strsplit(selectedYear, split="_")[[1]][2],
    nbBehaviour=2,
    sumNSD=sum(NSD),
    meanNSD=mean(NSD),
    maxNSD=max(NSD, na.rm = T),
    sdNSD=sd(NSD, na.rm=T),
    sumSL=sum(SL, na.rm=T),
    meanSL=mean(SL, na.rm=T),
    sdSL=sd(SL, na.rm=T),
    STRA=as.numeric(selectedTraj%>%amt::make_track(.,x,y,t, crs=2056)%>%straightness()),
    UI=as.numeric(selectedTraj%>%amt::make_track(.,x,y,t, crs=2056)%>%intensity_use() ),
    MCP95=suppressWarnings(mcp(selectedTraj, percent=95, unin = "m", unout="ha")$area),
    TAC=as.numeric(selectedTraj%>%amt::make_track(.,x,y,t, crs=2056)%>%tac2(),
                   MSD=as.numeric(selectedTraj%>%amt::make_track(.,x,y,t, crs=2056)%>%msd() )
    )
    
  )
  k<-k+1
} 

fullTrajStat_migrants<-do.call(rbind, fullTrajStatList)


#####  merge trajectory measures for residents (unimodal pattern) and migrants (bimodal pattern)

fullTrajStat<-rbind(fullTrajStat_migrants, fullTrajStat_residents)
fullTrajStat$key<-paste(fullTrajStat$prenom, fullTrajStat$deerYear)
rownames(fullTrajStat)<-fullTrajStat$key
fullTrajStat$nbBehaviour<-as.factor(fullTrajStat$nbBehaviour)

fviz_pca_biplot(acpTraj,
                habillage= fullTrajStat$nbBehaviour,
                col.var="black",
                addEllipses = T, ellipse.level=0.95, legend.title="nb Behaviour")+
  scale_color_manual(
    values= c("red","blue"))+theme(plot.title = element_blank())


###### Check relationship between chosen divider size (maxZscore) and trajectory measures #######

#read CSV with z-score from chosen step sizes after refinement
maxZscore<-read.csv("/media/loreto/Grande/ie-ofev-24-25/cerf_movement_patterns/plateau_selected/chosen_divider_sizes_chisq_zscores_upConf_plateau.csv",h=T)
maxZscore$key<-paste(maxZscore$prenom, unlist(lapply(maxZscore$deerYear, function(x) strsplit(x, "_")[[1]][2])), sep=" ")
names(maxZscore)
str(maxZscore)

#read CSV with best Migrate R behavior and support (manual classification between "weak" to "very_strong" support)
mR<-read.csv("/media/loreto/Grande/ie-ofev-24-25/cerf_movement_patterns/plateau_selected/results_migrateR_all_plateau_deers_add_manual_classification.csv", h=T)
mR$key<-paste(mR$prenom, unlist(lapply(mR$deerYear, function(x) strsplit(x, "_")[[1]][2])), sep=" ")

#add information on z-score and migrateR models

fullTrajStat_1<-merge(fullTrajStat, maxZscore, by="key", all.x=T)
fullTrajStat_2<-merge(fullTrajStat_1, mR, by="key", all.x=T)
fullTrajStat_2$nbBehaviour<-as.factor(fullTrajStat_2$nbBehaviour)
rownames(fullTrajStat_2)<-fullTrajStat_2$key

#compute only statistics of PHD Loreto
acpTraj<-PCA(fullTrajStat_2[,c("zscore", "maxNSD","STRA","UI","MCP95")], scale.unit = T, graph=F)

### PCA with nb Behaviour ####
  fviz_pca_biplot(acpTraj,
                  habillage=as.factor(fullTrajStat_2$nbBehaviour),
                  col.var="black", 
                  pointsize=2, 
                  pointshape=19,
                  addEllipses = T, repel = T, legend.title="nb Behaviour"
  )+scale_color_manual(values=c("red", "blue"))+theme(plot.title = element_blank())

### PCA with migrateR status ######
fviz_pca_biplot(acpTraj,
                habillage=as.factor(paste(fullTrajStat_2$status, fullTrajStat_2$support)),
                col.var="black",
                geom=c("point","text"),
                pointsize=4,
                pointshape=19,
                invisible="quali",
                legend.title="migrateR class"
  )+scale_color_manual(values=c("ambiguous weak" = "grey",
                              "disperser/nomad weak"="yellow",
                              "disperser/nomad strong"="orange",
                              "mixmig/migrant weak"="lightblue",
                              "mixmig/migrant strong"="blue",
                              "mixmig/migrant very_strong"="purple",
                              "mixmig/disperser very_strong"="lightgreen"))

### PCA with migrateR status only for bimodals ###### optionel ######
# fullTrajStat2_migrants<-subset(fullTrajStat_2, subset=nbBehaviour  %in% c(2))
# acpTraj_migrants<-PCA(fullTrajStat2_migrants[,c("zscore", "maxNSD","STRA","UI","MCP95")], scale.unit = T, graph=F)
# 
# fviz_pca_biplot(acpTraj_migrants,
#                 habillage=as.factor(paste(fullTrajStat2_migrants$status, fullTrajStat2_migrants$support)),
#                 col.var="black",
#                 geom=c("point","text"),
#                 pointsize=4,
#                 pointshape=19,
#                 invisible="quali",
#                 legend.title="migrateR class"
# )+scale_color_manual(values=c("ambiguous weak" = "grey",
#                               "disperser/nomad weak"="yellow",
#                               "disperser/nomad strong"="orange",
#                               "mixmig/migrant weak"="lightblue",
#                               "mixmig/migrant strong"="blue",
#                               "mixmig/migrant very_strong"="purple",
#                               "mixmig/disperser very_strong"="lightgreen"))
