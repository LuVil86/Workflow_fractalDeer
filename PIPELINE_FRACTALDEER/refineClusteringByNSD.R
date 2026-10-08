
refineClusteringByNSD<-function (behaviourFile,
                                 display.plot=FALSE,
                                 save.plot=FALSE, 
                                 save.data.frame=TRUE,
                                 save.BIC=FALSE,
                                 show.breakpoints=TRUE,
                                 parameter="ratio_meanNSD", return="Chisq" ){
  np<-parameter
  chosenParameter<-switch(parameter,
         ratio_cumulativeNSD="ratio_cumulativeNSD",
         speed="speed",
         ratio_endNSD="ratio_endNSD",
         ratio_maxNSD="ratio_maxNSD",
         ratio_medianNSD="ratio_medianNSD",
         ratio_meanNSD="ratio_meanNSD",
         ratio_varNSD="ratio_varNSD",
         sdNSD="sdNSD",
         STRA="STRA",
         IU="IU",
         MSD="MSD",
         covar="covar",
         optimum=c("ratio_endNSD", "ratio_cumulative")
         )
  if(is.null(chosenParameter)){stop("the parameter you specified does not exist : choose from list :
                                    ratio_cumulativeNSD,
                                    speed,
                                    ratio_endNSD,
                                    ratio_maxNSD,
                                    ratio_medianNSD,
                                    ratio_meanNSD,
                                    ratio_varNSD,
                                    sdNSD,
                                    STRA,
                                    IU,
                                    MSD,
                                    covar,
                                    optimum")}
  
  #################### load packages and functions ############################################
    suppressPackageStartupMessages({
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
  require(fs)
  })
  sourceDir <- function(path, trace = TRUE, ...) {
    op <- options(); on.exit(options(op)) # to reset after each 
    for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
      source(file.path(path, nm), ...)
      options(op)
    }
  }
  sourceDir("./migrateR_1.0.9/migrateR/R/")
  
  
  coordSplit<-function(x){
    tmp<-unlist(strsplit(x, split=" "))
    if(any(grepl("Z", tmp))){
      coordX<-as.numeric(substring(tmp[3],2,nchar(tmp[3])))
      coordY<-as.numeric(tmp[4])
    }
    else{
      coordX<-as.numeric(substring(tmp[2],2,nchar(tmp[2])))
      coordY<-as.numeric(substring(tmp[3],1,nchar(tmp[3])-1))
    }
    return(list("coordX"=coordX, "coordY"=coordY))
  }
  
  #####################################################################################################
  
  ########### load data ###############
  plotTitle<-strsplit(path_file(behaviourFile), split=".csv")[[1]][1]
  outDir=path_dir(behaviourFile)
  plotTitle<-gsub("_"," ",tail(strsplit(plotTitle, split="/")[[1]],n=1))
  animalTest<-strsplit(plotTitle, split=" ")[[1]][1]
  yearTest<-paste0("deerYear_",strsplit(plotTitle, split=" ")[[1]][3])
  stepSize<-strsplit(plotTitle, split=" ")[[1]][which(grepl("stepSize", strsplit(plotTitle, split=" ")[[1]]))+1]
  if(!grepl(".csv", behaviourFile)){stop("The file you provided does not have '.csv' extension")}
  if(grepl("REFINED", behaviourFile)){stop("The file you provided is already 'REFINED' : abort script")}
  
  datCerfs<-readr::read_csv(behaviourFile)%>%
    mutate(coord_X=apply(as.data.frame(geometry), 1,function(x) coordSplit(x)$coordX ))%>%
    mutate(coord_Y=apply(as.data.frame(geometry), 1,function(x) coordSplit(x)$coordY ))%>%
    mutate(dateTime=as.POSIXct(paste(UTC_DATE, UTC_TIME, sep=" "), origin="1970-01-01", tz="UTC"))%>%
    mutate(behaviour=factor(paste("b", behaviour, sep="_")))%>%
    mutate(deerYear=factor(deerYear))%>%
    mutate(prenom=factor(prenom))%>%
    mutate(saison=factor(saison))%>%
    mutate(jourNuit=factor(jourNuit))%>%
    mutate(path_no=factor(path_no))%>%
    dplyr::select(x="coord_X", y="coord_Y",t="dateTime", saison="saison",jourNuit="jourNuit",id="prenom",deerYear="deerYear","behaviour"=behaviour, "path_no"=path_no)%>%
  suppressWarnings(datCerfs<-datCerfs%>%transform(saison = fct_relevel(saison, c("Mars-Mai","Juin-Aout","Septembre-Novembre","Decembre-Fevrier"))))
  

  ################### COMPUTE SUB-PATH MEASURES FROM THE INPUT #################################33
  datCerfs$NSD_fullPath<-datCerfs%>%amt::make_track(.,x,y,t, crs=2056, all_cols = TRUE)%>%nsd()
  print("computing...")
  #####################################################################################################
  NSD<-c()
  sl_<-c()
  timeElapsed<-c()
  endNSD<-c()
  sumNSD<-c()
  maxNSD<-c()
  meanNSD<-c()
  medianNSD<-c()
  varNSD<-c()
  sdNSD<-c()
  
  duration<-c()
  startTime<-c()
  endTime<-c()
  STRA<-c()
  IU<-c()
  MSD<-c()
  
  sumSL<-c()
  meanSL<-c()
  maxSL<-c()
  varSL<-c()
  nFixes<-c()
  
  for(i in levels(datCerfs$path_no)) {
    sub<-subset(datCerfs, datCerfs$path_no==i)
    sub$NSD<-amt::make_track(sub,x,y,t, crs=2056)%>%nsd()
    sub$sl_<-amt::make_track(sub,x,y,t, crs=2056, all_cols = TRUE)%>%step_lengths()
    sub$NSD[which(is.na(sub$NSD))]<-0.00001
    nFixes<-c(nFixes,nrow(sub))
    startTime<-sub$t[1]
    endTime<-sub$t[nrow(sub)]
    timeElapsed<-c(timeElapsed,unlist(lapply(sub$t, function(x) difftime(x,startTime, units="hours"))))
    NSD<-c(NSD,sub$NSD)
    sl_<-c(sl_, sub$sl_)
    
    if(nrow(sub)>=3){
    endNSD<-c(endNSD,c(sub$NSD[nrow(sub)]-sub$NSD[1]))
    sumNSD<-c(sumNSD, sum(sub$NSD))
    maxNSD<-c(maxNSD,max(sub$NSD))
    medianNSD<-c(medianNSD,median(sub$NSD))
    meanNSD<-c(meanNSD, mean(sub$NSD))
    varNSD<-c(varNSD, var(sub$NSD, na.rm=T))
    sdNSD<-c(sdNSD, sd(sub$NSD))
    duration<-c(duration, difftime(sub$t[nrow(sub)], sub$t[1], units = "hours"))
    startTime<-c(startTime, sub$t[1])
    endTime<-c(endTime, sub$t[nrow(sub)])
    STRA<-c(STRA,amt::make_track(sub,x,y,t, crs=2056)%>%straightness())
    IU<-c(IU,amt::make_track(sub,x,y,t, crs=2056)%>%intensity_use())
    MSD<-c(MSD,amt::make_track(sub,x,y,t, crs=2056)%>%msd())
          
    
    
    sumSL<-c(sumSL, sum(sub$sl_, na.rm = TRUE))
    meanSL<-c(meanSL, mean(sub$sl_, na.rm=TRUE))
    varSL<-c(varSL, var(sub$sl_, na.rm = TRUE))
    maxSL<-c(maxSL, max(sub$sl_, na.rm = TRUE))
    }
    else{
      endNSD<-c(endNSD,0)
      sumNSD<-c(sumNSD,0)
      maxNSD<-c(maxNSD,0)
      meanNSD<-c(meanNSD,0)
      medianNSD<-c(medianNSD,0)
      varNSD<-c(varNSD,0)
      sdNSD<-c(sdNSD,0)
      meanSL<-c(meanSL,0)
      sumSL<-c(sumSL,0)
      maxSL<-c(maxSL, 0)
      varSL<-c(varSL, 0)
      duration<-c(duration,0)
      startTime<-c(startTime, 0)
      endTime<-c(endTime,0)
      STRA<-c(STRA,0)
      IU<-c(IU,0)
      MSD<-c(MSD,0)
    }
  }
  
  datCerfs$NSD<-NSD
  datCerfs$sl_<-sl_
  datCerfs$timeElapsed<-timeElapsed
  finDat<-data.frame(path_no=levels(datCerfs$path_no),
                     maxNSD,medianNSD,meanNSD, endNSD,sumNSD,varNSD,sdNSD,
                     meanSL,maxSL, sumSL, varSL,
                     duration,STRA,IU,MSD)
  
   finDat$ratio_cumulativeNSD<-finDat$sumNSD/(finDat$duration+1)
   finDat$ratio_endNSD<-finDat$endNSD/(finDat$duration+1)
   finDat$ratio_maxNSD<-finDat$maxNSD/(finDat$duration+1)
   finDat$ratio_medianNSD<-finDat$medianNSD/(finDat$duration+1)
   finDat$speed<-finDat$sumSL/(finDat$duration+1)
   finDat$ratio_meanNSD<-finDat$meanNSD/(finDat$duration+1)
   finDat$ratio_varNSD<-finDat$varNSD/(finDat$duration+1)
  
  
  
  finDat<-merge(finDat, datCerfs[,c("path_no", "behaviour")], by="path_no")
  finDat<-finDat[!duplicated(finDat),]
  finDat<-finDat[order(as.numeric(finDat$path_no)),]
  rownames(finDat)<-finDat$path_no
  
  distrib<-table(finDat$behaviour)
  
 ### I cannot obtain Chi-square statistic if there's only one path of a given behaviour....
    if(any(distrib==1)){
    Ttest_raw<-list(statistic=NA, p.value=NA)
    chisq_raw<-list(statistic=NA, p.value=NA)
  }else{
    sub<-subset(finDat, finDat$behaviour!="b_-1")
    sub<-drop.levels(sub)
    chisq_raw<-kruskal.test(sub[,parameter]~sub[,"behaviour"])  
    Ttest_raw<-t.test(sub[,parameter]~sub[,"behaviour"])
  }
  
  
  
  
    pl1<-ggplot(aes(x=t, y=NSD), data=datCerfs)+
    geom_point(aes(colour=behaviour))+
    ggtitle("before refining")+
    theme_bw()+
    ylab("Net-Squared Displacement [km²]")+
    scale_x_datetime(date_labels = "%b %y", date_breaks = "1 month")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.text.y=element_text(),
          axis.title.x=element_blank(),
          legend.title = element_blank(), legend.position = "top",plot.title = element_textbox_simple(halign = 0.5))
  
  if(show.breakpoints==TRUE){
    breakpoints<-c()
    tmp<-as.data.frame(datCerfs)
    
    for(i in 2:nrow(datCerfs)){
      if(tmp$behaviour[i-1]!=tmp$behaviour[i]){
        breakpoints<-c(breakpoints, i)
      }
     
    }
    tBreakPoints<-datCerfs$t[breakpoints]
    ymax <- ggplot_build(pl1)$layout$panel_params[[1]]$y.range[2]
    pl1<-pl1+geom_vline(xintercept = tBreakPoints, linetype="dotted")
    
    rm(tmp)
  }
  pl2<-ggplot(aes(x=as.factor(as.numeric(path_no)), y=!!sym(chosenParameter),colour= behaviour), data=finDat)+
    geom_point()+
    theme_bw()+
    scale_colour_manual(values=c("black","green", "red"))+theme(legend.position = "top", legend.title = element_blank())
  pl3<-ggplot(aes(x=behaviour, y=!!sym(chosenParameter)), data=finDat)+geom_boxplot()
  #################### ********* CLUSTERING ********* ######################################
  
  
 #subFinDat<-finDat%>%filter(behaviour%in%c("b_0", "b_1"))%>%droplevels()%>%select(-path_no)
  
#  plot(lda(behaviour~sl_ , data=datCerfs))
  
  
 #   ggplot(aes(x=maxSL, y=maxNSD, colour=behaviour), data=finDat)+geom_point()+geom_label(label=finDat$path_no)
  #
  
  #chosenParameter<-c("timeDiff","meanSL")
  # fit<-Mclust(finDat[,chosenParameter], G=2)
  # finDat$clustering<-apply(fit$z,1, function(x){
  #   ifelse(x[1]>0.52,"mc_1","mc_0") }
  # )
  
  BIC <- mclustBIC(finDat[,chosenParameter], G=2)
  mod1 <- Mclust(finDat[,chosenParameter], x = BIC)
  
  if (save.BIC==TRUE){
    if(!file.exists("summary_BIC.csv")){
      sink(file="summary_BIC.csv")
      cat("prenom,deerYear,stepSize,BIC.1,BIC.2\n")
      sink()
    }else{
    sink(file="summary_BIC.csv", append=T)
    cat(animalTest,",",yearTest,",",stepSize, ",",mod1$BIC[1,1],",",mod1$BIC[1,2],"\n" )
    sink()}
  }
  finDat$clustering<-paste0("mc_",mod1$classification)
  #hc<-hclust(dist(finDat[,c("maxSL", "sumSL")]), method="complete")
  #finDat$clustering<-cutree(hc, h=45000000)
  
  
  
  pl5<-ggplot(aes(x=as.factor(as.numeric(path_no)), y=!!sym(chosenParameter),colour= clustering), data=finDat)+
    geom_point()+
    theme_bw()+scale_colour_manual(values=c("green", "red"))+theme(legend.position = "top", legend.title = element_blank())
  pl6<-ggplot(aes(x=clustering, y=!!sym(chosenParameter)), data=finDat)+geom_boxplot()
  
  
  ############# CHANGE OLD BEHAVIOUR BY NEW BEHAVIOUR #####################
  
  
  datCerfs<-merge(datCerfs, finDat[,c("path_no", "clustering")], by="path_no")
  datCerfs<-datCerfs%>%dplyr::select(-c("behaviour"))%>%dplyr::rename("behaviour"="clustering")%>%dplyr::arrange(t)
  
    ###### change path_no according to new clustering ######
  path_merged<-c(0)
  newPathNo<-0
  breakpoints<-c()
  tmp<-as.data.frame(datCerfs)
  for(i in 2:nrow(tmp)){
    if(tmp$behaviour[i-1]!=tmp$behaviour[i]){
      newPathNo<-newPathNo+1
      breakpoints<-c(breakpoints, i)
    }
    path_merged<-c(path_merged, newPathNo)
  }
  rm(tmp)
  datCerfs<-datCerfs%>%mutate(new_path_no=path_merged)
  
  
  ########## RECOMPUTE SUB-PATH MEASURES FOR THE NEW SUB-PATHS (after clustering)  ##############
  NSD<-c()
  sl_<-c()
  timeElapsed<-c()
  endNSD<-c()
  sumNSD<-c()
  maxNSD<-c()
  meanNSD<-c()
  medianNSD<-c()
  varNSD<-c()
  sdNSD<-c()
  
  duration<-c()
  startTime<-c()
  endTime<-c()
  STRA<-c()
  IU<-c()
  MSD<-c()
  sumSL<-c()
  meanSL<-c()
  maxSL<-c()
  varSL<-c()
  nFixes<-c()
  
  for(i in levels(as.factor(datCerfs$new_path_no))) {
    sub<-subset(datCerfs, datCerfs$new_path_no==i)
    sub$NSD<-amt::make_track(sub,x,y,t, crs=2056)%>%nsd()
    sub$sl_<-amt::make_track(sub,x,y,t, crs=2056, all_cols = TRUE)%>%step_lengths()
    sub$NSD[which(is.na(sub$NSD))]<-0.00001
    nFixes<-c(nFixes,nrow(sub))
    startTime<-c(startTime,sub$t[1])
    endTime<-c(endTime,sub$t[nrow(sub)])
    timeElapsed<-c(timeElapsed,unlist(lapply(sub$t, function(x) difftime(x,sub$t[1], units="hours"))))
    
    NSD<-c(NSD,sub$NSD)
    sl_<-c(sl_, sub$sl_)
    if(nrow(sub)>=3){
      endNSD<-c(endNSD,c(sub$NSD[nrow(sub)]-sub$NSD[1]))
      sumNSD<-c(sumNSD, sum(sub$NSD))
      maxNSD<-c(maxNSD,max(sub$NSD))
      medianNSD<-c(medianNSD,median(sub$NSD))
      meanNSD<-c(meanNSD, mean(sub$NSD))
      varNSD<-c(varNSD, var(sub$NSD, na.rm=T))
      sdNSD<-c(sdNSD, sd(sub$NSD))
      duration<-c(duration, difftime(sub$t[nrow(sub)], sub$t[1], units = "hours"))   
      STRA<-c(STRA,amt::make_track(sub,x,y,t, crs=2056)%>%straightness())
      IU<-c(IU, amt::make_track(sub,x,y,t, crs=2056)%>%intensity_use())
      MSD<-c(MSD,amt::make_track(sub,x,y,t, crs=2056)%>%msd())
      
      meanSL<-c(meanSL, mean(sub$sl_, na.rm=TRUE))
      sumSL<-c(sumSL, sum(sub$sl_, na.rm = TRUE))
      maxSL<-c(maxSL, max(sub$sl_, na.rm = TRUE))
      varSL<-c(varSL, var(sub$sl_, na.rm = TRUE))
    }
    else{
      endNSD<-c(endNSD,0)
      sumNSD<-c(sumNSD,0)
      maxNSD<-c(maxNSD,0)
      medianNSD<-c(medianNSD,0)
      meanNSD<-c(meanNSD,0)
      varNSD<-c(varNSD,0)
      sdNSD<-c(sdNSD, 0)
      duration<-c(duration,0)
      STRA<-c(STRA,0)
      IU<-c(IU,0)
      MSD<-c(MSD,0)
      meanSL<-c(meanSL,0)
      sumSL<-c(sumSL,0)
      maxSL<-c(maxSL, 0)
      varSL<-c(varSL, 0)
    }
  }
  datCerfs$NSD<-NSD
  datCerfs$sl_<-sl_
  datCerfs$timeElapsed<-timeElapsed
  finDat<-data.frame(new_path_no=levels(as.factor(datCerfs$new_path_no)),
                     startTime=as.POSIXct(startTime), 
                     endTime=as.POSIXct(endTime), duration,STRA,IU,MSD,
                     maxNSD,medianNSD,meanNSD, endNSD,sumNSD,varNSD,sdNSD,meanSL,
                     maxSL, sumSL, varSL)
  
  finDat$speed<-finDat$sumSL/(finDat$duration+1)
  
  finDat$ratio_cumulativeNSD<-finDat$sumNSD/(finDat$duration+1)
  finDat$ratio_endNSD<-finDat$endNSD/(finDat$duration+1)
  finDat$ratio_maxNSD<-finDat$maxNSD/(finDat$duration+1)
  finDat$ratio_medianNSD<-finDat$medianNSD/(finDat$duration+1)
  finDat$ratio_meanNSD<-finDat$meanNSD/(finDat$duration+1)
  finDat$ratio_varNSD<-finDat$varNSD/(finDat$duration+1)
  
  
  finDat<-merge(finDat, datCerfs[,c("new_path_no", "behaviour")], by="new_path_no")
  finDat<-finDat[!duplicated(finDat),]
  finDat<-finDat[order(as.numeric(finDat$new_path_no)),]
  rownames(finDat)<-finDat$new_path_no
 
  
  distrib<-table(finDat$behaviour)
  if(any(distrib==1)){
    Ttest_refined<-list(statistic=NA, p.value=NA)
    chisq_refined<-list(statistic=NA, p.value=NA)
  }else{
    
    Ttest_refined<-t.test(finDat[,parameter]~finDat[,"behaviour"])
    chisq_refined<-kruskal.test(finDat[,parameter]~finDat[,"behaviour"])  
  }
  
  testStats<-data.frame(prenom=animalTest, 
             deerYear=yearTest,
             stepSize=as.numeric(stepSize),
             chisq.raw=as.numeric(chisq_raw$statistic),
             chisq.raw.pval=as.numeric(chisq_raw$p.value),
             Ttest.raw=as.numeric(Ttest_raw$statistic),
             Ttest.raw.pval=as.numeric(Ttest_raw$p.value),
             chisq.refined=as.numeric(chisq_refined$statistic),
             chisq.refined.pval=as.numeric(chisq_refined$p.value),
             Ttest.refined=as.numeric(Ttest_refined$statistic),
             Ttest.refined.pval=as.numeric(Ttest_refined$p.value)
  )
  
  
  #  corresp_pathNo<-datCerfs[,c("path_no", "new_path_no")]
  # corresp_pathNo<-corresp_pathNo[!duplicated(corresp_pathNo),]
  # 
  #  mapPathNo<-function(x){   
  #    res<-corresp_pathNo$new_path_no[corresp_pathNo$path_no==as.character(x)]
  #    if(length(res)==0){
  #    return(NA)
  #    }else{return(res)}
  #  }
  #  corresp_clustering<-datCerfs[,c("new_path_no", "behaviour")]
  #  corresp_clustering<-corresp_clustering[!duplicated(corresp_clustering),]
  # # print(head(corresp_clustering))
  #  mapClustering<-function(x){
  #    res<-corresp_clustering$behaviour[corresp_clustering$new_path_no==x]
  #    if(length(res)==0){return(NA)}else{
  #    return(res)}
  #  }
  # # 
  #  rebuildDF$new_path_no<-unlist(lapply(rebuildDF$path_no, mapPathNo))
  # # 
  #  rebuildDF$clustering<-unlist(lapply(rebuildDF$new_path_no, mapClustering))
  #  for(i in 1:length(levels(as.factor(rebuildDF$new_path_no)))){
  #  }
   # pl0<-ggplot(aes(x=timeElapsed, y=NSD/1e6, group=new_path_no), data=rebuildDF)+
   #   geom_line()+
   #   theme_bw()+
   #   xlab("Duration [hrs]")
   #   facet_wrap(~as.factor(clustering))
   #   
  
  ############## PLOT NSD vs Time elapsed for new path number ######################################################
     pl0_mc_1<-ggplot(aes(x=timeElapsed, y=NSD/1e6, group=as.factor(new_path_no)), data=datCerfs[datCerfs$behaviour=="mc_1",])+
       geom_line(linetype="dashed", linewidth=0.5, colour="grey")+
       geom_point(size=1, colour="black")+
       scale_x_continuous(minor_breaks = seq(0,3000,100))+
       theme_bw()+
       ylab("Net square displacement [km²]")+
       xlab("Duration [hrs]")+
       theme(panel.grid.major.x = element_line(linewidth = 2))

     pl0_mc_2<-ggplot(aes(x=timeElapsed, y=NSD/1e6, group=as.factor(new_path_no)), data=datCerfs[datCerfs$behaviour=="mc_2",])+
       geom_line(linetype="dashed", linewidth=0.5, colour="grey")+
       geom_point(size=1,colour="black")+
       scale_x_continuous(minor_breaks = seq(0,100,20))+
       theme_bw()+
       ylab("Net square displacement [km²]")+
       xlab("Duration [hrs]")+
       theme(panel.grid.major.x = element_line(linewidth = 2))
       
     
     #grid.arrange(pl0_mc_1,pl0_mc_2, nrow=2)
     #
  ################################################################
  
  pl4<-ggplot(aes(x=t, y=NSD), data=datCerfs)+
    geom_point(aes(colour=behaviour))+
    ggtitle("after refining")+
    theme_bw()+
    ylab("Net-Squared Displacement [km²]")+
    scale_x_datetime(date_labels = "%b %y", date_breaks = "1 month")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.text.y=element_text(),
          axis.title.x=element_blank(),
          legend.title = element_blank(), legend.position = "top",plot.title = element_textbox_simple(halign = 0.5))
  if(show.breakpoints==TRUE){

    tBreakPoints<-datCerfs$t[breakpoints]
    ymax <- ggplot_build(pl3)$layout$panel_params[[1]]$y.range[2]
    pl4<-pl4+geom_vline(xintercept = tBreakPoints, linetype="dotted")
    
  }
    if (save.plot==TRUE){
    print("trying to save plot...");
    fp<-file.path(getwd(),animalTest,paste("REFINEMENT_PLOTS ",toupper(np)," ", plotTitle,".png",sep=""), fsep = "/" )
    print(fp)
    ggsave(filename =fp,
           plot = pl4,
           device = "png")
    }
  if (save.data.frame==TRUE){
    datCerfs%>%dplyr::select(-c("NSD"))%>%write_csv(file=file.path(outDir,paste(plotTitle," REFINED ",toupper(np),".csv",sep=""),fsep = "/"))
    
  }
  if(display.plot==TRUE){
  gridExtra::grid.arrange(pl1,pl4,pl3,pl6,pl2,pl5,layout_matrix=rbind(c(1,1,3), c(2,2,4), c(5,6,NA)), 
                          top=textGrob(paste(animalTest, yearTest, stepSize," : ", toupper(np), sep=" "),gp=gpar(fontsize=14,font=3)))
    }
  #datCerfs%>%dplyr::select(-c("NSD"))%>%return()
  if(return=="Chisq"){
  return(testStats)
  }else{
    datCerfs%>%dplyr::select(-c("NSD"))%>%return()
  }
}
