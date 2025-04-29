setwd("/media/loreto/Grande/R_lore/")

##############################################################################
#### USAGE :
#     run the "mvmtclass" function from migrateR for a given animal and deer-year. print the results in a separate file
#    - here, the input is a list of red deers and deer years.
#     this script is a copy of the script get_values_migrateR_fit.R in FractalDeer elaborated by Lucas Villard
##############################################################################


require(multidplyr)
require(tidyverse)
require(lubridate)
require(adehabitatLT)
require(arulesViz)
require(gdata)
require(ggplot2)
require(ggtext)
sourceDir <- function(path, trace = TRUE, ...) {
  op <- options(); on.exit(options(op)) # to reset after each 
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    source(file.path(path, nm), ...)
    options(op)
  }
}
#sourceDir("./migrateR_1.0.9/migrateR/R/")
sourceDir("/media/loreto/Grande/R_lore/migrateR_1.0.9/migrateR/R/")
options(readr.show_col_types = FALSE)

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


compVector<-c("disperser","migrant", "mixmig", "nomad", "resident")


AIClist<-list()
is_refinment_file<-TRUE

folderPath="/media/loreto/Grande/ie-ofev-24-25/cerf_movement_patterns/plateau_selected/bimodal_cerf_plateau/"
for (file in list.files(path=folderPath)){


behaviourFile=paste0(folderPath,file )
plotTitle<-strsplit(behaviourFile, split=".csv")[[1]][1]
plotTitle<-gsub("_"," ",tail(strsplit(plotTitle, split="/")[[1]],n=1))
print(plotTitle)
selectedDeer<-strsplit(plotTitle, split=" ")[[1]][1]
selectedYear<-paste(strsplit(plotTitle, split=" ")[[1]][2], strsplit(plotTitle, split=" ")[[1]][3], sep="_")
selectedStepSize<-strsplit(plotTitle, split=" ")[[1]][5]
plotTitle_for_classical<-paste(strsplit(plotTitle, split=" ")[[1]][1:3], collapse = " ")
if(is_refinment_file==FALSE){
  datCerfs<-readr::read_csv(behaviourFile)%>%
    mutate(coord_X=apply(as.data.frame(geometry), 1,function(x) coordSplit(x)$coordX ))%>%
    mutate(coord_Y=apply(as.data.frame(geometry), 1,function(x) coordSplit(x)$coordY ))%>%
    mutate(dateTime=as.POSIXct(paste(UTC_DATE, UTC_TIME, sep=" "), origin="1970-01-01", tz="UTC"))%>%
    mutate(behaviour=factor(paste("b", behaviour, sep="_")))%>%
    mutate(deerYear=factor(deerYear))%>%
    mutate(prenom=factor(prenom))%>%
    mutate(saison=factor(saison))%>%
    mutate(path_no=factor(path_no))%>%
    dplyr::select(x="coord_X", y="coord_Y",t="dateTime", id="prenom",deerYear="deerYear",saison="saison","behaviour"=behaviour, "path_no"=path_no)
  suppressWarnings(datCerfs<-datCerfs%>%transform(saison = fct_relevel(saison, c("Mars-Mai","Juin-Aout","Septembre-Novembre","Decembre-Fevrier"))))
}else{
  datCerfs<-readr::read_csv(behaviourFile)%>%    
    mutate(path_no=factor(new_path_no))%>%
    mutate(behaviour=ifelse(behaviour=="mc_1","in-patch", "in-matrix"))%>%
    suppressWarnings(datCerfs<-datCerfs%>%transform(saison = fct_relevel(saison, c("Mars-Mai","Juin-Aout","Septembre-Novembre","Decembre-Fevrier"))))
  
}

RD_traj<-as.ltraj(xy=datCerfs[,c("x","y")], date=datCerfs$t, id=datCerfs$id)
suppressWarnings(RD.nsd1 <- mvmtClass(RD_traj))

aicres <- sapply(RD.nsd1[[selectedDeer]]@models, AIC)
missingComp<-compVector[which(!compVector%in%names(aicres))]
toAdd<-c(rep(NA, length(missingComp)))
names(toAdd)<-missingComp
fin<-c(aicres, toAdd)
AIClist[[file]]<- fin[order(names(fin))]
}
#create dataframe with AICs values 
finDat<-data.frame(do.call(rbind, AIClist))

#calculate delta AIC values
minAIC<-do.call(pmin, c(finDat, na.rm=TRUE))
DatDeltaAIC<- sweep(finDat, 1, minAIC, "-")  #1:row 2:columns
#prepare results to export
id <- data.frame(do.call(rbind,strsplit(as.character(rownames(DatDeltaAIC)),"\\s+")))
id <- id %>% unite("deerYear",X2:X3, na.rm=TRUE, sep="_") %>% mutate(deerYear)
finDatDeltaAIC <- DatDeltaAIC%>% mutate(prenom=id$X1, deerYear=id$deerYear, stepSize=id$X5) %>% relocate (prenom:stepSize, .before = 1)
#export results
write.csv(finDatDeltaAIC,file="/media/loreto/Grande/ie-ofev-24-25/cerf_movement_patterns/plateau_selected/results_migrateR_bi_plateau_deers.csv", row.names=FALSE)
