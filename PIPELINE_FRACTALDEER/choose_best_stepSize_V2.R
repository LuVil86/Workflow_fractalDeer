
###################### @@@@ choose_best_stepSize.r ########################

# Author : lucas villard
# Last revised : jan 28 2025
# Usage : The script computes parametric (t-test) and non-parametric ( kruskal-wallis) test on a trajectory metric defined by the user to estimate
# 	  how well a stepSize trajectory classification explains the observed trajectory 
#         

### ---Input---
# - the output table generated from the "get_Zscores_CRW.R" script
# - the metric you want to compute statistics on. by default, "ratio_meanNSD"

### --output ---
## a .csv table with the calculated statistics

#############################################################################



###################################################################
############### * PACKAGES, FUNCTIONS AND WORKDIR * #########################
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
  require(snow)
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

############ *** input parameters parsing ***** ################

cat("\n****** CHOOSE_BEST_STEPSIZES.R ********\n")

cerf_path<-commandArgs(trailingOnly = T)
inFolder=cerf_path[1]
tmpPath=unlist(strsplit(inFolder, "/"))
#print(tmpPath)
selectedDeer<-strsplit(tmpPath[length(tmpPath)], split="_")[[1]][2]
selectedYear<-paste(strsplit(tmpPath[length(tmpPath)], split="_")[[1]][3],strsplit(tmpPath[length(tmpPath)], split="_")[[1]][4], sep="_")


stepSizeListName=cerf_path[2]

###############################################################33


cat("- input file : ",stepSizeListName,"\n")

################# *************************************************************####################################################
################ ****** COMPUTE CHI-SQUARE FROM STEPSIZE LIST ****** ########
################# *************************************************************####################################################


options(scipen=999)


stepSizeList<-read.csv(stepSizeListName, h=T)


endNSDlist<-list()
meanNSDlist<-list()
cumulativeNSDlist<-list()


for(k in 1:nrow(stepSizeList)){
  
  selectedStepSize<-round(stepSizeList$stepSize[k])
  selectedDeer<-stepSizeList$prenom[k]
  selectedYear<-stepSizeList$deerYear[k]
  behaviourFile=paste0(inFolder,"/",selectedDeer,"_",selectedYear,"_stepSize_",selectedStepSize,"_autocorrelation_timeSeriesKmeans_2_classes.csv")
  cat("\n******* ", selectedDeer, " :: ", selectedYear, " --> ", selectedStepSize, "  ******\n")
  if(!file.exists(behaviourFile)){print(paste("the file with stepSize ",stepSizeList$stepSize[k],"has not been generated.. skipping it"));next}
  
  
  datCerfs<-readr::read_csv(behaviourFile,show_col_types = FALSE)%>%
    mutate(coord_X=apply(as.data.frame(geometry), 1,function(x) coordSplit(x)$coordX ))%>%
    mutate(coord_Y=apply(as.data.frame(geometry), 1,function(x) coordSplit(x)$coordY ))%>%
    mutate(dateTime=as.POSIXct(paste(UTC_DATE, UTC_TIME, sep=" "), origin="1970-01-01", tz="UTC"))%>%
    mutate(behaviour=factor(paste("b", behaviour, sep="_")))%>%
    mutate(deerYear=factor(deerYear))%>%
    mutate(prenom=factor(prenom))%>%
    mutate(path_no=factor(path_no))%>%
    dplyr::select(x="coord_X", y="coord_Y",t="dateTime",id="prenom",deerYear="deerYear","behaviour"=behaviour, "path_no"=path_no)
  
  
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
    chisq_raw_endNSD<-kruskal.test(sub[,"ratio_endNSD"]~sub[,"behaviour"])  
    chisq_raw_meanNSD<-kruskal.test(sub[,"ratio_meanNSD"]~sub[,"behaviour"])  
    chisq_raw_cumulativeNSD<-kruskal.test(sub[,"ratio_cumulativeNSD"]~sub[,"behaviour"])  

  }
  
  
  endNSDlist[[k]]<-data.frame(prenom=selectedDeer, 
             deerYear=selectedYear,
             stepSize=as.numeric(selectedStepSize),
             chisq.raw=as.numeric(chisq_raw_endNSD$statistic),
             chisq.raw.pval=as.numeric(chisq_raw_endNSD$p.value),
             fullPath=behaviourFile
 			 ) 
 			 
  meanNSDlist[[k]]<-data.frame(prenom=selectedDeer, 
             deerYear=selectedYear,
             stepSize=as.numeric(selectedStepSize),
             chisq.raw=as.numeric(chisq_raw_meanNSD$statistic),
             chisq.raw.pval=as.numeric(chisq_raw_meanNSD$p.value),
             fullPath=behaviourFile
 			 ) 			 

  cumulativeNSDlist[[k]]<-data.frame(prenom=selectedDeer, 
             deerYear=selectedYear,
             stepSize=as.numeric(selectedStepSize),
             chisq.raw=as.numeric(chisq_raw_cumulativeNSD$statistic),
             chisq.raw.pval=as.numeric(chisq_raw_cumulativeNSD$p.value),
             fullPath=behaviourFile
 			 )			 
  
  
  
  
  
  
  
  
  
  
  
  
}



finChisq_endNSD<-do.call(rbind,endNSDlist)
finChisq_meanNSD<-do.call(rbind,meanNSDlist)
finChisq_cumulativeNSD<-do.call(rbind,cumulativeNSDlist)



tot1<-cbind(finChisq_endNSD, zscore=stepSizeList$zscore)
tot1<-tot1[with(tot1, order(-zscore, -chisq.raw)), ]


tot2<-cbind(finChisq_meanNSD, zscore=stepSizeList$zscore)
tot2<-tot2[with(tot2, order(-zscore, -chisq.raw)), ]


tot3<-cbind(finChisq_cumulativeNSD, zscore=stepSizeList$zscore)
tot3<-tot3[with(tot3, order(-zscore, -chisq.raw)), ]

outNameSplit<-unlist(strsplit(path_file(stepSizeListName), split="_"))
outSuffix<-paste(selectedDeer, selectedYear, outNameSplit[6], outNameSplit[7]
  ,sep="_")
write.csv(tot1, file=paste0(inFolder,"/results_Chisq_endNSD_",outSuffix) , row.names=F)
write.csv(tot2, file=paste0(inFolder,"/results_Chisq_meanNSD_",outSuffix) , row.names=F)
write.csv(tot3, file=paste0(inFolder,"/results_Chisq_cumulativeNSD_",outSuffix) , row.names=F)





