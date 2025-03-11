
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

argmax <- function(x, y, w=1, ...) {
  require(zoo)
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, 
                     align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}

test <- function(w, span) {
  peaks <- argmax(x, y, w=w, span=span)
  
  plot(x, y, cex=0.75, col="Gray", main=paste("w = ", w, ", 
              span = ", span, sep=""))
  lines(x, peaks$y.hat,  lwd=2) #$
  y.min <- min(y)
  sapply(peaks$i, function(i) lines(c(x[i],x[i]), c(y.min, 
                                                    peaks$y.hat[i]),
                                    col="Red", lty=2))
  points(x[peaks$i], peaks$y.hat[peaks$i], col="Red", pch=19, 
         cex=1.25)
  
}

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

if(length(cerf_path)>2){
  parameter=cerf_path[3]
}else{
  parameter="ratio_meanNSD"
}
###############################################################33


cat("- input file : ",stepSizeListName,"\n")
cat("- measure to evaluate : ", parameter,"\n")

################# *************************************************************####################################################
################ ****** COMPUTE CHI-SQUARE FROM STEPSIZE LIST ****** ########
################# *************************************************************####################################################


source("./refineClusteringByNSD_workflow.R")


options(scipen=999)


chisqList<-list()
stepSizeList<-read.csv(stepSizeListName, h=T)

#### cluster for parallel processing 

clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 8), type = clusterType))


behaviourFileList<-list()
for(k in 1:nrow(stepSizeList)){
  
  selectedStepSize<-round(stepSizeList$stepSize[k])
  selectedDeer<-stepSizeList$prenom[k]
  selectedYear<-stepSizeList$deerYear[k]
  behaviourFile=paste0(inFolder,"/",selectedDeer,"_",selectedYear,"_stepSize_",selectedStepSize,"_autocorrelation_timeSeriesKmeans_2_classes.csv")
  cat("\n******* ", selectedDeer, " :: ", selectedYear, " --> ", selectedStepSize, "  ******\n")
  if(!file.exists(behaviourFile)){print(paste("the file with stepSize ",stepSizeList$stepSize[k],"has not been generated.. skipping it"));next}
  
  #cat("\n******* ",behaviourFile," ******\n")
  behaviourFileList[[k]]<-behaviourFile
}


chisqList<-clusterApply(clust,behaviourFileList,refineClusteringByNSD_workflow,parameter=parameter)


finChisq<-do.call(rbind,chisqList)



tot<-cbind(finChisq, zscore=stepSizeList$zscore)
tot<-tot[with(tot, order(-zscore, -chisq.raw)), ]

outNameSplit<-unlist(strsplit(path_file(stepSizeListName), split="_"))
outSuffix<-paste(selectedDeer, selectedYear, outNameSplit[6], outNameSplit[7]
  ,sep="_")
write.csv(tot, file=paste0(inFolder,"/results_Chisq_raw_refined_",parameter,"_",outSuffix) , row.names=F)


finChisq$key1<-paste(finChisq$prenom, finChisq$deerYear, sep="_")
ggplot(aes(x=stepSize, y=chisq.raw, colour=key1), data=finChisq)+geom_point()+geom_line()




