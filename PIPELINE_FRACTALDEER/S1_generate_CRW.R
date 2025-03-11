

###################### @@@@ generate_CRW.r ########################

# Author : lucas villard
# Last revised : jan 28 2025
# Usage : the script create a correlated random walk  (CRW) dataframe based on the caracteristics of the input file. 
# 	  It outputs a CRW based on mean location of the entire GPS point series and uses mean and standard-deviation of the step lengths as the center 
#         of the step length distribution for the CRW

### ---Input---
# .csv file by deerYear generated with the "format_GPS_data.py" script. (specify deer name and deerYear in the input parameters section)

### --output ---
## a CRW on the form "Simulated_CRW_[animalName]_[deerYear].csv

###################################################################
#####################################
suppressPackageStartupMessages({
  
require(trajr)
require(gdata)
require(tidyverse)
require(fs)
})
#############################

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

cerf_path=commandArgs(trailingOnly = T)
inFile=path_file(cerf_path[1])
selectedDeer<-strsplit(inFile[1], split="_")[[1]][1]
selectedYear<-path_ext_remove(paste(strsplit(inFile[1], split="_")[[1]][2],strsplit(inFile[1], split="_")[[1]][3], sep="_"))
cat("***** SIMULATE CRW TRAJECTORY ****\n")
cat("\n deer Name --> ",selectedDeer, ", deer Year --> ", selectedYear, "\n")

trajCerf<-readr::read_csv(cerf_path[1],show_col_types = FALSE)%>%
  mutate(coord_X=apply(as.data.frame(geometry), 1,function(x) coordSplit(x)$coordX ))%>%
  mutate(coord_Y=apply(as.data.frame(geometry), 1,function(x) coordSplit(x)$coordY ))%>%
  mutate(dateTime=as.POSIXct(paste(UTC_DATE, UTC_TIME, sep=" "), origin="1970-01-01", tz="UTC"))%>%
  mutate(deerYear=factor(deerYear))%>%
  mutate(saison=factor(saison))%>%
  mutate(prenom=factor(prenom))%>%
  mutate(jourNuit=factor(jourNuit))%>%
  #filter(deerYear==selectedYear)%>%
  droplevels()%>%
  TrajFromCoords(xCol="coord_X", yCol="coord_Y",spatialUnits = "m")

print(exp(mean(log(TrajStepLengths(trajCerf))[log(TrajStepLengths(trajCerf))!=-Inf])))
print( sd(TrajStepLengths(trajCerf)))
#hist(TrajStepLengths(trajCerf),breaks = 100,main = " histogram of step_length [m]")
## generate correlated random walk from stepLength mean :

### ORIGINAL DEFINITION #####
# trjSim<-TrajGenerate(n=nrow(trajCerf)-1, 
#                       random=T,
#                       stepLength = mean(TrajStepLengths(trajCerf)),
#                       linearErrorSd = sd(TrajStepLengths(trajCerf))
#                       )
###########################
trjSim <- TrajGenerate(n=nrow(trajCerf)-1,
                       random=T,
                       linearErrorDist = function(n) rlnorm(n)*median(TrajStepLengths(trajCerf))
)

# trjSim <- TrajGenerate(n=nrow(trajCerf)-1, 
#                        random=T,
#                        linearErrorDist = rlnorm)

#plot(trjSim)
#hist(TrajStepLengths(TrajFromCoords(trjSim,xCol="x", yCol="y",spatialUnits = "m")),breaks = 100,main = " histogram of simulated step length [m]")

trjSim_scaled<-data.frame(prenom=paste(selectedDeer,"_CRW",sep=""),
                          deerYear=selectedYear,
                          xcoord=trjSim$x+trajCerf$x[1],
                          ycoord=trjSim$y+trajCerf$y[1],
                          dateTime=trajCerf$dateTime)
outFile=paste0(cerf_path[2],"/Simulated_CRW_",selectedDeer,"_", selectedYear,".csv")

write.csv(trjSim_scaled, file=outFile)
