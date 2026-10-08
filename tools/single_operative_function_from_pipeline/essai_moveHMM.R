library(moveHMM)
library(raster)
library(lubridate)
library(tidyverse)
library(maptools)
library(sf)
library(stringr)
library(ggtext)
library(janitor)
library(landscapemetrics)
library(vegan)
library(FactoMineR)
library(factoextra)
library(gdata)
library(amt)

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
################ example from moveHMM-starting_values #######



hmmdataExample <- prepData(haggis_data, type = "UTM")
plot(hmmdataExample)


resultFolder="/media/luvil/T7_ROUGE/results_grisons/"

#### **** ICI  copier le contenu de la colonne "fullPath" de la ligne de stepSize à analyser ***** #####
fullPath="/home/luvil/results_grisons/OUT_20083_deerYear_2012-2013/20083_deerYear_2012-2013_stepSize_898_autocorrelation_timeSeriesKmeans_2_classes.csv"

#### *** ICI la mesure à prendre en compte : "ratio_endNSD", "ratio_meanNSD" ou "ratio_cumulativeNSD" ***** #######
parameter="ratio_meanNSD"
#parameter="ratio_endNSD"
#parameter="ratio_cumulativeNSD"

################################################################################################################################



inBehaviourFile=paste0(resultFolder,str_split_fixed(fullPath, '/', 5)[1,5])
dirPath=path_dir(inBehaviourFile)
fileName=path_file(inBehaviourFile)

stepCerfs<-readr::read_csv(stepFile)%>%
  mutate(x=apply(as.data.frame(geometry), 1,function(x) coordSplit(x)$coordX ))%>%
  mutate(y=apply(as.data.frame(geometry), 1,function(x) coordSplit(x)$coordY ))%>%
  mutate(t=c(0:c(nrow(.)-1)))%>%
  mutate(id=factor(selectedDeer))

trackStep<-amt::make_track(stepCerfs,x,y,t, crs=2056,all_cols = TRUE)%>%
  dplyr::select(ID="id",Easting="x_", Northing="y_")%>%as.data.frame()


datCerfs<-readr::read_csv(inBehaviourFile)%>%
  mutate(coord_X=apply(as.data.frame(geometry), 1,function(x) coordSplit(x)$coordX ))%>%
  mutate(coord_Y=apply(as.data.frame(geometry), 1,function(x) coordSplit(x)$coordY ))%>%
  mutate(dateTime=as.POSIXct(paste(UTC_DATE, UTC_TIME, sep=" "), origin="1970-01-01", tz="UTC"))%>%
  mutate(behaviour=factor(paste("b", behaviour, sep="_")))%>%
  #slice(stepCerfs$nearestPoint)%>%
  mutate(deerYear=factor(deerYear))%>%
  mutate(saison=factor(saison))%>%
  mutate(prenom=factor(prenom))%>%
  mutate(path_no=factor(path_no))%>%
  mutate(jourNuit=factor(jourNuit))%>%
  dplyr::select(x="coord_X", y="coord_Y",t="dateTime", id="prenom",deerYear="deerYear",saison="saison","behaviour"=behaviour,tod="jourNuit",  "path_no"=path_no)%>%
  arrange(t)%>%
  amt::make_track(.,x,y,t, crs=2056,all_cols = TRUE)%>%
  track_resample(rate = minutes(60), tolerance = minutes(4))%>%
  dplyr::select(ID="id",Easting="x_", Northing="y_")%>%as.data.frame()




#datCerfs$Easting <- datCerfs$Easting/1000
#datCerfs$Northing <- datCerfs$Northing/1000
hmmdata <- moveHMM::prepData(datCerfs,type="UTM",coordNames=c("Easting","Northing"))

hist(hmmdata$step, xlab="step length", breaks=50)
hist(hmmdata$angle, breaks = seq(-pi, pi, length = 15), xlab = "angle", main = "")

length(which(hmmdata==0))/nrow(hmmdata)

# For reproducibility
set.seed(12345)
# Number of tries with different starting values
niter <- 10
# Save list of fitted models
minP<-c(1000,1500)  ### x,y
maxP<-c(1500, 2500)
allm <- list()
for(i in 1:niter) {
  # Step length mean
  stepMean0 <- runif(2, min = minP,max = maxP) ### min(x,y) et max(x,y) ou x correspond au state 1 et y correspond au state 2
  #stepMean0<-c(0.1,0.5)
  # Step length standard deviation
   stepSD0 <- runif(2,  min = minP,max = maxP)
 # stepSD0<-c(0.1,0.5)
  stepPar0<-c(stepMean0,stepSD0)
  # Turning angle mean
  angleMean0 <- c(pi, 0)
  angleCon0 <- c(1, 10)
  #angleCon0 <- runif(2,min = c(0.5, 2), max = c(1, 4))
  
  anglePar0 <- c(angleMean0, angleCon0)
  
  # zeromass0 <- c(0.05,0.1) # step zero-mass
  
  # Turning angle concentration
  # angleCon0 <- runif(2,
  #                    min = c(0.5, 5),
  #                    max = c(2, 15))
  # Fit model
  allm[[i]] <- fitHMM(data = hmmdata, nbStates = 2, stepPar0 = stepPar0,
                      anglePar0 = anglePar0, stationary=T)
  print(i)
}

allnllk <- unlist(lapply(allm, function(m) m$mod$minimum))
allnllk
whichbest <- which.min(allnllk)
mbest <- allm[[whichbest]]
mbest
plot(mbest ,plotCI=TRUE)
#### output 
vit<-viterbi(mbest)
datCerfs%>%mutate(state=vit)%>%write_csv(file=paste("HMM_results_",unique(.$ID),sep=""))
