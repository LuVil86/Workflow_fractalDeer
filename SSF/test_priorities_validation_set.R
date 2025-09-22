

knitr::opts_chunk$set(warning = FALSE, message = FALSE, cache = TRUE)
#' ## Load libraries and read in data
#+ warning=FALSE, message=FALSE
library(survival)
library(TwoStepCLogit)
#library(INLA)
library(glmmTMB)
library(tidyverse)
library(dplyr)
library(vegan)
library(ggplot2)
library(ggtext)
library(gridExtra)
library(terra)
library(sp)
library(amt)
library(psych)
library(GGally)
library(janitor)
library(reshape2)
#option R- function where "width" control de maximum number of columns
options(width=150)

############ functions ############# 
corClass<-function(x, id_column){
  if(is.na(x)){return(NA)}
  else if(!x%in%corresp_pix_classes[,id_column]){return(NA)}
  else{
    return(corresp_pix_classes[which(corresp_pix_classes[,id_column]==x), "new_short_class"][1]) #new_short_class is the new name of the column of the habitat correspondences file
  }
}
printEquation<-function(params){
  f_fix<-paste("Loc ~ -1+",paste(params, collapse="+"),"+ sl_ + (1|new_step_id_)" )
  f_rnd<-unlist(lapply(params, FUN=function(x){paste("(0 + ",x,"| ANIMAL_ID)", sep="")}))
  f_final<-paste(f_fix,paste(f_rnd, collapse="+"), sep="+")
  cat("fixed only : \n\n")
  cat(f_fix, ",\n\n")
  cat("fixed and random :\n\n")
  cat(f_final, ",\n\n")
}

### "extractCovariates" function extract raster values to points
extractCovariates<-function(raster, dfLoc, covarExtractionType="begin-end"){
  x<-c(rep(NA, nrow(dfLoc)))
  if(covarExtractionType=="begin-end"){
    x[which(dfLoc$case_==TRUE)]<-terra::extract(raster,vect(SpatialPoints(dfLoc[which(dfLoc$case_==TRUE),c("x1_","y1_")])))[,2]
    x[which(dfLoc$case_==FALSE)]<-terra::extract(raster,vect(SpatialPoints(dfLoc[which(dfLoc$case_==FALSE),c("x2_","y2_")])))[,2]
  }else if(covarExtractionType=="end"){
    x<-terra::extract(raster,vect(SpatialPoints(dfLoc[,c("x2_","y2_")])))[,2]
  }else{stop("choose between 'end' and 'begin-end' in the parameter 'covarExtractionType")}
  return(x)
}


###### create the data frame of observed points ######
behaviourFolder="/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/valais_selected/Bimodal/"
### load GPS data #####
if(exists("output")){rm(output)}
output<-list()
for(f in list.files(path=behaviourFolder,full.names = TRUE,pattern=".csv")){
  cat("\n***********************************\n")
  print(f)
  cat("**************************************\n")
  output[[f]]<-readr::read_csv(f)%>%
    mutate(behaviour=ifelse(behaviour=="mc_1","in-patch", "in-matrix"))%>%
    filter(behaviour%in%c("in-matrix"))
}

obsPoints<-vect(do.call(rbind,output), geom=c("x", "y"),crs="epsg:2056")


##### create mask for Switzerland only #######

prio<-rast("/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Priorities/Prio_cerf_alpes20A_250915_10NB.tif")
CHMask<-ifel(prio>=0,1,NA )

qStep=0.1


#######################################################################################
###### **** EXTRACT CURRENT VALUES FROM GPS FIXES ***** ###############################
#######################################################################################


### NOTE : the current values are transformed in 1-100th quantiles to standardize the different models in order to compare them


### load and extract current values ####
current<-terra::rast("/media/luvil/NAS_DEVELOPPEMENT/IE_OFEV/csc_cerf_alps/mosaic_20A_test_vs/mosaic_20A_test_vs.tif")
current_masked<-mask(current, CHMask)
classVect<- global(current_masked, quantile, probs=seq(0, 1, by=qStep), na.rm=T)
currQuant<-classify(current_masked, t(as.matrix(classVect)))
obsPoints$current<-as.numeric(terra::extract(currQuant,obsPoints)[,2] )
obsDF<-melt(obsPoints$current, value.name="current_quantile", na.rm = T)
obsDF$model<-"20A"
g1<-ggplot(aes(x=model,y=as.factor(current_quantile)),data=obsDF)+
  geom_jitter(col="firebrick1")+
  theme_bw()+
  ggtitle("Model 20A")



########
current<-terra::rast("/media/luvil/NAS_DEVELOPPEMENT/IE_OFEV/csc_cerf_alps/mosaic_20B_test_vs/mosaic_20B_test_vs.tif")
current_masked<-mask(current, CHMask)
classVect<- global(current_masked, quantile, probs=seq(0, 1, by=qStep), na.rm=T)
currQuant<-classify(current_masked, t(as.matrix(classVect)))
obsPoints$current<-as.numeric(terra::extract(currQuant,obsPoints)[,2] )
obsDF<-melt(obsPoints$current, value.name="current_quantile", na.rm = T)
obsDF$model<-"20B"
g2<-ggplot(aes(x=model,y=as.factor(current_quantile)),data=obsDF)+
  geom_jitter(col="chartreuse2")+
  theme_bw()+
  ggtitle("Model 20B")

########
current<-terra::rast("/media/luvil/NAS_DEVELOPPEMENT/IE_OFEV/csc_cerf_alps/mosaic_21_test_vs/mosaic_21_test_vs.tif")
classVect<- global(current, quantile, probs=seq(0, 1, by=qStep), na.rm=T)
currQuant<-classify(current, t(as.matrix(classVect)))
obsPoints$current<-as.numeric(terra::extract(currQuant,obsPoints)[,2] )
obsDF<-melt(obsPoints$current, value.name="current_quantile", na.rm = T)
obsDF$model<-"21"
g3<-ggplot(aes(x=model,y=as.factor(current_quantile)),data=obsDF)+
  geom_jitter(col="lightblue")+
  theme_bw()+
  ggtitle("Model 21")


#### boxplots for the threee models side-by-side
grid.arrange(g1,g2,g3, nrow=1)



#######################################################################################
###### **** EXTRACT Priority VALUES FROM GPS FIXES ***** ###############################
#######################################################################################


qStep=0.1
### load and extract priority values ####
priority<-terra::rast("/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Priorities/Prio_cerf_alpes20A_250915.tif")
#priority_masked<-mask(priority, CHMask)
classVect<- global(priority, quantile, probs=seq(0, 1, by=qStep), na.rm=T)
quantPrio<-classify(priority, t(as.matrix(classVect)), include.lowest=TRUE)
obsPoints$priority<-as.numeric(as.factor(terra::extract(quantPrio,obsPoints)[,2]))
obsDF<-melt(obsPoints$priority, value.name="priority_values", na.rm = T)
obsDF$model<-"20A"
p1<-ggplot(aes(x=model,y=as.factor(priority_values)),data=obsDF)+
  geom_jitter(col="lightblue")+
  theme_bw()+
  scale_y_discrete(limits=factor(c(1:10)))+
  ggtitle("Model 20A")

### load and extract priority values ####
priority<-terra::rast("/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Priorities/Prio_cerf_alpes20B_250915.tif")
classVect<- global(priority, quantile, probs=seq(0, 1, by=qStep), na.rm=T)
quantPrio<-classify(priority, t(as.matrix(classVect)), include.lowest=TRUE)
obsPoints$priority<-as.numeric(as.factor(terra::extract(quantPrio,obsPoints)[,2]))
obsDF<-melt(obsPoints$priority, value.name="priority_values", na.rm = T)
obsDF$model<-"20B"
p2<-ggplot(aes(x=model,y=as.factor(priority_values)),data=obsDF)+
  geom_jitter(col="chartreuse2")+
  theme_bw()+
  scale_y_discrete(limits=factor(c(1:10)))+
  ggtitle("Model 20B")


### load and extract priority values ####
priority<-terra::rast("/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Priorities/Prio_cerf_alpes21_250915.tif")
classVect<- global(priority, quantile, probs=seq(0, 1, by=qStep), na.rm=T)
quantPrio<-classify(priority, t(as.matrix(classVect)), include.lowest=TRUE)
obsPoints$priority<-as.numeric(as.factor(terra::extract(quantPrio,obsPoints)[,2]))
obsDF<-melt(obsPoints$priority, value.name="priority_values", na.rm = T)
obsDF$model<-"21"
p3<-ggplot(aes(x=model,y=as.factor(priority_values)),data=obsDF)+
  geom_jitter(col="firebrick1")+
  theme_bw()+
  scale_y_discrete(limits=factor(c(1:10)))+
  ggtitle("Model 21")

grid.arrange(p1,p2,p3, nrow=1)



############## **** KERNEL TECHNIQUE ***** ############################
### -NOTE : the idea is to use the kernel density of GPS fixes for a given animal
# and use it as a "weighted probability" raster to sample points which would constitue the "available" space

########################################################################################


### this is an example for a single animal
library(move)
library(adehabitatHR)
library(raster)
library(spatialEco)


tmp<-readr::read_csv("/media/luvil/NAS_DEVELOPPEMENT/IE_OFEV/valais_selected/Bimodal/ID003 deerYear 2018-2019 stepSize 980 autocorrelation timeSeriesKmeans 2 classes REFINED RATIO_MEANNSD.csv")%>%
  mutate(behaviour=ifelse(behaviour=="mc_1","in-patch", "in-matrix"))
obsTraj<-move(x=tmp$x, y=tmp$y, time=tmp$t,proj = CRS("epsg:2056"))


### need to use the "sp" package for kernelUD 
coordinates(tmp) <- c("x", "y")
proj4string(tmp) <- CRS("epsg:2056")

### the raster of priorities has to be cropped to the extend and readable as "grid" for kernelUD
priorityRAST<-raster("/media/luvil/NAS_DEVELOPPEMENT/IE_OFEV/Priorities/Prio_cerf_alpes21_250915.tif")
trim<-crop(priorityRAST, extent(tmp)) ### crop to the extent
spDF<-as(trim, "SpatialPixels")
kern<-kernelUD(tmp, grid=spDF)

## retransform in "terra" format
kern2<-rast(as(kern, "SpatialPixelsDataFrame"))

#### classify in quantiles to create weights
qStep=0.01
qtVect<- global(kern2, quantile, probs=seq(0, 1, by=qStep), na.rm=T)
qt<-classify(kern2, t(as.matrix(qtVect)), include.lowest=TRUE)
poids<-raster.invert(qt/100)

### sample from weights
hs<-res(poids)/2
ptscell = sample(1:ncell(poids), 1000, prob=poids[], replace=TRUE)
centres = xyFromCell(poids,ptscell)
pts = cbind(runif(nrow(centres),centres[,1]-hs[1],centres[,1]+hs[1]),runif(nrow(centres),centres[,2]-hs[2],centres[,2]+hs[2]))

### display the results (for display purposes only)
plot(poids)
points(pts)
plot(obsTraj, add=T)


##### extract values

sampleValues<-extract(trim, pts)
tmp_inMatrix<-subset(tmp, tmp$behaviour=="in-matrix")
obsValues<-extract(trim,tmp_inMatrix )

#### plot the result (boxplots) #####

DF_for_plot<-data.frame(type=c(rep("obs",length(obsValues)), rep("rnd", length(sampleValues))), 
                               value = c(obsValues, sampleValues))

ggplot(aes(x=type, y=value),data=DF_for_plot)+geom_boxplot()+theme_bw()+ggtitle("priority values : observed VS random")+xlab("priority value")
t.test(value~type, data=DF_for_plot)

### write the example (for display purposes only)
writeRaster(poids, filename = "test_weighted_raster.tif", overwrite=T)
writeVector(vect(pts), filename = "test_sample_weighted.shp",overwrite=T)


############## **** RANDOM STEP VALIDATION TECHNIQUE ***** ############################
### -NOTE : the creation of random steps are almost always in the vicinity of observed data,
##          which "bias" the true availiability of current since the values obtained at random steps
##          still fall in the corridors
########################################################################################
source("./create_random_steps_along_paths_by_behaviour.R")
behaviourFolder="/media/luvil/NAS_DEVELOPPEMENT/IE_OFEV/valais_selected/Bimodal"

if(exists("rndSteps")){rm(rndSteps)}
if(exists("output")){rm(output)}
output<-list()
for(f in list.files(path=behaviourFolder,full.names = TRUE,pattern=".csv")){
  cat("\n***********************************\n")
  print(f)
  cat("**************************************\n")
  output[[f]]<-create_random_step_along_path_by_behaviour(behaviourFile = f, chosenBehaviour="in-matrix",
                                                          nbRandom=10)
  
}
rndSteps<-do.call(rbind,output)



if(exists("rndSteps")){rm(rndSteps)}
if(exists("output")){rm(output)}
output<-list()
for(f in list.files(path=behaviourFolder,full.names = TRUE,pattern=".csv")){
  cat("\n***********************************\n")
  print(f)
  cat("**************************************\n")
  output[[f]]<-read.csv(f, h=T)%>%filter(behaviour%in%"in-matrix")
}

obsPoints<-do.call(rbind,output)

### remove output after calculation because it is not needed after ###
rm(output)

### reformat the table of used (TRUE) /available (FALSE) points "rndSteps" to continue with the script ####
### - also reindex step_IDs to join every animal trajectories sequentially #####
rndSteps$Loc<-as.integer(rndSteps$case_)
rndSteps$step_id_<-as.integer(rndSteps$step_id_)

newStepID<-rep(NA, nrow(rndSteps))
stepIDNo<-1
newStepID[1]<-stepIDNo
for(i in 2:nrow(rndSteps)){
  if(rndSteps$step_id_[i-1]!=rndSteps$step_id_[i]){
    stepIDNo<-stepIDNo+1
  }
  newStepID[i]<-stepIDNo
}
rndSteps$new_step_id_<-newStepID
rm(newStepID)
rndSteps$numeric_animalID<-as.numeric(as.factor(rndSteps$animalID))




### scaling values of step length and turning angle
rndSteps$sl_<-scale(rndSteps$sl_)
rndSteps$ta_<-scale(rndSteps$ta_)
### eliminating values at final coordinate of steps
rndSteps<-rndSteps%>%filter(!is.na(rndSteps$x2_))

writeVector(vect(rndSteps[,c("case_", "x2_", "y2_")], geom=c("x2_", "y2_"), crs="epsg:2056"), filename = "./Random_steps_projection.geojson", overwrite=T)



############################# BIOMOD dev ###########################################################

require(biomod2)
bmData<-BIOMOD_FormatingData(resp.var=case_,
                             expl.var = allRast,
                             resp.xy=geom(allPoints)[,c("x", "y")],
                             resp.name="habitat suitability",
                             filter.raster = TRUE)

#writeVector(obsPoints,filename="model_20A_quantile_obsPoints.geojson",overwrite=T)


