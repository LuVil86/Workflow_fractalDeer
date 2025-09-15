

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

############# create random steps along paths by behaviour for many animals ############# 
source("./create_random_steps_along_paths_by_behaviour.R")
behaviourFolder="/media/luvil/NAS_DEVELOPPEMENT/IE_OFEV/valais_selected/Bimodal"


priority<-terra::rast("/media/luvil/NAS_DEVELOPPEMENT/IE_OFEV/csc_cerf_alps/mosaic_20A_test_vs/mosaic_20A_test_vs.tif")

### Create random steps along paths if chosenBehaviour="in-matrix" then long and directed displacements else if chosenBehaviour="in-habitat" then short and tortuous displacements ###

if(exists("rndSteps")){rm(rndSteps)}
if(exists("output")){rm(output)}
output<-list()
for(f in list.files(path=behaviourFolder,full.names = TRUE,pattern=".csv")){
  cat("\n***********************************\n")
  print(f)
  cat("**************************************\n")
  output[[f]]<-create_random_step_along_path_by_behaviour(behaviourFile = f, chosenBehaviour="in-matrix",
                                                          nbRandom=100)
  
}
rndSteps<-do.call(rbind,output)

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


rndSteps$priority<-extractCovariates(priority,rndSteps, covarExtractionType = "end") 

ggplot(aes(x=animalID,y=priority, fill=case_), data=rndSteps)+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1), legend.position="none")



issfPrio<-amt::fit_issf(Loc~-1+priority+strata(new_step_id_), data=rndSteps)
summary(issfPrio)
summISSFpriority<-as.data.frame(broom::tidy(issfPrio$model, conf.int=TRUE))
