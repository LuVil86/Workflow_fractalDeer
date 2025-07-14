############ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #####################
#########   RUN GLMM_TMB on behaviour pathsPATHS    ##########################
##### - By Loreto Urbina et Lucas Villard
#### --------------------------------------- ################
#+ include = FALSE
#'  
#'-This code replicates the analysis presented in Muff, Signer, Fieberg (2019) Section 4.2 "Habitat selection of otters: an SSF analysis".
#'
#### --------------------------------------- ################


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
behaviourFolder="/media/loreto/Grande/ie-ofev-24-25/cerf_movement_patterns/alpes_selected/bimodal_pattern"
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
newStepID<-c(1)
stepIDNo<-1
for(i in 2:nrow(rndSteps)){
  if(rndSteps$step_id_[i-1]!=rndSteps$step_id_[i]){
    stepIDNo<-stepIDNo+1
  }
  newStepID<-c(newStepID, stepIDNo)
}
rndSteps$new_step_id_<-newStepID
rm(newStepID)
rndSteps$numeric_animalID<-as.numeric(as.factor(rndSteps$animalID))


### scaling values of step length and turning angle
rndSteps$sl_<-scale(rndSteps$sl_)
rndSteps$ta_<-scale(rndSteps$ta_)
### eliminating values at final coordinate of steps
rndSteps<-rndSteps%>%filter(!is.na(rndSteps$x2_))
### selecting end coordinate of steps to do the srsf
exType<-"end"




##### step resource selection function using habitat categories and a reference class and continuous variables ##################################################

##select source folder
rasterFolder="/media/loreto/Grande/ie-ofev-24-25/variables"

##### habitat raster file + correspondence  ##################################################

##select source map
habitatRasterPath=paste0(rasterFolder,"/habitat_cerf_15avr25/HabitatMap_cerf_15avril25.tif")
##create habitat raster object
land_use<-terra::rast(habitatRasterPath)
##create a data frame of correspondences
habitatCorrespondence=paste0(rasterFolder,"/habitat_cerf_15avr25/habitat_classification_cerf_15avril25.csv")
corresp_pix_classes<-read.csv(habitatCorrespondence,h=T)

##extract land use values ##
rndSteps$landUse<-extractCovariates(land_use,rndSteps, covarExtractionType = exType) ###corClass is a function to add at the beginning of the script
rndSteps<-rndSteps%>%mutate(landUse=as.factor(purrr::map_chr(landUse, corClass, id_column="new_raster_value"))) ### i supposed id_column is the raster value in the correspondances files
tLandUse<-matrix(unlist(lapply(rndSteps$landUse, function(x) table(x))), ncol=length(levels(rndSteps$landUse)), byrow=T)
colnames(tLandUse)<-levels(rndSteps$landUse) ###tLandUse is a table of used=1/available=0 units by landuse 

## see used and available ##
t_use<-table(rndSteps$animalID[rndSteps$Loc==1], rndSteps$landUse[rndSteps$Loc==1])
t_avail<-table(rndSteps$animalID[rndSteps$Loc==0], rndSteps$landUse[rndSteps$Loc==0])
print(" ####### USED HABITATS ######")
print(t_use)
print(" \n####### AVAILABLE HABITATS ######")
print(t_avail)
## writing table of used units and available units by habitat class
##plateau
# write.csv(t_avail, file="/media/loreto/Grande/ie-ofev-24-25/ssf_plateau/available_habitatRL_nbRandom100.csv")
# write.csv(t_use, file="/media/loreto/Grande/ie-ofev-24-25/ssf_plateau/used_habitatRL_nbRandom100.csv")
##alps
write.csv(t_avail, file="/media/loreto/Grande/ie-ofev-24-25/ssf_alps/available_habitatRL_nbRandom100.csv")
write.csv(t_use, file="/media/loreto/Grande/ie-ofev-24-25/ssf_alps/used_habitatRL_nbRandom100.csv")

## eliminate used units in unlikely habitats such as in Auto, Bât, Goudr, PassInfM, Rivie classes ##
#toRemove<-rndSteps[which(rndSteps$Loc==1 & rndSteps$landUse%in%c("Auto", "Bât", "Goudr", "PassInfM", "Rivie")),]$new_step_id_#plateau
toRemove<-rndSteps[which(rndSteps$Loc==1 & rndSteps$landUse%in%c("Bât", "Lacs","Goudr", "Ouvr","PassInfM", "ParRo", "PassInfM", "Rivie")),]$new_step_id_#alps
rndStepsF<-subset(rndSteps, !rndSteps$new_step_id_%in%toRemove)
tF_use<-table(rndStepsF$animalID[rndStepsF$Loc==1], rndStepsF$landUse[rndStepsF$Loc==1])
tF_use
tF_avail<-table(rndStepsF$animalID[rndStepsF$Loc==0], rndStepsF$landUse[rndStepsF$Loc==0])
tF_avail
## writing table of used units and available units without unlikely units by habitat class
##plateau
#write.csv(tF_use, file="/media/loreto/Grande/ie-ofev-24-25/ssf_plateau/used_habitatRL_nbRandom10_F.csv")#to do if necessary by charging the correspondent file *.Rdata for plateau
#write.csv(tF_use, file="/media/loreto/Grande/ie-ofev-24-25/ssf_plateau/used_habitatRL_nbRandom10_F.csv")#to do if necessary by charging the correspondent file *.Rdata for plateau
##alps
write.csv(tF_use, file="/media/loreto/Grande/ie-ofev-24-25/ssf_alps/used_habitatRL_nbRandom100_F.csv")

##### adding log scaled variables of densities #####

##  build density files from scaled log values ##

buildDensPath<-c(paste0(rasterFolder, "/human/log_scaled_density/Density_Buildings_50_opt2_log+1_scaled.tif"))
buildDens50_2<-terra::rast(buildDensPath)
rndStepsF$buildDens50_2<-extractCovariates(buildDens50_2, rndStepsF, covarExtractionType = exType)

buildDensPath<-c(paste0(rasterFolder, "/human/log_scaled_density/Density_Buildings_100_opt2_log+1_scaled.tif"))
buildDens100_2<-terra::rast(buildDensPath)
rndStepsF$buildDens100_2<-extractCovariates(buildDens100_2, rndStepsF, covarExtractionType = exType)

buildDensPath<-c(paste0(rasterFolder, "/human/log_scaled_density/Density_Buildings_200_opt2_log+1_scaled.tif"))
buildDens200_2<-terra::rast(buildDensPath)
rndStepsF$buildDens200_2<-extractCovariates(buildDens200_2, rndStepsF, covarExtractionType = exType)

buildDensPath<-c(paste0(rasterFolder, "/human/log_scaled_density/Density_Buildings_400_opt2_log+1_scaled.tif"))
buildDens400_2<-terra::rast(buildDensPath)
rndStepsF$buildDens400_2<-extractCovariates(buildDens400_2, rndStepsF, covarExtractionType = exType)

## main road density file from scaled log values ##

MroadDistPath<-c(paste0(rasterFolder,"/human/log_scaled_density/Density_Merge_RoadPrimary__50_opt_log+1_scaled.tif"))
MroadDens50_2<-terra::rast(MroadDistPath)
rndStepsF$MroadDens50_2<-extractCovariates(MroadDens50_2, rndStepsF, covarExtractionType = exType)

MroadDistPath<-c(paste0(rasterFolder,"/human/log_scaled_density/Density_Merge_RoadPrimary__100_opt_log+1_scaled.tif")) #
MroadDens100_2<-terra::rast(MroadDistPath)
rndStepsF$MroadDens100_2<-extractCovariates(MroadDens100_2, rndStepsF, covarExtractionType = exType)

MroadDistPath<-c(paste0(rasterFolder,"/human/log_scaled_density/Density_Merge_RoadPrimary__200_opt_log+1_scaled.tif")) #
MroadDens200_2<-terra::rast(MroadDistPath)
rndStepsF$MroadDens200_2<-extractCovariates(MroadDens200_2, rndStepsF, covarExtractionType = exType)

MroadDistPath<-c(paste0(rasterFolder,"/human/log_scaled_density/Density_Merge_RoadPrimary__400_opt_log+1_scaled.tif")) #
MroadDens400_2<-terra::rast(MroadDistPath)
rndStepsF$MroadDens400_2<-extractCovariates(MroadDens400_2, rndStepsF, covarExtractionType = exType)

## secondary road density file from scaled log values ##

SroadDistPath<-c(paste0(rasterFolder,"/human/log_scaled_density/Density_Merge_RoadSecondary__50_opt_log+1_scaled.tif"))
SroadDens50_2<-terra::rast(SroadDistPath)
rndStepsF$SroadDens50_2<-extractCovariates(SroadDens50_2, rndStepsF,covarExtractionType = exType)

SroadDistPath<-c(paste0(rasterFolder,"/human/log_scaled_density/Density_Merge_RoadSecondary__100_opt_log+1_scaled.tif"))
SroadDens100_2<-terra::rast(SroadDistPath)
rndStepsF$SroadDens100_2<-extractCovariates(SroadDens100_2, rndStepsF, covarExtractionType = exType)

SroadDistPath<-c(paste0(rasterFolder,"/human/log_scaled_density/Density_Merge_RoadSecondary__200_opt_log+1_scaled.tif"))
SroadDens200_2<-terra::rast(SroadDistPath)
rndStepsF$SroadDens200_2<-extractCovariates(SroadDens200_2, rndStepsF, covarExtractionType = exType)

SroadDistPath<-c(paste0(rasterFolder,"/human/log_scaled_density/Density_Merge_RoadSecondary_400_opt_log+1_scaled.tif"))
SroadDens400_2<-terra::rast(SroadDistPath)
rndStepsF$SroadDens400_2<-extractCovariates(SroadDens400_2, rndStepsF, covarExtractionType = exType)

## forest scaled log density file ##

forDensPath<-c(paste0(rasterFolder, "/habitat/log_scaled_density/Density_Forest_50_opt2_log+1_scaled.tif"))
forDens50_2<-terra::rast(forDensPath)
rndStepsF$forDens50_2<-extractCovariates(forDens50_2, rndStepsF, covarExtractionType = exType)

forDensPath<-c(paste0(rasterFolder, "/habitat/log_scaled_density/Density_Forest_100_opt2_log+1_scaled.tif"))
forDens100_2<-terra::rast(forDensPath)
rndStepsF$forDens100_2<-extractCovariates(forDens100_2, rndStepsF, covarExtractionType = exType)

forDensPath<-c(paste0(rasterFolder, "/habitat/log_scaled_density/Density_Forest_200_opt2_log+1_scaled.tif"))
forDens200_2<-terra::rast(forDensPath)
rndStepsF$forDens200_2<-extractCovariates(forDens200_2, rndStepsF, covarExtractionType = exType)

forDensPath<-c(paste0(rasterFolder, "/habitat/log_scaled_density/Density_Forest_400_opt2_log+1_scaled.tif"))
forDens400_2<-terra::rast(forDensPath)
rndStepsF$forDens400_2<-extractCovariates(forDens400_2, rndStepsF, covarExtractionType = exType)

## highway scaled log density file ##

AutoDensPath<-c(paste0(rasterFolder, "/human/log_scaled_density/Density_Merge_Autobahn__50_opt_log+1_scaled.tif"))
AutoDens50_2<-terra::rast(AutoDensPath)
rndStepsF$AutoDens50_2<-extractCovariates(AutoDens50_2, rndStepsF, covarExtractionType = exType)

AutoDensPath<-c(paste0(rasterFolder, "/human/log_scaled_density/Density_Merge_Autobahn__100_opt_log+1_scaled.tif"))
AutoDens100_2<-terra::rast(AutoDensPath)
rndStepsF$AutoDens100_2<-extractCovariates(AutoDens100_2, rndStepsF, covarExtractionType = exType)

AutoDensPath<-c(paste0(rasterFolder, "/human/log_scaled_density/Density_Merge_Autobahn__200_opt_log+1_scaled.tif"))
AutoDens200_2<-terra::rast(AutoDensPath)
rndStepsF$AutoDens200_2<-extractCovariates(AutoDens200_2, rndStepsF, covarExtractionType = exType)

AutoDensPath<-c(paste0(rasterFolder, "/human/log_scaled_density/Density_Merge_Autobahn__400_opt_log+1_scaled.tif"))
AutoDens400_2<-terra::rast(AutoDensPath)
rndStepsF$AutoDens400_2<-extractCovariates(AutoDens400_2, rndStepsF, covarExtractionType = exType)


##### select the better resolution scaled log density variables #####

## select building density resolution ##

issfBati50_2<-amt::fit_issf(Loc~-1+buildDens50_2+strata(new_step_id_), data=rndStepsF)
AIC(issfBati50_2)
issfBati100_2<-amt::fit_issf(Loc~-1+buildDens100_2+strata(new_step_id_), data=rndStepsF)
AIC(issfBati100_2)
issfBati200_2<-amt::fit_issf(Loc~-1+buildDens200_2+strata(new_step_id_), data=rndStepsF)
AIC(issfBati200_2)
issfBati400_2<-amt::fit_issf(Loc~-1+buildDens400_2+strata(new_step_id_), data=rndStepsF)
AIC(issfBati400_2)

## select main road density resolution ##
issfMroad50_2<-amt::fit_issf(Loc~-1+MroadDens50_2+strata(new_step_id_), data=rndStepsF)
AIC(issfMroad50_2)
issfMroad100_2<-amt::fit_issf(Loc~-1+MroadDens100_2+strata(new_step_id_), data=rndStepsF)
AIC(issfMroad100_2)
issfMroad200_2<-amt::fit_issf(Loc~-1+MroadDens200_2+strata(new_step_id_), data=rndStepsF)
AIC(issfMroad200_2)
issfMroad400_2<-amt::fit_issf(Loc~-1+MroadDens400_2+strata(new_step_id_), data=rndStepsF)
AIC(issfMroad400_2)

## select secondary road density resolution ##
issfSroad50_2<-amt::fit_issf(Loc~-1+SroadDens50_2+strata(new_step_id_), data=rndStepsF)
AIC(issfSroad50_2)
issfSroad100_2<-amt::fit_issf(Loc~-1+SroadDens100_2+strata(new_step_id_), data=rndStepsF)
AIC(issfSroad100_2)
issfSroad200_2<-amt::fit_issf(Loc~-1+SroadDens200_2+strata(new_step_id_), data=rndStepsF)
AIC(issfSroad200_2)
issfSroad400_2<-amt::fit_issf(Loc~-1+SroadDens400_2+strata(new_step_id_), data=rndStepsF)
AIC(issfSroad400_2)

## select forest density resolution ##
issfFor50_2<-amt::fit_issf(Loc~-1+forDens50_2+strata(new_step_id_), data=rndStepsF)
AIC(issfFor50_2)
issfFor100_2<-amt::fit_issf(Loc~-1+forDens100_2+strata(new_step_id_), data=rndStepsF)
AIC(issfFor100_2)
issfFor200_2<-amt::fit_issf(Loc~-1+forDens200_2+strata(new_step_id_), data=rndStepsF)
AIC(issfFor200_2)
issfFor400_2<-amt::fit_issf(Loc~-1+forDens400_2+strata(new_step_id_), data=rndStepsF)
AIC(issfFor400_2)

## select highway density resolution ##
issfAuto50_2<-amt::fit_issf(Loc~-1+AutoDens50_2+strata(new_step_id_), data=rndStepsF)
AIC(issfAuto50_2)
issfAuto100_2<-amt::fit_issf(Loc~-1+AutoDens100_2+strata(new_step_id_), data=rndStepsF)
AIC(issfAuto100_2)
issfAuto200_2<-amt::fit_issf(Loc~-1+AutoDens200_2+strata(new_step_id_), data=rndStepsF)
AIC(issfAuto200_2)
issfAuto400_2<-amt::fit_issf(Loc~-1+AutoDens400_2+strata(new_step_id_), data=rndStepsF)
AIC(issfAuto400_2)

## plot histograms
par( mfrow=c(2,2), mar=c(4,4,1,0) ) 
hist(rndStepsF$forDens50_log1, breaks=10
     , xlim=c(0,10) , col=rgb(1,0,0,0.5) , xlab="forDens50_log1" , ylab="" ,
     main="" ) 

hist(rndStepsF$forDens50_2, breaks=10 , xlim=c(0,1.5) ,
     col=rgb(0,0,1,0.5) , xlab="forDens50_log2" , ylab="" , main="")

hist(rndStepsF$forDens200_log1, breaks=10 , xlim=c(0,10) , col=rgb(0,1,0,0.5)
     , xlab="forDens200_log1" , ylab="" , main="" )

hist(rndStepsF$forDens200_2, breaks=10 , xlim=c(0,1.5) , col=rgb(0,0,0,0.5)
     , xlab="forDens200_2" , ylab="" , main="")


#####  entry and extract distance file scaled log ##### 

## building distance file scaled log ##

buildDistPath<-c(paste0(rasterFolder,"/human/others/log_scaled_distances/Dist_Merge_Bati_16b_log+1_scaled.tif"))
distBati_scld2<-terra::rast(buildDistPath)
rndStepsF$distBati_scld2<-extractCovariates(distBati_scld2, rndStepsF, covarExtractionType = exType)

## main road distance file scaled log ##
roadDistPath<-c(paste0(rasterFolder, "/human/others/log_scaled_distances/Dist_Merge_RoadPrimary_16b_log+1_scaled.tif"))
distMroad_scld2<-terra::rast(roadDistPath)
rndStepsF$distMroad_scld2<-extractCovariates(distMroad_scld2, rndStepsF, covarExtractionType = exType)

## secondary road distance file scaled log ##
roadDistPath<-c(paste0(rasterFolder, "/human/others/log_scaled_distances/Dist_Merge_RoadSecondary_16b_log+1_scaled.tif"))
distSroad_scld2<-terra::rast(roadDistPath)
rndStepsF$distSroad_scld2<-extractCovariates(distSroad_scld2, rndStepsF, covarExtractionType = exType)

##highways distance file scaled log ##
roadDistPath<-c(paste0(rasterFolder, "/human/others/log_scaled_distances/Dist_Merge_Autobahn_16b_log+1_scaled.tif"))
distAuto_scld2<-terra::rast(roadDistPath)
rndStepsF$distAuto_scld2<-extractCovariates(distAuto_scld2, rndStepsF, covarExtractionType = exType)

## plot histograms
par( mfrow=c(2,2), mar=c(4,4,1,0) ) 
hist(rndStepsF$distBati_scld2, breaks=10
     , xlim=c(0,1) , col=rgb(1,0,0,0.5) , xlab="distBati_scld2" , ylab="" ,
     main="" ) 

hist(rndStepsF$distBati_scld, breaks=10 , xlim=c(0,1) ,
     col=rgb(0,0,1,0.5) , xlab="distBati_scld" , ylab="" , main="")

hist(rndStepsF$distSroad_scld2, breaks=10 , xlim=c(-5,1) , col=rgb(0,1,0,0.5)
     , xlab="distSroad_scld2" , ylab="" , main="" )

hist(rndStepsF$distSroad_scld, breaks=10 , xlim=c(-1,1) , col=rgb(0,0,0,0.5)
     , xlab="distSroad_scld" , ylab="" , main="")

##### select between scaled or scaled-log distance variables ####################################################
#### this was only doing for the plateau and for the other bio-regions we used scaled-log distance variables ####

## building distance
issfdBati_1<-amt::fit_issf(Loc~-1+distBati_scld+strata(new_step_id_), data=rndStepsF)
AIC(issfdBati_1)
issfdBati_2<-amt::fit_issf(Loc~-1+distBati_scld2+strata(new_step_id_), data=rndStepsF)
AIC(issfdBati_2)

## main road distance
issfdMroad_1<-amt::fit_issf(Loc~-1+distMroad_scld+strata(new_step_id_), data=rndStepsF)
AIC(issfdMroad_1)
issfdMroad_2<-amt::fit_issf(Loc~-1+distMroad_scld2+strata(new_step_id_), data=rndStepsF)
AIC(issfdMroad_2)

## secondary road distance
issfdSroad_1<-amt::fit_issf(Loc~-1+distSroad_scld+strata(new_step_id_), data=rndStepsF)
AIC(issfdSroad_1)
issfdSroad_2<-amt::fit_issf(Loc~-1+distSroad_scld2+strata(new_step_id_), data=rndStepsF)
AIC(issfdSroad_2)

#### entry and extract relief scaled files #####################################
#### this was only doing for the plateau and for the jura (not for plateau) ####

## scaled altitude
altPath<-c(paste0(rasterFolder,"/relief/scaled/Altitude_5m_16b_scaled.tif"))
alt_scld<-terra::rast(altPath)
rndStepsF$alt_scld<-extractCovariates(alt_scld, rndStepsF, covarExtractionType = exType)

## scaled exposition
expPath<-c(paste0(rasterFolder,"/relief/scaled/Exposition_5m_16b_scaled.tif"))
exp_scld<-terra::rast(expPath)
rndStepsF$exp_scld<-extractCovariates(exp_scld, rndStepsF, covarExtractionType = exType)

## scaled slope
slopePath<-c(paste0(rasterFolder,"/relief/scaled/Slope_5m_8b_scaled.tif"))
slope_scld<-terra::rast(slopePath)
rndStepsF$slope_scld<-extractCovariates(slope_scld, rndStepsF, covarExtractionType = exType)


##### Investigating the correlation between continuous variables pre-selected for the analysis

##plateau second model
# dat <- rndStepsF %>% dplyr::select(buildDens100_2, MroadDens100_2,SroadDens50_2,AutoDens400_2, forDens200_2,distBati_scld2,distMroad_scld2,distSroad_scld2, distAuto_scld2, slope_scld)
# library(corrplot)
# corrplot(cor(dat),
#          method = "number",
#          type = "upper") # show only upper side

## alps
dat <- rndStepsF %>% dplyr::select(buildDens100_2, MroadDens400_2,SroadDens50_2,AutoDens400_2, forDens100_2,distBati_scld2,distMroad_scld2, distSroad_scld2, distAuto_scld2 , alt_scld, exp_scld, slope_scld)
library(corrplot)
corrplot(cor(dat),
         method = "number",
         type = "upper") # show only upper side

####### Elaborate a step resource selection function using  habitat categories with a reference class and all categorical variables scaled log ###
## set reference class
finDFsubRL <-within(rndStepsF, landUse <-relevel(landUse, ref = "FoFe"))

## elaborate the model

##plateau first model
#issf.all2<-amt::fit_issf(Loc ~ -1+ buildDens100_2+MroadDens100_2+SroadDens50_2+forDens200_2+distBati_scld2+distMroad_scld2+distSroad_scld2+slope_scld+ landUse+ sl_  + ta_ + strata(new_step_id_), data=finDFsubRL)

##plateau second model
# issf.all3<-amt::fit_issf(Loc ~ -1+ buildDens100_2+MroadDens100_2+SroadDens50_2+AutoDens400_2+forDens200_2+distBati_scld2+distMroad_scld2+distSroad_scld2+distAuto_scld2+slope_scld+ landUse+ sl_  + ta_ + strata(new_step_id_), data=finDFsubRL)
# summary(issf.all3)
# summISSFAllRL3<-as.data.frame(broom::tidy(issf.all3$model, conf.int=TRUE))

##plateau third model
# issf.all4<-amt::fit_issf(Loc ~ -1+ buildDens100_2+MroadDens100_2+SroadDens50_2+AutoDens400_2 + I(AutoDens400_2^2)+forDens200_2+distBati_scld2+distMroad_scld2+distSroad_scld2+distAuto_scld2+slope_scld+ landUse+ sl_  + ta_ + strata(new_step_id_), data=finDFsubRL)
# summary(issf.all4)
# summISSFAllRL4<-as.data.frame(broom::tidy(issf.all4$model, conf.int=TRUE))

##alps first model
#issf.all2<-amt::fit_issf(Loc ~ -1+ buildDens100_2+MroadDens400_2+SroadDens50_2+forDens100_2+distBati_scld2+distMroad_scld2+distSroad_scld2+alt_scld+exp_scld+slope_scld+ landUse+ sl_  + ta_ + strata(new_step_id_), data=finDFsubRL)
#summary(issf.all2)
#summISSFAllRL2<-as.data.frame(broom::tidy(issf.all2$model, conf.int=TRUE))

##alps second model
issf.all5<-amt::fit_issf(Loc ~ -1+ buildDens100_2+MroadDens400_2+SroadDens50_2+forDens100_2+distBati_scld2+distMroad_scld2+distSroad_scld2+distAuto_scld2+alt_scld+exp_scld+slope_scld+ landUse+ sl_  + ta_ + strata(new_step_id_), data=finDFsubRL)
summary(issf.all5)
summISSFAllRL5<-as.data.frame(broom::tidy(issf.all5$model, conf.int=TRUE))

## write results
#plateau first model
#write.csv(summISSFAllRL2, file="/media/loreto/Grande/ie-ofev-24-25/ssf_plateau/issf_all_2_nbRandom100_F.csv")#model used for plateau

#plateau second model
#write.csv(summISSFAllRL3, file="/media/loreto/Grande/ie-ofev-24-25/ssf_plateau/issf_all_3_nbRandom100_F.csv")#model used for plateau
#plateau third model
#write.csv(summISSFAllRL4, file="/media/loreto/Grande/ie-ofev-24-25/ssf_plateau/issf_all_4_nbRandom100_F.csv")#model used for plateau

##alps
#write.csv(summISSFAllRL2, file="/media/loreto/Grande/ie-ofev-24-25/ssf_alps/issf_all_2_nbRandom100_F.csv")#model used for alps
##alps second model
write.csv(summISSFAllRL2, file="/media/loreto/Grande/ie-ofev-24-25/ssf_alps/issf_all_5_nbRandom100_F.csv")#model used for alps
##alps second model
####### Other step resource selection functions #######

### step resource selection function using only habitat categories and a reference class ### 
## built the function using the reference class "FoFe"
issf.habitat.RL<-amt::fit_issf(Loc ~ -1+landUse+ sl_ + ta_ + strata(new_step_id_), data=finDFsubRL, model=T)
##see results
summary(issf.habitat.RL)
summISSF.habitat.RL<-as.data.frame(broom::tidy(issf.all$model, conf.int=TRUE))
#summISSF.habitat.RL[order(summISSF.habitat.RL$estimate, decreasing=T),]
##export results
#write.csv(summISSFAllRL, file="/media/loreto/Grande/ie-ofev-24-25/ssf_plateau/issf_habitat.RL_nbRandom100_F.csv")#model for plateau but not used

### srsf with each class as a variable ### 
##prepare data
finDF2<-cbind(rndSteps, tLandUse)
##subset to remove used points that are aberrant e.g in Highways,etc.. 
##using "to remove" defined before
finDF2<-subset(finDF2, !finDF2$new_step_id_%in%toRemove)
##not: not in data set 'Landes', 'PassSupL'
issf.habitat<-amt::fit_issf(Loc ~ -1+ Alluv+ Auto+ Bât+ Buisso+ Chem+ Clai+ CultHerb+ CultLign+ Etangs+ FoCnf+ FoFe+ FoTourb+ GazoPrai + Goudr+ Lacs  + Ouvr + PassInfM+ PelPrai+ Rivie+ RouPrim+ RouSec+ Ruiss+ Vfer+ Zohum + sl_ + ta_ + strata(new_step_id_), data=finDF2, model=T)
##see results
summary(issf.habitat)
summISSFhabitat<-as.data.frame(broom::tidy(issf.habitat$model, conf.int=TRUE))
# rownames(summISSFhabitat)<-summISSFhabitat$term
# print(summISSFhabitat)
# summISSFhabitat[order(summISSFhabitat$estimate, decreasing = T),]
##export results
#write.csv(summISSFhabitat, file="/media/loreto/Grande/ie-ofev-24-25/ssf_plateau/issf_habitat_nbRandom100_F.csv")#model for plateau but not used

### srsf with only human perturbations ### 
issf.anthropo<-amt::fit_issf(Loc ~ -1 + buildDens100_log1 + MroadDens100_log1 + SroadDens50_log1 + distBati_scld + distMroad_scld + distSroad_scld + sl_  + ta_ + strata(new_step_id_), data=finDFsubRL)
summISSFanthropo<-as.data.frame(broom::tidy(issf.anthropo$model, conf.int=TRUE))
# summISSFanthropo[order(summISSFanthropo$estimate, decreasing=T),]
##export results
#write.csv(summISSFanthropo, file="/media/loreto/Grande/ie-ofev-24-25/ssf_plateau/issf_anthropo_nbRandom100_F.csv")#model for plateau but not used

######srsf with with each habitat class as a variable and other continuous and human perturbations variables
issf.all.all<-amt::fit_issf(Loc ~ -1 + Alluv+ Auto+ Bât+ Buisso+ Chem+ Clai+ CultHerb+ CultLign+ Etangs+ FoCnf+ FoFe+ FoTourb+ GazoPrai + Goudr+ Lacs  + Ouvr + PassInfM+ PelPrai+ Rivie+ RouPrim+ RouSec+ Ruiss+ Vfer+ Zohum+ forDens50_log1+slope_scld+ buildDens100_log1 + MroadDens100_log1 + SroadDens50_log1 + distBati_scld + distMroad_scld + distSroad_scld + sl_  + ta_ + strata(new_step_id_), data=finDF2)
##not: not in data set 'Landes', 'PassSupL'
summISSFHabAnthropo<-as.data.frame(broom::tidy(issf.all.all$model, conf.int=TRUE))
#summISSFanthropo[order(summISSFanthropo$estimate, decreasing=T),]
#print(summISSFanthropo)
##export results
write.csv(summISSFHabAnthropo, file="/media/loreto/Grande/ie-ofev-24-25/ssf_plateau/issf_habitat_anthropo_nbRandom100_F.csv")#model for plateau but not used

###save results in an image file .Rdata
#save.image("/media/loreto/Grande/R_lore/ssf_26juin_plateau_2_100avlble.Rdata")
save.image("/media/loreto/Grande/R_lore/ssf_30juin_alpes_100avlble.Rdata")
