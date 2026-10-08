############ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #####################
#########  Generate random steps along paths by behaviour    ##########################
##### - By Loreto Urbina et Lucas Villard
#### --------------------------------------- ################

library(lubridate)
library(amt)
library(tidyverse)
library(sf)
library(stringr)

create_random_step_along_path_by_behaviour<-function(behaviourFile=character(), chosenBehaviour=character(), nbRandom=numeric()){
  ###adding a only one "refined" trajectory, please verify that var1 correspond to "in-habitat" behaviour and is equal to "mc_1"
  datCerf<-readr::read_csv(behaviourFile)%>%
    mutate(id=factor(id))%>%
    mutate(path_no=factor(new_path_no))%>%
    mutate(deerYear=factor(deerYear))%>%
    mutate(saison=factor(saison))%>%
    mutate(behaviour=factor(behaviour))%>%
    mutate(jourNuit=factor(jourNuit))
  # mutate(jourNuit=factor(tod))#if cerfs from article FractalDeer
  suppressWarnings(datCerf<-datCerf%>%transform(saison = forcats::fct_relevel(saison, c("Mars-Mai","Juin-Aout","Septembre-Novembre","Decembre-Fevrier"))))
  var1<-"mc_1"
  datCerf<- datCerf%>%mutate(behaviour=ifelse(behaviour==var1,"in-patch", "in-matrix"))%>%
    mutate(behaviour=factor(behaviour))%>%
    dplyr::select(c("id","path_no", "x", "y","t","saison", "jourNuit","deerYear","behaviour"))
  print(datCerf%>%tabyl(behaviour))
  

  datCerfSub<-datCerf%>%
    filter(behaviour%in%chosenBehaviour)%>%
    droplevels() 
  
    datPath<-datCerfSub%>%
    nest_legacy(-path_no)

    
  trk_all<-datPath%>%mutate(trk=lapply(data, function(d){
    amt::make_track(d,x,y,t, crs=2056,all_cols = TRUE)
  }))
  

  sLsTa_all<-datPath%>%mutate(sLsTa=lapply(data, function(d){
    amt::make_track(d,x,y,t, crs=2056,all_cols = TRUE)%>%
      #summarize_sampling_rate()
      steps()
  }))
  datslsTa<-sLsTa_all%>%unnest_legacy(sLsTa)
  
  #built a Vonmises distribution from observed turning angles and a Gamma distribution from observed step lengths
  mod_Ta<-fit_distr(datslsTa$ta_, "vonmises", na.rm = TRUE)
  mod_sLs<-fit_distr(datslsTa$sl_, "gamma", na.rm = TRUE)
  
  #Filter data to keep only paths with 10 or more locations
  toFilter1<-trk_all%>%mutate(isNotEmpty=unlist(lapply(trk, function(x){ nrow(x)>=10})))%>%
    filter(isNotEmpty)%>%
    mutate(isNotEmpty=NULL)%>%
    droplevels()
  #Estimating random points
  rndSteps<-toFilter1%>%mutate(stps=lapply(trk, function(x){
    x %>% track_resample(rate = minutes(60), tolerance = minutes(4))%>%
      steps( keep_cols="end")%>%
      random_steps(n_control=nbRandom, sl_distr=mod_sLs, ta_distr=mod_Ta)
  }
  )
  )%>%unnest_legacy(stps)
rndSteps<-rndSteps%>%mutate(animalID=paste(levels(datCerf$id)[1], levels(datCerf$deerYear),sep="_"))  

rndSteps%>%return()  
}
