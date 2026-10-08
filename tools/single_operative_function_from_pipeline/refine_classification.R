

###################### @@@@ refine_classification.r ########################

# Author : lucas villard
# Last revised : jan 28 2025
# Usage : compute the refined classification from a classified GPS series output from the  "get_behaviour_vector_from_list.py"

### ---Input---
# - animalName, deerYear and integer stepSize you want to compute the refinment on. the script will look in the [animalName] folder for the dataframe and return
#   an error if it cannot find it
#   
# - the metric you want to refine on
#
### --output ---
## the GPS data series in input added with the refined classification and the new path number associated with it
## it will also output the summary plots showing the different effects of the refinement
### 
###  

###################################################################


### *** input parameters **** #############
  
selectedDeer<-"20074"
selectedYear<-"deerYear_2011-2012_denoised"
selectedStepSize<-1046
parameter="ratio_meanNSD"
##########################################





##### for RStudio : change to directory where the script file is ######
srcdir <- getSrcDirectory(function(){})[1]
if (srcdir=="") {
  # - run
  srcdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
}
setwd(srcdir)
rm(srcdir)



source("./getNSDValues.R")
source("./getNSDValuesFromList.R")
source("./refineClusteringByNSD.R")

behaviourFile=paste0("./",selectedDeer,"/",selectedDeer,"_",selectedYear,"_stepSize_",selectedStepSize,"_autocorrelation_timeSeriesKmeans_2_classes.csv")



x<-getNSDValues(behaviourFile, display="behaviour",save.plot=FALSE , mutate.df=TRUE, show.breakpoints = T)


newDat<-refineClusteringByNSD(behaviourFile,
                                parameter=parameter,
                                display.plot=TRUE,
                                save.plot = FALSE,
                                save.data.frame = TRUE,
                              show.breakpoints = TRUE,
                              save.BIC = FALSE)
#behaviourRefined=paste0("./20074_v1/20074 deerYear 2011-2012 stepSize 953 autocorrelation timeSeriesKmeans 2 classes REFINED RATIO_MEANNSD.csv")
behaviourRefined<-paste0("./",selectedDeer,"/",gsub("_", " " ,strsplit((strsplit(behaviourFile, split="\\.")[[1]][2]),split="/")[[1]][3])," REFINED ",toupper(parameter),".csv")
x<-getNSDValues(behaviourRefined, display="behaviour",save.plot=TRUE , mutate.df=FALSE, show.breakpoints = T)

###########################################





