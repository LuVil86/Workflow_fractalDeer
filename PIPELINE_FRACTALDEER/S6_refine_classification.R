

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
source("./getNSDValues.R")
source("./refineClusteringByNSD.R")
require(fs)

### *** input parameters **** #############

inBehaviourFile="/home/luvil/IE-OFEV/IE-OFEV SCRIPTS/Workflow_fractalDeer/PIPELINE_FRACTALDEER/OUT_20206_deerYear_2019-2020/20206_deerYear_2019-2020_stepSize_1863_autocorrelation_timeSeriesKmeans_2_classes.csv"
dirPath=path_dir(inBehaviourFile)
fileName=path_file(inBehaviourFile)
parameter="ratio_endNSD"
##########################################




x<-getNSDValues(inBehaviourFile, display="behaviour",save.plot=FALSE , mutate.df=TRUE, show.breakpoints = T)


newDat<-refineClusteringByNSD(inBehaviourFile,
                                parameter=parameter,
                                display.plot=TRUE,
                                save.plot = FALSE,
                                save.data.frame = TRUE,
                              show.breakpoints = TRUE,
                              save.BIC = FALSE)

behaviourRefined<-paste0(dirPath,"/",gsub("_", " ", strsplit(fileName, split=".csv")[[1]][1])," REFINED ",toupper(parameter),".csv")

x<-getNSDValues(behaviourRefined, display="behaviour",save.plot=TRUE , mutate.df=FALSE, show.breakpoints = T)

###########################################




