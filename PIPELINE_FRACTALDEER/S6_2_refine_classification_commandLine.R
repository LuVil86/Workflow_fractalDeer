

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
if(!require(stringr)) {
  install.packages("stringr"); require(stringr)}


### *** input parameters **** #############
### mettre le chemin au dossier des résultats 

#resultFolder="/media/loreto/Grande/"
#resultFolder="~/Documents/hepiaBy24/IE-OFEV-24-25/"

#### **** ICI  copier le contenu de la colonne "fullPath" de la ligne de stepSize à analyser ***** #####


inputParam<-commandArgs(trailingOnly = T)
fullPath=inputParam[1]




#### *** ICI la mesure à prendre en compte : "ratio_endNSD", "ratio_meanNSD" ou "ratio_cumulativeNSD" ***** #######
parameter=inputParam[2]
#parameter="ratio_endNSD"
#parameter="ratio_cumulativeNSD"

################################################################################################################################

#inBehaviourFile=paste0(resultFolder,str_split_fixed(fullPath, '/', 5)[1,5])
inBehaviourFile <- fullPath
dirPath=path_dir(inBehaviourFile)
fileName=path_file(inBehaviourFile)
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




