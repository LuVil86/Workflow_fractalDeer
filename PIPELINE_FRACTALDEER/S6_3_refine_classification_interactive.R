

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



getNSDValues_forGUI<-function(behaviourFile=character()){
  suppressPackageStartupMessages({
  require(multidplyr)
  require(tidyverse)
  require(dplyr)
  require(lubridate)
  require(adehabitatLT)
  require(arulesViz)
  require(gdata)
  require(ggplot2)
  require(ggtext)
  require(fs)
  })
  sourceDir <- function(path, trace = TRUE, ...) {
    op <- options(); on.exit(options(op)) # to reset after each 
    for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
      source(file.path(path, nm), ...)
      options(op)
    }
  }
  sourceDir("./migrateR_1.0.9/migrateR/R/")
  options(readr.show_col_types = FALSE)
  options(tidyverse.quiet = TRUE)
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
  
  
  plotTitle<-strsplit(path_file(behaviourFile), split=".csv")[[1]][1]
  plotTitle<-gsub("_"," ",tail(strsplit(plotTitle, split="/")[[1]],n=1))
  outDir=path_dir(behaviourFile)
  animalTest<-strsplit(plotTitle, split=" ")[[1]][1]
  plotTitle_for_classical<-paste(strsplit(plotTitle, split=" ")[[1]][1:3], collapse = " ")
  if(!grepl(".csv", behaviourFile)){stop("The file you provided does not have '.csv' extension")}
    datCerfs<-readr::read_csv(behaviourFile)%>%    
      mutate(path_no=factor(new_path_no))%>%
      mutate(behaviour=ifelse(behaviour=="mc_1","in-patch", "in-matrix"))%>%
      suppressWarnings(datCerfs<-datCerfs%>%transform(saison = fct_relevel(saison, c("Mars-Mai","Juin-Aout","Septembre-Novembre","Decembre-Fevrier"))))
    
  
  
  RD_traj<-as.ltraj(xy=datCerfs[,c("x","y")], date=datCerfs$t, id=datCerfs$id)
  suppressWarnings(RD.nsd1 <- mvmtClass(RD_traj))
  topmodel<-rownames(summary(topmvmt(RD.nsd1)))

  nsdDF<-data.frame(dateTime=datCerfs["t"],
                    NSD=RD.nsd1[[as.character(unique(datCerfs$id))]]@data[,2],
                    saison=datCerfs["saison"], 
                    jourNuit=datCerfs["jourNuit"],
                    behaviour=datCerfs["behaviour"],path_no=datCerfs["path_no"])

    pl<-ggplot(aes(x=t, y=NSD), data=nsdDF)+
      geom_point(aes(colour=behaviour),size=0.5)+
      #  ggtitle(plotTitle_for_classical)+
      theme_bw()+
      ylab("NSD [km²]")+
      scale_x_datetime(date_labels = "%b %y", date_breaks = "1 month")+
      scale_colour_manual(values = c("red","black","blue","chartreuse3","gold1"))+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=8),
            axis.text.y=element_text(size=8),
            axis.title.x=element_blank(),
            axis.title.y=element_text(size=8),
            legend.title = element_blank(),
            legend.position = "none",
            legend.text = element_text(size=8),
            legend.margin = margin(0,0,0,0),
            legend.box.margin = margin(0,0,-10,0),
            #plot.title = element_textbox_simple(halign = 0.5, size=8),
            panel.background = element_rect(fill='transparent'), #transparent panel bg
            plot.background = element_rect(fill='transparent', color=NA)
      )

      tmp<-as.data.frame(datCerfs)
      breakpoints<-c()
      for(i in 2:nrow(tmp)){
        if(tmp$behaviour[i-1]!=tmp$behaviour[i]){
          breakpoints<-c(breakpoints, i)
        }
      }
      rm(tmp)
      tBreakPoints<-datCerfs$t[breakpoints]
      ymax <- ggplot_build(pl)$layout$panel_params[[1]]$y.range[2]
      pl<-pl+geom_vline(xintercept = tBreakPoints, linetype="dotted")
      

      fp<-file.path(outDir,"tmp_S6_3.png", fsep = "/" )
      ggsave(filename =fp,
             plot = pl,
             device = "png",
             width = 174, height = 116, units = "mm",
             dpi=600)
   
}







###################################################################
source("./getNSDValues.R")
source("./refineClusteringByNSD.R")
  suppressPackageStartupMessages({

require(fs)
if(!require(stringr)) {
  install.packages("stringr"); require(stringr)}

})
inputParam<-commandArgs(trailingOnly = T)
fullPath=inputParam[1]




parameter=inputParam[2]

################################################################################################################################

#inBehaviourFile=paste0(resultFolder,str_split_fixed(fullPath, '/', 5)[1,5])
inBehaviourFile <- fullPath
dirPath=path_dir(inBehaviourFile)
fileName=path_file(inBehaviourFile)
x<-getNSDValues(inBehaviourFile, display="behaviour",save.plot=FALSE , mutate.df=TRUE, show.breakpoints = T,print_output=FALSE)


newDat<-refineClusteringByNSD(inBehaviourFile,
                                parameter=parameter,
                                display.plot=TRUE,
                                save.plot = FALSE,
                                save.data.frame = TRUE,
                              show.breakpoints = TRUE,
                              save.BIC = FALSE)

behaviourRefined<-paste0(dirPath,"/",gsub("_", " ", strsplit(fileName, split=".csv")[[1]][1])," REFINED ",toupper(parameter),".csv")

x<-getNSDValues_forGUI(behaviourRefined)

###########################################




