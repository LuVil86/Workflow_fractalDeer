
plotMigrateR<-function(behaviourFile=character(), outDir=character(),mutate.df=TRUE){
  suppressPackageStartupMessages({
  require(multidplyr)
  require(tidyverse)
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
  #outDir=path_dir(behaviourFile)
  outDir=outDir
  animalTest<-strsplit(plotTitle, split=" ")[[1]][1]
  plotTitle_for_classical<-paste(strsplit(plotTitle, split=" ")[[1]][1:3], collapse = " ")
  if(!grepl(".csv", behaviourFile)){stop("The file you provided does not have '.csv' extension")}
  if(mutate.df==TRUE){
    datCerfs<-readr::read_csv(behaviourFile)%>%
      mutate(coord_X=apply(as.data.frame(geometry), 1,function(x) coordSplit(x)$coordX ))%>%
      mutate(coord_Y=apply(as.data.frame(geometry), 1,function(x) coordSplit(x)$coordY ))%>%
      mutate(dateTime=as.POSIXct(paste(UTC_DATE, UTC_TIME, sep=" "), origin="1970-01-01", tz="UTC"))%>%
      mutate(deerYear=factor(deerYear))%>%
      mutate(prenom=factor(prenom))%>%
      mutate(saison=factor(saison))%>%
      dplyr::select(x="coord_X", y="coord_Y",t="dateTime", id="prenom",deerYear="deerYear", jourNuit="jourNuit",saison="saison")
    suppressWarnings(datCerfs<-datCerfs%>%transform(saison = fct_relevel(saison, c("Mars-Mai","Juin-Aout","Septembre-Novembre","Decembre-Fevrier"))))
  }else{
    datCerfs<-readr::read_csv(behaviourFile)%>%    
      mutate(path_no=factor(new_path_no))%>%
      mutate(behaviour=ifelse(behaviour=="mc_1","in-patch", "in-matrix"))%>%
      suppressWarnings(datCerfs<-datCerfs%>%transform(saison = fct_relevel(saison, c("Mars-Mai","Juin-Aout","Septembre-Novembre","Decembre-Fevrier"))))
    
  }
  
  RD_traj<-as.ltraj(xy=datCerfs[,c("x","y")], date=datCerfs$t, id=datCerfs$id)
  suppressWarnings(RD.nsd1 <- mvmtClass(RD_traj))
  cat("**** BEST MODEL ****\n")
  topmodel<-rownames(summary(topmvmt(RD.nsd1)))
  print(topmodel)
  
  fp<-file.path(outDir,paste("MIGRATE_R PLOT ", plotTitle_for_classical,".png",sep=""), fsep = "/" )
  png(file=fp,
      width=800, height=650)
    plot(RD.nsd1, main=plotTitle_for_classical)
    if (topmodel == "mixmig"){
      
      breakpoints<-mvmt2dt(RD.nsd1,mod="mixmig" )
      abline(v=breakpoints[[1]]["str1","dday"], col="red", lty="dashed")
      abline(v=breakpoints[[1]]["end1","dday"], col="red", lty="dashed")
      abline(v=breakpoints[[1]]["str2","dday"], col="orange", lty="dashed")
      abline(v=breakpoints[[1]]["end2","dday"], col="orange", lty="dashed")
      
    }
    if(topmodel=="migrant"){
      breakpoints<-mvmt2dt(RD.nsd1,mod="migrant" )
      abline(v=breakpoints[[1]]["str1","dday"], col="red", lty="dashed")
      abline(v=breakpoints[[1]]["end1","dday"], col="red", lty="dashed")
      abline(v=breakpoints[[1]]["str2","dday"], col="orange", lty="dashed")
      abline(v=breakpoints[[1]]["end2","dday"], col="orange", lty="dashed")
    }
  

    dev.off()
}
  


arg<-commandArgs(trailingOnly = T)

plotMigrateR(arg[1], arg[2], arg[3])