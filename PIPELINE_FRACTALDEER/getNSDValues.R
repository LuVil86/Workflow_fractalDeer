
getNSDValues<-function(behaviourFile=character(), display="saison", save.plot=FALSE, mutate.df=TRUE, show.breakpoints=TRUE,print_output=TRUE){
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
cat(" input behaviour file -->  ",behaviourFile,"\n")
if(!grepl(".csv", behaviourFile)){stop("The file you provided does not have '.csv' extension")}
if(mutate.df==TRUE){
datCerfs<-readr::read_csv(behaviourFile, name_repair = "minimal")%>%
  mutate(coord_X=apply(as.data.frame(geometry), 1,function(x) coordSplit(x)$coordX ))%>%
  mutate(coord_Y=apply(as.data.frame(geometry), 1,function(x) coordSplit(x)$coordY ))%>%
  mutate(dateTime=as.POSIXct(paste(UTC_DATE, UTC_TIME, sep=" "), origin="1970-01-01", tz="UTC"))%>%
  mutate(behaviour=factor(paste("b", behaviour, sep="_")))%>%
  mutate(deerYear=factor(deerYear))%>%
  mutate(prenom=factor(prenom))%>%
  mutate(saison=factor(saison))%>%
  mutate(path_no=factor(path_no))%>%
  dplyr::select(x="coord_X", y="coord_Y",t="dateTime", id="prenom",deerYear="deerYear", jourNuit="jourNuit",saison="saison","behaviour"=behaviour, "path_no"=path_no)
  suppressWarnings(datCerfs<-datCerfs%>%transform(saison = fct_relevel(saison, c("Mars-Mai","Juin-Aout","Septembre-Novembre","Decembre-Fevrier"))))
}else{
  datCerfs<-readr::read_csv(behaviourFile, id_repair = "unique")%>%    
  mutate(path_no=factor(new_path_no))%>%
  mutate(behaviour=ifelse(behaviour=="mc_1","in-patch", "in-matrix"))%>%
  suppressWarnings(datCerfs<-datCerfs%>%transform(saison = fct_relevel(saison, c("Mars-Mai","Juin-Aout","Septembre-Novembre","Decembre-Fevrier"))))
  
}

RD_traj<-as.ltraj(xy=datCerfs[,c("x","y")], date=datCerfs$t, id=datCerfs$id)
suppressWarnings(RD.nsd1 <- mvmtClass(RD_traj))
topmodel<-rownames(summary(topmvmt(RD.nsd1)))
if(print_output == TRUE){
cat("**** BEST MODEL ****\n")
  print(topmodel)
}



nsdDF<-data.frame(dateTime=datCerfs["t"],
                  NSD=RD.nsd1[[as.character(unique(datCerfs$id))]]@data[,2],
                  saison=datCerfs["saison"], 
                  jourNuit=datCerfs["jourNuit"],
                  behaviour=datCerfs["behaviour"],path_no=datCerfs["path_no"])
if(display=="saison"){
pl<-ggplot(aes(x=t, y=NSD), data=nsdDF)+
  geom_point(aes(colour=saison))+
  ggtitle(plotTitle)+
  theme_bw()+
  ylab("Net-Squared Displacement [km²]")+
  scale_x_datetime(date_labels = "%b %y", date_breaks = "1 month")+
  geom_vline(xintercept = tBreakPoints, linetype="dotted")+
  scale_colour_manual(values = c("chartreuse",
                                             "goldenrod1",
                                             "tan3",
                                             "deepskyblue"))+
                                               theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                                                     axis.text.y=element_text(),
                                                     axis.title.x=element_blank(),
                                                     legend.title = element_blank(), legend.position = "top",plot.title = element_textbox_simple(halign = 0.5))
plot(pl)
if (save.plot==TRUE){
  print("trying to save plot...");
  fp<-file.path(getwd(),animalTest,paste("NSD WITH SEASONS ", plotTitle,".png",sep=""), fsep = "/" )
  ggsave(filename =fp,
         plot = pl,
         device = "png")
}

}else if(display=="behaviour"){
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
  if(show.breakpoints==TRUE){
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
    
    #  annotate("text", x=tBreakPoints, y=rep(ymax,length(tBreakPoints)), label=finDat$path_no[-1], size=4)
  
  }
  plot(pl)
  if (save.plot==TRUE){
    print("trying to save plot...");
    fp<-file.path(outDir,paste("NSD WITH BEHAVIOUR ", plotTitle,".png",sep=""), fsep = "/" )
    print(fp)
    ggsave(filename =fp,
           plot = pl,
           device = "png",
           #device = "tiff",
           width = 87, height = 58, units = "mm",
           #width = 174, height = 116, units = "mm",
           dpi=600
    )
  }
}else if (display=="migrateR"){
  plot(RD.nsd1, main=plotTitle_for_classical)
  if (topmodel == "mixmig"){
    
    breakpoints<-mvmt2dt(RD.nsd1,mod="mixmig" )
    abline(v=breakpoints[[1]]["str1","dday"], col="red", lty="dashed")
    abline(v=breakpoints[[1]]["end1","dday"], col="red", lty="dashed")
    abline(v=breakpoints[[1]]["str2","dday"], col="orange", lty="dashed")
    abline(v=breakpoints[[1]]["end2","dday"], col="orange", lty="dashed")
    
  }
  if (save.plot=="TRUE"){
    print("trying to save plot...");
    fp<-file.path(outDir,paste("MIGRATE_R PLOT ", plotTitle_for_classical,".png",sep=""), fsep = "/" )
    png(file=fp,
        width=800, height=650)
    plot(RD.nsd1, main=plotTitle_for_classical)
    dev.off()
  }
  }else if  (display=="raw"){
    pl<-ggplot(aes(x=t, y=NSD), data=nsdDF)+
      geom_point(col="black")+
      ggtitle(plotTitle_for_classical)+
      theme_bw()+
      ylab("Net-Squared Displacement [km²]")+
      scale_x_datetime(date_labels = "%b %y", date_breaks = "1 month")
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            axis.text.y=element_text(),
            axis.title.x=element_blank(),
            legend.title = element_blank(), legend.position = "top",plot.title = element_textbox_simple(halign = 0.5))
    
      if (save.plot=="TRUE"){
        print("trying to save plot...");
        fp<-file.path(outDir,paste("RAW NSD PLOT ", plotTitle_for_classical,".png",sep=""), fsep = "/" )
        png(file=fp,
            width=800, height=650)
        plot(RD.nsd1, main=plotTitle_for_classical)
        dev.off()
      }
}  else {print("no plot was required. You can specifiy 'saison' ,'behaviour' or 'migrateR' ")}
return(nsdDF)
}


