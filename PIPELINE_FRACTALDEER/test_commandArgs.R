
cerf_path<-commandArgs(trailingOnly = T)
inFolder=cerf_path[1]
tmpPath=unlist(strsplit(inFolder, "/"))
print(tmpPath)
selectedDeer<-strsplit(tmpPath[length(tmpPath)], split="_")[[1]][2]
selectedYear<-paste(strsplit(tmpPath[length(tmpPath)], split="_")[[1]][3],strsplit(tmpPath[length(tmpPath)], split="_")[[1]][4], sep="_")

cat("selectedDeer : ", selectedDeer, "  selectedYear : ", selectedYear)
