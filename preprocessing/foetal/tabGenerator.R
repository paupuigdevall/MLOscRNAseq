
library(readxl)

tt <- read_excel("metadata/metadata_mid_organoids_10x_sample.xlsx", sheet = "foetal")
tt <- as.data.frame(tt)
list_perFile <- split(tt$Donors, tt$Supplier_Sample_Name)
list_perFile <- sapply(list_perFile, function(x) x[!duplicated(x)], simplify=F)

tt2 <- read_excel("metadata/metadata_mid_organoids_10x_sample.xlsx", sheet = "chipInfo")
tt2 <- as.data.frame(tt2)
tt2 <- subset(tt2, sampleOrigin=="foetal")

vecIds <- tt2$chipId
names(vecIds) <- tt2$sample


for (i in 1:length(list_perFile)){
  
  print(i)
  tmp <- as.data.frame(unname(vecIds[unlist(strsplit(list_perFile[[i]], ","))]))
  colnames(tmp) <- NULL
  fileName=names(list_perFile)[i]
  pathToDemuxlet <- "preprocessing/foetal/"
  write.table(tmp, file = paste0(pathToDemuxlet,fileName,"_sample_list.txt"),
              quote=F,
              sep="\t",
              col.names=F,
              row.names=F)
}






  