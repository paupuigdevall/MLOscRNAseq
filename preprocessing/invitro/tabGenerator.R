

library(readxl)

tt <- read_excel("metadata/metadata_mid_organoids_10x_sample.xlsx", sheet = "invitro")
tt <- as.data.frame(tt)
list_perFile <- split(tt$Donors, tt$Sample)


tt2 <- read_excel("metadata/metadata_mid_organoids_10x_sample.xlsx", sheet = "chipInfo")
tt2 <- as.data.frame(tt2)
tt2 <- subset(tt2, sampleOrigin=="invitro")

vecIds <- tt2$chipId
names(vecIds) <- tt2$sample


for (i in 1:length(list_perFile)){
  
  tmp <- as.data.frame(unname(vecIds[list_perFile[[x]]]))
  colnames(tmp) <- NULL
  fileName=names(list_perFile)[x]
  pathToDemuxlet <- "preprocessing/invitro/"
  write.table(tmp, file = paste0(pathToDemuxlet,fileName,"_sample_list.txt"),
              quote=F,
              sep="\t",
              col.names=F,
              row.names=F)
}





  