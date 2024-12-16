
library(Biobase)
library(GEOquery)
library(Seurat)
library(readxl)
library(ggplot2)
library(dplyr)
library(harmony)
library(GenomicRanges)
library(Seurat)
library(patchwork)
library(cowplot)
library(data.table)
library(scales)
library(org.Hs.eg.db)
library(rtracklayer)



getGEOSuppFiles('GSE76381', baseDir="downloads/datasetsToIntegrate/LeManno/")
dirProcessing <- "downloads/datasetsToIntegrate/LeManno/GSE76381/"
leManno <- read.table(paste0(dirProcessing,'GSE76381_EmbryoMoleculeCounts.cef.txt.gz'), sep="\t",
                      header = T, skip=1)
matrixExp <- leManno[4:dim(leManno)[1], 3:dim(leManno)[2]]
colnames(matrixExp) <- gsub("^X","", colnames(matrixExp))
rownames(matrixExp) <- leManno[,1][4:dim(leManno)[1]]

metaData <- data.frame(CellType=as.character(leManno[1,])[-c(1,2)],
                       Timepoint=as.character(leManno[2,])[-c(1,2)])
rownames(metaData) <- colnames(matrixExp)

## Create Seurat object
leMannoExp <- CreateSeuratObject(counts = matrixExp,
                                 meta.data = metaData,
                                 project = "LeMannoMidBrain2016")

saveRDS(leMannoExp, "referenceDatasets/LeMannoSeurat_RawExp.RDS")




