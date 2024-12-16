library(Seurat)
library(rhdf5)
library(Matrix)

lapply(c("dplyr","Seurat","patchwork","ggplot2","tidyr","openxlsx","harmony", "miloR", "SingleCellExperiment", "scater","SeuratWrappers"), library, character.only = T)


pathToFolder <- "downloads/datasetsToIntegrate/Linnearson/"
h5ls(paste0(pathToFolder, "HumanFetalBrainPool.h5"))
tissue <- h5read(paste0(pathToFolder,"HumanFetalBrainPool.h5"), "/shoji/Tissue")
ValidGenes <- h5read(paste0(pathToFolder,"HumanFetalBrainPool.h5"), "/shoji/ValidGenes")
mask <- tissue=="Ventral midbrain"

data <- h5read(paste0(pathToFolder,"HumanFetalBrainPool.h5"), "/shoji/Expression", index=list(NULL,which(mask)))
data <- data[which(ValidGenes==TRUE),]

genes <- h5read(paste0(pathToFolder,"HumanFetalBrainPool.h5"), "/shoji/Gene")[which(ValidGenes==TRUE)]
cellID <-  h5read(paste0(pathToFolder,"HumanFetalBrainPool.h5"), "/shoji/CellID")[mask]
rownames(data) <- genes
colnames(data) <- cellID

## metaData
age <-  round(h5read(paste0(pathToFolder,"HumanFetalBrainPool.h5"), "/shoji/Age"),1)[mask]
CellClass <- h5read(paste0(pathToFolder,"HumanFetalBrainPool.h5"), "/shoji/CellClass")[mask]
CellCycleFraction <- h5read(paste0(pathToFolder,"HumanFetalBrainPool.h5"), "/shoji/CellCycleFraction")[mask]
Chemistry <- h5read(paste0(pathToFolder,"HumanFetalBrainPool.h5"), "/shoji/Chemistry")[mask]
Clusters <- h5read(paste0(pathToFolder,"HumanFetalBrainPool.h5"), "/shoji/Clusters")[mask]
Donor <- h5read(paste0(pathToFolder,"HumanFetalBrainPool.h5"), "/shoji/Donor")[mask]

DoubletFlag <- h5read(paste0(pathToFolder,"HumanFetalBrainPool.h5"), "/shoji/DoubletFlag")[mask]
MitoFraction <- h5read(paste0(pathToFolder,"HumanFetalBrainPool.h5"), "/shoji/MitoFraction")[mask]
Region <- h5read(paste0(pathToFolder,"HumanFetalBrainPool.h5"), "/shoji/Region")[mask]
SampleID <- h5read(paste0(pathToFolder,"HumanFetalBrainPool.h5"), "/shoji/SampleID")[mask]

metadata <- list(age, CellClass, CellCycleFraction, Chemistry, Clusters, Donor,
                 DoubletFlag, MitoFraction, Region, SampleID)
metadata <- do.call("cbind", metadata)

colnames(metadata) <- c("age","CellClass","CellCycleFraction","Chemistry","Clusters","Donor",
                        "DoubletFlag","MitoFraction","Region","SampleID")
metadata <- as.data.frame(metadata)

metadata$age <- as.numeric(metadata$age)
metadata$CellCycleFraction <- as.numeric(metadata$CellCycleFraction)
metadata$MitoFraction <- as.numeric(metadata$MitoFraction)
rownames(metadata) <- cellID

dupGenes <- rownames(data)[which(duplicated(rownames(data)))]
data2 <- data[is.na(match(rownames(data), dupGenes)),]
early_ref <- CreateSeuratObject(counts = data2, meta.data=metadata, project="Linnarsson_VM")
early_ref$age <- paste0("PCW", early_ref$age )

annotDefinition <- h5read(paste0(pathToFolder,"HumanFetalBrainPool.h5"), "/shoji/AnnotationDefinition")
annotDescription <- h5read(paste0(pathToFolder,"HumanFetalBrainPool.h5"), "/shoji/AnnotationDescription")
annotName <- h5read(paste0(pathToFolder,"HumanFetalBrainPool.h5"), "/shoji/AnnotationName")
annotPosterior <- h5read(paste0(pathToFolder,"HumanFetalBrainPool.h5"), "/shoji/AnnotationPosterior")

saveRDS(early_ref, file="referenceDatasets/ventralMidbrain_Linnarson.RDS")



