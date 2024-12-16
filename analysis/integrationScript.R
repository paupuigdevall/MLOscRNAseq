
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
library(gghighlight)
library(dplyr)
library(Seurat)
library(patchwork)
library(openxlsx)
library(RColorBrewer)
library(tidyr)

##############################
### LeManno et al (Foetal) ###
##############################

leManno <- readRDS("referenceDatasets/LeMannoSeurat_RawExp.RDS")
leManno$system <- "Foetal"

############################
### Braun et al (Foetal) ###
############################

braun <- readRDS("referenceDatasets/ventralMidbrain_Linnarson.RDS")
## Removal of genes expressed in less than 0.1% total cells
selected_f <- names(which(rowSums(braun$RNA["counts"]>0)>0.001*dim(braun)[1]))
braun <- subset(braun, features = selected_f)
braun$system <- "Foetal"
braun$timePoint <- braun$age
braun$age <- NULL

########################
### Smits et al (3D) ###
########################

pathToFiles <- "downloads/datasetsToIntegrate/schwamborn/GSE133894_RAW/"
wtSamples <- paste0(pathToFiles, dir(pathToFiles)[grepl("wt", dir(pathToFiles))])
mtCountsWT <- sapply(wtSamples, function(x){
  tmp <- read.table(x, header = TRUE, row.names = 1)
}, simplify=F)


toMerge <- sapply(1:length(mtCountsWT), function(x){
  
  rawMatrix <- mtCountsWT[[x]]
  nameFile <- names(mtCountsWT)[x]
  nameFile <- sapply(strsplit(nameFile, "/"), function(x) x[length(x)])
  metaData <- data.frame(timePoint=rep(gsub(".+wt","",sapply(strsplit(nameFile, "_"), function(x) x[2])),length(colnames(rawMatrix))),
                         condition=rep(ifelse(grepl("wt", gsub(".txt.gz","",nameFile)),"WT"),length(colnames(rawMatrix))),
                         type=rep("organoid", length(colnames(rawMatrix))))
  metaData$timePoint <- gsub("[A-Z].+","",metaData$timePoint)
  rownames(metaData) <- colnames(rawMatrix)
  stopifnot(!any(duplicated(rownames(metaData))))
  rawMatrix <- CreateSeuratObject(counts = rawMatrix, meta.data = metaData, project = "SmitsLetAl2022")
  rawMatrix
  
}, simplify=F)

# Initialize the Seurat object with the raw (non-normalized data).
smitsEtAl <- merge(toMerge[[1]], y=do.call("c",toMerge[2:length(toMerge)]))
smitsEtAl <- JoinLayers(smitsEtAl)

## genes expressed in less than 0.1% total cells are removed
selected_f <- names(which(rowSums(smitsEtAl$RNA["counts"]>0)>0.001*dim(smitsEtAl)[1]))
smitsEtAl <- subset(smitsEtAl, features = selected_f)
smitsEtAl$system <- "Organoid"



###################################
### Agarwal et al (Post mortem) ###
###################################

##Post-mortem (single nuclei)

pathToDir <- "downloads/datasetsToIntegrate/Agarwal/GSE140231/"
files <- paste0(pathToDir, dir(pathToDir)[!grepl("tar", dir(pathToDir))],"/")
files <- files[grepl("_N",sapply(strsplit(files, "/"), function(x) x[12]))]

toMerge <- sapply(files, function(x){
  
  rawMatrix <- Read10X(data.dir = x)
  fileSpec <- gsub("/","",sapply(strsplit(x,"_"), function(y) y[length(y)]))
  
  ##metaData was downloaded from Supplementary Data 2 of the original publication. 
  metaData <- read.xlsx("downloads/datasetsToIntegrate/Agarwal/ctypeAnnot_Agarwal.xlsx")
  metaData <- subset(metaData, Library==fileSpec)
  colnames(rawMatrix) <- gsub("-.+","",paste0(metaData$Library,"_",colnames(rawMatrix)))
  rawMatrix <- rawMatrix[,colnames(rawMatrix) %in% metaData$Library_Barcode]
  metaData <- metaData[match(colnames(rawMatrix), metaData$Library_Barcode),]
  stopifnot(all(metaData$Library_Barcode==colnames(rawMatrix)))
  rownames(metaData) <- metaData$Library_Barcode
  rawMatrix <- CreateSeuratObject(counts = rawMatrix, meta.data = metaData, project = "Agarwal2020")
  rawMatrix
  
}, simplify=F)


# Initialize the Seurat object with the raw (non-normalized data).
names(toMerge) <- sapply(strsplit(files, "/"), function(x) x[12])
Agarwal2020 <- merge(toMerge[[1]], y=do.call("c",toMerge[2:length(toMerge)]))
Agarwal2020 <- JoinLayers(Agarwal2020)
stopifnot(all(colnames(Agarwal2020)==Agarwal2020@meta.data$Library_Barcode))

## genes expressed in less than 0.1% total cells are removed
selected_f <- names(which(rowSums(Agarwal2020$RNA["counts"]>0)>0.001*dim(Agarwal2020)[1]))
Agarwal2020 <- subset(Agarwal2020, features = selected_f)
Agarwal2020$system <- "PostMortem"


##############################
### Birtele et al (Foetal) ###
##############################

## Metadata pulled from Parmar Lab ".rds" Seurat object sent upon request.
tt <- readRDS("downloads/datasetsToIntegrate/Parmar/metadata/2_fetal_vm_week6_12.rds")
sampleList <- unique(tt$Sample)

parmarRef <- "downloads/datasetsToIntegrate/Parmar/GSE192405/"

parmarReference <- sapply(sampleList, function(x){
  
  fileToRead <- dir(parmarRef)[grepl(x, dir(parmarRef))]
  print(paste0("Processing file ", paste0(parmarRef, fileToRead)))
  tmp <- read.csv(paste0(parmarRef,fileToRead))
  #tmp <- as.data.frame(tmp)
  colnames(tmp) <- gsub("\\.","-", colnames(tmp))
  rownames(tmp) <- tmp[,1]
  tmp[,1] <- NULL
  
  meta <- tt[,colnames(tt) %in% colnames(tmp)]@meta.data
  
  parmarReference <- CreateSeuratObject(counts = tmp, meta.data=meta)
  stopifnot(all(colnames(tt)==rownames(tt@meta.data)))
  rm(meta)
  rm(tmp)
  gc()
  
  return(parmarReference)
  
}, simplify=F)


## merge
parmarReference <- merge(parmarReference[[1]], y=do.call("c",parmarReference[2:length(parmarReference)]),
                   project="ParmarFoetal")

parmarReference$PCW <- as.numeric(gsub("wks", "", sapply(strsplit(parmarReference$Sample, "-"), function(x) x[length(x)])))
parmarReference$PCW <- gsub(5,11, parmarReference$PCW)
## Label is wrong (11-5wks and 10-5wks) is aparently 11 weeks in the paper: pending explanation from authors

parmarReference <- JoinLayers(parmarReference)
## genes expressed in less than 0.1% total cells are removed
selected_f <- names(which(rowSums(parmarReference$RNA["counts"]>0)>0.001*dim(parmarReference)[1]))
parmarReference <- subset(parmarReference, features = selected_f)

parmarReference$AnnotType <- parmarReference$NamedClusters
parmarReference$AnnotType <- gsub(" ","", parmarReference$AnnotType)
parmarReference$AnnotType <- gsub("\\/","-", parmarReference$AnnotType)
parmarReference$AnnotType <- gsub("RG-1", "RG1-NPC", parmarReference$AnnotType)
parmarReference$AnnotType <- gsub("RG-2", "RG2", parmarReference$AnnotType)
parmarReference$AnnotType <- gsub("FPP", "RG3-FPP", parmarReference$AnnotType)
parmarReference$system <- "Foetal"



#############################
### Fiorenzano et al (3D) ###
#############################

fiorenzano <- readRDS("downloads/datasetsToIntegrate/Parmar/4_vm_organoids_d20_d30_d30_d60_d120.rds")
fiorenzano[["RNA5"]] <- as(object = fiorenzano[["RNA"]], Class = "Assay5")
fiorenzano[["RNA"]] <- as(object = fiorenzano[["RNA5"]], Class = "Assay5")
fiorenzano[["RNA5"]] <- NULL

selected_f <- names(which(rowSums(fiorenzano$RNA["counts"]>0)>0.001*dim(fiorenzano)[1]))
fiorenzano <- subset(fiorenzano, features = selected_f)
fiorenzano$system <- "Organoid"


#################
### Our query ###
#################


query <- readRDS("saved/toZenodo/mlo_resolution075_Annot.RDS")
clustAnnot <- c(0:23)
names(clustAnnot) <- c("hRgl2/immAstro","hNbDA","hProgFPM","OPC_1","VascLepto","hDA1b","hRgl1","hDA1a","hRgl3_caudal","hDA2","hProgM",
                       "hPreDA","hMidPre","hMgl","hEndo","hNbGaba","hNPro","hDA3/hGABA/hSer","Unk","hRgl4/MultiEpend","Astro","hPeric","Eryth","OPC_2")
query$seurat_clusters_24_Annot <- names(clustAnnot[query$seurat_clusters])


#################
#################
#################


## Create a list of datasets to integrate ##

obj <- list(query, leManno, smitsEtAl, Agarwal2020, parmarReference, fiorenzano, braun)
names(obj) <- c("MLOQuery","LeManno","Schwamborn","Agarwal","Birtele","Fiorenzano", "Braun")
gc()

obj <- merge(obj[[1]], y=do.call("c",obj[2:length(obj)]),
                         project="Integration")


## Release memory usage by removing the individual datasets ## 
rm(query)
rm(leManno)
rm(smitsEtAl)
rm(Agarwal2020)
rm(parmarReference)
rm(fiorenzano)
rm(braun)
gc()

#### Apply scRNA-seq pipeline steps to all datasets (Test for unintegrated approach) #### 
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
gc()
obj <- RunPCA(obj)
gc()

obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 0.75, cluster.name = "unintegrated_clusters")
gc()
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")


# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth

##defineVariablePerProject

obj$Dataset <- NA
obj$Dataset[obj$orig.ident=="MLO"] <- "MLO"
obj$Dataset[obj$orig.ident=="SmitsLetAl2022"] <- "Schwamborn"
obj$Dataset[grepl("standardorg", obj$orig.ident)] <- "Fiorenzano-ParmarOrg"
obj$Dataset[grepl("hVM", obj$orig.ident)] <- "Birtele-ParmarFoetal"
obj$Dataset[grepl("177", obj$orig.ident)] <- "LeManno"
obj$Dataset[grepl("N.+", obj$orig.ident)] <- "Agarwal"
obj$Dataset[grepl("10X", obj$orig.ident)] <- "Braun"



toPlot <- DimPlot(obj, reduction = "umap.unintegrated", group.by = c("Dataset"))
pdf(file="/home/jovyan/seurat_v5/testIntegration.pdf", width=8, height=6)
plot(toPlot)
dev.off()


## We first tried with CCA integration.

######################
### CCAIntegration ###
######################

obj <- IntegrateLayers(
  object = obj, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

gc()

obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.75, cluster.name = "cca_clusters")

obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

toPlot2 <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by = c("Dataset"))

pdf(file="/home/jovyan/seurat_v5/ccaIntegration.pdf", width=8, height=6)
plot(toPlot2)
dev.off()

toPlot3 <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by = c("Dataset","cca_clusters", "seurat_clusters_24_Annot"),
  combine = TRUE
)

pdf(file="/home/jovyan/seurat_v5/ccaIntegration.pdf", width=20, height=6)
plot(toPlot3)
dev.off()



##############################
## Harmonisation of columns ##
##############################

## After this initial integration, we renamed and cleaned the metadata from the multiple datasets to create an harmonised version.

## orig.ident (10x samples)
mask <- obj$orig.ident=="MLO"
obj$orig.ident[mask] <- obj$midBrainId[mask]

mask2 <- obj$orig.ident=="SmitsLetAl2022" & obj$timePoint==35
obj$orig.ident[mask2] <- "GSM3929395_org3wt35_S1"

mask3 <- obj$orig.ident=="SmitsLetAl2022" & obj$timePoint==70
obj$orig.ident[mask3] <- "GSM3929396_org1wt70II_S1"

  
## Model system
mask <- obj$old.ident=="MLO" & obj$originDimensions=="2D"
obj$system[mask] <- "2D"

mask2 <- obj$old.ident=="MLO" & obj$originDimensions=="3D"
obj$system[mask2] <- "Organoid"

mask3 <- obj$old.ident=="MLO" & obj$originDimensions=="foetal"
obj$system[mask3] <- "Foetal"

obj$origin <- NULL

##dropColumns
obj$midBrainId <- NULL
obj$midBrainIndex <- NULL


## timePoint (mixed)
mask3 <- obj$old.ident=="MLO" & obj$originDimensions=="foetal" & !is.na(obj$old.ident)
obj$timePoint[mask3] <- paste0("PCW", obj$timePoint[mask3])


#copymetaData <- obj@meta.data
mask4 <- obj$old.ident=="MLO" & ( obj$originDimensions=="2D" | obj$originDimensions=="3D" ) & !is.na(obj$old.ident)
obj$timePoint[mask4] <- paste0("day", obj$timePoint[mask4])

mask5 <- obj$Dataset=="Schwamborn" 
obj$timePoint[mask5] <- paste0("day", obj$timePoint[mask5])

mask6 <- obj$Dataset=="LeManno"
obj$timePoint[mask6] <- gsub("week_","GW", obj$Timepoint[mask6])

mask7 <- obj$Dataset=="Birtele-ParmarFoetal"
obj$timePoint[mask7] <- paste0("PCW",obj$PCW[mask7])

mask8 <- obj$Dataset=="Fiorenzano-ParmarOrg"
obj$timePoint[mask8] <- gsub("day0","day",obj$groups[mask8])

mask9 <- obj$Dataset=="Agarwal"
obj$timePoint[mask9] <- "PostMortem"


##dropColumns
obj$old.ident <- NULL
obj$RNA_snn_res.0.3 <- NULL
obj$RNA_snn_res.0.5 <- NULL
obj$predicted.celltype <- NULL
obj$predicted.celltype.score <- NULL
obj$seurat_clustersAnnot_SB <- NULL
obj$scTypeIntegrated <- NULL

##renameAnnotations

obj$annotation_MLO <- obj$seurat_clusters_24_Annot
obj$annotation_LeManno <- obj$CellType
obj$CellType <- NULL
obj$Timepoint <- NULL
obj$RNA_snn_res.2 <- NULL
obj$condition <- NULL
obj$type <- NULL
obj$Library <- NULL
obj$Library_Barcode <- NULL
obj$Brain_Region <- NULL
obj$batch_Agarwal <- obj$Batch
obj$Batch <- NULL
obj$annotation_Agarwal_level1_cellType <- obj$Level_1_cell_type
obj$annotation_Agarwal_level2_cellType <- obj$Level_2_cell_type
obj$Level_1_cell_type <- NULL
obj$Level_2_cell_type <- NULL
obj$Sample <- NULL

obj$cell <- NULL
obj$Files <- NULL
obj$Group <- NULL

obj$Cells.before.QC <- NULL
obj$Cells.after.QC <- NULL
obj$stress <- NULL
obj$RNA_snn_res.0.1 <- NULL

mask <- obj$orig.ident=="Birtele-ParmarFoetal"
stopifnot(all(obj$cell_id[mask]==colnames(obj)[mask]))
obj$cell_id <- NULL
obj$integrated_snn_res.0.1 <- NULL
obj$integrated_snn_res.0.4 <- NULL

maskBirtele <- obj$Dataset=="Birtele-ParmarFoetal"
obj$annotation_Birtele <- NA
obj$annotation_Birtele[maskBirtele] <- obj$NamedClusters[maskBirtele]
obj$annotation_Birtele2 <- obj$AnnotType

maskFiorenzano <- obj$Dataset=="Fiorenzano-ParmarOrg"
obj$annotation_Fiorenzano <- NA
obj$annotation_Fiorenzano[maskFiorenzano] <- obj$NamedClusters[maskFiorenzano]

maskBraun <- obj$Dataset=="Braun"
obj$annotation_Braun <- NA
obj$annotation_Braun[maskBraun] <- obj$CellClass[maskBraun]

obj$celltype <- NULL
obj$seurat.unc.rds <- NULL

mask <- obj$orig.ident=="Birtele-ParmarFoetal"
stopifnot(all(obj$shortIdent[mask]==colnames(obj)[mask]))
obj$shortIdent <- NULL
obj$first.labels <- NULL
obj$MajorLabel <- NULL
obj$NamedClusters <- NULL
obj$PCW <- NULL
obj$AnnotType <- NULL

obj$groups <- NULL
obj$RNA_snn_res.0.07 <- NULL
obj$CC.Difference<- NULL
obj$S.Score.log<- NULL
obj$G2M.Score.log<- NULL
obj$Cycling<- NULL



## create single Annot mix column at the end
obj$annotation_mixed <- NA
maskMLO <- obj$Dataset=="MLO"
maskLeManno <- obj$Dataset=="LeManno"
maskBirtele <- obj$Dataset=="Birtele-ParmarFoetal"
maskFiorenzano <- obj$Dataset=="Fiorenzano-ParmarOrg"
maskAgarwal <- obj$Dataset=="Agarwal"
maskBraun <- obj$Dataset=="Braun"


#maskSchwamborn (no cell annotation provided)

obj$annotation_mixed[maskMLO] <- obj$annotation_MLO[maskMLO]
obj$annotation_mixed[maskLeManno] <- obj$annotation_LeManno[maskLeManno]
obj$annotation_mixed[maskBirtele] <- obj$annotation_Birtele2[maskBirtele]
obj$annotation_mixed[maskFiorenzano] <- obj$annotation_Fiorenzano[maskFiorenzano]
obj$annotation_mixed[maskAgarwal] <- obj$annotation_Agarwal_level2_cellType[maskAgarwal]
obj$annotation_mixed[maskBraun] <- obj$CellClass[maskBraun]




## After the harmonisation, we run the harmony and RPCA integrations

######################
### HarmonyIntegra ###
######################


obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = TRUE
)

gc()

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.75, cluster.name = "harmony_clusters")

obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

gc()



##################
### RPCAIntegr ###
##################


obj <- IntegrateLayers(
  object = obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = TRUE
)

gc()
obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.75, cluster.name = "rpca_clusters")
obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
gc()



## Finally, we run the FASTMNN integration.

###############
### FastMNN ###
###############

library(SeuratWrappers)

obj <- IntegrateLayers(
  object = obj, method = FastMNNIntegration,
  orig.reduction = "pca", new.reduction = "integrated.mnn",
  verbose = TRUE
)

gc()
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.75, cluster.name = "mnn_clusters")
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")
gc()



## We visualised the dimensional reduction with UMAP for each integration.

ccaIntegration <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by = "Dataset",
  combine=FALSE)

harmonyIntegration <- DimPlot(
  obj,
  reduction = "umap.harmony",
  group.by = "Dataset",
  combine=FALSE)

rpcaIntegration <- DimPlot(
  obj,
  reduction = "umap.rpca",
  group.by = "Dataset",
  combine=FALSE)

fastmnnIntegration <- DimPlot(
  obj,
  reduction = "umap.mnn",
  group.by = "Dataset",
  combine=FALSE)

pdf(file="/home/jovyan/seurat_v5/allIntegrations.pdf", width=24, height=6)
wrap_plots(c(ccaIntegration, harmonyIntegration, rpcaIntegration, fastmnnIntegration), ncol=4, nrow=1)
dev.off()

saveRDS(obj, file="saved/toZenodo/midBrainIntegration.RDS")




