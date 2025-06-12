
library(Seurat)
library(ggplot2)
library(gghighlight)
library(ggbeeswarm)
library(ggpubr)
library(RColorBrewer)
library(clustree)


querySeurat <- readRDS("saved/toZenodo/mlo_resolution075_Annot.RDS")

clustAnnot <- c(0:23)
names(clustAnnot) <- c("hRgl2/immAstro","hNbDA","hProgFPM","OPC_1","VLMC","hDA1b","hRgl1","hDA1a","hRgl3_caudal","hDA2","hProgM",
                       "hPreDA","hMidPre","hMgl","hEndo","hNbGaba","hNPro","hDA3/hGABA/hSer","Unk","hRgl4/MultiEpend","Astro","hPeric","Eryth","OPC_2")

### 

querySeurat$seurat_clusters_24_Annot <- names(clustAnnot[querySeurat$seurat_clusters])
querySeurat$toPlotAnnot <- querySeurat$seurat_clusters_24_Annot
querySeurat$toPlotAnnot <- gsub("_","", querySeurat$toPlotAnnot)


clustAnnot_simplified <- c("hProgFPM"="Early Midbrain Prog", "hProgM"="Highly prolif NSC", "hMidPre"="Early Midbrain Prog", "hRgl1"="Early Midbrain Prog", "hRgl4/MultiEpend"="Early Midbrain Prog",
                           "hPreDA"="Neuron Prog", "hNPro"="Neuron Prog",
                           "OPC1"="Late Midbrain Prog", "OPC2"="Late Midbrain Prog", "hRgl2/immAstro"="Late Midbrain Prog", "hRgl3caudal"="Late Midbrain Prog",
                           "hNbDA"="Immature Neurons",
                           "hNbGaba"="Immature Neurons",
                           "hDA1a"="Mature Neurons", "hDA1b"="Mature Neurons", "hDA2"="Mature Neurons", "hDA3/hGABA/hSer"="Mature Neurons",
                           "hMgl"="Microglia",
                           "Astro"="Astrocytes",
                           "VLMC"="Perivascular cells", "hPeric"="Perivascular cells", "hEndo"="Perivascular cells",
                           "Unk"="Others", "Eryth"="Others")

querySeurat$lowRes <- unname(clustAnnot_simplified[match(querySeurat$toPlotAnnot, names(clustAnnot_simplified))])
querySeurat$lowRes <- factor(querySeurat$lowRes,
                               levels=c(sort(unique(querySeurat$lowRes))[-match("Others", sort(unique(querySeurat$lowRes)))],"Others"))

querySeurat <- querySeurat[,!(grepl("Patient", querySeurat$donorDimensions) | grepl("Crispr", querySeurat$donorDimensions))]

colVec_final <- readRDS("saved/others/colVec_lowRes_final.RDS")
colVec <- readRDS("/scratch/project_2007686/midbrainOrganoids/saved/colVec_ctypesOriginal.RDS")




##### Creation of an h5ad object with log-normalised counts ##### 

library(SeuratDisk)

querySeurat@assays$prediction.score.celltype <- NULL
DefaultAssay(querySeurat) <- "RNA"
querySeurat <- querySeurat[,querySeurat$originDimensions=="foetal"]

VariableFeatures <- VariableFeatures(querySeurat)
tmp <- as.data.frame(VariableFeatures)

write.csv(tmp, "saved/scanpy/HVG.csv")
write.csv(Embeddings(querySeurat[["pca"]]), "saved/scanpy/pca_embeddings.csv")
write.csv(Loadings(querySeurat[["pca"]]), "saved/scanpy/pca_loadings.csv")

write.csv(Embeddings(querySeurat[["harmony"]]), "saved/scanpy/harmony_embeddings.csv")
write.csv(Loadings(querySeurat[["harmony"]]), "saved/scanpy/harmony_loadings.csv")

write.csv(Embeddings(querySeurat[["umap"]]), "saved/scanpy/umap_embeddings.csv")

querySeurat$toPlotAnnot <- as.character(querySeurat$toPlotAnnot)
querySeurat$lowRes <- as.character(querySeurat$lowRes)
querySeurat$old.ident <- as.character(querySeurat$midBrainId)

querySeurat_lognorm <- DietSeurat(querySeurat, scale.data=FALSE)

SaveH5Seurat(querySeurat_lognorm, file="saved/scanpy/mlo_with_new_annot_foetal.h5Seurat")
Convert("saved/scanpy/mlo_with_new_annot_foetal.h5Seurat", dest = "h5ad")



##### Creation of an h5ad object with raw counts #####


querySeurat@assays$RNA@data <- querySeurat@assays$RNA@counts

querySeurat_rawCounts <- DietSeurat(querySeurat, scale.data=FALSE)

SaveH5Seurat(querySeurat_rawCounts, file="saved/scanpy/mlo_with_new_annot_rawdata_foetal.h5Seurat")
Convert("saved/scanpy/mlo_with_new_annot_rawdata_foetal.h5Seurat", dest = "h5ad")

















