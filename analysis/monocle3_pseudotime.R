################
### Monocle3 ###
################


library(Seurat)
library(monocle3)

rootMain <- "figures/main/"
rootSupp <- "figures/supp/"
rootDir <- "otherPlots/monocle3/"

obj <- readRDS("saved/toZenodo/midBrainIntegration.RDS")


obj$Dataset2 <- obj$Dataset
obj$Dataset2[grepl("Agarwal", obj$Dataset2 )] <- "Agarwal et al. 2020 (PostMortem)"
obj$Dataset2[grepl("Birtele-ParmarFoetal", obj$Dataset2 )] <- "Birtele et al. 2022 (Foetal)"
obj$Dataset2[grepl("Fiorenzano-ParmarOrg", obj$Dataset2 )] <- "Fiorenzano et al. 2021 (Organoids)"
obj$Dataset2[grepl("LeManno", obj$Dataset2 )] <- "La Manno et al. 2016 (Foetal)"
obj$Dataset2[grepl("MLO", obj$Dataset2 )] <- "This work (2D, Organoids, Foetal)"
obj$Dataset2[grepl("Schwamborn", obj$Dataset2 )] <- "Zagare et al. 2022 (Organoids)"
obj$Dataset2[grepl("Braun", obj$Dataset2 )] <- "Braun et al. 2023 (Foetal)"


obj$Dataset2 <- factor(obj$Dataset2, levels=c("This work (2D, Organoids, Foetal)",
                                              "La Manno et al. 2016 (Foetal)",
                                              "Birtele et al. 2022 (Foetal)",
                                              "Braun et al. 2023 (Foetal)",
                                              "Fiorenzano et al. 2021 (Organoids)",
                                              "Agarwal et al. 2020 (PostMortem)",
                                              "Zagare et al. 2022 (Organoids)"))


allctypes <- names(table(obj$annotation_mixed))
namesUnified_allctypes <- c("Astrocytes", "Astrocytes","Astrocytes", "Astrocytes",
                            "Neurons","Neurons","Neurons","Endothelial/Pericytes",
                            "Erythrocytes","Erythrocytes","Fibroblasts","FPP",
                            "FPP","FPP","FPP","Neurons",
                            "Glioblasts","Neurons","Neurons","Neurons",
                            "Neurons","Neurons","Neurons","Neurons",
                            "Endothelial/Pericytes","Neurons","Microglia","Precursors",
                            "Neuroblasts","Neuroblasts","Neuroblasts","Neuroblasts",
                            "Neuroblasts","Progenitors","Progenitors","Neurons",
                            "OPC","Endothelial/Pericytes","Precursors","Progenitors",
                            "Progenitors","Progenitors","Progenitors","RadialGlia",
                            "RadialGlia","RadialGlia","RadialGlia","RadialGlia",
                            "RadialGlia","RadialGlia","RadialGlia","Neurons",
                            "Neurons","Immune","Microglia","Unknown",
                            "Neuroblasts","Neurons","IPC","ODC",
                            "ODC","ODC","ODC","OPC",
                            "OPC","OPC","Neurons","Endothelial/Pericytes",
                            "Progenitors","RadialGlia","Erythrocytes","RadialGlia",
                            "RadialGlia","RadialGlia","Unknown","VLMC",
                            "VLMC","VLMC")

allctypesUnified <- setNames(namesUnified_allctypes, allctypes)

obj$annotation_unified <- "NA"
obj$annotation_unified <- unname(allctypesUnified[obj$annotation_mixed])
obj$annotation_unified[is.na(obj$annotation_unified)] <- "Unknown"

obj$samplesToPseudobulk <- as.character(obj$Dataset2)
maskMLO <- obj$Dataset2=="This work (2D, Organoids, Foetal)"
mask2d <- obj$system=="2D"
mask3d <- obj$system=="Organoid"
maskfoet <- obj$system=="Foetal"

obj$samplesToPseudobulk[maskMLO & mask2d] <- "This work (2D)"
obj$samplesToPseudobulk[maskMLO & mask3d] <- "This work (3D)"
obj$samplesToPseudobulk[maskMLO & maskfoet] <- "This work (Foetal)"

obj <- obj[,!(grepl("Patient", obj$donorDimensions) | grepl("Crispr", obj$donorDimensions))]
gc()

obj$samplesToPseudobulk <- paste0(obj$timePoint, " - ", obj$samplesToPseudobulk)
obj$samplesToPseudobulk <- gsub("^(.)", "\\U\\1", obj$samplesToPseudobulk, perl = TRUE)

obj$samplesToPseudobulk <- factor(obj$samplesToPseudobulk,
                                  levels=c("Day40 - This work (2D)", "Day70 - This work (2D)",
                                           "Day40 - This work (3D)", "Day70 - This work (3D)", "Day120 - This work (3D)",
                                           "PCW10 - This work (Foetal)","PCW12 - This work (Foetal)","PCW16 - This work (Foetal)","PCW20 - This work (Foetal)",
                                           "GW6 - La Manno et al. 2016 (Foetal)",
                                           "GW7 - La Manno et al. 2016 (Foetal)",
                                           "GW8 - La Manno et al. 2016 (Foetal)",
                                           "GW9 - La Manno et al. 2016 (Foetal)",
                                           "GW10 - La Manno et al. 2016 (Foetal)",
                                           "GW11 - La Manno et al. 2016 (Foetal)",
                                           "PCW8 - Braun et al. 2023 (Foetal)", "PCW14 - Braun et al. 2023 (Foetal)",
                                           "PCW6 - Birtele et al. 2022 (Foetal)", "PCW8 - Birtele et al. 2022 (Foetal)", "PCW11 - Birtele et al. 2022 (Foetal)",
                                           "Day35 - Zagare et al. 2022 (Organoids)", "Day70 - Zagare et al. 2022 (Organoids)",
                                           "Day15 - Fiorenzano et al. 2021 (Organoids)", "Day30 - Fiorenzano et al. 2021 (Organoids)", "Day60 - Fiorenzano et al. 2021 (Organoids)",
                                           "Day90 - Fiorenzano et al. 2021 (Organoids)", "Day120 - Fiorenzano et al. 2021 (Organoids)",
                                           "PostMortem - Agarwal et al. 2020 (PostMortem)"))


#data <- as(as.matrix(GetAssayData(obj, assay = "mnn.reconstructed")), 'sparseMatrix')
pd <- data.frame(obj@meta.data)
#keep only the columns that are relevant
#pData <- pd %>% select(orig.ident, nCount_RNA, nFeature_RNA)
fData <- data.frame(gene_short_name = rownames(GetAssayData(obj, assay = "RNA")),
                    row.names = rownames(GetAssayData(obj, assay = "RNA")))

#Construct monocle cds
mn.obj <- new_cell_data_set(expression_data = GetAssayData(obj, assay = "RNA"), cell_metadata = pd, gene_metadata = fData)
#mn.sct <- monocle3::estimate_size_factors(mn.sct)
#rm(data)
mn.obj@int_colData@listData[["reducedDims"]][["UMAP"]] <- obj@reductions[["umap.mnn"]]@cell.embeddings
mn.obj <- monocle3::cluster_cells(mn.obj, resolution=1e-4, reduction_method = "UMAP")
p0 <- plot_cells(mn.obj, color_cells_by = "annotation_unified", show_trajectory_graph = FALSE)
p2 <- plot_cells(mn.obj, color_cells_by = "partition", show_trajectory_graph = FALSE)


pdf(file=paste0(rootDir,"clusterPartition.pdf"), width=6, height = 4)
p0+p2
dev.off()

mn.obj <- learn_graph(mn.obj, use_partition = TRUE, verbose = FALSE)


p1 <- plot_cells(mn.obj,
                 color_cells_by = "annotation_unified",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE)
pdf(file=paste0(rootDir,"umap_annotation_unified.pdf"), width=8, height = 6)
p1
dev.off()


p4 <- plot_cells(mn.obj,
                 color_cells_by = "timePoint",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 label_roots = FALSE,
                 rasterize = TRUE,
                 alpha=0.5)

pdf(file=paste0(rootDir,"DevTime.pdf"), width=12, height = 6)
p4
dev.off()


#FPP

get_earliest_principal_node <- function(cds, time_bin="FPP"){
  cell_ids <- which(colData(cds)[, "annotation_unified"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
mn.obj <- order_cells(mn.obj, root_pr_nodes=get_earliest_principal_node(mn.obj))


####################
## Supp Figure 8A ##
####################


figS8A <- plot_cells(mn.obj,
                 color_cells_by = "pseudotime",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 label_roots = FALSE)

pdf(file=paste0(rootSupp,"suppFigure8A.pdf"), width=8, height = 6)
plot(figS8A)
dev.off()


pdf(file=paste0(rootDir,"pseudoTimeCellTypes_unified.pdf"), width=12, height = 6)
p1+p3
dev.off()

pdf(file=paste0(rootDir,"pseudoTimeDevTime_unified.pdf"), width=12, height = 6)
p3+p4
dev.off()


## Calculate size factors using built-in function in monocle3
## Add gene names into CDS
mn.obj@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(obj[["RNA"]])
rowData(mn.obj)$gene_short_name <- rownames(mn.obj)
rownames(rowData(mn.obj)) <- NULL


mn.obj <- estimate_size_factors(mn.obj)
mn.obj@clusters@listData[["UMAP"]][["clusters"]] <- obj$annotation_mixed


set.seed(123)
geneChangePT <- graph_test(mn.obj, neighbor_graph="principal_graph")
saveRDS(geneChangePT, "saved/monocle3/geneChangePT_integratedlogNormCountsAllGenes.RDS")
deg_ids <- row.names(subset(geneChangePT, q_value < 0.05))

set.seed(123)
mn.obj <- preprocess_cds(mn.obj, num_dim = 50)

#set.seed(123)
gene_module_df <- find_gene_modules(mn.obj[deg_ids,], resolution=c(10^seq(-6,-1)), random_seed=123)
table(gene_module_df$module)

saveRDS(mn.obj, "saved/monocle3/monocle_integratedlogNormCountsAllGenes.RDS")
saveRDS(gene_module_df, "saved/monocle3/monocle_geneModules_integratedlogNormCountsAllGenes.RDS")


## Export Table S8
suppTable8 <- as.data.frame(gene_module_df)
colnames(suppTable8)[1] <- "cellId"

write.table(suppTable8, "saved/suppTables/TableS8.txt",
            quote=F, col.names=T, row.names=F, sep="\t")








## QC
stopifnot(all(names(pseudotime(mn.obj)) %in% colnames(mn.obj)))
colData(mn.obj)$pseudotime_df <- unname(pseudotime(mn.obj))



## add samplesToPseudobulk2
library(stringr)
colData(mn.obj)$samplesToPseudobulk <- str_replace_all(colData(mn.obj)$samplesToPseudobulk , "Organoids", "3D")


# Replace "3D" with "organoids" using stringr
colData(mn.obj)$samplesToPseudobulk2 <- as.character(colData(mn.obj)$samplesToPseudobulk)
colData(mn.obj)$samplesToPseudobulk2[grepl("PostMortem", colData(mn.obj)$samplesToPseudobulk2)] <- "Agarwal et al. 2020 (PostMortem)"
colData(mn.obj)$samplesToPseudobulk2 <- factor(colData(mn.obj)$samplesToPseudobulk2,
                                               levels=c("Day40 - This work (2D)", "Day70 - This work (2D)",
                                                        "Day40 - This work (3D)", "Day70 - This work (3D)", "Day120 - This work (3D)",
                                                        "PCW10 - This work (Foetal)","PCW12 - This work (Foetal)","PCW16 - This work (Foetal)","PCW20 - This work (Foetal)",
                                                        "GW6 - La Manno et al. 2016 (Foetal)",
                                                        "GW7 - La Manno et al. 2016 (Foetal)",
                                                        "GW8 - La Manno et al. 2016 (Foetal)",
                                                        "GW9 - La Manno et al. 2016 (Foetal)",
                                                        "GW10 - La Manno et al. 2016 (Foetal)",
                                                        "GW11 - La Manno et al. 2016 (Foetal)",
                                                        "PCW8 - Braun et al. 2023 (Foetal)", "PCW14 - Braun et al. 2023 (Foetal)",
                                                        "PCW6 - Birtele et al. 2022 (Foetal)", "PCW8 - Birtele et al. 2022 (Foetal)", "PCW11 - Birtele et al. 2022 (Foetal)",
                                                        "Day35 - Zagare et al. 2022 (3D)", "Day70 - Zagare et al. 2022 (3D)",
                                                        "Day15 - Fiorenzano et al. 2021 (3D)", "Day30 - Fiorenzano et al. 2021 (3D)", "Day60 - Fiorenzano et al. 2021 (3D)",
                                                        "Day90 - Fiorenzano et al. 2021 (3D)", "Day120 - Fiorenzano et al. 2021 (3D)",
                                                        "Agarwal et al. 2020 (PostMortem)"))

allctypes <- names(table(colData(mn.obj)$annotation_mixed))
namesUnified_allctypes <- c("Astrocytes", "Astrocytes","Astrocytes", "Astrocytes",
                            "Neurons","Neurons","Neurons","Endothelial/Pericytes",
                            "Erythrocytes","Erythrocytes","Fibroblasts","FPP",
                            "FPP","FPP","FPP","Neurons",
                            "Glioblasts","Neurons","Neurons","Neurons",
                            "Neurons","Neurons","Neurons","Neurons",
                            "Endothelial/Pericytes","Neurons","Microglia","Precursors",
                            "Neuroblasts","Neuroblasts","Neuroblasts","Neuroblasts",
                            "Neuroblasts","Progenitors","Progenitors","Neurons",
                            "OPC","Endothelial/Pericytes","Precursors","Progenitors",
                            "Progenitors","Progenitors","Progenitors","RadialGlia",
                            "RadialGlia","RadialGlia","RadialGlia","RadialGlia",
                            "RadialGlia","RadialGlia","RadialGlia","Neurons",
                            "Neurons","Immune","Microglia","Unknown",
                            "Neuroblasts","Neurons","IPC","ODC",
                            "ODC","ODC","ODC","OPC",
                            "OPC","OPC","Neurons","Endothelial/Pericytes",
                            "Progenitors","RadialGlia","Erythrocytes","RadialGlia",
                            "RadialGlia","RadialGlia","Unknown","VLMC",
                            "VLMC","VLMC")

allctypesUnified <- setNames(namesUnified_allctypes, allctypes)

colData(mn.obj)$annotation_unified  <- "NA"
colData(mn.obj)$annotation_unified <- unname(allctypesUnified[colData(mn.obj)$annotation_mixed])
colData(mn.obj)$annotation_unified[is.na(mn.obj$annotation_unified)] <- "Unknown"


pseudotime_df <- colData(mn.obj)[,c("samplesToPseudobulk2","annotation_unified","annotation_mixed","pseudotime_df")]

#rownames(pseudotime_df) <- colnames(mn.obj)
pseudotime_df$samplesToPseudobulk <- pseudotime_df$samplesToPseudobulk2

pseudotime_df$samplesToPseudobulk <- as.character(pseudotime_df$samplesToPseudobulk)
pseudotime_df$samplesToPseudobulk <- gsub(".+ - ","", pseudotime_df$samplesToPseudobulk)
pseudotime_df$pseudotime_ranked <- rank(pseudotime_df$pseudotime_df)
pseudotime_df <- pseudotime_df[!is.infinite(pseudotime_df$pseudotime_df),]


pseudotime_df$samplesToPseudobulk2 <- factor(pseudotime_df$samplesToPseudobulk2,
                                             levels=c("Day40 - This work (2D)", "Day70 - This work (2D)",
                                                      "Day40 - This work (3D)", "Day70 - This work (3D)", "Day120 - This work (3D)",
                                                      "PCW10 - This work (Foetal)","PCW12 - This work (Foetal)","PCW16 - This work (Foetal)","PCW20 - This work (Foetal)",
                                                      "GW6 - La Manno et al. 2016 (Foetal)",
                                                      "GW7 - La Manno et al. 2016 (Foetal)",
                                                      "GW8 - La Manno et al. 2016 (Foetal)",
                                                      "GW9 - La Manno et al. 2016 (Foetal)",
                                                      "GW10 - La Manno et al. 2016 (Foetal)",
                                                      "GW11 - La Manno et al. 2016 (Foetal)",
                                                      "PCW8 - Braun et al. 2023 (Foetal)", "PCW14 - Braun et al. 2023 (Foetal)",
                                                      "PCW6 - Birtele et al. 2022 (Foetal)", "PCW8 - Birtele et al. 2022 (Foetal)", "PCW11 - Birtele et al. 2022 (Foetal)",
                                                      "Day35 - Zagare et al. 2022 (3D)", "Day70 - Zagare et al. 2022 (3D)",
                                                      "Day15 - Fiorenzano et al. 2021 (3D)", "Day30 - Fiorenzano et al. 2021 (3D)", "Day60 - Fiorenzano et al. 2021 (3D)",
                                                      "Day90 - Fiorenzano et al. 2021 (3D)", "Day120 - Fiorenzano et al. 2021 (3D)",
                                                      "Agarwal et al. 2020 (PostMortem)"))



## plot pseudotime per cell type and sample


library(RColorBrewer)
library(ggthemes)
library(ggrastr)
library(base64enc)
library(ggplot2)  
library(readxl)
library(Biobase)
library(ggbeeswarm)
library(cowplot)
library(stringr)

pseudotime_df <- as.data.frame(pseudotime_df)
colVec2 <- readRDS(file="saved/colVec_ctypesColours.RDS")


setLast <- function(ctypesVec, lastCtype="Unk"){
  
  idLast <- match(lastCtype, ctypesVec)
  ordVec <- c(sort(ctypesVec[-idLast]), ctypesVec[idLast])
  return(ordVec)
}

pseudotime_df$annotation_unified <- factor(pseudotime_df$annotation_unified, levels=rev(setLast(unique(pseudotime_df$annotation_unified), lastCtype="Unknown")))


toPlotUnified_mn <- ggplot(pseudotime_df, 
                           aes(x = pseudotime_ranked, 
                               y = annotation_unified, colour = annotation_unified)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values=colVec) + theme_bw() +
  xlab("Monocle 3 ranked pseudotime [Num.cells]") +
  ylab("") +
  ggtitle("Cells ordered by pseudotime")+
  facet_wrap(~samplesToPseudobulk, scales="free_y")+
  theme(legend.position="none")


pdf(file=paste0(rootDir,paste0("diffusionPseudotimeRanked_MNN_samples.pdf")), width=16, height = 10)
plot(toPlotUnified_mn)
dev.off()


## The figure below is a panel with all the combinations for dataset, timepoint and model (as indicated by column "samplesToPseudobulk2").


####################
## Main Figure 4G ##
####################

pseudotime_df_subset <- subset(pseudotime_df, samplesToPseudobulk2=="PCW8 - Braun et al. 2023 (Foetal)" | samplesToPseudobulk2=="PCW14 - Braun et al. 2023 (Foetal)")

fig4G <- ggplot(pseudotime_df_subset, 
                              aes(x = pseudotime_ranked, 
                                  y = annotation_unified, colour = annotation_unified)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values=colVec2) + theme_bw() +
  xlab("Monocle 3 ranked pseudotime [Num.cells]") +
  ylab("") +
  ggtitle("Cells ordered by pseudotime")+
  facet_wrap(~samplesToPseudobulk2, scales="free")+
  theme(legend.position="none",
        strip.text.x=element_text(size=6),
        axis.text.x=element_text(size=6))+
  scale_x_continuous(limits=c(0,180000), breaks=seq(0,180000, 50000))


pdf(file=paste0(rootMain,"mainFigure4G.pdf"), width=8, height = 4.5)
plot(fig4G)
dev.off()


###############################
## Supplemental Figure 8C-8D ##
###############################

pseudotime_df_subset2 <- subset(pseudotime_df, samplesToPseudobulk2=="Day120 - This work (3D)" | samplesToPseudobulk2=="PCW10 - This work (Foetal)")

suppfig8CD <- ggplot(pseudotime_df_subset2, 
                aes(x = pseudotime_ranked, 
                    y = annotation_unified, colour = annotation_unified)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values=colVec2) + theme_bw() +
  xlab("Monocle 3 ranked pseudotime [Num.cells]") +
  ylab("") +
  ggtitle("Cells ordered by pseudotime")+
  facet_wrap(~samplesToPseudobulk2, scales="free")+
  theme(legend.position="none",
        strip.text.x=element_text(size=6),
        axis.text.x=element_text(size=6))+
  scale_x_continuous(limits=c(0,180000), breaks=seq(0,180000, 50000))


pdf(file=paste0(rootMain,"suppFigure8C8D.pdf"), width=8, height = 4.5)
plot(suppfig8CD)
dev.off()


## All the combinations (TP-model-Dataset)

toPlotUnified_mn_tp <- ggplot(pseudotime_df, 
                              aes(x = pseudotime_ranked, 
                                  y = annotation_unified, colour = annotation_unified)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values=colVec2) + theme_bw() +
  xlab("Monocle 3 ranked pseudotime [Num.cells]") +
  ylab("") +
  ggtitle("Cells ordered by pseudotime")+
  facet_wrap(~samplesToPseudobulk2, scales="free")+
  theme(legend.position="none",
        strip.text.x=element_text(size=6),
        axis.text.x=element_text(size=6))+
  scale_x_continuous(limits=c(0,180000), breaks=seq(0,180000, 50000))


pdf(file=paste0(rootDir,paste0("diffusionPseudotimeRanked_MNN_samplesTP.pdf")), width=18, height = 10)
plot(toPlotUnified_mn_tp)
dev.off()

## calculations in the paper

mean(subset(pseudotime_df, annotation_unified=="Neurons" & samplesToPseudobulk2=="PCW14 - Braun et al. 2023 (Foetal)")$pseudotime_ranked)
#[1] 130214.7
mean(subset(pseudotime_df, annotation_unified=="Neurons" & samplesToPseudobulk2=="PCW8 - Braun et al. 2023 (Foetal)")$pseudotime_ranked)
#[1] 109852.1



mean(subset(pseudotime_df, annotation_unified=="Progenitors" & samplesToPseudobulk2=="Day120 - This work (3D)")$pseudotime_ranked)
# [1] 15561.82
mean(subset(pseudotime_df, annotation_unified=="Progenitors" & samplesToPseudobulk2=="PCW10 - This work (Foetal)")$pseudotime_ranked)
# [1] 107952.2


mean(subset(pseudotime_df, annotation_unified=="Neurons" & samplesToPseudobulk2=="Day120 - This work (3D)")$pseudotime_ranked)
# [1] 145459.4
mean(subset(pseudotime_df, annotation_unified=="Neurons" & samplesToPseudobulk2=="PCW10 - This work (Foetal)")$pseudotime_ranked)
# [1] 108737.6


pseudotime_df[is.na(pseudotime_df$annotation_mixed),]$annotation_mixed <- "Unknown"


## Save table with pseudotime values
saveRDS(pseudotime_df, file="saved/pseudotime/pseudotimePerCellMonocle3.RDS")













