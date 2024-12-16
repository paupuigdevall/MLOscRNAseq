library(Seurat)
library(destiny)
library(slingshot)
library(conflicted)
library(scran)
library(purrr)
library(SingleCellExperiment)
library(ggthemes)
library(ggrastr)
library(base64enc)
library(ggplot2)  
library(readxl)
library(Biobase)
library(ggbeeswarm)
library(cowplot)
library(stringr)


obj <- readRDS("saved/toZenodo/midBrainIntegration.RDS")

allctypes <- names(table(obj$annotation_mixed))
namesUnified_allctypes <- c("Astrocytes", "Astrocytes","Astrocytes", "Astrocytes",
                            "DopaminergicN","DopaminergicN","DopaminergicN","Endothelial",
                            "RBC","RBC","Fibroblasts","FPP",
                            "FPP","FPP","FPP","GabaN",
                            "Glioblasts","SerotonergicN","DopaminergicN","DopaminergicN",
                            "DopaminergicN","DopaminergicN","DopaminergicN","hDA3/hGABA/hSer",
                            "Endothelial","GabaN","Microglia","Precursors",
                            "Neuroblasts","Neuroblasts","Neuroblasts","Neuroblasts",
                            "Neuroblasts","Progenitors","Progenitors","OMTN",
                            "OPC","Pericytes","Precursors","Progenitors",
                            "Progenitors","Progenitors","Progenitors","RadialGlia",
                            "RadialGlia/immAstro","RadialGlia","RadialGlia","RadialGlia",
                            "RadialGlia","RadialGlia","RadialGlia","hRN",
                            "SerotonergicN","Immune","Microglia","Unknown",
                            "Neuroblasts","Neurons","IPC","ODC",
                            "ODC","ODC","ODC","OPC",
                            "OPC","OPC","OMTN","Pericytes",
                            "Progenitors","RadialGlia","RBC","RadialGlia",
                            "RadialGlia","RadialGlia","Unknown","VLMC",
                            "VLMC","VLMC")

allctypesUnified <- setNames(namesUnified_allctypes, allctypes)

obj$annotation_unified <- "NA"
obj$annotation_unified <- base::unname(allctypesUnified[obj$annotation_mixed])
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

obj$samplesToPseudobulk2 <- paste0(obj$timePoint, " - ", obj$samplesToPseudobulk)
obj$samplesToPseudobulk2 <- gsub("^(.)", "\\U\\1", obj$samplesToPseudobulk2, perl = TRUE)

obj$samplesToPseudobulk2 <- factor(obj$samplesToPseudobulk2,
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

obj$samplesToPseudobulk <- str_replace_all(obj$samplesToPseudobulk, "Organoids", "3D")


# Replace "3D" with "organoids" using stringr
obj$samplesToPseudobulk2 <- as.character(obj$samplesToPseudobulk2)
obj$samplesToPseudobulk2[grepl("PostMortem", obj$samplesToPseudobulk2)] <- "Agarwal et al. 2020 (PostMortem)"
obj$samplesToPseudobulk2 <- str_replace_all(obj$samplesToPseudobulk2, "Organoids", "3D")
obj$samplesToPseudobulk2 <- factor(obj$samplesToPseudobulk2,
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



## Convert to SingleCellExperiment
## Need to diet those matrix of counts that do not fit the dimensions of the lognorm counts matrix

integratedMatrix <- obj@assays$mnn.reconstructed["data"]
obj@assays$mnn.reconstructed <- NULL
obj@assays$prediction.score.celltype <- NULL

obj.sce <- as.SingleCellExperiment(obj)
integrated.mnn <- reducedDim(obj.sce, "INTEGRATED.MNN")

obj.sce$MNN1 <- integrated.mnn[, 1]
obj.sce$MNN2 <- integrated.mnn[, 2]

library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

colourCount = length(unique(colData(obj.sce)$annotation_unified))
colVec <- setNames(getPalette(colourCount),
                   unique(colData(obj.sce)$annotation_unified))


## Plot the cells taking into account the two main MNN components (MNN1, MNN2)

toPlot <- ggplot(as.data.frame(colData(obj.sce)), aes(x = MNN1, y = MNN2, color = annotation_unified)) + geom_quasirandom_rast(groupOnX = FALSE) +
  scale_color_manual(values=colVec) + theme_classic() +
  xlab("MNN1") + ylab("MNN2") + ggtitle("MNN biplot")

rootDir <- "otherPlots/destiny/"

pdf(file=paste0(rootDir,paste0("mnn1_mnn2_plot.pdf")), width=8, height = 8)
plot(toPlot)
dev.off()

## Cells ordered by first MNN component (split by cell type - unified annotation - and dataset)

obj.sce$pseudotime_MNN1 <- rank(obj.sce$MNN1)  # rank cells by their MNN1 score
toPlot2 <- ggplot(as.data.frame(colData(obj.sce)), aes(x = pseudotime_MNN1, y = annotation_unified,
                                                       colour = annotation_unified)) +
  geom_quasirandom_rast(groupOnX = FALSE) +
  scale_color_manual(values=colVec) + theme_classic() +
  xlab("MNN11") + ylab("") +
  ggtitle("Cells ordered by first MNN component")+
  theme(legend.position = "none")+
  facet_wrap(~samplesToPseudobulk, scales="free_y")

pdf(file=paste0(rootDir,paste0("pseudotime_mnn1.pdf")), width=20, height = 14)
plot(toPlot2)
dev.off()


#####
#####

maskNA <- is.na(colData(obj.sce)$annotation_mixed)
colData(obj.sce)$annotation_mixed[maskNA] <- "Unknown"

colourCount = length(unique(colData(obj.sce)$annotation_mixed))
colVec <- setNames(getPalette(colourCount),
                   unique(colData(obj.sce)$annotation_mixed))

## Cells ordered by first MNN component (split by cell type - dataset specific annotation - and dataset)

obj.sce$pseudotime_MNN1 <- rank(obj.sce$MNN1)  # rank cells by their MNN1 score
toPlot3 <- ggplot(as.data.frame(colData(obj.sce)),
                  aes(x = pseudotime_MNN1, y = annotation_mixed,
                      colour = annotation_mixed)) +
  geom_quasirandom_rast(groupOnX = FALSE) +
  scale_color_manual(values=colVec) + theme_classic() +
  xlab("MNN11") + ylab("") +
  ggtitle("Cells ordered by first MNN component")+
  theme(legend.position = "none")+
  facet_wrap(~samplesToPseudobulk, scales="free_y")

pdf(file=paste0(rootDir,paste0("pseudotime_mnn1_ctypes_mixed.pdf")), width=20, height = 14)
plot(toPlot3)
dev.off()


### integrated-mnn (T3)
dm <- DiffusionMap(reducedDims(obj.sce)$INTEGRATED.MNN)
dpt <- DPT(dm)

## "dm" contains the results for the eigenvectors (also referred as diffusion components) and the 
## eigenvalues of the diffusion distance matrix (diffusion components importance)

## "dpt" is used to calculate the diffusion pseudo time


colourCount = length(unique(colData(obj.sce)$annotation_unified))
colVec <- setNames(getPalette(colourCount),
                   unique(colData(obj.sce)$annotation_unified))

## Data are samples from a diffusion process
## Get the first two diffusion components

## Destiny: Diffusion component 1 vs diffusion component 2 (shown as by cell type)


tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  Ctype = obj.sce$annotation_unified)
test <- ggplot(tmp, aes(x = DC1, y = DC2, colour = Ctype)) +
  geom_point() + scale_color_manual(values=colVec) +
  xlab("Diffusion component 1") +
  ylab("Diffusion component 2") +
  theme_bw()+
  facet_wrap(~Ctype)+
  theme(legend.position="none")


pdf(file=paste0(rootDir,"diffusionMap_MNN_MLOFoetal_ctypes.pdf"), width=10, height = 8)
plot(test)
dev.off()


## Destiny: Diffusion component 1 vs diffusion component 2 (shown as by dataset)

tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  Ctype = obj.sce$annotation_unified,
                  samplesToPseudobulk=obj.sce$samplesToPseudobulk)
test <- ggplot(tmp, aes(x = DC1, y = DC2, colour = Ctype)) +
  geom_point() + scale_color_manual(values=colVec) +
  xlab("Diffusion component 1") +
  ylab("Diffusion component 2") +
  theme_bw()+
  facet_wrap(~samplesToPseudobulk)


pdf(file=paste0(rootDir,paste0("diffusionMap_MNN_MLOFoetal_samples.pdf")), width=10, height = 8)
plot(test)
dev.off()


tmp2 <- data.frame(DC1 = eigenvectors(dm)[, 1],
                   DC2 = eigenvectors(dm)[, 2],
                   DC3 = eigenvectors(dm)[, 3],
                   Ctype = obj.sce$annotation_unified,
                   samplesToPseudobulk=obj.sce$samplesToPseudobulk,
                   samplesToPseudobulk2=obj.sce$samplesToPseudobulk2)

diffcompPlot <- ggplot(tmp2, aes(x = DC1, y = DC2, colour = Ctype)) +
  geom_point() + scale_color_manual(values=colVec) +
  xlab("Diffusion component 1") +
  ylab("Diffusion component 2") +
  theme_bw()+theme(legend.position="none")+
  facet_grid(Ctype~samplesToPseudobulk)

pdf(file=paste0(rootDir,paste0("diffusionMap_MNN_MLOFoetal_ctypeAndsamples.pdf")), width=16, height = 20)
plot(diffcompPlot)
dev.off()

diffcompPlot_dc3 <- ggplot(tmp2, aes(x = DC1, y = DC3, colour = Ctype)) +
  geom_point() + scale_color_manual(values=colVec) +
  xlab("Diffusion component 1") +
  ylab("Diffusion component 3") +
  theme_bw()+theme(legend.position="none")+
  facet_grid(Ctype~samplesToPseudobulk)

pdf(file=paste0(rootDir,paste0("diffusionMap_MNN_MLOFoetal_ctypeAndsamples_DC3.pdf")), width=16, height = 20)
plot(diffcompPlot_dc3)
dev.off()


## Calculate first diffusion component (DC1)

obj.sce$pseudotime_MNN1 <- rank(eigenvectors(dm)[,1])  # rank cells by their MNN1 score
obj.sce$pseudotime_MNN1_noranked <- eigenvectors(dm)[,1]  


## Plot cells ordered by DC1 by unified annot
colourCount = length(unique(colData(obj.sce)$annotation_unified))
colVec <- setNames(getPalette(colourCount),
                   unique(colData(obj.sce)$annotation_unified))

toPlot2 <- ggplot(as.data.frame(colData(obj.sce)), aes(x = pseudotime_MNN1, y = annotation_unified,
                                                       colour = annotation_unified)) +
  geom_quasirandom_rast(groupOnX = FALSE) +
  scale_color_manual(values=colVec) +
  xlab("DC1") + ylab("") + theme_bw() +
  ggtitle("Cells ordered by first diffusion component")+
  theme(legend.position = "none",
        strip.text.x = element_text(size=7))+
  facet_wrap(~samplesToPseudobulk2, scales="free_y")

pdf(file=paste0(rootDir,paste0("pseudotime_mnn1_samplesTP.pdf")), width=20, height = 14)
plot(toPlot2)
dev.off()


toSub <- as.data.frame(colData(obj.sce))
toSumToAllData <- -min(toSub$pseudotime_MNN1_noranked)+0.0001
toSub$transformed_pt <- toSub$pseudotime_MNN1_noranked+toSumToAllData

toPlot2b <- ggplot(toSub, aes(x = transformed_pt, y = annotation_unified,
                              colour = annotation_unified)) +
  geom_quasirandom_rast(groupOnX = FALSE) +
  scale_color_manual(values=colVec) +
  xlab("DC1") + ylab("") + theme_bw() +
  ggtitle("Cells ordered by first diffusion component")+
  theme(legend.position = "none",
        strip.text.x = element_text(size=7))+
  facet_wrap(~samplesToPseudobulk2, scales="free_y")+
  scale_x_log10()

pdf(file=paste0(rootDir,paste0("pseudotime_mnn1_noRanked_samplesTP2.pdf")), width=20, height = 14)
plot(toPlot2b)
dev.off()


### Plot cells ordered by DC1 by mixed annot
colourCount = length(unique(colData(obj.sce)$annotation_mixed))
colVec <- setNames(getPalette(colourCount),
                   unique(colData(obj.sce)$annotation_mixed))

toPlot2_mixed <- ggplot(as.data.frame(colData(obj.sce)), aes(x = pseudotime_MNN1, y = annotation_mixed,
                                                             colour = annotation_mixed)) +
  geom_quasirandom_rast(groupOnX = FALSE) +
  scale_color_manual(values=colVec) +
  xlab("DC1") + ylab("") + theme_bw() + 
  ggtitle("Cells ordered by first diffusion component")+
  theme(legend.position = "none",
        strip.text.x = element_text(size=7))+
  facet_wrap(~samplesToPseudobulk2, scales="free_y")

pdf(file=paste0(rootDir,paste0("pseudotime_mnn1_samplesTP_mixedAnnot.pdf")), width=20, height = 14)
plot(toPlot2_mixed)
dev.off()


## Visualization of >2 diffusion components
plotQplot <- qplot(y = eigenvalues(dm)) + theme_minimal() +
  labs(x = 'Diffusion component (DC)', y = 'Eigenvalue')

pdf(file=paste0(rootDir,paste0("elbowPlot_diffComponents.pdf")), width=6, height = 4)
plot(plotQplot)
dev.off()

colourCount = length(unique(colData(obj.sce)$annotation_unified))
colVec <- setNames(getPalette(colourCount),
                   unique(colData(obj.sce)$annotation_unified))


annot <- base::unname(colData(obj.sce)$annotation_unified)

library(Biobase)
kk <- dataset(dm)
kk <- cbind(dataset(dm), annotation_unified=annot)
kk <- as.ExpressionSet(as.data.frame(kk))
dataset(dm) <- kk

phenoData(dataset(dm))$annotation_unified <- factor(phenoData(dataset(dm))$annotation_unified,
                                                    levels=unique(phenoData(dataset(dm))$annotation_unified))
#palette(cube_helix(colourCount))
par(mar = c(5, 5, 5, 20)) 
pdf(file=paste0(rootDir,paste0("severalBranches_diffComponents.pdf")), width=10, height = 6)
plot(dm, pch=20, pal=colVec, col_by= 'annotation_unified', legend_main = '')
dev.off()



# Plot DC1 vs DC2 and color the cells by their inferred diffusion pseudotime.
# We can accesss diffusion pseudotime via dpt$dpt.


df <- data.frame(DC1 = eigenvectors(dm)[, 1], DC2 = eigenvectors(dm)[, 2], 
                 dptval = dpt$dpt, cell_type2 = obj.sce$annotation_unified)
p1 <- ggplot(df) + geom_point(aes(x = DC1, y = DC2, color = dptval))+theme_bw()
p2 <- ggplot(df) + geom_point(aes(x = DC1, y = DC2, color = cell_type2))+theme_bw()
p <- plot_grid(p1, p2, rel_widths = c(0.4,0.6))

pdf(file=paste0(rootDir,paste0("diffusionComponents_diffusionPST_MNN.pdf")), width=12, height = 6)
plot(p)
dev.off()



### Calculate now the diffusion pseudotime

obj.sce$pseudotime_dpt <- rank(dpt$dpt) 
toPlotUnified <- ggplot(as.data.frame(colData(obj.sce)), 
                        aes(x = pseudotime_dpt, 
                            y = annotation_unified, colour = annotation_unified)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values=colVec) + theme_bw() +
  xlab("Diffusion map pseudotime (dpt)") +
  ylab("Timepoint") +
  ggtitle("Cells ordered by diffusion map pseudotime")+
  theme(legend.position="none")


pdf(file=paste0(rootDir,paste0("diffusionPseudotimeRanked_MNN.pdf")), width=12, height = 8)
plot(toPlotUnified)
dev.off()

colourCount = length(unique(colData(obj.sce)$annotation_mixed))
colVec <- setNames(getPalette(colourCount),
                   unique(colData(obj.sce)$annotation_mixed))


toPlotUnified2 <- ggplot(as.data.frame(colData(obj.sce)), 
                         aes(x = pseudotime_dpt, 
                             y = annotation_mixed, colour = annotation_mixed)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values=colVec) + theme_bw() +
  xlab("Diffusion map pseudotime (dpt)") +
  ylab("Timepoint") +
  ggtitle("Cells ordered by diffusion map pseudotime")+
  facet_wrap(~samplesToPseudobulk, scales="free_y")+
  theme(legend.position="none")


pdf(file=paste0(rootDir,paste0("diffusionPseudotimeRanked_MNN_samples.pdf")), width=16, height = 10)
plot(toPlotUnified2)
dev.off()



#####
#####


obj.sce$pseudotime_MNN1_ordCells <- rank(eigenvectors(dm)[,1])  # rank cells by their MNN1 score
obj.sce$pseudotime_MNN1_raw <- eigenvectors(dm)[,1]  

obj.sce$pseudotime_dpt_ordCells <- rank(dpt$dpt) 
obj.sce$pseudotime_dpt_raw <- dpt$dpt



pseudotimePerCellDestiny <- as.data.frame(colData(obj.sce))
pseudotimePerCellDestiny <- pseudotimePerCellDestiny[,c("MNN1","pseudotime_MNN1_raw","pseudotime_MNN1_ordCells","pseudotime_dpt_raw","pseudotime_dpt_ordCells")]
saveRDS(pseudotimePerCellDestiny, file="saved/pseudotime/pseudotimePerCellDestiny.RDS")







































