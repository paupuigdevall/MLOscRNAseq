library(slingshot)
library(uwot)
library(Seurat)
library(SingleCellExperiment)
library(RColorBrewer)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)


obj <- readRDS("saved/toZenodo/midBrainIntegration.RDS")

# for "midBrainDatasets_ccaIntegration_v3.RDS"
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


rareCtypes <- names(which(table(obj$annotation_unified)<50))
obj@assays$prediction.score.celltype <- NULL

obj <- obj[,is.na(match(obj$annotation_unified, rareCtypes))]
gc()

integratedMatrix <- obj@assays$mnn.reconstructed[]
obj@assays$mnn.reconstructed <- NULL
obj@assays$prediction.score.celltype <- NULL

obj.sce <- as.SingleCellExperiment(obj)
integrated.mnn <- reducedDim(obj.sce, "INTEGRATED.MNN")




### to consider: https://github.com/kstreet13/slingshot/issues/87
## Run slingshot on big datasets and specifying which annotation system and which starting cluster (FPP)

stopifnot(!any(is.nan(reducedDim(obj.sce, "PCA"))))
stopifnot(!any(is.nan(reducedDim(obj.sce, "PCA"))))
stopifnot(!any(is.na(obj.sce$annotation_unified)))

dim(cov(reducedDim(obj.sce,'INTEGRATED.MNN')))
integrated.mnn <- integrated.mnn[,1:49]
reducedDim(obj.sce, "INTEGRATED.MNN2") <- integrated.mnn


###############################################
### MNN reduced dim using FPP as start.clus ###
###############################################

sce <- slingshot(obj.sce, reducedDim = 'INTEGRATED.MNN2', clusterLabels = 'annotation_unified', start.clus = 'FPP')

####
####

slingshot_df <- colData(sce)$slingshot
maxWeightLineages <- colnames(pathStats(slingshot_df)$weights)[max.col(pathStats(slingshot_df)$weights)]

pseudotime_df <- pathStats(slingshot_df)$pseudotime
extracted_values <- pseudotime_df[cbind(seq_len(nrow(pseudotime_df)), match(maxWeightLineages, colnames(pseudotime_df)))]


pseudotimePerCell <- as.data.frame(setNames(extracted_values, rownames(pathStats(slingshot_df)$weights)))
colnames(pseudotimePerCell) <- "pseudotimeMaxWeight"
pseudotimePerCell$lineageMaxWeight <- maxWeightLineages

## QC
stopifnot(all(rownames(pseudotimePerCell)==rownames(colData(sce))))

pseudotimePerCell$annotation_unified <- colData(sce)$annotation_unified
pseudotimePerCell$annotation_mixed <- colData(sce)$annotation_mixed
pseudotimePerCell$samplesToPseudobulk <- colData(sce)$samplesToPseudobulk
pseudotimePerCell$samplesToPseudobulk2 <- colData(sce)$samplesToPseudobulk2
pseudotimePerCell$lineageMaxWeight <- gsub("Lineage","", pseudotimePerCell$lineageMaxWeight)
pseudotimePerCell$donorId <- colData(sce)$donorId


saveRDS(pseudotimePerCell, file="saved/pseudotime/pseudotimePerCellSlingshot.RDS")


####
####

rootDir <- "otherPlots/slingshot/"

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colourCount = length(unique(pseudotimePerCell$annotation_unified))
colVec <- setNames(getPalette(colourCount),
                   unique(pseudotimePerCell$annotation_unified))


plotUnified <- ggplot(pseudotimePerCell, aes(x = pseudotimeMaxWeight, y = annotation_unified, 
                                             colour = annotation_unified)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values=colVec) + theme_bw() +
  xlab("Slingshot pseudotime (corr. to Max weights)") + ylab("Cell types") +
  ggtitle("Cells ordered by Slingshot pseudotime")+
  theme(legend.position="none")

pdf(file=paste0(rootDir,paste0("slingshot_maxWeightPseudotime_annotUnified.pdf")), width=8, height = 8)
plot(plotUnified)
dev.off()

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colourCount2 = length(unique(pseudotimePerCell$annotation_mixed))
colVec2 <- setNames(getPalette(colourCount2),
                    unique(pseudotimePerCell$annotation_mixed))


plotMixed <- ggplot(pseudotimePerCell, aes(x = pseudotimeMaxWeight, y = annotation_mixed, 
                                           colour = annotation_mixed)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values=colVec2) + theme_bw() +
  xlab("Slingshot pseudotime (corr. to Max weights)") + ylab("Cell types") +
  ggtitle("Cells ordered by Slingshot pseudotime")+
  theme(legend.position="none")

pdf(file=paste0(rootDir,paste0("slingshot_maxWeightPseudotime_annotUnified_mixed.pdf")), width=8, height = 8)
plot(plotMixed)
dev.off()

library(stringr)

# Replace "3D" with "organoids" using stringr
pseudotimePerCell$samplesToPseudobulk2 <- as.character(pseudotimePerCell$samplesToPseudobulk2)
pseudotimePerCell$samplesToPseudobulk2[grepl("PostMortem", pseudotimePerCell$samplesToPseudobulk2)] <- "Agarwal et al. 2020 (PostMortem)"
pseudotimePerCell$samplesToPseudobulk2 <- str_replace_all(pseudotimePerCell$samplesToPseudobulk2, "Organoids", "3D")
pseudotimePerCell$samplesToPseudobulk2 <- factor(pseudotimePerCell$samplesToPseudobulk2,
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



plotUnified_Lineages <- ggplot(pseudotimePerCell, aes(x = pseudotimeMaxWeight, y = annotation_unified, 
                                                      colour = annotation_unified)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values=colVec) + theme_bw() +
  xlab("Slingshot pseudotime (corr. to Max weights)") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime per lineage")+
  facet_wrap(~lineageMaxWeight)+
  theme(legend.position="none")

pdf(file=paste0(rootDir,paste0("slingshot_maxWeightPseudotime_annotUnifiedLineages.pdf")), width=14, height = 10)
plot(plotUnified_Lineages)
dev.off()


## produce barplot of cell types per lineage + number of cells per lineage

library(dplyr)

# Compute counts of each "annotation_unified" cell type per "lineageMaxWeight" category
counts <- pseudotimePerCell %>%
  group_by(lineageMaxWeight, annotation_unified) %>%
  summarise(count = n())

# Compute total counts per "lineageMaxWeight" category
total_counts <- counts %>%
  group_by(lineageMaxWeight) %>%
  summarise(total_count = sum(count))

# Calculate the percentage of each "annotation_unified" cell type per "lineageMaxWeight" category
percentage <- as.data.frame(counts %>%
                              left_join(total_counts, by = "lineageMaxWeight") %>%
                              mutate(percentage = count / total_count * 100))

# Print the result

percentage$lineageMaxWeight <- factor(percentage$lineageMaxWeight, levels=sort(unique(as.numeric(percentage$lineageMaxWeight))))
percentage$annotation_unified <- factor(percentage$annotation_unified, levels=names(colVec))

barplotStckd <- ggplot(percentage, aes(fill=annotation_unified, x=lineageMaxWeight, y=percentage))+
  geom_bar(position="fill", stat="identity")+
  ggtitle("")+
  theme_bw()+
  theme(plot.title=element_text(size=18, face="bold", hjust=0.5),
        legend.title = element_text(size = 12, face="bold", hjust=0.5),
        legend.text  = element_text(size = 11),
        legend.key.size = unit(1, "lines"),
        axis.text.x = element_text(size=11, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=11),
        axis.title=element_text(size=12)) +
  guides(shape = guide_legend(override.aes = list(size = 8)),
         fill = guide_legend(override.aes = list(size = 8), ncol=1))+
  scale_fill_manual(name="Cell types",
                    values = colVec)+
  scale_y_continuous(labels = scales::percent, breaks=seq(0,1,0.2))+
  xlab("")+
  ylab("Cell type abundance [%]")


barplotStckdNoLegend <- barplotStckd + theme(legend.position="none")
barplotStckdOnlyLegend <- cowplot::get_legend(barplotStckd)


total_counts <- as.data.frame(total_counts)
total_counts$lineageMaxWeight <- factor(total_counts$lineageMaxWeight, levels=sort(unique(as.numeric(total_counts$lineageMaxWeight))))


barplotLines <- ggplot(total_counts, aes(x=lineageMaxWeight, y=total_count))+
  geom_bar(stat="identity", fill="black")+
  geom_text(aes(label=total_count), vjust=-0.3, size=3)+
  theme_bw()+
  theme(axis.text.x = element_text(size=9, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=8),
        axis.title=element_text(size=11),
        strip.text.x = element_text(size = 10))+
  ylab("Number of cells")+
  xlab("")+
  scale_y_continuous(labels = scales::comma, breaks=seq(0,30000,5000), limits=c(0,30000))


allplotsWithOutLegend <- ggarrange(barplotStckdNoLegend, barplotLines, nrow =2, heights=c(0.8,0.2))+bgcolor("white") 


allplotsWithLegend <- ggarrange(allplotsWithOutLegend, barplotStckdOnlyLegend, ncol=2, widths=c(0.85,0.15))+bgcolor("white")  
pdf(file=paste0(rootDir,"combinedPlot_barplot_nCells.pdf"), width=12, height = 10)
plot(allplotsWithLegend)
dev.off()

ggsave(plot = allplotsWithLegend,
       filename = paste0(rootDir,"combinedPlot_barplot_nCells.png"),
       height = 12, width =12, units = "in", dpi = 300,
       device = "png",limitsize = FALSE,bg="white")



pseudotimePerCell$lineageMaxWeight <- factor(pseudotimePerCell$lineageMaxWeight, levels=sort(as.numeric(unique(pseudotimePerCell$lineageMaxWeight))))

plotUnified_Lineages_pseudo <- ggplot(pseudotimePerCell, aes(x = pseudotimeMaxWeight, y = lineageMaxWeight, 
                                                             colour = annotation_unified)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values=colVec) + theme_bw() +
  xlab("Slingshot pseudotime (lineage with max weight)") + ylab("Lineages") +
  ggtitle("Cells ordered by Slingshot pseudotime per lineage")+
  facet_wrap(~samplesToPseudobulk2)+
  guides(colour = guide_legend(override.aes = list(size = 8)))+
  theme(strip.text = element_text(size=8))


pdf(file=paste0(rootDir,paste0("slingshot_unified_lineages_pseudo.pdf")), width=17, height = 10)
plot(plotUnified_Lineages_pseudo)
dev.off()

### produce barplot per lineage

# Compute counts of each "annotation_unified" cell type per "lineageMaxWeight" category
counts2 <- pseudotimePerCell %>%
  group_by(lineageMaxWeight, annotation_unified, samplesToPseudobulk2) %>%
  summarise(count = n())

# Compute total counts per "lineageMaxWeight" category
total_counts2 <- counts2 %>%
  group_by(lineageMaxWeight, samplesToPseudobulk2) %>%
  summarise(total_count = sum(count))

# Calculate the percentage of each "annotation_unified" cell type per "lineageMaxWeight" category
percentage2 <- as.data.frame(counts2 %>%
                               left_join(total_counts2, by = c("samplesToPseudobulk2","lineageMaxWeight")) %>%
                               mutate(percentage = count / total_count * 100))


percentage2$lineageMaxWeight <- factor(percentage2$lineageMaxWeight, levels=sort(as.numeric(unique(percentage2$lineageMaxWeight))))
percentage2$annotation_unified <- factor(percentage2$annotation_unified, levels=sort(unique(percentage2$annotation_unified)))


plotUnified_Lineages_pseudo_barplot <- ggplot(percentage2, aes(x = percentage, y = lineageMaxWeight, 
                                                               fill = annotation_unified)) +
  geom_bar(stat="identity", position="fill")+
  theme_bw()+
  ggtitle("Cells ordered by Slingshot pseudotime per lineage")+
  theme(plot.title=element_text(size=18, face="bold", hjust=0.5),
        legend.title = element_text(size = 12, face="bold", hjust=0.5),
        legend.text  = element_text(size = 11),
        legend.key.size = unit(1, "lines"),
        axis.text.x = element_text(size=11, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=11),
        axis.title=element_text(size=12)) +
  guides(shape = guide_legend(override.aes = list(size = 8)),
         fill = guide_legend(override.aes = list(size = 8), ncol=1))+
  scale_fill_manual(name="Cell types",
                    values = colVec)+
  scale_x_continuous(labels = scales::percent, breaks=seq(0,1,0.2))+
  xlab("Cell type composition per lineage") + ylab("Lineages") +
  facet_wrap(~samplesToPseudobulk2)+
  theme(strip.text = element_text(size=8))


pdf(file=paste0(rootDir,paste0("slingshot_unified_lineages_pseudo_barplot.pdf")), width=17, height = 10)
plot(plotUnified_Lineages_pseudo_barplot)
dev.off()


total_counts2 <- as.data.frame(total_counts2)

plotUnified_Lineages_pseudo_barplot_numCells <- ggplot(total_counts2, aes(x = total_count, y = lineageMaxWeight)) +
  geom_bar(stat="identity", fill="black")+
  theme_bw()+
  ggtitle("Cells ordered by Slingshot pseudotime per lineage")+
  theme(plot.title=element_text(size=18, face="bold", hjust=0.5),
        legend.title = element_text(size = 12, face="bold", hjust=0.5),
        legend.text  = element_text(size = 11),
        legend.key.size = unit(1, "lines"),
        axis.text.x = element_text(size=11, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=11),
        axis.title=element_text(size=12)) +
  #scale_x_continuous(labels = scales::percent, breaks=seq(0,1,0.2))+
  xlab("Number of cells per lineage") + ylab("Lineages") +
  facet_wrap(~samplesToPseudobulk2)+
  theme(strip.text = element_text(size=8))


pdf(file=paste0(rootDir,paste0("slingshot_unified_lineages_pseudo_barplot_numCells.pdf")), width=17, height = 10)
plot(plotUnified_Lineages_pseudo_barplot_numCells)
dev.off()

### produce pseudotime plot ignoring lineages (per sample)



getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colourCount = length(unique(pseudotimePerCell$annotation_unified))
colVec <- setNames(getPalette(colourCount),
                   unique(pseudotimePerCell$annotation_unified))

plotUnified_pseudo_samples <- ggplot(pseudotimePerCell, aes(x = pseudotimeMaxWeight, y = annotation_unified, 
                                                            colour = annotation_unified)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values=colVec) + theme_bw() +
  xlab("Slingshot pseudotime") + ylab("Cell types") +
  ggtitle("Cells ordered by Slingshot pseudotime")+
  facet_wrap(~samplesToPseudobulk2, scales="free_y")+
  theme(legend.position="none",
        strip.text = element_text(size=8))


pdf(file=paste0(rootDir,paste0("slingshot_pseudotime_samples.pdf")), width=20, height = 14)
plot(plotUnified_pseudo_samples)
dev.off()


### produce pseudotime plot ignoring lineages (per sample), barplot NCells

# Compute counts of each "annotation_unified" cell type per "samplesToPseudobulk2" category
counts3 <- as.data.frame(pseudotimePerCell %>%
                           group_by(annotation_unified, samplesToPseudobulk2) %>%
                           summarise(count = n()))

counts3$annotation_unified <- factor(counts3$annotation_unified, levels=sort(unique(counts3$annotation_unified)))

plotUnified_annot_barplot_numCells <- ggplot(counts3, aes(x = count, y = annotation_unified)) +
  geom_bar(stat="identity", fill="black")+
  theme_bw()+
  ggtitle("Number of cells per sample and annotation")+
  theme(plot.title=element_text(size=18, face="bold", hjust=0.5),
        legend.title = element_text(size = 12, face="bold", hjust=0.5),
        legend.text  = element_text(size = 11),
        legend.key.size = unit(1, "lines"),
        axis.text.x = element_text(size=9, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=9),
        axis.title=element_text(size=11)) +
  #scale_x_continuous(labels = scales::percent, breaks=seq(0,1,0.2))+
  xlab("Number of cells per cell type") + ylab("") +
  facet_wrap(~samplesToPseudobulk2, scales="free")+
  theme(strip.text = element_text(size=7.5))


pdf(file=paste0(rootDir,paste0("slingshot_unified_annot_barplot_numCells.pdf")), width=20, height = 12)
plot(plotUnified_annot_barplot_numCells)
dev.off()




### Plot per Sample (is this driving the bimodal projection in our pseudotime?)

maskOursInVitro <- grepl("This work \\(2D\\)", pseudotimePerCell$samplesToPseudobulk2) | grepl("This work \\(3D\\)", pseudotimePerCell$samplesToPseudobulk2)
pseudotimePerCell_subset <- pseudotimePerCell[maskOursInVitro,]

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colourCount = length(unique(pseudotimePerCell_subset$annotation_unified))
colVec <- setNames(getPalette(colourCount),
                   unique(pseudotimePerCell_subset$annotation_unified))

plotUnified_pseudo_samples2 <- ggplot(pseudotimePerCell_subset, aes(x = pseudotimeMaxWeight, y = annotation_unified, 
                                                                    colour = annotation_unified)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values=colVec) + theme_bw() +
  xlab("Slingshot pseudotime") + ylab("Cell types") +
  ggtitle("Cells ordered by Slingshot pseudotime")+
  facet_grid(donorId~samplesToPseudobulk2, scales="free_y")+
  theme(legend.position="none",
        strip.text = element_text(size=8))


pdf(file=paste0(rootDir,paste0("slingshot_pseudotime_samples_donorID.pdf")), width=20, height = 14)
plot(plotUnified_pseudo_samples2)
dev.off()


plotUnified_pseudo_samples3 <- ggplot(pseudotimePerCell_subset, aes(x = pseudotimeMaxWeight, y = annotation_unified, 
                                                                    colour = annotation_unified)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values=colVec) + theme_bw() +
  xlab("Slingshot pseudotime") + ylab("Cell types") +
  ggtitle("Cells ordered by Slingshot pseudotime")+
  facet_grid(lineageMaxWeight~samplesToPseudobulk2, scales="free_y")+
  theme(legend.position="none",
        strip.text = element_text(size=8))


pdf(file=paste0(rootDir,paste0("slingshot_pseudotime_samples_lineages.pdf")), width=20, height = 14)
plot(plotUnified_pseudo_samples3)
dev.off()





