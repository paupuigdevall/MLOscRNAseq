
library(Seurat)
library(monocle3)
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
library(ggridges)
library(tidyverse)


rootMain <- "figures/main/"
rootSupp <- "figures/supp/"
rootDir <- "otherPlots/monocle3/"


setLast <- function(ctypesVec, lastCtype="Unk"){
  
  idLast <- match(lastCtype, ctypesVec)
  ordVec <- c(sort(ctypesVec[-idLast]), ctypesVec[idLast])
  return(ordVec)
}

# mn.obj <- readRDS("/lustre/scratch126/cellgen/kilpinen/pp9/scRNA-seq/midbrainOrganoid/monocle3/saved/monocle_integratedlogNormCountsAllGenes.RDS")
# gene_module_df <- readRDS("/lustre/scratch126/cellgen/kilpinen/pp9/scRNA-seq/midbrainOrganoid/monocle3/saved/monocle_geneModules_integratedlogNormCountsAllGenes.RDS")
mn.obj <- readRDS("saved/monocle3/monocle_integratedlogNormCountsAllGenes.RDS")
gene_module_df <- readRDS("saved/monocle3/monocle_geneModules_integratedlogNormCountsAllGenes.RDS")


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





pseudotime_df <- as.data.frame(pseudotime_df)

colVec <- readRDS(file="saved/colVec_datasetColours.RDS")
colVec2 <- readRDS(file="saved/colVec_ctypesColours.RDS")

pseudotime_df$annotation_unified <- factor(pseudotime_df$annotation_unified, levels=rev(setLast(unique(pseudotime_df$annotation_unified), lastCtype="Unknown")))
pseudotime_df[is.na(pseudotime_df$annotation_mixed),]$annotation_mixed <- "Unknown"


###########

####################
## Main Figure 4H ##
####################


pseudotime_df$pseudotime_ranked_norm <- pseudotime_df$pseudotime_ranked*100/(max(pseudotime_df$pseudotime_ranked))
pseudotime_df_unified_DA <- subset(pseudotime_df, annotation_unified=="Neurons")
pseudotime_df_unified_DA_ours <- pseudotime_df[grepl("This work", pseudotime_df$samplesToPseudobulk2) & pseudotime_df$annotation_unified=="Neurons",]


### DA neurons for all
pseudotime_df_unified_DA$colorDatasets <- pseudotime_df_unified_DA$samplesToPseudobulk
pseudotime_df_unified_DA$colorDatasets[grepl("This work", pseudotime_df_unified_DA$colorDatasets)] <- "This work (2D, 3D, Foetal)"
pseudotime_df_unified_DA$tpPoint <- paste0(gsub(".*\\(([^)]+)\\).*", "\\1", pseudotime_df_unified_DA$samplesToPseudobulk2),"-",
                                           gsub(" ", "", sapply(strsplit(as.character(pseudotime_df_unified_DA$samplesToPseudobulk2), "-"), function(x) x[1])))
pseudotime_df_unified_DA$tpPoint <- gsub("PostMortem-.+","PostMortem", pseudotime_df_unified_DA$tpPoint)
pseudotime_df_unified_DA$tpPoint <- gsub("Day","day", pseudotime_df_unified_DA$tpPoint)

maskThisStudy <- grepl("This work", pseudotime_df_unified_DA$samplesToPseudobulk2)
pseudotime_df_unified_DA[!maskThisStudy,]$tpPoint <- gsub("Foetal-","",pseudotime_df_unified_DA[!maskThisStudy,]$tpPoint)

##remove LeManno (few number of cells)
pseudotime_df_unified_DA_nLM <- pseudotime_df_unified_DA[!grepl("La Manno", pseudotime_df_unified_DA$samplesToPseudobulk),]

pseudotime_df_unified_DA_nLM$tpPoint <- paste0(gsub("This","This work",gsub(" .+", "", pseudotime_df_unified_DA_nLM$colorDatasets)),"-", pseudotime_df_unified_DA_nLM$tpPoint)

pseudotime_df_unified_DA_nLM$tpPoint <- factor(pseudotime_df_unified_DA_nLM$tpPoint,
                                           levels=rev(c("This work-2D-day40","This work-2D-day70","This work-3D-day40","This work-3D-day70","This work-3D-day120",
                                                        "This work-Foetal-PCW10","This work-Foetal-PCW12","This work-Foetal-PCW16","This work-Foetal-PCW20",
                                                        "Birtele-PCW6", "Birtele-PCW8","Birtele-PCW11",
                                                        "Braun-PCW8","Braun-PCW14",
                                                        "Fiorenzano-3D-day15","Fiorenzano-3D-day30","Fiorenzano-3D-day60","Fiorenzano-3D-day90","Fiorenzano-3D-day120",
                                                        "Agarwal-PostMortem")))


tt <- ggplot(pseudotime_df_unified_DA_nLM, aes(x = pseudotime_ranked, y = tpPoint, fill = colorDatasets)) +
  theme_bw()+
  ggridges::geom_density_ridges(scale = 2, show.legend = FALSE, alpha=0.9) +
  scale_x_continuous(name = "Monocle3 ranked cell position",
                     limits = c(0, 180000)) +
  # control space at top and bottom of plot
  scale_y_discrete(name = "", expand = c(0.02, 0, .08, 0)) +
  scale_fill_manual(values=colVec) # colourblind-safe colours


distributionsToTest <- list(c("This work-2D-day40","This work-2D-day70"),
                            c("This work-3D-day40","This work-3D-day120"),
                            c("This work-Foetal-PCW10","This work-Foetal-PCW20"),
                            c("Birtele-PCW6","Birtele-PCW11"),
                            c("Braun-PCW8","Braun-PCW14"),
                            c("Fiorenzano-3D-day15","Fiorenzano-3D-day120"))



library(ggsignif)
library(ggpubr)

fig4H <- tt + stat_compare_means(
  comparisons = distributionsToTest, aes(label = ..p.signif..),
  step.increase = 0, size=4, label.x = 50, label.y = 50, color="blue")

pdf(file=paste0(rootMain,"mainFigure4H.pdf"), width=4.5, height = 6)
plot(fig4H)
dev.off()



### DA neurons subtypes from ours (ridge plot), vertical display

pseudotime_df_unified_DA_ours$tpPoint <- paste0(gsub(".*\\(([^)]+)\\).*", "\\1", pseudotime_df_unified_DA_ours$samplesToPseudobulk2),"-",
                                                gsub(" ", "", sapply(strsplit(as.character(pseudotime_df_unified_DA_ours$samplesToPseudobulk2), "-"), function(x) x[1])))

pseudotime_df_unified_DA_ours$tpPointAnnot <- paste0(pseudotime_df_unified_DA_ours$annotation_mixed,"-",pseudotime_df_unified_DA_ours$tpPoint)


toDiscard <- names(which(table(pseudotime_df_unified_DA_ours$tpPointAnnot)<50))
pseudotime_df_unified_DA_ours <- pseudotime_df_unified_DA_ours[!pseudotime_df_unified_DA_ours$tpPointAnnot %in% toDiscard,]


## breaks by ctype (recover colVec from figure 3)


####################
## Main Figure 4I ##
####################


colVec0 <- readRDS(file="saved/colVec_ctypesOriginal.RDS")

pseudotime_df_unified_DA_ours$tpPointAnnot <- factor(pseudotime_df_unified_DA_ours$tpPointAnnot,
                                               levels=rev(c("hDA1a-2D-Day70",
                                                            "hDA1a-3D-Day40","hDA1a-3D-Day70","hDA1a-3D-Day120","hDA1a-Foetal-PCW10",
                                                            "hDA1b-2D-Day40","hDA1b-2D-Day70","hDA1b-3D-Day40","hDA1b-3D-Day70","hDA1b-3D-Day120",
                                                            "hDA1b-Foetal-PCW10","hDA1b-Foetal-PCW12",
                                                            "hDA2-2D-Day40","hDA2-2D-Day70","hDA2-3D-Day40","hDA2-3D-Day70","hDA2-3D-Day120",
                                                            "hDA2-Foetal-PCW10","hDA2-Foetal-PCW12",
                                                            "hDA3/hGABA/hSer-Foetal-PCW10","hDA3/hGABA/hSer-Foetal-PCW12")))


tt2 <- ggplot(pseudotime_df_unified_DA_ours, aes(x = pseudotime_ranked, y = tpPointAnnot, fill = annotation_mixed)) +
  theme_bw()+
  ggridges::geom_density_ridges(scale = 2, show.legend = FALSE, alpha=0.9) +
  scale_x_continuous(name = "Monocle3 ranked cell position",
                     limits = c(0, 180000)) +
  # control space at top and bottom of plot
  scale_y_discrete(name = "", expand = c(0.02, 0, .08, 0)) +
  scale_fill_manual(values=colVec0) # colourblind-safe colours

pdf(file=paste0(rootDir,paste0("pseudotime_Neurons_originalOurs.pdf")), width=5, height = 6)
plot(tt2)
dev.off()

distributionsToTest <- list(c("hDA1b-2D-Day40","hDA1b-2D-Day70"),
                            c("hDA1b-3D-Day40","hDA1b-3D-Day120"),
                            c("hDA2-2D-Day40","hDA2-2D-Day70"))

library(ggsignif)
library(ggpubr)

fig4I <- tt2 + stat_compare_means(
  comparisons = distributionsToTest, aes(label = ..p.signif..),
  step.increase = 0, size=4, label.x = 50, label.y = 50, color="blue")

pdf(file=paste0(rootMain,"mainFigure4I.pdf")), width=5, height = 6)
plot(fig4I)
dev.off()




############################
## Supplemental Figure 8E ##
############################



cell_group_df<- tibble::tibble(cell=row.names(colData(mn.obj)),
                               cell_group=colData(mn.obj)$annotation_unified)
agg_mat <- aggregate_gene_expression(mn.obj, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
suppFig8E <- pheatmap::pheatmap(agg_mat,
                         scale="column", clustering_method="ward.D2")
pdf(file=paste0(rootSupp,"suppFigure8E.pdf"), width=12, height = 8)
suppFig8E
dev.off()


