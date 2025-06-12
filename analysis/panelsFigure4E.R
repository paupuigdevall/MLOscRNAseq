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
library(ggpubr)
library(ggh4x)


rootMain <- "figures/main/"
rootSupp <- "figures/supp/"

obj <- readRDS("saved/toZenodo/midBrainIntegration.RDS")
obj <- obj[,(obj$Dataset=="MLO" &  grepl("Control", obj$donorDimensions)) | obj$system=="Foetal"]
gc()

obj$Dataset2 <- obj$Dataset
obj$Dataset2[grepl("Agarwal", obj$Dataset2 )] <- "Agarwal et al. 2020 (PostMortem)"
obj$Dataset2[grepl("Birtele-ParmarFoetal", obj$Dataset2 )] <- "Birtele et al. 2022 (Foetal)"
obj$Dataset2[grepl("Fiorenzano-ParmarOrg", obj$Dataset2 )] <- "Fiorenzano et al. 2021 (3D)"
obj$Dataset2[grepl("LeManno", obj$Dataset2 )] <- "La Manno et al. 2016 (Foetal)"
obj$Dataset2[grepl("MLO", obj$Dataset2 )] <- "This work (2D, 3D, Foetal)"
obj$Dataset2[grepl("Schwamborn", obj$Dataset2 )] <- "Zagare et al. 2022 (3D)"
obj$Dataset2[grepl("Braun", obj$Dataset2 )] <- "Braun et al. 2023 (Foetal)"
obj$Dataset2 <- factor(obj$Dataset2, levels=c("This work (2D, 3D, Foetal)",
                                              "La Manno et al. 2016 (Foetal)",
                                              "Birtele et al. 2022 (Foetal)",
                                              "Braun et al. 2023 (Foetal)",
                                              "Fiorenzano et al. 2021 (3D)",
                                              "Agarwal et al. 2020 (PostMortem)",
                                              "Zagare et al. 2022 (3D)"))

integrationList  <- "umap.mnn"

library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(6, "Set1"))

colourCount = length(levels(obj$Dataset2))
colVec <- setNames(getPalette(colourCount),
                   levels(obj$Dataset2))

obj$facetUnit <- NA
maskMLO <- obj$Dataset=="MLO"
obj$facetUnit[maskMLO] <- gsub("foetal-","PCW-", obj$dimensionsTime[maskMLO])
obj$facetUnit[!maskMLO] <- gsub("GW","GW-", obj$timePoint[!maskMLO])
obj$facetUnit[!maskMLO] <- gsub("PCW","PCW-", obj$facetUnit[!maskMLO])

obj$facetUnit <- factor(obj$facetUnit,
                        levels=c("2D-40","2D-70","3D-40","3D-70","3D-120",
                                 "PCW-6","PCW-8","GW-6","GW-7","GW-8","PCW-10",
                                 "GW-9","PCW-11","GW-10","PCW-12","GW-11",
                                 "PCW-14","PCW-16","PCW-20"))

####################
## Main Figure 4E ##
####################

plot2 <- DimPlot(obj, reduction = "umap.mnn", label = FALSE, raster=NULL, pt.size = 0.1, group.by="Dataset2", split.by="facetUnit",cols=colVec, ncol=10)+theme_bw()+
  theme(plot.title=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        strip.text.x=element_text(size=17),
        legend.position="right",
        legend.text=element_text(size=15))+
  ggtitle("")+gghighlight(keep_scales = T)


pdf(file=paste0(rootMain,"mainFigure4E.pdf"), width=34, height = 8)
plot(plot2)
dev.off()


