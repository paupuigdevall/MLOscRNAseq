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


pdf(file=paste0(rootMain,"mainFigure4E.pdf"), width=14, height = 4)
plot(plot2)
dev.off()


