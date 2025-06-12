
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


allctypes <- names(table(obj$annotation_mixed))
namesUnified_allctypes_level1 <- c("Astrocytes", "Astrocytes","Astrocytes", "Astrocytes",
                                   "Neurons","Neurons","Neurons","Endothelial/Pericytes",
                                   "Erythrocytes","Erythrocytes","Fibroblasts","FPP",
                                   "FPP","FPP","FPP","Neurons",
                                   "Glioblasts","Neurons","Neurons","Neurons",
                                   "Neurons","Neurons","Neurons","Neurons",
                                   "Endothelial/Pericytes","Neurons","Microglia","Precursors",
                                   "Immature Neurons","Immature Neurons","Immature Neurons","Immature Neurons",
                                   "Immature Neurons","Progenitors","Progenitors","Neurons",
                                   "OPC","Endothelial/Pericytes","Precursors","Progenitors",
                                   "Progenitors","Progenitors","Progenitors","Progenitors",
                                   "Progenitors","Progenitors","Progenitors","Progenitors",
                                   "Progenitors","Progenitors","Progenitors","Neurons",
                                   "Neurons","Immune","Microglia","Unknown",
                                   "Immature Neurons","Neurons","IPC","ODC",
                                   "ODC","ODC","ODC","OPC",
                                   "OPC","OPC","Neurons","Endothelial/Pericytes",
                                   "Progenitors","Progenitors","Erythrocytes","Progenitors",
                                   "Progenitors","Progenitors","Unknown","VLMC",
                                   "VLMC","VLMC")

namesUnified_allctypes_level2 <- c("Others", "Others","Others", "Others",
                                   "DopaminergicN","DopaminergicN","DopaminergicN","Others",
                                   "Others","Others","Others","Others",
                                   "Others","Others","Others","Others",
                                   "Others","Others","DopaminergicN","DopaminergicN",
                                   "DopaminergicN","DopaminergicN","DopaminergicN","DopaminergicN",
                                   "Others","Others","Others","Others",
                                   "Others","Others","Others","Others",
                                   "Others","Others","Others","Others",
                                   "Others","Others","Others","Others",
                                   "Others","Others","Others","Others",
                                   "Others","Others","Others","Others",
                                   "Others","Others","Others","Others",
                                   "Others","Others","Others","Others",
                                   "Others","Others","Others","Others",
                                   "Others","Others","Others","Others",
                                   "Others","Others","Others","Others",
                                   "Others","Others","Others","Others",
                                   "Others","Others","Others","Others",
                                   "Others","Others")

allctypesUnified <- setNames(namesUnified_allctypes_level1, allctypes)
allctypesUnified2 <- setNames(namesUnified_allctypes_level2, allctypes)


obj$annotation_unified <- "NA"
obj$annotation_unified <- unname(allctypesUnified[obj$annotation_mixed])
obj$annotation_unified[is.na(obj$annotation_unified)] <- "Unknown"


obj$annotation_unified_DAN <- "NA"
obj$annotation_unified_DAN <- unname(allctypesUnified2[obj$annotation_mixed])
obj$annotation_unified_DAN[is.na(obj$annotation_unified_DAN)] <- "Others"



## Export Table S6
suppTab6 <- obj@meta.data[,c("annotation_mixed","annotation_unified","Dataset2")]
rownames(suppTab6) <- NULL
suppTab6 <- distinct(suppTab6)
colnames(suppTab6)[3] <- "Dataset"
suppTab6 <- suppTab6[suppTab6$Dataset!="Zagare et al. 2022 (3D)",]
write.table(suppTab6, file="saved/suppTables/TableS6.txt",
            sep="\t", col.names=TRUE, row.names=F, quote=F)




####################
## Main Figure 4A ##
####################


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

saveRDS(colVec, file="saved/colVec_datasetColours.RDS")

sapply(integrationList, function(y){
  
  fig4A <- DimPlot(obj, reduction = y, label = FALSE, raster=FALSE, group.by="Dataset2", cols=colVec)+theme_bw()+
    theme(plot.title=element_blank(),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title=element_blank(),
          strip.text.x=element_text(size=15))+
    ggtitle("")+facet_wrap(~Dataset2, nrow=2)+gghighlight()
  s
  
  fig4A_2 <- cowplot::get_legend(plotToPdf)
  fig4A_1 <- plotToPdf + theme(legend.position = "none")
  
  
  pdf(file=paste0(rootMain,"mainFigure4A.pdf"), width=14, height = 8)
  plot(fig4A_1)
  dev.off()
  
  pdf(file=paste0(rootMain,"mainFigure4A_onlyLegend.pdf"))
  plot(fig4A_2)
  dev.off()
  
})



####################
## Main Figure 4C ##
####################


setLast <- function(ctypesVec, lastCtype="Unk"){
  
  idLast <- match(lastCtype, ctypesVec)
  ordVec <- c(sort(ctypesVec[-idLast]), ctypesVec[idLast])
  return(ordVec)
}

obj$toPlotAnnot <- obj$annotation_unified
getPalette = colorRampPalette(brewer.pal(8, "Accent"))

colourCount = length(unique(obj$toPlotAnnot))
colVec2 <- setNames(getPalette(colourCount),
                    setLast(unique(obj$toPlotAnnot), lastCtype="Unknown"))

obj$toPlotAnnot <- factor(obj$toPlotAnnot, levels=names(colVec2))
saveRDS(colVec2, file="saved/colVec_ctypesColours.RDS")


sapply(integrationList, function(y){
  
  p2 <- DimPlot(obj, reduction = y, group.by = "toPlotAnnot", label = FALSE, cols=alpha(colVec2,0.9),
                repel = TRUE)+ggtitle("")+
    theme_bw()+
    guides(col=guide_legend(ncol=1, override.aes = list(size=4)))+
    theme(legend.title = element_text(size = 9, face="bold", hjust=0.5),
          legend.text  = element_text(size = 9),
          legend.key.size = unit(0.12, "cm"),
          plot.title=element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank())
  
  pdf(file=paste0(rootMain,"mainFigure4C.pdf"), width=10, height = 8)
  plot(p2)
  dev.off()

})


####################
## Main Figure 4D ##
####################


obj$toPlotAnnot2 <- obj$annotation_unified_DAN
obj$toPlotAnnot2 <- factor(obj$toPlotAnnot2, levels=sort(unique(obj$toPlotAnnot2)))

sapply(integrationList, function(y){
  
  p2 <- DimPlot(obj, reduction = y, group.by = "toPlotAnnot2", label = FALSE, cols=c("red","grey"),
                repel = TRUE)+ggtitle("")+
    theme_bw()+
    guides(col=guide_legend(ncol=1, override.aes = list(size=4)))+
    theme(legend.title = element_text(size = 9, face="bold", hjust=0.5),
          legend.text  = element_text(size = 9),
          legend.key.size = unit(0.12, "cm"),
          plot.title=element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank())
  

  pdf(file=paste0(rootMain,"mainFigure4D.pdf"), width=10, height = 8)
  plot(p2)
  dev.off()
})



####################
## Main Figure 4B ##
####################

## Plot consenus annotation for each dataset

obj$plotUnit <- paste0(obj$Dataset,"-",obj$system,"-",obj$timePoint)
obj$plotUnit <- gsub("Agarwal-PostMortem-PostMortem","Agarwal-PostMortem", obj$plotUnit)
obj$plotUnit <- gsub("Birtele-ParmarFoetal-Foetal","Birtele-ParmarFoetal", obj$plotUnit)
obj$plotUnit <- gsub("Fiorenzano-ParmarOrg-Organoid","Fiorenzano-ParmarOrganoid", obj$plotUnit)

resultMixed <- obj@meta.data %>%
  group_by(plotUnit, annotation_unified) %>%
  summarise(cell_type_fraction = n()) %>%
  mutate(cell_type_fraction = cell_type_fraction*100 / sum(cell_type_fraction)) %>%
  as.data.frame()

resultMixed$Dataset <- sapply(strsplit(resultMixed$plotUnit,"-"), function(x) x[1])
resultMixed <- subset(resultMixed, Dataset!="Schwamborn")

maskMLO <- grepl("MLO", resultMixed$plotUnit)
resultMixed$plotUnit[!maskMLO] <- sapply(strsplit(resultMixed$plotUnit[!maskMLO], "-"), function(x) x[length(x)])
resultMixed$plotUnit[maskMLO] <- sapply(strsplit(resultMixed$plotUnit[maskMLO], "-"), function(x) paste(x[-1], collapse="-"))

resultMixed$plotUnit <- gsub("Organoid-","3D-", resultMixed$plotUnit)
resultMixed[resultMixed$Dataset=="Fiorenzano",]$plotUnit <- paste0("3D-", resultMixed[resultMixed$Dataset=="Fiorenzano",]$plotUnit)

resultMixed$plotUnit <- factor(resultMixed$plotUnit, levels=c("PostMortem",
                                                              "PCW6","PCW8","PCW11","PCW14",
                                                              "3D-day15","3D-day30","3D-day60","3D-day90",
                                                              "GW6","GW7","GW8","GW9","GW10","GW11",
                                                              "2D-day40","2D-day70",
                                                              "3D-day40","3D-day70","3D-day120",
                                                              "Foetal-PCW10","Foetal-PCW12","Foetal-PCW16","Foetal-PCW20"))

resultMixed$Dataset2 <- resultMixed$Dataset
resultMixed$Dataset2[grepl("Agarwal", resultMixed$Dataset2 )] <- "Agarwal et al. 2020 (PostMortem)"
resultMixed$Dataset2[grepl("Birtele", resultMixed$Dataset2 )] <- "Birtele et al. 2022 (Foetal)"
resultMixed$Dataset2[grepl("Fiorenzano", resultMixed$Dataset2 )] <- "Fiorenzano et al. 2021 (3D)"
resultMixed$Dataset2[grepl("LeManno", resultMixed$Dataset2 )] <- "La Manno et al. 2016 (Foetal)"
resultMixed$Dataset2[grepl("MLO", resultMixed$Dataset2 )] <- "This work"
resultMixed$Dataset2[grepl("Braun", resultMixed$Dataset2 )] <- "Braun et al. 2023 (Foetal)"

resultMixed$Dataset2 <- factor(resultMixed$Dataset2, levels=c("This work",
                                                              "La Manno et al. 2016 (Foetal)",
                                                              "Birtele et al. 2022 (Foetal)",
                                                              "Braun et al. 2023 (Foetal)",
                                                              "Fiorenzano et al. 2021 (3D)",
                                                              "Agarwal et al. 2020 (PostMortem)"))

resultMixed$annotation_unified <- factor(resultMixed$annotation_unified, levels=names(colVec2))

barplotStckd <- ggplot(resultMixed, aes(fill=annotation_unified, x=plotUnit, y=cell_type_fraction))+
  geom_bar(position="fill", stat="identity")+
  ggtitle("")+
  theme_bw()+
  facet_grid(.~Dataset2, scales="free_x", space="free")+
  theme(plot.title=element_text(size=18, face="bold", hjust=0.5),
        legend.title = element_text(size = 11, face="bold", hjust=0.5),
        legend.text  = element_text(size = 9),
        legend.key.size = unit(0.5, "lines"),
        legend.position="none",
        axis.text.x = element_text(size=9, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=9),
        axis.title=element_text(size=11),
        strip.text.x = element_text(size = 10)) +
  guides(shape = guide_legend(override.aes = list(size = 5)),
         fill = guide_legend(override.aes = list(size = 5), ncol=1))+
  scale_fill_manual(name="Cell types",
                    values = colVec2)+
  scale_y_continuous(labels = scales::percent, breaks=seq(0,1,0.2))+
  xlab("")+
  ylab("Cell type abundance [%]")


pdf(file=paste0(rootMain,"mainFigure4B.pdf"), width=14, height = 4)
barplotStckd
dev.off()

rm(list=ls())
gc()











