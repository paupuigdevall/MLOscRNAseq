library(Seurat)
library(ggplot2)
library(gghighlight)
library(ggbeeswarm)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(purrr)
library(ggplot2)
library(cowplot)


rootMain <- "figures/main/"
rootSupp <- "figures/supp/"

querySeurat <- readRDS("saved/toZenodo/mlo_resolution075_Annot.RDS")



clustAnnot <- c(0:23)
names(clustAnnot) <- c("hRgl2/immAstro","hNbDA","hProgFPM","OPC_1","VascLepto","hDA1b","hRgl1","hDA1a","hRgl3_caudal","hDA2","hProgM",
                       "hPreDA","hMidPre","hMgl","hEndo","hNbGaba","hNPro","hDA3/hGABA/hSer","Unk","hRgl4/MultiEpend","Astro","hPeric","Eryth","OPC_2")

### 

querySeurat$seurat_clusters_24_Annot <- names(clustAnnot[querySeurat$seurat_clusters])
querySeurat$toPlotAnnot <- querySeurat$seurat_clusters_24_Annot
querySeurat$toPlotAnnot <- gsub("_","", querySeurat$toPlotAnnot)
querySeurat <- querySeurat[,!(grepl("Patient", querySeurat$donorDimensions) | grepl("Crispr", querySeurat$donorDimensions))]


setLast <- function(ctypesVec, lastCtype="Unk"){
  
  idLast <- match(lastCtype, ctypesVec)
  ordVec <- c(sort(ctypesVec[-idLast]), ctypesVec[idLast])
  return(ordVec)
}

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colourCount = length(unique(querySeurat$toPlotAnnot))
colVec <- setNames(getPalette(colourCount),
                   setLast(unique(querySeurat$toPlotAnnot), lastCtype="Unk"))

saveRDS(colVec, file="saved/others/colVec_ctypesOriginal.RDS")


####################
## Main Figure 3C ##
####################

fig3C <- DimPlot(querySeurat, reduction = "umap", group.by = "toPlotAnnot", label = FALSE, label.size = 2.1,cols=alpha(colVec,0.6),
              repel = TRUE, )+ggtitle("")+
  theme_bw()+
  #legend.spacing.y = unit(0.05, 'cm'))+
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
        axis.title=element_blank())+
  theme(legend.position="none")

pdf(file=paste0(rootMain,"mainFigure3C.pdf"), width=4, height = 4)
fig3C
dev.off()

####################
## Supp Figure 3B ##
####################

figS3B <- DimPlot(querySeurat, reduction =  "umap", label = FALSE, group.by="toPlotAnnot", cols=colVec, pt.size=0.1)+theme_bw()+
  theme(plot.title=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        strip.text.x=element_text(size=7.5))+
  ggtitle("")+facet_wrap(~toPlotAnnot)+gghighlight()

pdf(file=paste0(rootSupp,"suppFigure3B.pdf"))
figS3B
dev.off()


####################
## Main Figure 3E ##
####################


querySeurat$originDimensions <- gsub("foetal","Foetal",querySeurat$originDimensions)

getPalette = colorRampPalette(brewer.pal(3, "Set1"))
colourCount = length(unique(querySeurat$originDimensions))
colVec <- setNames(getPalette(colourCount),
                   sort(unique(querySeurat$originDimensions)))

fig3E <- DimPlot(querySeurat, reduction =  "umap", label = FALSE, group.by="originDimensions", cols=colVec, pt.size=0.1)+theme_bw()+
  theme(plot.title=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        strip.text.x=element_text(size=10))+
  ggtitle("")+facet_wrap(~originDimensions)+gghighlight()

pdf(file=paste0(rootMain,"mainFigure3E.pdf"), width=8, height = 3)
fig3E
dev.off()




# Markers used for cell type annotation (available in Table S3)

markers <- c("SOX18", "SOX17", "ERG", "BCL6B", "EPAS1", "FOXF2",
             "FOXS1", "TBX2", "FOXD1", "FOXD2", "FOXC1", "JUNB",
             "HBG1", "HBG2", "HBA1", "HBA2", "HBB", "IGFBP7",
             "COL1A1", "COL1A2", "PDGFRA", "LUM", "FBLN2", "VCAN",
             "SPI1", "IRF8", "IKZF1", "REL", "MEF2C",
             "S100B", "GFAP", "AQP4", "ALDH1A1", "SOX9",
             "OLIG1", "OLIG2", "ETV1", "ETV5", "ZMAT3","SOX10",
             "MSX1", "GABPB2", "SAMD13", "CORIN", "SOX6", "PBX1", "FOXA2",
             "ALDH1L1", 
             "GLIS3", "DLK1", "TTR", "TFPI2", "SDC2",
             "ZMYND10", "CCNO", "RSPH1", "PIFO", "TUBB4B",
             "ENO1", "GTF2F1","GTF2F2", "TOX", "JADE1", "HMGB2",
             "HES6", "TFDP2", "NEUROG1", "OTX2",
             "EON1", "ASLC1", "CALB1", "DCC",
             "TH", "NETO2",
             "ARF4", "MBD4", "VGF",
             "DCX", "NEUROD6", "NEUROD1", "ZNF33A", "NHLH1", "EBF2", "LMO3", "DDC",
             "GATA3", "GAD1", "GAD2", "SLC32A1", "SLC17A6", "MEIS2",
             "KCNJ6",
             "APP",
             "SLC6A3", "PITX3",
             "TPH2",
             "COL2A1", "COL9A3", "COL9A1", "COL9A2", "COL11A1")



####################
## Supp Figure 3A ##
####################

querySeurat$toPlotAnnot <- factor(querySeurat$toPlotAnnot,
                                  levels=rev(c("hEndo","hPeric","Eryth","VascLepto",
                                           "hMgl","Astro","OPC1","OPC2",
                                           "hRgl1","hRgl2/immAstro","hRgl3caudal","hRgl4/MultiEpend",
                                           "hProgFPM","hProgM","hNPro","hMidPre","hPreDA",
                                           "hNbDA","hNbGaba","hDA1a","hDA1b","hDA2","hDA3/hGABA/hSer","Unk")))

Idents(querySeurat) <- "toPlotAnnot"

figS3A <- DotPlot(querySeurat, features = markers) + 
  RotatedAxis()+
  theme(axis.text.x=element_text(size=7),
        legend.position="top")+
  xlab("")+ylab("")


pdf(file=paste0(rootSupp,"suppFigure3A.pdf"), width=20, height = 6)
figS3A
dev.off()


# Gene modules with gene markers
modulesToCompute <- list("hEndo"=c("SOX18", "SOX17", "ERG", "BCL6B", "EPAS1", "FOXF2"),
                         "hPeric"=c("FOXS1", "TBX2", "FOXD1", "FOXD2", "FOXC1", "JUNB"),
                         "Eryth"=c("HBG1", "HBG2", "HBA1", "HBA2", "HBB", "IGFBP7"),
                         "VascLepto"=c("COL1A1", "COL1A2", "PDGFRA", "LUM", "FBLN2", "VCAN"),
                         "hMgl"=c("SPI1", "IRF8", "IKZF1", "REL", "MEF2C"),
                         "Astro"=c("S100B", "GFAP", "AQP4", "ALDH1A1", "SOX9"),
                         "OPC1"=c("OLIG1", "OLIG2", "ETV1", "ETV5", "ZMAT3"),
                         "OPC2"=c("OLIG1", "OLIG2", "ETV1", "ETV5", "ZMAT3", "SOX10"),
                         "hRgl1"=c("MSX1", "GABPB2", "SAMD13", "CORIN", "SOX6", "PBX1", "FOXA2"),
                         "hRgl2/immAstro"=c("MSX1", "GABPB2", "SAMD13", "CORIN", "SOX6", "PBX1", "AQP4", "ALDH1L1", "S100B"),
                         "hRgl3caudal"=c("GLIS3","DLK1","TTR","TFPI2", "SDC2"),
                         "hRgl4/MultiEpend"=c("ZMYND10", "CCNO", "RSPH1", "PIFO", "TUBB4B"),
                         "hProgFPM"=c("FOXA2", "ENO1", "GTF2F1", "GTF2F2", "TOX", "JADE1", "HMGB2"),
                         "hProgM"=c("ENO1", "TOX", "JADE1", "HMGB2"),
                         "hNPro"=c("HES6", "TFDP2", "NEUROG1", "OTX2"),
                         "hMidPre"=c("DLK1", "ASCL1", "CALB1", "DCC"),
                         "hPreDA"=c("ARF4", "MBD4", "VGF"),
                         "hNbDA"=c("DCX", "NEUROD6", "NEUROD1", "ZNF33A", "NHLH1", "EBF2", "TH", "LMO3", "CALB1", "DDC", "PBX1"),
                         "hNbGaba"=c("DCX", "NEUROD6", "NEUROD1", "ZNF33A", "NHLH1", "EBF2", "GATA3", "GAD1", "GAD2", "SLC32A1", "SLC17A6", "MEIS2"),
                         "hDA1b"=c("TH", "LMO3", "KCNJ6", "DCX"),
                         "hDA2"=c("TH", "KCNJ6", "APP", "DDC", "GAD1", "GAD2"),
                         "hDA3/hGABA/hSer"=c("ALDH1A1", "LMO3", "SLC6A3", "PITX3", "GAD1", "GAD2", "TPH2"),
                         "Unk"=c("COL2A1", "COL9A3", "COL9A1", "COL9A2", "COL11A1"))


moduleValues <- sapply(modulesToCompute, function(x){
  
  moduleScore <- AddModuleScore(querySeurat, list(x))
  return(moduleScore$Cluster1)
  
}, simplify=F)


moduleValues <- do.call("cbind", moduleValues)
querySeurat@meta.data <- cbind(querySeurat@meta.data, moduleValues[match(rownames(querySeurat@meta.data), rownames(moduleValues)),])


moduleDim <- cbind(querySeurat@reductions$umap[[]], moduleValues)
moduleDim_long <- as.data.frame(as.data.frame(moduleDim) %>% pivot_longer(-c("UMAP_1","UMAP_2"), names_to="CellType", values_to="ActivationScore"))

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colourCount = length(unique(moduleDim_long$CellType))
colVec <- setNames(getPalette(colourCount),
                   setLast(unique(moduleDim_long$CellType), lastCtype="Unk"))


moduleDim_long$CellType <- factor(moduleDim_long$CellType, levels=names(colVec))

# facetModules2 <- moduleDim_long %>% group_split(CellType) %>%
#   map(
#     ~ggplot(., aes(x=UMAP_1, y=UMAP_2, col=ActivationScore))+
#         geom_point(size=0.05)+
#         theme_bw()+
#         scale_colour_gradient(name="", low = "grey80", high = "red")+
#         #scale_color_gradientn(name="",colours = terrain.colors(10))+
#         theme(plot.title=element_blank(),
#               panel.border = element_blank(), panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               axis.text.x=element_blank(),
#               axis.ticks.x=element_blank(),
#               axis.text.y=element_blank(),
#               axis.ticks.y=element_blank(),
#               axis.title=element_blank(),
#               strip.text.x=element_text(size=6.5),
#               legend.title=element_text(size=7),
#               legend.text=element_text(size=7),
#               legend.key.size=unit(0.5,"lines"))+
#         ggtitle("")+facet_wrap(~CellType)) %>% 
#   plot_grid(plotlist = ., align = 'hv', ncol = 4)
#   
# 
# 
# pdf(file=paste0(rootDir,"facetModules_perCellType2.pdf"))
# plot(facetModules2)
# dev.off()


####################
## Main Figure 3B ##
####################

moduleDim2 <- cbind(querySeurat@meta.data[,c("toPlotAnnot","seurat_clusters_24_Annot")], moduleValues)
moduleDim2$seurat_clusters_24_Annot <- NULL
moduleDim2_long <- as.data.frame(as.data.frame(moduleDim2) %>% pivot_longer(-c("toPlotAnnot"), names_to="ModulesCellType", values_to="ActivationScore"))

stopifnot(all(unique(moduleDim2_long$ModulesCellType) %in% levels(moduleDim2_long$toPlotAnnot)))

moduleDim2_long$ModulesCellType <- factor(moduleDim2_long$ModulesCellType, levels=setLast(unique(moduleDim2_long$ModulesCellType), lastCtype="Unk"))

allctypes <- sapply(levels(moduleDim2_long$ModulesCellType), function(x){
  
  tmp <- subset(moduleDim2_long, ModulesCellType==x)
  tmp$toPlotAnnot <- as.character(tmp$toPlotAnnot)
  tmp[!tmp$toPlotAnnot %in% x,]$toPlotAnnot <- "OC"
  return(tmp)
  
}, simplify=F)

allctypes <- do.call("rbind", allctypes)
rownames(allctypes) <- NULL

allctypes$toPlotAnnot <- factor(allctypes$toPlotAnnot , levels=setLast(unique(allctypes$toPlotAnnot), lastCtype="OC"))

fig3B <- ggplot(data=allctypes, aes(x=toPlotAnnot, y=ActivationScore,))+
  theme_bw()+
  facet_wrap(~ModulesCellType, scales="free")+
  geom_violin(fill="grey70")+
  geom_boxplot(width=0.05, outlier.shape=NA)+
  xlab("Module activation score in Table S2")

pdf(file=paste0(rootMain,"mainFigure3B.pdf"), width=10, height = 6)
plot(fig3B)
dev.off()




#################################
## Cell type composition plots ##
#################################

querySeurat$dimensionsTime2 <- gsub("foetal-","foetal-PCW", querySeurat$dimensionsTime)


computeCtypeProp <- function(querySeurat, ctype="toPlotAnnot"){
  
  allres <- sapply(unique(querySeurat$dimensionsTime2), function(x){
    
    tmp <- subset(querySeurat@meta.data, dimensionsTime2==x)
    col_to_Pull <- match(ctype, colnames(tmp))
    tmp2 <- as.data.frame(table(tmp[,col_to_Pull]))
    colnames(tmp2) <- c(ctype,"nCells")
    tmp2$Percentage <- tmp2$nCells*100/sum(tmp2$nCells)
    tmp2$dimensionsTime2 <- x
    return(tmp2)

  }, simplify=F)
  
  allres <- do.call("rbind", allres)
  rownames(allres) <- NULL
  return(allres)
  
}

ctypeProp_annot <- computeCtypeProp(querySeurat, ctype="toPlotAnnot")
#ctypeProp_annot2 <- computeCtypeProp(querySeurat, ctype="toPlotAnnot2")


sampleCellsNum <- as.data.frame(table(querySeurat$dimensionsTime2))
colnames(sampleCellsNum) <- c("dimensionsTime2","nCells")
sampleCellsNum$dimensionsTime2 <- factor(sampleCellsNum$dimensionsTime2,
                                           levels=c("2D-40","2D-70","3D-40","3D-70","3D-120",
                                                    "foetal-PCW10","foetal-PCW12","foetal-PCW16","foetal-PCW20"))

## number of cells for each 10x sample


####################
## Main Figure 3D ##
####################

fig3D_b <- ggplot(sampleCellsNum, aes(x=dimensionsTime2, y=nCells))+
  geom_bar(stat="identity", fill="black")+
  geom_text(aes(label=nCells), vjust=-0.3, size=3)+
  theme_bw()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
        axis.title=element_text(size=9))+
  ylab("")+xlab("")+
scale_y_continuous(labels = scales::comma, breaks=seq(0,10000,2500), limits=c(0,11000))

pdf(file=paste0(rootMain,"mainFigure3D_bot.pdf"), width=4, height=2)
plot(fig3D_b)
dev.off()


## Specific-level (barplot + cells)


getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colourCount = length(sort(levels(ctypeProp_annot$toPlotAnnot)))
colVec <- setNames(getPalette(colourCount),
                   setLast(sort(levels(ctypeProp_annot$toPlotAnnot)), lastCtype=c("Unk")))

ctypeProp_annot$toPlotAnnot <- factor(ctypeProp_annot$toPlotAnnot, levels=names(colVec))
ctypeProp_annot$dimensionsTime2 <- factor(ctypeProp_annot$dimensionsTime2,
                                           levels=c("2D-40","2D-70","3D-40","3D-70","3D-120",
                                                    "foetal-PCW10","foetal-PCW12","foetal-PCW16","foetal-PCW20"))

fig3D_u <- ggplot(ctypeProp_annot, aes(fill=toPlotAnnot, x=dimensionsTime2, y=Percentage))+
  geom_bar(position="fill", stat="identity")+
  theme_bw()+
  theme(plot.title=element_text(size=18, face="bold", hjust=0.5),
        legend.title = element_text(size = 9, face="bold", hjust=0.5),
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "lines"),
        axis.text.x = element_text(size=14, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=14),
        axis.title=element_text(size=16),
        strip.text.x = element_text(size = 10)) +
  guides(shape = guide_legend(override.aes = list(size = 5)),
         fill = guide_legend(override.aes = list(size = 5), ncol=1))+
  scale_fill_manual(name="Cell types",
                    values = colVec)+
  scale_y_continuous(labels = scales::percent, breaks=seq(0,1,0.2))+
  xlab("")+
  ylab("Percentage")


pdf(file=paste0(rootMain,"mainFigure3D_upp.pdf"), width=6, height=8)
plot(fig3D_u)
dev.off()

