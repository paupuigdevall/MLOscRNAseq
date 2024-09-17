library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(scuttle)
library(ggrepel)
library(Seurat)
library(ggplot2)
library(gghighlight)
library(ggbeeswarm)
library(ggpubr)
library(RColorBrewer)
library(knitr)
library(cowplot)

##############################
#### Early vs Late Foetal ####
##############################


otherFigs <- "figures/others/"
rootSupp <- "figures/supp/"

querySeurat <- readRDS("saved/toZenodo/mlo_resolution075_Annot.RDS")

clustAnnot <- c(0:23)
names(clustAnnot) <- c("hRgl2/immAstro","hNbDA","hProgFPM","OPC_1","VascLepto","hDA1b","hRgl1","hDA1a","hRgl3_caudal","hDA2","hProgM",
                       "hPreDA","hMidPre","hMgl","hEndo","hNbGaba","hNPro","hDA3/hGABA/hSer","Unk","hRgl4/MultiEpend","Astro","hPeric","Eryth","OPC_2")

querySeurat$seurat_clusters_24_Annot <- names(clustAnnot[querySeurat$seurat_clusters])
querySeurat$toPlotAnnot <- querySeurat$seurat_clusters_24_Annot
querySeurat$toPlotAnnot <- gsub("_","", querySeurat$toPlotAnnot)

querySeurat <- querySeurat[,querySeurat$originDimensions=="foetal"]

library(RColorBrewer)

setLast <- function(ctypesVec, lastCtype="Unk"){
  
  idLast <- match(lastCtype, ctypesVec)
  ordVec <- c(sort(ctypesVec[-idLast]), ctypesVec[idLast])
  return(ordVec)
}

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colourCount = length(unique(querySeurat$toPlotAnnot))
colVec <- setNames(getPalette(colourCount),
                   setLast(unique(querySeurat$toPlotAnnot), lastCtype="Unk"))


querySeurat.sce <- as.SingleCellExperiment(querySeurat)
milo.obj <- Milo(querySeurat.sce)

tmp <- as.data.frame(table(querySeurat.sce$toPlotAnnot))
colnames(tmp) <- c("seurat_clusters_24_Annot", "Number of cells")
kable(tmp)


## plot by cell type
UMAP_plot <- plotReducedDim(milo.obj, colour_by="toPlotAnnot",
                            dimred = "UMAP", point_size=0.1)+
  guides(colour = guide_legend(override.aes = list(size=2),
                               title="Cell types"))

pdf(file=paste0(otherFigs,"UMAP_milor.pdf"), width=8, height = 4)
UMAP_plot
dev.off()

## plot by foetal age
UMAP_timepoint_plot <- plotReducedDim(milo.obj, colour_by="dimensionsTime",
                                      dimred = "UMAP", point_size=0.3)+
  guides(colour = guide_legend(override.aes = list(size=3),
                               title="Foetal ages"))

pdf(file=paste0(otherFigs,"UMAP_tp_milor.pdf"), width=8, height = 4)
UMAP_timepoint_plot
dev.off()



milo.obj <- buildGraph(milo.obj, k = 30, d = 30, reduced.dim = "HARMONY")
milo.obj <- makeNhoods(milo.obj, prop = 0.2, k = 30, d=30, refined = TRUE, reduced_dims = "HARMONY")


histPlot <- plotNhoodSizeHist(milo.obj)
pdf(file=paste0(otherFigs,"NhoodSizeHist.pdf"), width=8, height = 4)
histPlot
dev.off()

## Test linear effect of time on cell type composition abundance
milo.obj$Sample <- paste0(milo.obj$midBrainId,"-",milo.obj$timePoint)
tmp <- as.data.frame(table(milo.obj$Sample))
colnames(tmp) <- c("Sample","Number of cells")
kable(tmp)

milo.obj <- countCells(milo.obj, meta.data = data.frame(colData(milo.obj)), sample="Sample")
head(nhoodCounts(milo.obj))


milo_design <- data.frame(colData(milo.obj))[,c("Sample", "timePoint")]
milo_design <- distinct(milo_design)
rownames(milo_design) <- milo_design$Sample
kable(head(milo_design))

milo.obj <- calcNhoodDistance(milo.obj, d=30, reduced.dim = "HARMONY")
da_results <- testNhoods(milo.obj, design = ~ timePoint, design.df = milo_design)

tmp <- da_results %>%
  arrange(SpatialFDR) %>%
  head()
rownames(tmp) <- NULL
kable(tmp)


distrUncorrPval <- ggplot(da_results, aes(PValue)) +
  geom_histogram(bins=50, col="white", fill="black")+
  theme_bw()

pdf(file=paste0(otherFigs,"distrUncorrPval.pdf"), width=8, height = 4)
distrUncorrPval
dev.off()

spatialfdrPlot <- ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point(size=0.5) +
  geom_hline(yintercept = 1)+
  theme_bw()

pdf(file=paste0(otherFigs,"spatialFDR.pdf"), width=8, height = 4)
spatialfdrPlot
dev.off()

milo.obj <- buildNhoodGraph(milo.obj)

umap_pl <- plotReducedDim(milo.obj, dimred = "UMAP", colour_by="timePoint", text_by = "seurat_clusters_24_Annot", 
                          text_size = 3, point_size=0.5) +
  guides(colour = guide_legend(override.aes = list(size=5),
                               title="Time"))
nh_graph_pl <- plotNhoodGraphDA(milo.obj, da_results, layout="UMAP",alpha=0.1)+
  theme(legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.75, "lines"))

pdf(file=paste0(otherFigs,"nhoodAndUMAPClusters.pdf"), width=8, height = 4)
umap_pl + nh_graph_pl +
  plot_layout(guides="keep")
dev.off()



da_results <- annotateNhoods(milo.obj, da_results, coldata_col = "seurat_clusters_24_Annot")
kable(head(da_results))


hist_fractions <- ggplot(da_results, aes(seurat_clusters_24_Annot_fraction)) +
  geom_histogram(bins=50, col="white", fill="black")+
  theme_bw()+coord_cartesian(xlim=c(0,1))

pdf(file=paste0(otherFigs,"histFractions_clust.pdf"), width=8, height = 4)
hist_fractions
dev.off()


da_results$celltype <- ifelse(da_results$seurat_clusters_24_Annot_fraction < 0.75, "Mixed", da_results$seurat_clusters_24_Annot)


mask_late <- da_results$FDR<0.1 & da_results$logFC>0
mask_early <- da_results$FDR<0.1 & da_results$logFC<0
da_results$dotGroup <- "Non signif."
da_results$dotGroup[mask_late] <- "Late Foetal"
da_results$dotGroup[mask_early] <- "Early Foetal"

da_results2 <- subset(da_results, celltype!="Mixed")

## ctypes with minor DA are not plotted
ctypesToRemove <- c("hPeric","hRgl4/MultiEpend","Eryth", "Unk", "hNPro","Astro")
da_results2 <- da_results2[!da_results2$celltype %in% ctypesToRemove,]

da_results2$celltype <- factor(da_results2$celltype,
                               levels=rev(c("OPC_1","hRgl3_caudal","hRgl2/immAstro",
                                        "hRgl1","hProgM","hProgFPM",
                                        "hNbGaba","hNbDA","hMgl",
                                        "hEndo","hDA3/hGABA/hSer","hDA2",
                                        "hDA1b","VascLepto")))



####################
## Supp Figure 3E ##
####################

DAbeeswarmCustomPlot <- ggplot(da_results2, aes(x=logFC, y=celltype, col=dotGroup))+
  theme_bw()+
  geom_quasirandom(alpha=0.7, groupOnX = F)+
  scale_color_manual(name="FDR<0.1",values=c("#663300","#000099","grey80"))+
  xlab("Log Fold-Change (linear PCW effect)")+
  ylab("")+
  theme(axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=11, angle=90, hjust=1, vjust=0.5),
        axis.title=element_text(size=15),
        legend.position="top")+
  scale_x_continuous(limits=c(-1.25,1), breaks=seq(-1.25,1,0.25))+
  guides(color = guide_legend(override.aes = list(size = 9)))

DAbeeswarmCustomPlotNoLegend <- DAbeeswarmCustomPlot + theme(legend.position="none")
DAbeeswarmCustomPlotOnlyLegend <- cowplot::get_plot_component(DAbeeswarmCustomPlot, "guide-box", return_all = TRUE)[[4]]
figS3E <- ggarrange(DAbeeswarmCustomPlotOnlyLegend, DAbeeswarmCustomPlotNoLegend, ncol = 1, heights = c(0.1,0.9))+bgcolor("white")

pdf(file=paste0(rootSupp,"suppFigure3E.pdf"), width=5.5, height = 6)
plot(figS3E)
dev.off()





