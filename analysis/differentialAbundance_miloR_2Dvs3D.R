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


##################
#### 2D vs 3D ####
##################

otherFigs <- "figures/others/"
rootMain <- "figures/main/"

querySeurat <- readRDS("saved/toZenodo/mlo_resolution075_Annot.RDS")
clustAnnot <- c(0:23)
names(clustAnnot) <- c("hRgl2/immAstro","hNbDA","hProgFPM","OPC_1","VascLepto","hDA1b","hRgl1","hDA1a","hRgl3_caudal","hDA2","hProgM",
                       "hPreDA","hMidPre","hMgl","hEndo","hNbGaba","hNPro","hDA3/hGABA/hSer","Unk","hRgl4/MultiEpend","Astro","hPeric","Eryth","OPC_2")

### 

querySeurat$seurat_clusters_24_Annot <- names(clustAnnot[querySeurat$seurat_clusters])
querySeurat$toPlotAnnot <- querySeurat$seurat_clusters_24_Annot
querySeurat$toPlotAnnot <- gsub("_","", querySeurat$toPlotAnnot)


querySeurat <- querySeurat[,!(grepl("Patient", querySeurat$donorDimensions) | grepl("Crispr", querySeurat$donorDimensions))]



querySeurat.sce <- as.SingleCellExperiment(querySeurat)
querySeurat.sce <- querySeurat.sce[,colData(querySeurat.sce)$origin=="iPSC"]
milo.obj <- Milo(querySeurat.sce)
dim(milo.obj)


tmp <- as.data.frame(table(querySeurat.sce$toPlotAnnot))
colnames(tmp) <- c("Cluster", "Number of cells")
kable(tmp)



## Construct KNN graph

embryo_milo <- buildGraph(milo.obj, k = 30, d = 30, reduced.dim = "HARMONY")
embryo_milo <- makeNhoods(embryo_milo, prop = 0.3, k = 30, d=30, refined = TRUE, reduced_dims = "HARMONY")


nhoodPlot <- plotNhoodSizeHist(embryo_milo)
pdf(file=paste0(otherFigs, "nhood_Plot.pdf"))
plot(nhoodPlot)
dev.off()


embryo_milo$Sample <- paste0(embryo_milo$donorDimensions,"-",embryo_milo$timePoint)
tmp <- as.data.frame(table(embryo_milo$Sample))
colnames(tmp) <- c("Sample","Number of cells")
kable(tmp)

embryo_milo <- countCells(embryo_milo, meta.data = data.frame(colData(embryo_milo)), sample="Sample")
head(nhoodCounts(embryo_milo))


##Defining experimental design
embryo_milo$Time <- gsub(".+-","",embryo_milo$dimensionsTime)
embryo_design <- data.frame(colData(embryo_milo))[,c("Sample", "Time","originDimensions")]
head(embryo_design)

embryo_design <- distinct(embryo_design)
rownames(embryo_design) <- embryo_design$Sample

embryo_design <- embryo_design[colnames(nhoodCounts(embryo_milo)), , drop=FALSE]
table(embryo_design$Time)

## Computing neighbourhood connectivity
embryo_milo <- calcNhoodDistance(embryo_milo, d=30, reduced.dim = "HARMONY")

#da_results <- testNhoods(embryo_milo, design = ~ originDimensions, design.df = embryo_design)
da_results_cov <- testNhoods(embryo_milo, design = ~ Time + originDimensions, design.df = embryo_design)


## with covariates

tmp <- da_results_cov %>%
  arrange(SpatialFDR) %>%
  head()
rownames(tmp) <- NULL
kable(tmp)

pvaluehist1 <- ggplot(da_results_cov, aes(PValue)) +
  geom_histogram(bins=50, col="white", fill="black")+
  theme_bw()

volcanoplot1 <- ggplot(da_results_cov, aes(logFC, -log10(SpatialFDR))) + 
  geom_point(size=0.5) +
  geom_hline(yintercept = 1)+ ## Mark significance threshold (10% FDR)
  theme_bw()

pdf(file=paste0(otherFigs, "histPval_cov.pdf"))
plot(pvaluehist1)
dev.off()

pdf(file=paste0(otherFigs, "volcanoPlot_cov.pdf"))
plot(volcanoplot1)
dev.off()


embryo_milo <- buildNhoodGraph(embryo_milo)

umap_pl <- plotReducedDim(embryo_milo, dimred = "UMAP", colour_by="originDimensions", text_by = "toPlotAnnot", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")
nh_graph_pl <- plotNhoodGraphDA(embryo_milo, da_results_cov, layout="UMAP",alpha=0.1) 
library(patchwork)
combined_pl <- umap_pl + nh_graph_pl +
  plot_layout(guides="collect")

pdf(file=paste0(otherFigs, "umapAndNhoodGraph.pdf"))
plot(combined_pl)
dev.off()


da_results_cov <- annotateNhoods(embryo_milo, da_results_cov, coldata_col = "toPlotAnnot")
#kable(head(da_results_cov))

histFrac <- ggplot(da_results_cov, aes(toPlotAnnot_fraction)) +
  geom_histogram(bins=50, col="white", fill="black")+
  theme_bw()+coord_cartesian(xlim=c(0,1))

pdf(file=paste0(otherFigs, "hist_annotateNhoods.pdf"))
plot(histFrac)
dev.off()

da_results_cov$celltype <- ifelse(da_results_cov$toPlotAnnot_fraction < 0.7, "Mixed", da_results_cov$toPlotAnnot)



####################
## Main Figure 3F ##
####################

mask_3d <- da_results_cov$FDR<0.1 & da_results_cov$logFC>0
mask_2d <- da_results_cov$FDR<0.1 & da_results_cov$logFC<0
da_results_cov$dotGroup <- "Non signif."
da_results_cov$dotGroup[mask_2d] <- "2D"
da_results_cov$dotGroup[mask_3d] <- "3D"
da_results_cov2 <- subset(da_results_cov, celltype!="Mixed")

ctypesToRemove <- names(which(rowSums(table(da_results_cov2$celltype, da_results_cov2$SpatialFDR<0.1))<20))

if (length(ctypesToRemove)){
  da_results_cov2 <- da_results_cov2[!da_results_cov2$celltype %in% ctypesToRemove,]
}

da_results_cov2$celltype <- factor(da_results_cov2$celltype,
                               levels=rev(sort(unique(da_results_cov2$celltype))))

fig3F <- ggplot(da_results_cov2, aes(x=logFC, y=celltype, col=dotGroup))+
  theme_bw()+
  geom_quasirandom(alpha=0.7, groupOnX = F)+
  scale_color_manual(name="FDR<0.1",values=c("#E41A1C","#377EB8","grey80"))+
  xlab("Log Fold-Change (~ Day + Model)")+
  ylab("")+
  theme(axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=11, angle=90, hjust=1, vjust=0.5),
        axis.title=element_text(size=15),
        legend.position="top")+
  scale_x_continuous(limits=c(-8,8), breaks=seq(-8,8,2))+
  guides(color = guide_legend(override.aes = list(size = 9)))

pdf(file=paste0(rootMain,"mainFigure3F.pdf"), width=5, height = 6)
plot(fig3F)
dev.off()



## Transcriptomic differences within dopaminergic populations and between in vitro models

## Run buildNhoodGraph to store nhood adjacency matrix
embryo_milo <- buildNhoodGraph(embryo_milo)

## Find groups
da_results_cov <- groupNhoods(embryo_milo, da_results_cov, max.lfc.delta = 2)
head(da_results_cov)
table(da_results_cov$NhoodGroup, da_results_cov$toPlotAnnot)/rowSums(table(da_results_cov$NhoodGroup, da_results_cov$toPlotAnnot))*100

nhood_afterDA <- plotNhoodGroups(embryo_milo, da_results_cov, layout="UMAP")
DAbeeswarm_afterDA <- plotDAbeeswarm(da_results_cov, "NhoodGroup")

pdf(file=paste0(otherFigs,"nhood_afterDA.pdf"), width=5.5, height = 6)
plot(nhood_afterDA)
dev.off()

pdf(file=paste0(otherFigs,"DAbeeswarm_afterDA.pdf"), width=5.5, height = 6)
plot(DAbeeswarm_afterDA)
dev.off()


## define Nhoodgroup to be a cell type

vecNames <- setNames(1:length(unique(da_results_cov$celltype)),
                     sort(unique(da_results_cov$celltype)))

da_results_cov$NhoodGroup <- unname(vecNames[da_results_cov$celltype])


## hDA1a
dge_hDA1a <- testDiffExp(embryo_milo, da_results_cov, design = ~ Time + originDimensions, meta.data = data.frame(colData(embryo_milo)),
                     subset.nhoods=da_results_cov$NhoodGroup==which(names(vecNames)=="hDA1a"))

dge_hDA1b <- testDiffExp(embryo_milo, da_results_cov, design = ~ Time + originDimensions, meta.data = data.frame(colData(embryo_milo)),
                     subset.nhoods=da_results_cov$NhoodGroup==which(names(vecNames)=="hDA1b"))

dge_hDA2 <- testDiffExp(embryo_milo, da_results_cov, design = ~ Time + originDimensions, meta.data = data.frame(colData(embryo_milo)),
                        subset.nhoods=da_results_cov$NhoodGroup==which(names(vecNames)=="hDA2"))

resultsDE_nhood <- c(dge_hDA1a, dge_hDA1b, dge_hDA2)
names(resultsDE_nhood) <- c("dge_hDA1a","dge_hDA1b","dge_hDA2")

## freeze results
#saveRDS(resultsDE_nhood, file="saved/DEtab/dgeList_dopNctypes_nhood.RDS")
resultsDE_nhood <- readRDS("saved/DEtab/dgeList_dopNctypes_nhood.RDS") 




hDA1a_de <- resultsDE_nhood[[1]][resultsDE_nhood[[1]]$adj.P.Val<0.05 & resultsDE_nhood[[1]]$logFC>log(1.75)/log(10),]
module_hDA1a_de <- rownames(hDA1a_de)
  
hDA1b_de <- resultsDE_nhood[[2]][resultsDE_nhood[[2]]$adj.P.Val<0.05 & resultsDE_nhood[[2]]$logFC>log(1.75)/log(10),]
module_hDA1b_de  <- rownames(hDA1b_de)

hDA2_de <- resultsDE_nhood[[3]][resultsDE_nhood[[3]]$adj.P.Val<0.05 & resultsDE_nhood[[3]]$logFC>log(1.75)/log(10),]
module_hDA2_de  <- rownames(hDA2_de)


list_modules_neurons <- list(module_hDA1a_de, module_hDA1b_de, module_hDA2_de)
names(list_modules_neurons) <- c("module_hDA1a_de","module_hDA1b_de","module_hDA2_de")


upsetPlots <- sapply(1:length(list_modules_neurons), function(x){
  
  tmp <- as.data.frame(list_modules_neurons[[x]])
  colnames(tmp) <- "gene"
  tmp$module <- gsub("module_","",names(list_modules_neurons)[x])
  tmp$module <- gsub("_de","", tmp$module)
  return(tmp)
  
}, simplify=F)

upsetPlots <- do.call("rbind", upsetPlots)


df_upset <- data.frame(genes=names(split(upsetPlots$module, upsetPlots$gene)),
                       modules=NA)

df_upset$modules <- unname(split(upsetPlots$module, upsetPlots$gene))


## Level of intersection of DE genes (2D vs 3D) between dopaminergic cell types

library(ggupset)
plotUpset <- ggplot(data=df_upset, aes(x = modules)) +
  geom_bar() +
  scale_x_upset()+
  theme_bw()+
  ylab("Number of DE genes")


pdf(file=paste0(otherFigs,"upsetPlot_2d_vs_3d.pdf"), width=5.5, height = 6)
plot(plotUpset)
dev.off()


###################
## GO enrichment ##
###################


library(gprofiler2)
library(dplyr)
library(ggplot2)


# geneUniverseTab <- as.data.frame(rownames(querySeurat))
# colnames(geneUniverseTab) <- NULL
# write.table(geneUniverseTab,
#             file="saved/DEgeneUniverse.txt",
#             quote=F,
#             col.names=FALSE,
#             row.names=FALSE,
#             sep="\t")

geneUniverse <- read.table("saved/DEgeneUniverse.txt")

stopifnot(all(sapply(list_modules_neurons, function(y){
  all(y %in% geneUniverse$V1)
}, simplify=T)))


resultsGO <- sapply(1:length(list_modules_neurons), function(x){
  
  geneList_DE <- list_modules_neurons[[x]]
  
  gostres <- gost(query = geneList_DE, 
                  organism = "hsapiens", ordered_query = FALSE, 
                  multi_query = FALSE, significant = TRUE, exclude_iea = TRUE, 
                  measure_underrepresentation = FALSE, evcodes = TRUE, 
                  user_threshold = 0.05, correction_method = "g_SCS", 
                  domain_scope = "custom_annotated", custom_bg = geneUniverse$V1, 
                  numeric_ns = "", sources = c("GO:MF,GO:BP","GO:CC","KEGG"), as_short_link = FALSE, highlight = TRUE)
  
  dfResults <- gostres$result
  dfResults$test <- names(list_modules_neurons)[x]
  return(dfResults)
  
}, simplify=F)


resultsGO_clean <- resultsGO[sapply(resultsGO, function(x) is.data.frame(x))]
resultsGO_clean <- sapply(resultsGO_clean, function(x) {x[x$source=="GO:CC",] }, simplify=F )

allres <- sapply(resultsGO_clean, function(y){
  
  resPerRow <- sapply(1:dim(y)[1], function(z){
    
    ballsFunction <- function(numW, numB, numDrawn, numWdrawn){
      n21 <- numW - numWdrawn
      n12 <- numDrawn - numWdrawn
      n22 <- numB - n12
      odds_ratio <- (numWdrawn * n22)/(n12 * n21)
      expected <- (numWdrawn + n12) * (numWdrawn + n21)
      expected <- expected/(numWdrawn + n12 + n21 + n22)
      
      return(data.frame(OR=odds_ratio, Expected=expected))
    }
    
    #term_size                                            numW: Number of white balls in the urn (representing the number of genes in the reference set associated with the GO term).
    #effective_domain_size - term_size                    numB: Number of black balls in the urn (representing the number of genes not associated with the GO term).
    #query_size                                           numDrawn: Number of balls drawn from the urn (representing the total number of genes in the input list).
    #intersection_size                                    numWdrawn: Number of white balls drawn (representing the number of genes in the input list associated with the GO term).
    
    ballsFunction(numW=y[z,]$term_size,
                  numB=y[z,]$effective_domain_size-y[z,]$term_size,
                  numDrawn=y[z,]$query_size,
                  numWdrawn=y[z,]$intersection_size)
    
    
  }, simplify=F)
  
  resPerRow <- do.call("rbind", resPerRow)
  
  y$OddsRatio <- resPerRow$OR
  y$Expected <- resPerRow$Expected
  
  y <- y %>%
    arrange(desc(OddsRatio), p_value)
  
  return(y)
  
}, simplify=F)


allres <- do.call("rbind", allres)

allresToExport <- allres[,c("term_id","term_name","p_value","OddsRatio",
                            "Expected","intersection_size","query_size",
                            "term_size","precision","recall", "source",
                            "intersection","test","highlighted")]

## Export Supp. Table 4
write.table(allresToExport, file="saved/suppTables/TableS4.txt",
            row.names = F, col.names=T, quote = F, sep="\t")


shared_output <- intersect(intersect(subset(allres, test=="module_hDA1a_de")$term_name,
                                     subset(allres, test=="module_hDA1b_de")$term_name),
                           subset(allres, test=="module_hDA2_de")$term_name)

## plot heatmap with the common GO terms
allres_intersect <- allres[allres$term_name %in% shared_output,]
allres_intersect$testNew <- sapply(strsplit(allres_intersect$test, "_"), function(x) x[2])

allres_intersect$short <- allres_intersect$term_name
mask_long<- nchar(allres_intersect$term_name)>30

allres_intersect$short[mask_long] <- sapply(allres_intersect[mask_long,]$term_name, function(x){
  mask <- unlist(sapply(strsplit(x, " "), function(y) cumsum(nchar(y))>30, simplify=F))
  sapply(strsplit(x, " "), function(y) paste0(paste0(y[!mask], collapse=" "),"\n",paste0(y[mask], collapse=" ")))
})

allres_intersect$short <- factor(allres_intersect$short, levels=unique(allres_intersect$short))


####################
## Main Figure 3H ##
####################

fig3H <- ggplot(data=allres_intersect, aes(x=testNew, y=short, fill=OddsRatio))+
  theme_bw()+
  geom_tile()+
  geom_point(data=allres_intersect, aes(size=-log10(p_value)), inherit.aes=T)+
  ylab("")+xlab("")+
  ggtitle("Shared upregulated BP")+
  scale_fill_viridis_c()+
  theme(plot.title=element_text(hjust=0.5, face="bold"),
        axis.text=element_text(size=10))+
  labs(size=expression("-log"[10]*"pval"))


pdf(file=paste0(rootMain,"mainFigure3H.pdf"), width=6, height = 4)
plot(fig3H)
dev.off()























