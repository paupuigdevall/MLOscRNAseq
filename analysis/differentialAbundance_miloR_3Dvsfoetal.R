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
library(gprofiler2)

######################
#### 3D vs FOETAL ####
######################

otherFigs <- "figures/others/"
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


######################
#### 3D vs foetal ####
######################

querySeurat.sce <- as.SingleCellExperiment(querySeurat)
querySeurat.sce <- querySeurat.sce[,grepl("3D", colData(querySeurat.sce)$dimensionsTime) | grepl("foetal", colData(querySeurat.sce)$dimensionsTime)]
milo.obj <- Milo(querySeurat.sce)
dim(milo.obj)


tmp <- as.data.frame(table(querySeurat.sce$toPlotAnnot))
colnames(tmp) <- c("Cluster", "Number of cells")
kable(tmp)

#otherFigs <- "/home/jovyan/plots/midBrainOrganoid/paperMLO/figuresSingleCell/monocle/foetalvs3d/retest/"

plotUMAP <- plotReducedDim(milo.obj, colour_by="toPlotAnnot",
                           dimred = "UMAP", point_size=0.1)+
  guides(colour = guide_legend(override.aes = list(size=10),
                               title="Clusters"))

pdf(file=paste0(otherFigs, "umap_clusters.pdf"))
plot(plotUMAP)
dev.off()

plotUMAP2 <- plotReducedDim(milo.obj, colour_by="originDimensions",
                            dimred = "UMAP", point_size=0.3)+
  guides(colour = guide_legend(override.aes = list(size=10),
                               title="Clusters"))

pdf(file=paste0(otherFigs, "umap_donorDim.pdf"))
plot(plotUMAP2)
dev.off()

## Construct KNN graph

embryo_milo <- buildGraph(milo.obj, k = 30, d = 30, reduced.dim = "HARMONY")
embryo_milo <- makeNhoods(embryo_milo, prop = 0.3, k = 30, d=30, refined = TRUE, reduced_dims = "HARMONY")


nhoodPlot <- plotNhoodSizeHist(embryo_milo)
pdf(file=paste0(otherFigs, "nhood_Plot.pdf"))
plot(nhoodPlot)
dev.off()



embryo_milo$Sample <- paste0(embryo_milo$donorDimensions,"-",embryo_milo$timePoint)
embryo_milo <- countCells(embryo_milo, meta.data  = data.frame(colData(embryo_milo)), samples="Sample")

tmp <- as.data.frame(table(embryo_milo$Sample))
colnames(tmp) <- c("Sample","Number of cells")
kable(tmp)

embryo_milo <- countCells(embryo_milo, meta.data = data.frame(colData(embryo_milo)), sample="Sample")
head(nhoodCounts(embryo_milo))


##Defining experimental design
embryo_milo$Time <- gsub(".+-","",embryo_milo$dimensionsTime)
embryo_design <- data.frame(colData(embryo_milo))[,c("Sample", "Time","originDimensions")]
kable(head(embryo_design))

embryo_design <- distinct(embryo_design)
rownames(embryo_design) <- embryo_design$Sample

embryo_design <- embryo_design[colnames(nhoodCounts(embryo_milo)), , drop=FALSE]

table(embryo_design$Time)


## Computing neighbourhood connectivity
embryo_milo <- calcNhoodDistance(embryo_milo, d=30, reduced.dim = "HARMONY")
da_results <- testNhoods(embryo_milo, design = ~ originDimensions, design.df = embryo_design)


## without covariates

tmp2 <- da_results %>%
  arrange(SpatialFDR) %>%
  head()
rownames(tmp2) <- NULL
kable(tmp2)

pvaluehist2 <- ggplot(da_results, aes(PValue)) +
  geom_histogram(bins=50, col="white", fill="black")+
  theme_bw()

volcanoplot2 <- ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point(size=0.5) +
  geom_hline(yintercept = 1)+ ## Mark significance threshold (10% FDR)
  theme_bw()

pdf(file=paste0(otherFigs, "histPval_NoCov.pdf"))
plot(pvaluehist2)
dev.off()

pdf(file=paste0(otherFigs, "volcanoPlot_NoCov.pdf"))
plot(volcanoplot2)
dev.off()


embryo_milo <- buildNhoodGraph(embryo_milo)

umap_pl <- plotReducedDim(embryo_milo, dimred = "UMAP", colour_by="originDimensions", text_by = "toPlotAnnot", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")
nh_graph_pl <- plotNhoodGraphDA(embryo_milo, da_results, layout="UMAP",alpha=0.1) 
combined_pl <- umap_pl + nh_graph_pl +
  plot_layout(guides="collect")

pdf(file=paste0(otherFigs, "umapAndNhoodGraph.pdf"))
plot(combined_pl)
dev.off()


da_results <- annotateNhoods(embryo_milo, da_results, coldata_col = "toPlotAnnot")
kable(head(da_results))

histFrac <- ggplot(da_results, aes(toPlotAnnot_fraction)) +
  geom_histogram(bins=50, col="white", fill="black")+
  theme_bw()+coord_cartesian(xlim=c(0,1))

pdf(file=paste0(otherFigs, "hist_annotateNhoods.pdf"))
plot(histFrac)
dev.off()

da_results$celltype <- ifelse(da_results$toPlotAnnot_fraction < 0.7, "Mixed", da_results$toPlotAnnot)


## from foetal part

mask_foetal <- da_results$FDR<0.1 & da_results$logFC>0
mask_3d <- da_results$FDR<0.1 & da_results$logFC<0
da_results$dotGroup <- "Non signif."
da_results$dotGroup[mask_foetal] <- "Foetal"
da_results$dotGroup[mask_3d] <- "3D"
da_results2 <- subset(da_results, celltype!="Mixed")

ctypesToRemove <- names(which(rowSums(table(da_results2$celltype, da_results2$SpatialFDR<0.1))<20))
da_results2 <- da_results2[!da_results2$celltype %in% ctypesToRemove,]

da_results2$celltype <- factor(da_results2$celltype,
                               levels=rev(sort(unique(da_results2$celltype))))



####################
## Main Figure 3G ##
####################


fig3G <- ggplot(da_results2, aes(x=logFC, y=celltype, col=dotGroup))+
  theme_bw()+
  geom_quasirandom(alpha=0.7, groupOnX = F)+
  #scale_color_manual(name="FDR<0.1",values=c("#663300","#000099","grey80"))+
  scale_color_manual(name="FDR<0.1",values=c("#377EB8","#4DAF4A","grey80"))+
  xlab("Log Fold-Change (~Model)")+
  ylab("")+
  theme(axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=11, angle=90, hjust=1, vjust=0.5),
        axis.title=element_text(size=15),
        legend.position="top")+
  scale_x_continuous(limits=c(-8,8), breaks=seq(-8,8,2))+
  guides(color = guide_legend(override.aes = list(size = 9)))


pdf(file=paste0(rootMain,"mainFigure3G.pdf"), width=5.5, height = 6)
plot(fig3G)
dev.off()


#############################
## Differential expression ##
#############################

querySeurat <- readRDS("saved/toZenodo/mlo_resolution075_Annot.RDS")

clustAnnot <- c(0:23)
names(clustAnnot) <- c("hRgl2/immAstro","hNbDA","hProgFPM","OPC_1","VascLepto","hDA1b","hRgl1","hDA1a","hRgl3_caudal","hDA2","hProgM",
                       "hPreDA","hMidPre","hMgl","hEndo","hNbGaba","hNPro","hDA3/hGABA/hSer","Unk","hRgl4/MultiEpend","Astro","hPeric","Eryth","OPC_2")

### 

querySeurat$seurat_clusters_24_Annot <- names(clustAnnot[querySeurat$seurat_clusters])
querySeurat$toPlotAnnot <- querySeurat$seurat_clusters_24_Annot
querySeurat$toPlotAnnot <- gsub("_","", querySeurat$toPlotAnnot)


##################################
### Upregulated (Foetal vs 3D) ###
##################################

## DE and GO ontology
foldChange_base10 <- log(1.5)/log(2)

## DE expression: Tissue (hDA2) vs Organoids (hDA2)
obj.onlyMLOall <- querySeurat[,querySeurat$originDimensions!="2D"]
obj.onlyMLOall.comp3 <- obj.onlyMLOall[,obj.onlyMLOall$seurat_clusters_24_Annot=="hDA2"]
Idents(obj.onlyMLOall.comp3) <- "originDimensions"
resDE.comp3 <- FindMarkers(obj.onlyMLOall.comp3, ident.1="foetal", ident.2="3D")
degenes.comp3 <- resDE.comp3[which(resDE.comp3$p_val_adj<0.05 & resDE.comp3$avg_log2FC>foldChange_base10),]
rm(obj.onlyMLOall.comp3)
gc()

source("analysis/functionsToImport.R")

gtf <- rtracklayer::import('refgenomes/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf')
gtf_df=as.data.frame(gtf)
genesGTF <- subset(gtf_df, type=="gene")

geneUniverse <- data.frame(symbol=rownames(obj.onlyMLOall))
geneUniverse$symbol_amend <- NA
geneUniverse$symbol_amend[geneUniverse$symbol %in% gtf_df$gene_name] <- geneUniverse$symbol[geneUniverse$symbol %in% gtf_df$gene_name]
geneUniverse$symbol_amend[!geneUniverse$symbol %in% gtf_df$gene_name] <- gsub("\\..+","",geneUniverse$symbol[!geneUniverse$symbol %in% gtf_df$gene_name])
geneUniverse$ensembl <- genesGTF[match(geneUniverse$symbol_amend, genesGTF$gene_name),]$gene_id


tabCorr <- AnnotationDbi::select(org.Hs.eg.db, geneUniverse$ensembl, "ENTREZID", "ENSEMBL")
tabCorr <- tabCorr[!is.na(tabCorr$ENTREZID),]
tabCorr <- tabCorr[!duplicated(tabCorr$ENSEMBL),]
geneUniverse$entrezid <- tabCorr[match(geneUniverse$ensembl, tabCorr$ENSEMBL),]$ENTREZID

list_de_genes <- list(rownames(degenes.comp3))

list_de_genes_entrezid <- sapply(list_de_genes, function(x){
  moduleGenes <- geneUniverse[match(x, geneUniverse$symbol),]$entrezid
  moduleGenes <- unique(moduleGenes[!is.na(moduleGenes)])
  moduleGenes
}, simplify=F)


names(list_de_genes_entrezid) <- c("hDA2 (Foetal) vs hDA2 (Organoids)")

geneUniverse_entrezid <- unique(geneUniverse[!is.na(geneUniverse$entrezid),]$entrezid)

library(GOstats)
report_testClust <- sapply(1:length(list_de_genes_entrezid), function(x){
  
  labelClust= names(list_de_genes_entrezid)[x]
  report_testClust <- GOenrichmentAndReport(list_de_genes_entrezid[[x]], geneUniverse_entrezid, minSize=25, maxSize=500, minCount=10, p.value=0.05, label=labelClust)
  report_testClust$GeneSyms <- NULL
  report_testClust
}, simplify=F)

names(report_testClust) <- names(list_de_genes_entrezid)
saveRDS(report_testClust, file="saved/GOresults_DEcomparisons_upregulated.RDS")


####################
## Main Figure 3I ##
####################

## Top-5 enriched GO terms

nameObject <- "hDA2 (Foetal) vs hDA2 (Organoids)"
modName <- c("hDA1a (All) vs hDA1b (All)"="hDA1a_all_VS_hDA1b_all",
             "hDA1a/b (Foetal+Organoids) vs hDA2 (Foetal+Organoids)"="hDA1a-b_foetOrg_VS_hDA2_foetOrg",
             "hDA2 (Foetal) vs hDA2 (Organoids)"="hDA2_foet_VS_hDA2_Org")
nameObject2 <- unname(modName[nameObject])

kk <- report_testClust[[1]]

if (dim(kk)[1]>=5){
  kk <- kk[1:5,]
} else {
  kk <- kk[1:dim(kk)[1],]
}

kk$short <- kk$Term
mask_long<- nchar(kk$Term)>30

kk$short[mask_long] <- sapply(kk[mask_long,]$Term, function(x){
  
  mask <- unlist(sapply(strsplit(x, " "), function(y) cumsum(nchar(y))>30, simplify=F))
  #sapply(strsplit(x, " "), function(y) paste0(paste0(y[!mask], collapse=" "),"\n", y[mask]))
  sapply(strsplit(x, " "), function(y) paste0(paste0(y[!mask], collapse=" "),"\n",paste0(y[mask], collapse=" ")))
  
})
kk <- kk[order(kk$OddsRatio, decreasing=F),]

kk$short <- factor(kk$short, levels=kk$short)

fig3I <- ggplot(kk, aes(x=OddsRatio, y=short, size=Count, col=-log10(Pvalue)))+
  geom_point()+
  ggtitle("")+
  theme_bw()+
  scale_colour_gradient(name=expression("-log"[10]*"pval"), low = "red", high = "blue", na.value = NA)+
  xlab("Odds Ratio")+
  #ylab("GO term :: Biological process")+
  ylab("")+
  theme(axis.text.y=element_text(size=10),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(hjust=0.5, size=14, face="bold"))

pdf(file=paste0(rootMain,"mainFigure3I.pdf"), width=6, height = 3.5)
plot(fig3I)
dev.off()


####################
## Supp Figure 3C ##
####################

## Top-10 enriched GO terms

kk <- report_testClust[[1]]

if (dim(kk)[1]>=10){
  kk <- kk[1:10,]
} else {
  kk <- kk[1:dim(kk)[1],]
}

kk$short <- kk$Term
mask_long<- nchar(kk$Term)>30

kk$short[mask_long] <- sapply(kk[mask_long,]$Term, function(x){
  
  mask <- unlist(sapply(strsplit(x, " "), function(y) cumsum(nchar(y))>30, simplify=F))
  #sapply(strsplit(x, " "), function(y) paste0(paste0(y[!mask], collapse=" "),"\n", y[mask]))
  sapply(strsplit(x, " "), function(y) paste0(paste0(y[!mask], collapse=" "),"\n",paste0(y[mask], collapse=" ")))
  
})
kk <- kk[order(kk$OddsRatio, decreasing=F),]

kk$short <- factor(kk$short, levels=kk$short)

figS3C <- ggplot(kk, aes(x=OddsRatio, y=short, size=Count, col=-log10(Pvalue)))+
  geom_point()+
  ggtitle("")+
  theme_bw()+
  scale_colour_gradient(name=expression("-log"[10]*"pval"), low = "red", high = "blue", na.value = NA)+
  xlab("Odds Ratio")+
  #ylab("GO term :: Biological process")+
  ylab("")+
  theme(axis.text.y=element_text(size=10),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(hjust=0.5, size=14, face="bold"))

pdf(file=paste0(rootSupp,"suppFigure3C.pdf"), width=6, height = 5)
plot(figS3C)
dev.off()



####################################
### Downregulated (Foetal vs 3D) ###
####################################


## DE and GO ontology
foldChange_base10 <- -log(1.5)/log(2)

## DE expression: Tissue (hDA2) vs Organoids (hDA2)
obj.onlyMLOall <- querySeurat[,querySeurat$originDimensions!="2D"]
obj.onlyMLOall.comp3 <- obj.onlyMLOall[,obj.onlyMLOall$seurat_clusters_24_Annot=="hDA2"]
Idents(obj.onlyMLOall.comp3) <- "originDimensions"
resDE.comp3 <- FindMarkers(obj.onlyMLOall.comp3, ident.1="foetal", ident.2="3D")
degenes.comp3 <- resDE.comp3[which(resDE.comp3$p_val_adj<0.05 & resDE.comp3$avg_log2FC<foldChange_base10),]
rm(obj.onlyMLOall.comp3)
gc()

source("analysis/functionsToImport.R")

gtf <- rtracklayer::import('refgenomes/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf')
gtf_df=as.data.frame(gtf)
genesGTF <- subset(gtf_df, type=="gene")

geneUniverse <- data.frame(symbol=rownames(obj.onlyMLOall))
geneUniverse$symbol_amend <- NA
geneUniverse$symbol_amend[geneUniverse$symbol %in% gtf_df$gene_name] <- geneUniverse$symbol[geneUniverse$symbol %in% gtf_df$gene_name]
geneUniverse$symbol_amend[!geneUniverse$symbol %in% gtf_df$gene_name] <- gsub("\\..+","",geneUniverse$symbol[!geneUniverse$symbol %in% gtf_df$gene_name])
geneUniverse$ensembl <- genesGTF[match(geneUniverse$symbol_amend, genesGTF$gene_name),]$gene_id


tabCorr <- AnnotationDbi::select(org.Hs.eg.db, geneUniverse$ensembl, "ENTREZID", "ENSEMBL")
tabCorr <- tabCorr[!is.na(tabCorr$ENTREZID),]
tabCorr <- tabCorr[!duplicated(tabCorr$ENSEMBL),]
geneUniverse$entrezid <- tabCorr[match(geneUniverse$ensembl, tabCorr$ENSEMBL),]$ENTREZID

list_de_genes <- list(rownames(degenes.comp3))

list_de_genes_entrezid <- sapply(list_de_genes, function(x){
  moduleGenes <- geneUniverse[match(x, geneUniverse$symbol),]$entrezid
  moduleGenes <- unique(moduleGenes[!is.na(moduleGenes)])
  moduleGenes
}, simplify=F)


names(list_de_genes_entrezid) <- c("hDA2 (Foetal) vs hDA2 (Organoids)")


geneUniverse_entrezid <- unique(geneUniverse[!is.na(geneUniverse$entrezid),]$entrezid)

library(GOstats)
report_testClust <- sapply(1:length(list_de_genes_entrezid), function(x){
  
  labelClust= names(list_de_genes_entrezid)[x]
  report_testClust <- GOenrichmentAndReport(list_de_genes_entrezid[[x]], geneUniverse_entrezid, minSize=25, maxSize=500, minCount=10, p.value=0.05, label=labelClust)
  report_testClust$GeneSyms <- NULL
  report_testClust
}, simplify=F)

names(report_testClust) <- names(list_de_genes_entrezid)
saveRDS(report_testClust, file="saved/GOresults_DEcomparisons_downregulated.RDS")


####################
## Supp Figure 3D ##
####################


nameObject <- "hDA2 (Foetal) vs hDA2 (Organoids)"
modName <- c("hDA1a (All) vs hDA1b (All)"="hDA1a_all_VS_hDA1b_all",
             "hDA1a/b (Foetal+Organoids) vs hDA2 (Foetal+Organoids)"="hDA1a-b_foetOrg_VS_hDA2_foetOrg",
             "hDA2 (Foetal) vs hDA2 (Organoids)"="hDA2_foet_VS_hDA2_Org")
nameObject2 <- unname(modName[nameObject])

kk <- report_testClust[[1]]

if (dim(kk)[1]>=10){
  kk <- kk[1:10,]
} else {
  kk <- kk[1:dim(kk)[1],]
}

kk$short <- kk$Term
mask_long<- nchar(kk$Term)>30

kk$short[mask_long] <- sapply(kk[mask_long,]$Term, function(x){
  
  mask <- unlist(sapply(strsplit(x, " "), function(y) cumsum(nchar(y))>30, simplify=F))
  #sapply(strsplit(x, " "), function(y) paste0(paste0(y[!mask], collapse=" "),"\n", y[mask]))
  sapply(strsplit(x, " "), function(y) paste0(paste0(y[!mask], collapse=" "),"\n",paste0(y[mask], collapse=" ")))
  
})

kk <- kk[order(kk$OddsRatio, decreasing=F),]

kk$short <- factor(kk$short, levels=kk$short)

figS3D <- ggplot(kk, aes(x=OddsRatio, y=short, size=Count, col=-log10(Pvalue)))+
  geom_point()+
  ggtitle("")+
  theme_bw()+
  scale_colour_gradient(name=expression("-log"[10]*"pval"), low = "red", high = "blue", na.value = NA)+
  xlab("Odds Ratio")+
  ylab("")+
  theme(axis.text.y=element_text(size=10),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(hjust=0.5, size=14, face="bold"))

pdf(file=paste0(rootSupp,"suppFigure3D.pdf"))
plot(figS3D)
dev.off()


suppTab_up <- readRDS("saved/GOresults_DEcomparisons_upregulated.RDS")[[1]]
suppTab_down <- readRDS("saved/GOresults_DEcomparisons_downregulated.RDS")[[1]]

suppTab_up$direction <- "Upregulated"
suppTab_down$direction <- "Downregulated"

suppTab <- rbind(suppTab_up, suppTab_down)

## Export Table S5
write.table(suppTab, file="saved/suppTables/TableS5.txt",
            row.names = F, col.names=T, quote = F, sep="\t")






