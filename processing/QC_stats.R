library(Seurat)
library(readxl)
library(ggplot2)
library(dplyr)
library(harmony)

##QC read

rootPath <- "foetal/cellRangerOutput/"
summfiles <- paste0(rootPath, dir(rootPath),"outs/metrics_summary.csv")

rootPath <- "invitro/run36169/"
summfiles <- c(summfiles,(paste0(rootPath, dir(rootPath),"cellranger-hg38/outs/metrics_summary.csv")))


qcTab <- sapply(summfiles, function(x){
  
  tt <- read.csv(x)
  tt$analysis <- sapply(strsplit(x,"/"), function(x) x[9])
  tt
}, simplify=F)

qcTab <- do.call("rbind", qcTab)
rownames(qcTab) <- NULL
qcTab$file <- summfiles
qcTab$sample <- sapply(strsplit(qcTab$file,"/"), function(x) x[11])
qcTab$timePoint <- NA


## tp!
metadata <- read_excel(path = "metadata/metadata_mid_organoids_10x_sample.xlsx",
                       sheet = "foetal")
metadata <- as.data.frame(metadata)
mask_foetal <- qcTab$analysis=="foetal"


#qcTab$Supplier_Sample_Name <- metadata[match(qcTab$tenxRun, metadata$Sample),]$Supplier_Sample_Name

qcTab$donors <- NA
qcTab$donors[mask_foetal] <- metadata[match(qcTab$sample[mask_foetal], metadata$Supplier_Sample_Name),]$Donors
qcTab$timePoint[mask_foetal] <- gsub("PCW","",qcTab$donors[mask_foetal])


metadata <- read_excel(path = "metadata/metadata_mid_organoids_10x_sample.xlsx",
                       sheet = "invitro")
metadata <- as.data.frame(metadata)

donor_list <- split(metadata$DonorsDimensions, metadata$Sample)
donor_list <- sapply(donor_list[match(qcTab$sample[!mask_foetal], names(donor_list))], function(x) paste(x,collapse=","))
qcTab$donors[!mask_foetal] <- donor_list[match(qcTab$sample[!mask_foetal], names(donor_list))]


vectorSample <- gsub("_.+","",metadata$Supplier_Sample_Name)
names(vectorSample) <- metadata$Sample
vectorSample <- vectorSample[!duplicated(names(vectorSample))]

qcTab$timePoint[!mask_foetal] <- gsub("D","",unname(vectorSample[match(qcTab$sample[!mask_foetal], names(vectorSample))]))

list_10xfiles <- gsub("outs/metrics_summary.csv","outs/filtered_feature_bc_matrix/", qcTab$file)



## The sapply function below computes the percentage of singletons per 10x library. It assumes that demuxlet has been run previously.

perc <- sapply(list_10xfiles, function(x){
  
  mask <- grepl("foetal", x)
  
  if (mask){
    demuxlet_file <- dir(gsub("/filtered_feature_bc_matrix/","", x), full.names=T)[grepl(".best",dir(gsub("/filtered_feature_bc_matrix/","", x)))]
  } else {
    demuxlet_file <- dir(gsub("/filtered_feature_bc_matrix/","", x), full.names=T)[grepl("noFilter.best",dir(gsub("/filtered_feature_bc_matrix/","", x)))]
  }
  demuxlet <- read.table(demuxlet_file, header=T, sep="\t")
  singletons_perc <- unname(round(table(grepl("SNG-", demuxlet$BEST))*100/sum(table(grepl("SNG-", demuxlet$BEST))),2)["TRUE"])
  
  tmp <- data.frame(sample=sapply(strsplit(demuxlet_file, "/"), function(x) x[11]),
              percSingletons=singletons_perc)
  
  return(tmp)
  
}, simplify=F)

perc <- do.call("rbind", perc)
rownames(perc) <- NULL

qcTab <- merge(qcTab, perc, by="sample")
qcTab$sampleIndex <- 1:length(qcTab$sample)


## It saves a list of Seurat objects corresponding to each 10x library.

allobj <- sapply(list_10xfiles, function(x){
  
  tenXrun <- Read10X(data.dir = x)
  tenXrun <- CreateSeuratObject(counts = tenXrun, project = "MLO")
  sample <- sapply(strsplit(x,"/"), function(x) x[11])
  
  mask <- grepl("foetal", x)
  
  if (mask){
    demuxlet_file <- dir(gsub("/filtered_feature_bc_matrix/","", x), full.names=T)[grepl(".best",dir(gsub("/filtered_feature_bc_matrix/","", x)))]
    tenXrun@meta.data$origin <-  "foetal"
    
    
  } else {
    demuxlet_file <- dir(gsub("/filtered_feature_bc_matrix/","", x), full.names=T)[grepl("noFilter.best",dir(gsub("/filtered_feature_bc_matrix/","", x)))]
    tenXrun@meta.data$origin <-  "iPSC"
  
    metadata2 <- read_excel(path = "metadata/metadata_mid_organoids_10x_sample.xlsx",
                            sheet = "invitro")
    metadata2 <- as.data.frame(metadata2)
    metadata2 <- subset(metadata2, Sample==sample)
    
  }
  
  demuxlet <- read.table(demuxlet_file, header=T, sep="\t")
  
  tenXrun <- tenXrun[,which(!is.na(match(rownames(tenXrun@meta.data), demuxlet$BARCODE)))]
  tenXrun <- tenXrun[,which(grepl("SNG-", demuxlet$BEST))]
  demuxlet <- demuxlet[grepl("SNG-", demuxlet$BEST),]
  
  stopifnot(dim(demuxlet)[1]==dim(tenXrun)[2])
  
  tenXrun@meta.data$vcfId <- gsub("SNG-","",demuxlet[match(rownames(tenXrun@meta.data), demuxlet$BARCODE),]$BEST)
  tenXrun@meta.data$midBrainId <- sapply(strsplit(x,"/"), function(y) y[11])
  tenXrun@meta.data$midBrainIndex <- qcTab[match(tenXrun@meta.data$midBrainId, qcTab$sample),]$sampleIndex
  
  metadata <- read_excel(path = "metadata/metadata_mid_organoids_10x_sample.xlsx",
                         sheet = "chipInfo")
  metadata <- as.data.frame(metadata)
  
  vec <- metadata$sample
  names(vec) <- metadata$chipId
  
  tenXrun@meta.data$donorId <- unname(vec[tenXrun@meta.data$vcfId])
  tenXrun@meta.data$donorDimensions <- NA
  tenXrun@meta.data$timePoint <- NA
  
  if (!mask){
    vecMetadata2 <- metadata2$DonorsDimensions
    names(vecMetadata2) <- metadata2$Donors
    tenXrun@meta.data$donorDimensions <- unname(vecMetadata2[tenXrun@meta.data$donorId])
    tenXrun@meta.data$timePoint <- unique(as.numeric(gsub("_.+","",gsub("D","",metadata2$Supplier_Sample_Name))))
    limsNFeatures <- c(1500,6000)
    limsMito <- 10
    
  } else {
    
    tenXrun@meta.data$timePoint <- as.numeric(gsub("PCW","", tenXrun@meta.data$donorId))
    limsNFeatures <- c(1500,5000)
    limsMito <- 25
  }
  
  idDef <- unname(sapply(strsplit(x,"/"), function(y) y[11]))
  tenXrun[["percent.mt"]] <- PercentageFeatureSet(tenXrun, pattern = "^MT-")
  
  plot1 <- FeatureScatter(tenXrun, feature1 = "nCount_RNA", feature2 = "percent.mt")+
    ggtitle(idDef)+theme(plot.title=element_text(hjust=0.5, size=11),
                         legend.position="none")
  plot2 <- FeatureScatter(tenXrun, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
    ggtitle(idDef)+
    theme(plot.title=element_text(hjust=0.5, size=11),
          legend.position="none")
  qcplot1 <- plot1 + plot2
  
  tenXrun[["percent.ribo"]]<- PercentageFeatureSet(tenXrun, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
  ##QC-ribosomal content
  stopifnot(all(mean(tenXrun@meta.data$percent.ribo)>10 | mean(tenXrun@meta.data$percent.ribo)<20))
  
  dirQCplots <- "QC/QCplots/"
  pdf(file=paste0(dirQCplots,"QCplot1_",idDef, ".pdf"))
  plot(qcplot1)
  dev.off()
  
  tmpTabQC <- data.frame(nFeature_RNA=sort(tenXrun$nFeature_RNA),
                         index=1:length(sort(tenXrun$nFeature_RNA)))
  
  qcplot2 <- ggplot(tmpTabQC, aes(x=index, y=nFeature_RNA))+
    geom_point(alpha=0.5, colour="black", size=3, shape=1)+
    theme_bw()+
    ggtitle(paste0("nFeature_RNA distribution: ", idDef))+
    theme(plot.title=element_text(hjust=0.5, size=12, face="bold"))+
    geom_hline(yintercept=limsNFeatures, linetype="dashed", col="red")+
    ylab("Number of expressed genes per cell")+
    xlab("Cells sorted by nFeature_RNA")
  
  pdf(file=paste0(dirQCplots,"QCplot2_",idDef, ".pdf"))
  plot(qcplot2)
  dev.off()
  
  tmpTabQC <- data.frame(percent.mt=sort(tenXrun$percent.mt),
                         index=1:length(sort(tenXrun$percent.mt)))
  
  qcplot3<- ggplot(tmpTabQC, aes(x=index, y=percent.mt))+
    geom_point(alpha=0.5, colour="black", size=3, shape=1)+
    theme_bw()+
    ggtitle(paste0("Mitochondrial content: ", idDef))+
    theme(plot.title=element_text(hjust=0.5, size=14, face="bold"))+
    geom_hline(yintercept=limsMito, linetype="dashed", col="red")+
    ylab("% of Mitochondrial counts")+
    xlab("Cells sorted by mitochondrial content")
  
  pdf(file=paste0(dirQCplots,"QCplot3_",idDef, ".pdf"))
  plot(qcplot3)
  dev.off()
  
  tenXrun <- subset(tenXrun, subset = nFeature_RNA > limsNFeatures[1] & nFeature_RNA < limsNFeatures[2] & percent.mt < limsMito)
  
  tenXrun
  
})


## Produce plots of QC metrics in one final image.

library(scales)

##percSingletons
qcTab$ids <- NA
qcTab$ids[grepl(",", qcTab$timePoint)] <- paste0("Foetal-",qcTab$timePoint[grepl(",", qcTab$timePoint)])
qcTab$ids[!grepl(",", qcTab$timePoint)] <- paste0("Cultures-",qcTab$timePoint[!grepl(",", qcTab$timePoint)])

qcTab$ids <- factor(qcTab$ids, levels=c("Cultures-40","Cultures-70","Cultures-120","Foetal-10,12","Foetal-16,20"))

set.seed(123)
singletons <- ggplot(qcTab, aes(x=ids, y=percSingletons))+
  geom_jitter(width=0.15,height=0,alpha=0.5, size=1.5)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5, size=10),
        axis.text.y=element_text(size=10),
        plot.title=element_text(hjust=0.5))+
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,20))+
  geom_hline(yintercept=80, color="red", linetype="dashed")+
  xlab("")+ylab("")+
  ggtitle("DEMUXLET: % of singl.")


## Mean.Reads.per.Cell
qcTab$Mean.Reads.per.Cell <- as.numeric(gsub(",","", qcTab$Mean.Reads.per.Cell))
meanReads <- ggplot(qcTab, aes(x=ids, y=Mean.Reads.per.Cell))+
  geom_jitter(width=0.15,height=0,alpha=0.5, size=1.5)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5, size=10),
        axis.text.y=element_text(size=10),
        plot.title=element_text(hjust=0.5))+
  geom_hline(yintercept=c(25000,50000), color="red", linetype="dashed")+
  xlab("")+ylab("")+
  ggtitle("Mean reads per cell")+
  scale_y_continuous(labels = comma)


## Sequencing.Saturation
# Sequencing saturation is a measure of the fraction of library complexity that was sequenced in a given experiment.
# The inverse of the sequencing saturation can be interpreted as the number of additional reads it would take to detect a new transcript.

# Sequencing Saturation = 1 - (n_deduped_reads / n_reads)
# 
# where
# 
# n_deduped_reads = Number of unique (valid cell-barcode, valid UMI, gene) combinations among confidently mapped reads. 
# 
# n_reads = Total number of confidently mapped, valid cell-barcode, valid UMI reads.
# 
# Note that the numerator of the fraction is n_deduped_reads, not the non-unique reads that are mentioned in the definition. 
# n_deduped_reads is a degree of uniqueness, not a degree of duplication/saturation.
# Therefore we take the complement of (n_deduped_reads / n_reads) to measure saturation.


## Sequencing.Saturation
qcTab$Sequencing.Saturation <- as.numeric(gsub("%","", qcTab$Sequencing.Saturation))
seqSaturation <- ggplot(qcTab, aes(x=ids, y=Sequencing.Saturation))+
  geom_jitter(width=0.15,height=0,alpha=0.5, size=1.5)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5, size=10),
        axis.text.y=element_text(size=10),
        plot.title=element_text(hjust=0.5))+
  geom_hline(yintercept=25, color="red", linetype="dashed")+
  xlab("")+ylab("")+
  ggtitle("% Seq.Saturation")+
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,20))

## Fraction.Reads.in.Cells
qcTab$Fraction.Reads.in.Cells <- as.numeric(gsub("%","", qcTab$Fraction.Reads.in.Cells))
fracReads <- ggplot(qcTab, aes(x=ids, y=Fraction.Reads.in.Cells))+
  geom_jitter(width=0.15,height=0,alpha=0.5, size=1.5)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5, size=10),
        axis.text.y=element_text(size=10),
        plot.title=element_text(hjust=0.5))+
  geom_hline(yintercept=70, color="red", linetype="dashed")+
  xlab("")+ylab("")+
  ggtitle("% Frac Reads In Cells")+
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,20))

## Total.Genes.Detected
qcTab$Total.Genes.Detected <- as.numeric(gsub(",","", qcTab$Total.Genes.Detected))
totGenesDet <- ggplot(qcTab, aes(x=ids, y=Total.Genes.Detected))+
  geom_jitter(width=0.15,height=0,alpha=0.5, size=1.5)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5, size=10),
        axis.text.y=element_text(size=10),
        plot.title=element_text(hjust=0.5))+
  geom_hline(yintercept=20000, color="red", linetype="dashed")+
  xlab("")+ylab("")+
  ggtitle("Total genes detected")+
  scale_y_continuous(labels = comma, limits=c(0,30000),breaks=seq(0,30000,5000))

## Reads.Mapped.Confidently.to.Transcriptome

qcTab$Reads.Mapped.Confidently.to.Transcriptome <- as.numeric(gsub("%","", qcTab$Reads.Mapped.Confidently.to.Transcriptome))
readsMapToTrans <- ggplot(qcTab, aes(x=ids, y=Reads.Mapped.Confidently.to.Transcriptome))+
  geom_jitter(width=0.15,height=0,alpha=0.5, size=1.5)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5, size=10),
        axis.text.y=element_text(size=10),
        plot.title=element_text(hjust=0.5))+
  geom_hline(yintercept=40, color="red", linetype="dashed")+
  xlab("")+ylab("")+
  ggtitle("% Reads M.Transcr")+
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,20))


figureQC1 <- ggarrange(singletons, meanReads,seqSaturation,fracReads,totGenesDet,readsMapToTrans,
                      ncol=3, nrow=2, common.legend = F)

figureQC1 <- annotate_figure(figureQC1, top=text_grob("QC - Batch2+Foetal",size=15, face="bold"))
 
pathToQC <- "QC/QCplots/"
pdf(file=paste0(pathToQC,"metrics_plot.pdf"))
plot(figureQC1)
dev.off()

write.table(qcTab, file="QC/QCsummary_batch2AndFoetal.txt",
            col.names=T,
            row.names=F,
            quote=F,
            sep="\t")




## The sapply below orders the cells (yet unfiltered) by "nFeature_RNA" and "percent.mt"

dfInfoQC <- sapply(list_10xfiles, function(x){
  
  tenXrun <- Read10X(data.dir = x)
  tenXrun <- CreateSeuratObject(counts = tenXrun, project = "MLO")
  sample <- sapply(strsplit(x,"/"), function(x) x[11])
  
  mask <- grepl("foetal", x)
  
  if (mask){
    demuxlet_file <- dir(gsub("/filtered_feature_bc_matrix/","", x), full.names=T)[grepl(".best",dir(gsub("/filtered_feature_bc_matrix/","", x)))]
    tenXrun@meta.data$origin <-  "foetal"
    
  } else {
    demuxlet_file <- dir(gsub("/filtered_feature_bc_matrix/","", x), full.names=T)[grepl("noFilter.best",dir(gsub("/filtered_feature_bc_matrix/","", x)))]
    tenXrun@meta.data$origin <-  "iPSC"
    
    metadata2 <- read_excel(path = "metadata/metadata_mid_organoids_10x_sample.xlsx",
                            sheet = "invitro")
    metadata2 <- as.data.frame(metadata2)
    metadata2 <- subset(metadata2, Sample==sample)
    
  }
  
  demuxlet <- read.table(demuxlet_file, header=T, sep="\t")
  
  tenXrun <- tenXrun[,which(!is.na(match(rownames(tenXrun@meta.data), demuxlet$BARCODE)))]
  tenXrun <- tenXrun[,which(grepl("SNG-", demuxlet$BEST))]
  demuxlet <- demuxlet[grepl("SNG-", demuxlet$BEST),]
  
  stopifnot(dim(demuxlet)[1]==dim(tenXrun)[2])
  
  tenXrun@meta.data$vcfId <- gsub("SNG-","",demuxlet[match(rownames(tenXrun@meta.data), demuxlet$BARCODE),]$BEST)
  tenXrun@meta.data$midBrainId <- sapply(strsplit(x,"/"), function(y) y[11])
  tenXrun@meta.data$midBrainIndex <- qcTab[match(tenXrun@meta.data$midBrainId, qcTab$sample),]$sampleIndex
  
  metadata <- read_excel(path = "metadata/metadata_mid_organoids_10x_sample.xlsx",
                         sheet = "chipInfo")
  metadata <- as.data.frame(metadata)
  
  vec <- metadata$sample
  names(vec) <- metadata$chipId
  
  tenXrun@meta.data$donorId <- unname(vec[tenXrun@meta.data$vcfId])
  tenXrun@meta.data$donorDimensions <- NA
  tenXrun@meta.data$timePoint <- NA
  
  if (!mask){
    vecMetadata2 <- metadata2$DonorsDimensions
    names(vecMetadata2) <- metadata2$Donors
    tenXrun@meta.data$donorDimensions <- unname(vecMetadata2[tenXrun@meta.data$donorId])
    tenXrun@meta.data$timePoint <- unique(as.numeric(gsub("_.+","",gsub("D","",metadata2$Supplier_Sample_Name))))
    limsNFeatures <- c(1500,6000)
    limsMito <- 10
    
  } else {
    
    tenXrun@meta.data$timePoint <- as.numeric(gsub("PCW","", tenXrun@meta.data$donorId))
    limsNFeatures <- c(1500,5000)
    limsMito <- 25
  }
  
  idDef <- unname(sapply(strsplit(x,"/"), function(y) y[11]))
  tenXrun[["percent.mt"]] <- PercentageFeatureSet(tenXrun, pattern = "^MT-")
  tenXrun[["percent.ribo"]]<- PercentageFeatureSet(tenXrun, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
  stopifnot(all(mean(tenXrun@meta.data$percent.ribo)>10 | mean(tenXrun@meta.data$percent.ribo)<20))
  

  
  vec_ids <- qcTab$ids
  names(vec_ids) <- qcTab$sample
  tmpTabQC <- data.frame(nFeature_RNA=sort(tenXrun$nFeature_RNA),
                         index=1:length(sort(tenXrun$nFeature_RNA)),
                         sample=idDef,
                         sample2=unname(vec_ids[idDef]))
  rownames(tmpTabQC) <- NULL
  
  
  
  tmpTabQC2 <- data.frame(percent.mt=sort(tenXrun$percent.mt),
                         index=1:length(sort(tenXrun$percent.mt)),
                         sample=idDef,
                         sample2=unname(vec_ids[idDef]))
  rownames(tmpTabQC2) <- NULL
  
  

  list(tmpTabQC, tmpTabQC2)
  
  
}, simplify=F)


featureRNA <- do.call("rbind", sapply(dfInfoQC, function(x) x[1]))
rownames(featureRNA) <- NULL

mitoPercent <- do.call("rbind", sapply(dfInfoQC, function(x) x[2]))
rownames(mitoPercent) <- NULL


featuresPlot <- ggplot(featureRNA, aes(x=index, y=nFeature_RNA, group=sample, col=sample))+
  geom_line()+
  theme_bw()+
  theme(legend.position="top",
        plot.title=element_text(hjust=0.5, face="bold", size=13))+
  guides(col=guide_legend(nrow=4,byrow=TRUE))+
  theme(legend.title = element_text( size=8), legend.text=element_text(size=8))+
  scale_color_viridis_d(name="")+
  ggtitle("Mean number of genes per cell")+
  geom_hline(yintercept=c(1500,5000,6000), col="red", linetype="dashed")+
  scale_y_continuous(labels = comma, limits=c(0,12500),breaks=seq(0,12500,2500))+
  xlab("Cells ordered by number of genes expressed (in each 10x Sample)")+
  ylab("")

mitoPlot <- ggplot(mitoPercent, aes(x=index, y=percent.mt, group=sample, col=sample))+
  geom_line()+
  theme_bw()+
  theme(legend.position="top",
        plot.title=element_text(hjust=0.5, face="bold", size=13))+
  guides(col=guide_legend(nrow=4,byrow=TRUE))+
  theme(legend.title = element_text( size=8), legend.text=element_text(size=8))+
  scale_color_viridis_d(name="")+
  ggtitle("% of Mitochondrial Reads")+
  geom_hline(yintercept=c(10,25), col="red", linetype="dashed")+
  scale_y_continuous(labels = comma, limits=c(0,100),breaks=seq(0,100,10))+
  xlab("Cells ordered by % of mito reads (in each 10x Sample)")+
  ylab("")

dirQCplots <- "QC/QCplots/"
pdf(file=paste0(dirQCplots,"QCplot2_numGenesPerCell.pdf"))
plot(featuresPlot)
dev.off()

pdf(file=paste0(dirQCplots,"QCplot3_mitoPerc.pdf"))
plot(mitoPlot)
dev.off()



## The following code produces plots before and after QC, that is comparing CellRanger output with further filtered QC cells. 

##numCells
numCellsTab <- qcTab[,c("sample","Estimated.Number.of.Cells")]
colnames(numCellsTab)[2] <- "numCells" 
numCellsTab$numCells <- as.numeric(gsub(",","", numCellsTab$numCells))
numCellsTab$QC <- "CellRanger"

tmp <- as.data.frame(sapply(allobj, function(x) dim(x)[2]))
tmp$sample <- sapply(strsplit(rownames(tmp),"/"), function(x) x[11])
rownames(tmp) <- NULL
colnames(tmp)[1] <- "numCells"
tmp$QC <- "AfterQC"

numCellsTab <- rbind(numCellsTab, tmp)

vec_qcTab <- qcTab$ids
names(vec_qcTab) <- qcTab$sample
numCellsTab$sample2 <- unname(vec_qcTab[numCellsTab$sample])

numCellsTab$QC <- factor(numCellsTab$QC, levels=c("CellRanger","AfterQC"))
numCellsPlot <- ggplot(numCellsTab, aes(x=sample2, y=numCells, col=QC))+
  geom_jitter(width=0.15,height=0,alpha=0.5, size=1.5)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5, size=10),
        axis.text.y=element_text(size=10),
        plot.title=element_text(hjust=0.5, face="bold", size=6.5))+
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  xlab("")+ylab("")+
  ggtitle("Number of estimated cells")+
  scale_y_continuous(labels = comma, limits=c(0,25000),breaks=seq(0,25000,5000))
  

#Median.UMI.Counts.per.Cell//nCount_RNA
medianUMITab <- qcTab[,c("sample","Median.UMI.Counts.per.Cell")]
colnames(medianUMITab)[2] <- "medianUMI" 
medianUMITab$medianUMI <- as.numeric(gsub(",","", medianUMITab$medianUMI))
medianUMITab$QC <- "CellRanger"

tmp <- as.data.frame(sapply(allobj, function(x) median(x$nCount_RNA)))
tmp$sample <- sapply(strsplit(rownames(tmp),"/"), function(x) x[11])
rownames(tmp) <- NULL
colnames(tmp)[1] <- "medianUMI"
tmp$QC <- "AfterQC"

medianUMITab <- rbind(medianUMITab, tmp)

vec_qcTab <- qcTab$ids
names(vec_qcTab) <- qcTab$sample
medianUMITab$sample2 <- unname(vec_qcTab[medianUMITab$sample])

medianUMITab$QC <- factor(medianUMITab$QC, levels=c("CellRanger","AfterQC"))
medianUMIPlot <- ggplot(medianUMITab, aes(x=sample2, y=medianUMI, col=QC))+
  geom_jitter(width=0.15,height=0,alpha=0.5, size=1.5)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5, size=10),
        axis.text.y=element_text(size=10),
        plot.title=element_text(hjust=0.5, face="bold", size=6.5))+
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  xlab("")+ylab("")+
  ggtitle("Median UMI counts per cell")+
  scale_y_continuous(labels = comma, limits=c(0,12500),breaks=seq(0,12500,2500))


##Median.Genes.per.Cell//nFeature_RNA

medianGenesTab <- qcTab[,c("sample","Median.Genes.per.Cell")]
colnames(medianGenesTab)[2] <- "medianGenes" 
medianGenesTab$medianGenes <- as.numeric(gsub(",","", medianGenesTab$medianGenes))
medianGenesTab$QC <- "CellRanger"

tmp <- as.data.frame(sapply(allobj, function(x) median(x$nFeature_RNA)))
tmp$sample <- sapply(strsplit(rownames(tmp),"/"), function(x) x[11])
rownames(tmp) <- NULL
colnames(tmp)[1] <- "medianGenes"
tmp$QC <- "AfterQC"

medianGenesTab <- rbind(medianGenesTab, tmp)

vec_qcTab <- qcTab$ids
names(vec_qcTab) <- qcTab$sample
medianGenesTab$sample2 <- unname(vec_qcTab[medianGenesTab$sample])

medianGenesTab$QC <- factor(medianGenesTab$QC, levels=c("CellRanger","AfterQC"))
medianGenesPlot <- ggplot(medianGenesTab, aes(x=sample2, y=medianGenes, col=QC))+
  geom_jitter(width=0.15,height=0,alpha=0.5, size=1.5)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5, size=10),
        axis.text.y=element_text(size=10),
        plot.title=element_text(hjust=0.5, face="bold", size=9))+
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  xlab("")+ylab("")+
  ggtitle("Median genes per cell")+
  scale_y_continuous(labels = comma, limits=c(0,3500),breaks=seq(0,3500,500))



figureQC4 <- ggarrange(numCellsPlot,medianUMIPlot,medianGenesPlot,
                       ncol=3, nrow=1, common.legend = T)

figureQC4 <- annotate_figure(figureQC4, top=text_grob("QC - Batch2+Foetal",size=15, face="bold"))

pathToQC <- "QC/QCplots/"
pdf(file=paste0(pathToQC,"metricsBeforeAfterQC.pdf"), width=6, height = 4)
plot(figureQC4)
dev.off()




#####################################
## After QC: Pipeline continuation ##
#####################################

## merge
names(allobj) <- sapply(strsplit(names(allobj), "/"), function(x) x[11])

mlo.big <- merge(allobj[[1]], y=do.call("c",allobj[2:length(allobj)]),
                 add.cell.ids = qcTab[match(names(allobj), qcTab$sample),]$sampleIndex,
                 project="MLO")



## genes expressed in less than 0.1% total cells are removed
selected_f <- names(which(rowSums(mlo.big@assays$RNA[]>0)>0.001*dim(mlo.big)[1]))
mlo.big <- subset(mlo.big, features = selected_f)
mlo.big

#mlo.big@assays$raw <- mlo.big@assays$RNA

## normalise
mlo.big <- NormalizeData(mlo.big, normalization.method = "LogNormalize", scale.factor = 10000)

## Find highly variable features
mlo.big <- FindVariableFeatures(mlo.big, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mlo.big), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mlo.big)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2


# Read in the expression matrix The first row is a header row, the first column is rownames
pathToFile="metadata/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt"
exp.mat <- read.table(file = pathToFile,
                      header = TRUE,
                      as.is = TRUE, row.names = 1)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

mlo.big <- CellCycleScoring(mlo.big, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
mlo.big@meta.data$old.ident <- NULL
# mlo.big <- ScaleData(mlo.big,
#                      vars.to.regress = c("S.Score", "G2M.Score"),
#                      features = rownames(mlo.big))
mlo.big <- ScaleData(mlo.big, features = rownames(mlo.big))
mlo.big <- RunPCA(mlo.big, features = VariableFeatures(object = mlo.big))

# Examine and visualize PCA results a few different ways
print(mlo.big[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(mlo.big, dims = 1:2, reduction = "pca")

## PCA plot with cells colored by 10x library
dirProcessing <- "QC/QCplots/dimRed/"
plotToPdf <- DimPlot(mlo.big, reduction = "pca", pt.size = .1, group.by = 'midBrainIndex')+ggtitle("Before Harmony")
pdf(file=paste0(dirProcessing,"pooled_beforeHarmony.pdf"))
plot(plotToPdf)
dev.off()

## Batch key correction for harmony is the 10x library "midBrainIndex" (cells processed in each GEM well)

options(repr.plot.height = 2.5, repr.plot.width = 6)
mloHarmony <- mlo.big %>% 
  RunHarmony(c("midBrainIndex"), plot_convergence = TRUE)
harmony_embeddings <- Embeddings(mloHarmony, 'harmony')
harmony_embeddings[1:5, 1:5]

## PCA-corrected plot with cells colored by 10x library
plotToPdf <- DimPlot(mloHarmony, reduction = "harmony", pt.size = .1, group.by = 'midBrainIndex')+ggtitle("After Harmony")
pdf(file=paste0(dirProcessing,"pooled_AfterHarmony.pdf"))
plot(plotToPdf)
dev.off()

## PCA-corrected plot with cells colored by time point
plotToPdf <- DimPlot(mloHarmony, reduction = "harmony", pt.size = .1, group.by = 'timePoint')+ggtitle("After Harmony")
pdf(file=paste0(dirProcessing,"pooled_AfterHarmony_timePoint.pdf"))
plot(plotToPdf)
dev.off()

# Use pca-corrected embeddings from harmony to produce UMAP and run clustering
mloHarmony <- mloHarmony %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:10) %>%
  FindClusters(resolution = 0.75) %>%
  identity()

## UMAP plot: Cells colored by unannotated  clusters
plotToPdf <- DimPlot(mloHarmony, reduction = "umap", label = TRUE, pt.size = .1)+theme_gray()+
  theme(plot.title=element_text(hjust=0.5, size=13, face="bold"))+
  ggtitle("Clusters")+theme_bw()
pdf(file=paste0(dirProcessing,"cluteringUMAP_notAnnotated.pdf"))
plot(plotToPdf)
dev.off()

## UMAP plot: Cells colored by time point
plotToPdf <- DimPlot(mloHarmony, reduction = "umap", pt.size = .05, group.by="timePoint")+theme_gray()+
  theme(plot.title=element_text(hjust=0.5, size=13, face="bold"))+theme_bw()
pdf(file=paste0(dirProcessing,"cluteringUMAP_timePoint.pdf"))
plot(plotToPdf)
dev.off()

## UMAP plot: Cells colored by donor
plotToPdf <- DimPlot(mloHarmony, reduction = "umap", pt.size = .05, group.by="donorId")+theme_gray()+
  theme(plot.title=element_text(hjust=0.5, size=13, face="bold"))+theme_bw()
pdf(file=paste0(dirProcessing,"cluteringUMAP_donor.pdf"))
plot(plotToPdf)
dev.off()

## add value to donor dimensions for fetal cells
mloHarmony$donorDimensions[is.na(mloHarmony$donorDimensions)] <- mloHarmony$donorId[is.na(mloHarmony$donorDimensions)]


## UMAP plot: Cells colored by donor and model (2D,3D,foetal)
plotToPdf <- DimPlot(mloHarmony, reduction = "umap", pt.size = .05, group.by="donorDimensions")+theme_gray()+
  theme(plot.title=element_text(hjust=0.5, size=13, face="bold"))+theme_bw()
pdf(file=paste0(dirProcessing,"cluteringUMAP_dimensions.pdf"))
plot(plotToPdf)
dev.off()

## UMAP plot: Cells colored by 10x library
plotToPdf <- DimPlot(mloHarmony, reduction = "umap", pt.size = .05, group.by="midBrainId")+theme_gray()+
  theme(plot.title=element_text(hjust=0.5, size=13, face="bold"))+theme_bw()
pdf(file=paste0(dirProcessing,"cluteringUMAP_midBrainIndex.pdf"))
plot(plotToPdf)
dev.off()

## UMAP plot: Cells colored by cell cycling phase
plotToPdf <- DimPlot(mloHarmony, reduction = "umap", pt.size = .05, group.by="Phase")+theme_gray()+
  theme(plot.title=element_text(hjust=0.5, size=13, face="bold"))+theme_bw()
pdf(file=paste0(dirProcessing,"cluteringUMAP_Phase.pdf"))
plot(plotToPdf)
dev.off()

## UMAP plot: Cells colored by either in vitro (iPSC-derived) or in vivo model (foetal)
plotToPdf <- DimPlot(mloHarmony, reduction = "umap", pt.size = .05, group.by="origin")+theme_gray()+
  theme(plot.title=element_text(hjust=0.5, size=13, face="bold"))+theme_bw()
pdf(file=paste0(dirProcessing,"cluteringUMAP_origin.pdf"))
plot(plotToPdf)
dev.off()

mloHarmony$originDimensions <- mloHarmony$origin
mloHarmony$originDimensions[mloHarmony$originDimensions=="iPSC"] <- unname(gsub(".+_","",mloHarmony$donorDimensions[mloHarmony$originDimensions=="iPSC"]))


## UMAP plot: Cells colored by model, either in vitro 2D, in vitro 3D or foetal.
plotToPdf <- DimPlot(mloHarmony, reduction = "umap", pt.size = .05, group.by="originDimensions")+theme_gray()+
  theme(plot.title=element_text(hjust=0.5, size=13, face="bold"))+theme_bw()
pdf(file=paste0(dirProcessing,"cluteringUMAP_originDimensions.pdf"))
plot(plotToPdf)
dev.off()

mloHarmony$dimensionsTime <- paste0(mloHarmony$originDimensions,"-", mloHarmony$timePoint)

## UMAP plot: Cells colored by model, either in vitro 2D, in vitro 3D or foetal, and the corresponding timepoint.
plotToPdf <- DimPlot(mloHarmony, reduction = "umap", pt.size = .05, group.by="dimensionsTime")+theme_gray()+
  theme(plot.title=element_text(hjust=0.5, size=13, face="bold"))+theme_bw()
pdf(file=paste0(dirProcessing,"cluteringUMAP_dimensionsTime.pdf"))
plot(plotToPdf)
dev.off()


####
####
####


### run FindMarkers: show top-10 DE markers for each unnanotated celltype
mlo.markers <- FindAllMarkers(mloHarmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tabDiff <- mlo.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

tabDiff <- as.data.frame(tabDiff)

dirObjects <- "QC/intermediateTabs/"

write.table(tabDiff, file=paste0(dirObjects, "tableDiff_res75_n10.txt"),
            col.names=T,
            row.names=F,
            quote=F,
            sep="\t")


#scType annotation not used, and therefore not included

### add here the manual annotation from Serena Barral and Dimitri Budinger
clust <- c(0:23)
names(clust) <- c("hRgl2/immAstro","hNbDA","hProgFPM","OPC_1","VascLepto","hDA1b","hRgl1","hDA1a","hRgl3_caudal","hDA2","hProgM",
                  "hPreDA","hMidPre","hMgl","hEndo","hNbGaba","hNPro","hDA3/hGABA/hSer","Unk","hRgl4/MultiEpend","Astro","hPeric","Eryth","OPC_2")
mloHarmony$seurat_clusters_24_Annot <- names(clust[mloHarmony@meta.data$seurat_clusters])


## to keep in GitHub
saveRDS(mloHarmony, "saved/toZenodo/mlo_resolution075_Annot.RDS")











