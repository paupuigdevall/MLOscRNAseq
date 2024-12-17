library(Seurat)
library(miloR)
library(scater)
library(patchwork)
library(dplyr)
library(scran)
library(knitr)
library(ggrepel)


querySeurat <- readRDS("saved/toZenodo/mlo_resolution075_Annot.RDS")
querySeurat$condition <- gsub("[0-9]+","", querySeurat$donorId)


suppTab14 <- querySeurat@meta.data[,c("midBrainId", "donorId", "originDimensions","timePoint", "vcfId", "condition")] 
rownames(suppTab14) <- NULL
suppTab14[suppTab14$condition=="PCW",]$condition <- "Control"
suppTab14 <- distinct(suppTab14)
write.table(suppTab14, "saved/suppTables/TableS14.txt",
            quote=F, col.names=T, row.names=F, sep="\t")

querySeurat <- querySeurat[,querySeurat$originDimensions=="3D"]

clustAnnot <- c(0:23)
names(clustAnnot) <- c("hRgl2/immAstro","hNbDA","hProgFPM","OPC_1","VascLepto","hDA1b","hRgl1","hDA1a","hRgl3_caudal","hDA2","hProgM",
                       "hPreDA","hMidPre","hMgl","hEndo","hNbGaba","hNPro","hDA3/hGABA/hSer","Unk","hRgl4/MultiEpend","Astro","hPeric","Eryth","OPC_2")

querySeurat$seurat_clusters_24_Annot <- names(clustAnnot[querySeurat$seurat_clusters])


## Comparison controls vs patients
table(querySeurat$condition)
querySeurat.sub <- querySeurat[,querySeurat$condition!="Crispr"]

#######################
### Seurat approach ###
#######################

tpoints <- sort(unique(querySeurat.sub$timePoint))

outputPerTP <- sapply(tpoints, function(x){
  
  print(paste0("Processing timepoint day ", x))
  test <- querySeurat.sub[,querySeurat.sub$timePoint==x]
  valid_clusters <- names(which(rowSums(table(test@meta.data$seurat_clusters_24_Annot, test@meta.data$condition)>10)==2))
  
  outputCluster <- sapply(valid_clusters, function(y){
    
    print(paste0("Processing cluster ", y))
    testCluster <-  test[,test$seurat_clusters_24_Annot==y]
    Idents(testCluster) <- "condition"
    de.response <- FindMarkers(testCluster, ident.1 = "Patient", ident.2 = "Control",
                               verbose = FALSE, base = exp(1),
                               logfc.threshold = log10(1.5)/log10(10))
    de.response$timePoint <- x
    de.response$cluster <- y
    de.response$nCellsControl <- unname(table(testCluster$condition)["Control"])
    de.response$nCellsPatient <- unname(table(testCluster$condition)["Patient"])
    
    de.response <- de.response[de.response$p_val_adj<0.05,]
    
    if (dim(de.response)[1]>0){
      
      de.response$direction <- "NA"
      mask_upreg <- de.response$avg_logFC>0
      mask_down <- de.response$avg_logFC<0
      
      if (sum(mask_upreg)){
        de.response[mask_upreg,]$direction <- "up"
      }
      
      if (sum(mask_down)){
        de.response[mask_down,]$direction <- "down"
      }
      
      de.response$geneSymbol <- rownames(de.response)
      rownames(de.response) <- NULL
      de.response <- de.response[,c("geneSymbol",colnames(de.response)[-length(colnames(de.response))])]
      
      return(de.response)
      
    } else {
      NA
    }

    
  }, simplify=F)
  
  
  outputCluster <- outputCluster[!sapply(outputCluster, function(x) all(is.na(x)))]
  outputCluster <- do.call("rbind", outputCluster)
  rownames(outputCluster) <- NULL
  
  return(outputCluster)
  
  
}, simplify=F)


suppTab11 <- do.call("rbind", outputPerTP)
rownames(suppTab11) <- NULL

write.table(suppTab11, "saved/suppTables/TableS11.txt",
            quote=F, col.names=T, row.names=F, sep="\t")


geneUniverseTab <- as.data.frame(rownames(querySeurat))
colnames(geneUniverseTab) <- NULL
write.table(geneUniverseTab,
             file="saved/seurat/DEgeneUniverse.txt",
             quote=F,
             col.names=FALSE,
             row.names=FALSE,
             sep="\t")


