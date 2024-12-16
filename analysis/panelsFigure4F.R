library(Biobase)
library(SeuratWrappers)
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


rootMain <- "figures/main/"
rootSupp <- "figures/supp/"

obj <- readRDS("saved/toZenodo/midBrainIntegration.RDS")

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


obj$samplesToPseudobulk <- as.character(obj$Dataset2)
maskMLO <- obj$Dataset2=="This work (2D, 3D, Foetal)"
mask2d <- obj$system=="2D"
mask3d <- obj$system=="Organoid"
maskfoet <- obj$system=="Foetal"

obj$samplesToPseudobulk[maskMLO & mask2d] <- "This work (2D)"
obj$samplesToPseudobulk[maskMLO & mask3d] <- "This work (3D)"
obj$samplesToPseudobulk[maskMLO & maskfoet] <- "This work (Foetal)"


obj <- obj[,!(grepl("Patient", obj$donorDimensions) | grepl("Crispr", obj$donorDimensions))]
gc()

obj$samplesToPseudobulk <- paste0(obj$timePoint, " - ", obj$samplesToPseudobulk)
obj$samplesToPseudobulk <- gsub("^(.)", "\\U\\1", obj$samplesToPseudobulk, perl = TRUE)

obj$samplesToPseudobulk <- factor(obj$samplesToPseudobulk,
                                  levels=c("Day40 - This work (2D)", "Day70 - This work (2D)",
                                           "Day40 - This work (3D)", "Day70 - This work (3D)", "Day120 - This work (3D)",
                                           "PCW10 - This work (Foetal)","PCW12 - This work (Foetal)","PCW16 - This work (Foetal)","PCW20 - This work (Foetal)",
                                           "GW6 - La Manno et al. 2016 (Foetal)",
                                           "GW7 - La Manno et al. 2016 (Foetal)",
                                           "GW8 - La Manno et al. 2016 (Foetal)",
                                           "GW9 - La Manno et al. 2016 (Foetal)",
                                           "GW10 - La Manno et al. 2016 (Foetal)",
                                           "GW11 - La Manno et al. 2016 (Foetal)",
                                           "PCW8 - Braun et al. 2023 (Foetal)","PCW14 - Braun et al. 2023 (Foetal)",
                                           "PCW6 - Birtele et al. 2022 (Foetal)", "PCW8 - Birtele et al. 2022 (Foetal)", "PCW11 - Birtele et al. 2022 (Foetal)",
                                           "Day15 - Fiorenzano et al. 2021 (3D)", "Day30 - Fiorenzano et al. 2021 (3D)", "Day60 - Fiorenzano et al. 2021 (3D)",
                                           "Day90 - Fiorenzano et al. 2021 (3D)", "Day120 - Fiorenzano et al. 2021 (3D)",
                                           "PostMortem - Agarwal et al. 2020 (PostMortem)",
                                           "Day35 - Zagare et al. 2022 (3D)", "Day70 - Zagare et al. 2022 (3D)"))


########################################
##### Ours vs other datasets ###########
########################################

computePseudoBulk_Sample <- function(obj, samples=samples_mlo, annotCol="samplesToPseudobulk"){
  
  ctypesExpr <- sapply(samples, function(x){
    
    print(x)
    matchAnnot <- match(annotCol, colnames(obj@meta.data))
    tmp <- obj[,obj@meta.data[,matchAnnot]==x]
    gc()
    tmp <- as.matrix(tmp@assays$mnn.reconstructed["data"])
    gc()
    tmp <- data.frame(genes=names(rowMeans(tmp)),
                      avExpr=unname(rowMeans(tmp)),
                      group=x)
    return(tmp)
    
  }, simplify=F)
  
  return(ctypesExpr)
  
}



computeCorrelation <- function(obj1=foetalPseudoTP, obj2=othersTP){

  iterMain <- sapply(obj1, function(x){

    iterObj <- sapply(obj2, function(y){

      stopifnot(all(x$genes==y$genes))

      data.frame(group.x=unique(x$group),
                 group.y=unique(y$group),
                 corrPearson=cor(x$avExpr, y$avExpr, method="pearson"))
    }, simplify=F)

    iterObj <- do.call("rbind", iterObj); rownames(iterObj) <- NULL
    return(iterObj)


  }, simplify=F)

  iterMain <- do.call("rbind", iterMain); rownames(iterMain) <- NULL
  return(iterMain)

}


#####################################
##### In vivo VS In vitro ###########
#####################################


samples_tpFoetal <- levels(obj$samplesToPseudobulk)
samples_inVitro <- samples_tpFoetal[grepl("2D", samples_tpFoetal) | grepl("3D", samples_tpFoetal)]
samples_inVivo <- samples_tpFoetal[grepl("Foetal", samples_tpFoetal) | grepl("PostMortem", samples_tpFoetal) ]

inVitro <- computePseudoBulk_Sample(obj, samples=samples_inVitro, annotCol="samplesToPseudobulk")
inVivo <- computePseudoBulk_Sample(obj, samples=samples_inVivo, annotCol="samplesToPseudobulk")

inVivo_Vs_inVitro_TP <- computeCorrelation(inVivo, inVitro)

library(decoupleR)
library(ComplexHeatmap)

mtrixSampleTP <- pivot_wider_profile(
  data = inVivo_Vs_inVitro_TP,
  id_cols = group.x,
  names_from = group.y,
  values_from = corrPearson,
  to_matrix = TRUE
)


## compute the Ncells per group of comparison

nCellsPerTPFoetal <- obj@meta.data %>%
  group_by(samplesToPseudobulk) %>%
  summarise(cellCounts = n()) %>% as.data.frame()
colnames(nCellsPerTPFoetal)[1] <- "group"

nCellsPerTPFoetal <- nCellsPerTPFoetal[!is.na(nCellsPerTPFoetal$group),]

nCellsPerTPFoetal$group <- gsub("Organoids","3D", nCellsPerTPFoetal$group)
colnames(mtrixSampleTP) <- gsub("Organoids","3D",colnames(mtrixSampleTP))

nCellsInVivo <- nCellsPerTPFoetal[match(rownames(mtrixSampleTP), nCellsPerTPFoetal$group),]
nCellsInVitro <- nCellsPerTPFoetal[match(colnames(mtrixSampleTP), nCellsPerTPFoetal$group),]

row_ha = rowAnnotation(nCells = anno_barplot(setNames(nCellsInVivo$cellCounts, nCellsInVivo$group)))
column_ha = HeatmapAnnotation(nCells = anno_barplot(setNames(nCellsInVitro$cellCounts, nCellsInVitro$group)))


## fix N samples
barplotCorr2 <- Heatmap(mtrixSampleTP, name="Pearson",
                        top_annotation = column_ha, right_annotation = row_ha, column_names_rot = 90,
                        column_names_gp = gpar(fontsize = 8),
                        row_names_gp = gpar(fontsize = 9),
                        cluster_rows=TRUE,
                        cluster_columns=TRUE)



pdf(file=paste0(rootMain,"mainFigure4F.pdf"), width=7, height = 8)
plot(barplotCorr2)
dev.off()



## values in the paper
# Days 15-70 in MLO strongly correlated with gestational weeks (GW6-GW8) fetuses from La Manno
cols_comparison1 <- c("Day40 - This work (3D)","Day70 - This work (3D)", "Day15 - Fiorenzano et al. 2021 (3D)", "Day30 - Fiorenzano et al. 2021 (3D)", "Day60 - Fiorenzano et al. 2021 (3D)")
rows_comparison1 <- c("GW6 - La Manno et al. 2016 (Foetal)","GW7 - La Manno et al. 2016 (Foetal)","GW8 - La Manno et al. 2016 (Foetal)","PCW6 - Birtele et al. 2022 (Foetal)")
round(mean(as.numeric(mtrixSampleTP[match(rows_comparison1, rownames(mtrixSampleTP)),match(cols_comparison1, colnames(mtrixSampleTP))])),2)
#[1] 0.68

# Days 90-120 in MLO correlated with 12 PCW fetus from this study 
cols_comparison2 <- c("Day90 - Fiorenzano et al. 2021 (3D)","Day120 - Fiorenzano et al. 2021 (3D)","Day120 - This work (3D)")
rows_comparison2 <- c("PCW12 - This work (Foetal)")
round(mean(as.numeric(mtrixSampleTP[match(rows_comparison2, rownames(mtrixSampleTP)),match(cols_comparison2, colnames(mtrixSampleTP))])),2)
# [1] 0.47





