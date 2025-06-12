library(tradeSeq)
library(ggplot2)
library(tidyverse)
library(scales)
library(scico)

set.seed(8)

pathToDir <- "saved/scanpy/"
rootMain <- "figures/main/"
rootSupp <- "figures/supp/"
rootOthers <- "figures/others/"

w <- as.matrix(read.csv(paste0(pathToDir,"cellWeights_astrocytes.csv"), header = FALSE))
dpt <- as.matrix(read.csv(paste0(pathToDir,"pseudotime_astrocytes.csv"), header = FALSE))
cMatrix <- as.matrix(read.csv(paste0(pathToDir,"counts_astrocytes.csv"), header = FALSE, row.names = 1))

print("Starting fitGAM")
gamObj <- fitGAM(cMatrix, verbose = TRUE, pseudotime = dpt, cellWeights = w, nknots = 8, sce=FALSE)
names(gamObj) <- rownames(cMatrix)
print("fitGAM finalised")

dptseq <- seq(min(dpt), max(dpt), length.out = 5)

startRes <- startVsEndTest(gamObj, pseudotimeValues = c(min(dpt)+.01,max(dpt)-.01))

startResFilt <- startRes[startRes$pvalue <=0.001 & abs(startRes$logFClineage1)>=2,]
startResFilt$Gene <- rownames(startResFilt)
startResFilt$test <- "startVsEndTest"

nGenes <- 3
branchTag <- "Astro"

# Top-decreasing (positive to negative)
startResGenesPositive <- startResFilt[startResFilt$logFClineage1<0,]
startResGenesPositive <- startResGenesPositive[order(startResGenesPositive$pvalue,
                                                     startResGenesPositive$logFClineage1), ]
startResGenesPositive_plot <- startResGenesPositive[1:nGenes,]$Gene
startResGenesPositive$PatternType <- "decreasing"
colnames(startResGenesPositive)[colnames(startResGenesPositive)=="logFClineage1"] <- "logFC"

## Top-increasing (negative to positive)
startResGenesNegatives <- startResFilt[startResFilt$logFClineage1>0,]
startResGenesNegatives <- startResGenesNegatives[order(startResGenesNegatives$pvalue,
                                                       -startResGenesNegatives$logFClineage1), ]
startResGenesNegatives_plot <- startResGenesNegatives[1:nGenes,]$Gene
startResGenesNegatives$PatternType <- "increasing"
colnames(startResGenesNegatives)[colnames(startResGenesNegatives)=="logFClineage1"] <- "logFC"

rownames(startResGenesPositive) <- NULL
rownames(startResGenesNegatives) <- NULL


######################
## Export Table S16 ##
######################

CombinedDF <- rbind(startResGenesPositive, startResGenesNegatives)
CombinedDF$Branch <- branchTag
write.table(CombinedDF, file="saved/suppTables/TableS16.txt", row.names = F, col.names = T, sep="\t", quote=F)

ListOfInterest <- c("AQP4","GFAP","SLC1A2","SLC6A11","IL33","S100B","ALDH1L1","GJA1","HOPX","ZBTB20")

.getPredictRangeDf <- function(dm, lineageId, conditionId = NULL, nPoints = 100){
  vars <- dm[1, ]
  if ("y" %in% colnames(vars)) {
    vars <- vars[!colnames(vars) %in% "y"]
    off <- 1
  } else {
    off <- 0
  }
  offsetId <- grep(x = colnames(vars), pattern = "offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
  names(vars)[offsetId] <- offsetName
  # set all times on 0
  vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
  # set all lineages on 0
  vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
  # duplicate to nPoints
  vars <- rbind(vars, vars[rep(1, nPoints - 1), ])
  # set range of pseudotime for lineage of interest
  if (is.null(conditionId)) {
    lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId))  
  } else {
    lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId, conditionId))
  }
  if (length(lineageIds) == 1){
    lineageData <- dm[dm[, lineageIds + off] == 1,
                      paste0("t", lineageId)]
  } else {
    lineageData <- dm[rowSums(dm[, lineageIds + off]) == 1,
                      paste0("t", lineageId)]
  }
  # make sure lineage starts at zero
  if(min(lineageData) / max(lineageData) < .01) {
    lineageData[which.min(lineageData)] <- 0
  }
  vars[, lineageIds] <- 1 / length(lineageIds)
  # set lineage
  vars[, paste0("t", lineageId)] <- seq(min(lineageData),
                                        max(lineageData),
                                        length = nPoints)
  # set offset
  vars[, offsetName] <- mean(dm[, grep(x = colnames(dm),
                                       pattern = "offset")])
  return(vars)
}


datalist <- list()
for (g in ListOfInterest){

  localModel <- gamObj[[g]]
  data <- localModel$model
  y <- data$y
  nCurves <- length(localModel$smooth)

  for (jj in seq_len(nCurves)) {
    df <- .getPredictRangeDf(localModel$model, jj, nPoints = 5)
    yhat <- predict(localModel, newdata = df, type = "response")

    Newframe <- data.frame("fittedCounts" = yhat)
    colnames(Newframe) <- c(paste0("fittedCounts.",g))

    datalist[[paste(g,jj)]] <- Newframe

  }

}

assoResSSPandas_alt2 =  do.call(cbind, datalist)

print("Starting writing tables")
write.table(assoResSSPandas_alt2, file=paste0(pathToDir,"assoResSSPandas.alt.bin5.",branchTag,".tsv"), row.names = T, col.names = T, sep="\t", quote=F)



####################
## Supp Figure 9P ##
####################

tt <- read.table(paste0(pathToDir,"assoResSSPandas.alt.bin5.",branchTag,".tsv"))
colnames(tt) <- gsub("fittedCounts\\.","", colnames(tt))
tt <- as.data.frame(tt)
tt$bucket <- paste0("b",1:dim(tt)[1])

vecNames <- c("Transporter" = "AQP4",
	      "Transporter" = "SLC1A2",
	      "Transporter" = "SLC6A11",
	      "Transporter" = "GJA1",
	      "Structural" = "GFAP",
	      "Structural" = "S100B",
	      "Enzymes" = "ALDH1L1",
	      "Signalling" = "IL33",
	      "TFs" = "HOPX",
	      "TFs" = "ZBTB20")

genesOfInterest <- as.data.frame(tt %>% pivot_longer(-c("bucket"), names_to="genes", values_to="fittedExpression"))
genesOfInterest$type <- names(vecNames[match(genesOfInterest$genes, unname(vecNames))])

genesOfInterest$significant <- ifelse(genesOfInterest$genes %in% CombinedDF$Gene, "p<0.05","p>0.05")
genesOfInterest$logFC_trajectory <- CombinedDF[match(genesOfInterest$gene, CombinedDF$Gene),]$logFC
genesOfInterest$direction <- ifelse(CombinedDF[match(genesOfInterest$gene, CombinedDF$Gene),]$logFC>0, "Trajectory-up", "Trajectory-down")
genesOfInterest$direction[is.na(genesOfInterest$direction)] <- "No change"


figS9P <- ggplot(data = genesOfInterest, aes(x = bucket, y = fittedExpression, fill=type, col=type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~genes, scales = "free", nrow = 1) +  # one row
  theme_void() +  # remove everything
  theme(
    strip.background = element_blank(),       # remove box around facet title
    strip.text = element_text(size = 12),     # keep facet title
    panel.spacing = unit(1, "lines"),         # optional spacing between panels
    plot.margin = margin(10, 10, 10, 10)      # optional outer margin
  )+
  scale_fill_scico_d(name="",palette = "batlow")+
  scale_color_scico_d(name="", palette = "batlow")+
  theme(legend.position="top")


pdf(file = paste0(rootSupp, "suppFigure9P.pdf"), width=10, height=1.5)
plot(figS9P)
dev.off()
