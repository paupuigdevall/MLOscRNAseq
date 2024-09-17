

GOenrichmentAndReport <- function(geneIds, universeGeneIds,
                                  minSize=3, maxSize=300, minCount=3,
                                  minOddsRatio=1.5, p.value=0.05, highlightGenes=NULL, highlightStr="*%s*",
                                  label="allDE"){
  
  
  GOparams <- new("GOHyperGParams", geneIds=geneIds,
                  universeGeneIds=universeGeneIds,
                  annotation="org.Hs.eg.db",
                  ontology="BP", pvalueCutoff=0.05, conditional=TRUE,
                  minSizeCutoff=3, maxSizeCutoff=300, orCutoff=1.5,
                  testDirection="over")
  
  hgOverGOBP <- hyperGTest(GOparams)
  
  report <- data.frame(GOBPID=as.vector(names(geneIdUniverse(hgOverGOBP))),
                       Pvalue=pvalues(hgOverGOBP),
                       OddsRatio=oddsRatios(hgOverGOBP),
                       ExpCount=expectedCounts(hgOverGOBP),
                       Count=geneCounts(hgOverGOBP),
                       Size=universeCounts(hgOverGOBP),
                       stringsAsFactors=FALSE)
  
  ## discard gene sets that do not meet a minimum and maximum number of genes
  report <- report[report$Size >= minSize & report$Size <= maxSize, , drop=FALSE]
  
  ## discard gene sets that show a p.value>0.05
  report <- report[report$Pvalue < p.value, , drop=FALSE]
  
  ## discard gene sets that do not satisfy the OR cutoff
  report <- report[report$OddsRatio >= minOddsRatio & report$Count >= minCount, , drop=FALSE]
  
  ## apply the maximum-number-of-GO-terms-reported cutoff and sort by odds ratio
  maxReported <- min(nrow(report))
  report <- report[sort(report$OddsRatio, decreasing=TRUE, index.return=TRUE)$ix[1:maxReported], ]
  
  if (dim(report[complete.cases(report),])[1]==0){
    message <- "No GO terms enriched"
    print(message)
    return(message)
  }
  
  ## add the symbol and GO term description, information
  reportGenes <- geneIdsByCategory(hgOverGOBP, report$GOBPID)
  reportGeneSyms <- lapply(reportGenes, annotate::getSYMBOL, "org.Hs.eg.db")
  highlightGeneMask <- lapply(reportGenes, function(x, hgenes, fmt) x %in% hgenes, highlightGenes)
  reportGeneSyms <- mapply(function(genes, mask, fmt) ifelse(mask, sprintf(fmt, genes), genes),
                           reportGeneSyms, highlightGeneMask, MoreArgs=list(fmt=highlightStr), SIMPLIFY=FALSE)
  reportGeneSyms <- sapply(reportGeneSyms, paste, collapse=", ")
  reportGenes <- sapply(reportGenes, paste, collapse=", ")
  report <- cbind(report,
                  Term=sapply(mget(report$GOBPID, GOstats:::GOenv("TERM")), Term),
                  GeneSyms=reportGeneSyms,
                  label=label)
  rownames(report) <- NULL
  
  return(report)
  
}



