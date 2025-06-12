library(gprofiler2)
library(dplyr)
library(ggplot2)

rootMain <- "figures/main/"
rootSupp <- "figures/supp/"
rootDir <- "otherPlots/monocle3/"

##ctrl vs patients

geneUniverse <- read.table("saved/seurat/DEgeneUniverse.txt")
resDE <- read.table("saved/suppTables/TableS11.txt", header=TRUE)
stopifnot(all(resDE$geneSymbol %in% geneUniverse$V1))


##controls vs patients in D40 organoids, looking at hDA1, MidPre, and Rgl genes

listGOtests <- list(c(40,"hDA1a","up"),c(40,"hMidPre","up"),c(40,"hRgl1","up"),
                    c(40,"hDA1a","down"),c(40,"hMidPre","down"),c(40,"hRgl1","down"),
                    c(70,"hDA1a","up"),c(70,"hMidPre","up"),c(70,"hRgl1","up"),
                    c(70,"hDA1a","down"),c(70,"hMidPre","down"),c(70,"hRgl1","down"),
                    c(120,"hDA1a","up"),c(120,"hMidPre","up"),c(120,"hRgl1","up"),
                    c(120,"hDA1a","down"),c(120,"hMidPre","down"),c(120,"hRgl1","down"))



names(listGOtests) <- sapply(listGOtests, function(x) paste(x, collapse="/"))

resultsGO <- sapply(1:length(listGOtests), function(x){
  
  geneList_DE <- subset(resDE, timePoint==listGOtests[[x]][1] & cluster==listGOtests[[x]][2] & direction==listGOtests[[x]][3])$geneSymbol
  
  gostres <- gost(query = geneList_DE, 
                  organism = "hsapiens", ordered_query = FALSE, 
                  multi_query = FALSE, significant = TRUE, exclude_iea = TRUE, 
                  measure_underrepresentation = FALSE, evcodes = TRUE, 
                  user_threshold = 0.05, correction_method = "g_SCS", 
                  domain_scope = "custom_annotated", custom_bg = geneUniverse$V1, 
                  numeric_ns = "", sources = c("GO:MF,GO:BP","GO:CC","KEGG"), as_short_link = FALSE, highlight = TRUE)
  
  dfResults <- gostres$result
  dfResults$test <- names(listGOtests)[x]
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

allres$ctype <- sapply(strsplit(allres$test,"/"), function(x) x[2])
allres$timePoint <- as.numeric(sapply(strsplit(allres$test,"/"), function(x) x[1]))
allres$diffExpDir <- sapply(strsplit(allres$test,"/"), function(x) x[3])
allres$diffExpDir <- gsub("down","Downregulated",allres$diffExpDir)
allres$diffExpDir <- gsub("up","Upregulated", allres$diffExpDir)
allres$timePoint <- as.character(allres$timePoint)
allres$timePoint <- factor(allres$timePoint, levels=c("40","70","120"))


####################
## Main Figure 7G ##
####################

fig7G <- ggplot(data=subset(allres, ctype=="hRgl1"), aes(x=timePoint, y=term_name, fill=OddsRatio))+
  theme_bw()+
  geom_tile()+
  facet_grid(~diffExpDir)+
  xlab("Timepoint [in days]")+
  ylab("GO enrichment: Cellular component terms")+
  ggtitle("hRgl1 (DE genes: patients vs ctrl)")+
  scale_fill_viridis_c()+
  theme(plot.title=element_text(hjust=0.5, face="bold"))


pdf(file=paste0(rootMain,"mainFigure7G.pdf"), height =8, width=8)
plot(fig7G)
dev.off()


## hDA1a ##

genes_intersection <- subset(allres, ctype=="hDA1a")
genes_intersection$evidence_codes <- NULL

##40/70
genes_intersection_initial <- subset(genes_intersection, timePoint!="120")
relevantGenes_d40d70 <- unique(unlist(strsplit(genes_intersection_initial$intersection,",")))
res_initial <- resDE[resDE$geneSymbol %in% relevantGenes_d40d70,]
res_initial <- subset(res_initial, cluster=="hDA1a")
res_initial <- subset(res_initial, timePoint!="120")
res_initial$timePoint <- as.character(res_initial$timePoint)
res_initial$timePoint <- factor(res_initial$timePoint, levels=c("40","70","120"))

##120
genes_intersection_final <- subset(genes_intersection, timePoint=="120")
relevantGenes_d120 <- unique(unlist(strsplit(genes_intersection_final$intersection,",")))
res_final <- resDE[resDE$geneSymbol %in% relevantGenes_d120,]
res_final <- subset(res_final, cluster=="hDA1a")
res_final <- subset(res_final, timePoint=="120")
res_final$timePoint <- as.character(res_final$timePoint)
res_final$timePoint <- factor(res_final$timePoint, levels=c("40","70","120"))

allTP <- rbind(res_initial, res_final)
allTP$direction <- gsub("down","Downregulated",allTP$direction)
allTP$direction <- gsub("up","Upregulated", allTP$direction)



library(tidytext)

allTP <- allTP %>%
  mutate(group = tidytext::reorder_within(geneSymbol, avg_logFC, within=timePoint))


####################
## Main Figure 7H ##
####################

fig7H <- ggplot(allTP, aes(x=avg_logFC, y=group, fill=direction))+
  theme_bw()+
  ggtitle("Genes contributing to GO enrichment (hDA1a)")+
  geom_bar(position="dodge", stat="identity")+
  xlab("Average logFC in hDA1a (3D organoids) [Patient/Control]")+
  ylab("")+
  scale_fill_manual(name="DE genes",values=c("cornflowerblue","brown1"))+
  theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust = 0.5))+
  facet_wrap(vars(timePoint), scales="free")+
  tidytext::scale_y_reordered() 

pdf(file=paste0(rootMain,"mainFigure7H.pdf"), width=8, height = 4)
plot(fig7H)
dev.off()



## hRgl1 ##

genes_intersection <- subset(allres, ctype=="hRgl1")
genes_intersection$evidence_codes <- NULL

## upregulated
genes_intersection_up <- subset(genes_intersection, diffExpDir=="Upregulated")
term_consistent_up <- names(which(table(genes_intersection_up$term_name)==3))
genes_intersection_up <- genes_intersection_up[match(term_consistent_up, genes_intersection_up$term_name),]
genes_intersection_up <- unique(unlist(strsplit(genes_intersection_up$intersection,",")))
res_upreg <- resDE[resDE$geneSymbol %in% genes_intersection_up,]
res_upreg <- subset(res_upreg, cluster=="hRgl1" & direction=="up")
res_upreg$timePoint <- as.character(res_upreg$timePoint)
res_upreg$timePoint <- factor(res_upreg$timePoint, levels=c("40","70","120"))

consistent_up <- names(which(table(res_upreg$geneSymbol)==3))
res_upreg <- res_upreg[res_upreg$geneSymbol %in% consistent_up,]

## downregulated
genes_intersection_down <- subset(genes_intersection, diffExpDir=="Downregulated")
term_consistent_down <- names(which(table(genes_intersection_down$term_name)==3))
genes_intersection_down <- genes_intersection_down[match(term_consistent_down, genes_intersection_down$term_name),]
genes_intersection_down <- unique(unlist(strsplit(genes_intersection_down$intersection,",")))
res_downreg <- resDE[resDE$geneSymbol %in% genes_intersection_down,]
res_downreg <- subset(res_downreg, cluster=="hRgl1" & direction=="down")
res_downreg$timePoint <- as.character(res_downreg$timePoint)
res_downreg$timePoint <- factor(res_downreg$timePoint, levels=c("40","70","120"))

consistent_down <- names(which(table(res_downreg$geneSymbol)==3))
res_downreg <- res_downreg[res_downreg$geneSymbol %in% consistent_down,]

allTP <- rbind(res_upreg, res_downreg)
allTP$direction <- gsub("down","Downregulated",allTP$direction)
allTP$direction <- gsub("up","Upregulated", allTP$direction)


# pick legend

plotBarInitial0 <- ggplot(allTP, aes(x=avg_logFC, y=geneSymbol, fill=direction))+
  theme_bw()+
  ggtitle("Genes contributing to GO enrichment in all TP (hRgl1)")+
  geom_bar(position="dodge", stat="identity")+
  xlab("Average logFC in hRgl1 (3D organoids) [Patient/Control]")+
  ylab("")+
  scale_fill_manual(name="DE genes",values=c("cornflowerblue","brown1"))+
  theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust = 0.5))+
  facet_grid(direction~timePoint, scales="free")

legend <- get_legend(plotBarInitial0)


# row1 (upregulated)

allTP_up <- subset(allTP, direction=="Upregulated")
allTP_up <- allTP_up %>% mutate(group = tidytext::reorder_within(geneSymbol, avg_logFC, within=timePoint))


plotBarInitial <- ggplot(allTP_up, aes(x=avg_logFC, y=group, fill=direction))+
  theme_bw()+
  ggtitle("Upregulated")+
  geom_bar(position="dodge", stat="identity")+
  xlab("")+
  ylab("")+
  scale_fill_manual(name="DE genes",values=c("brown1"))+
  theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust = 0.5))+
  facet_wrap(vars(timePoint), scales="free", nrow=1)+
  tidytext::scale_y_reordered()+
  theme(legend.position="none")


# row2 (downregulated) 

allTP_down <- subset(allTP, direction=="Downregulated")
allTP_down <- allTP_down %>% mutate(group = tidytext::reorder_within(geneSymbol, avg_logFC, within=timePoint))


plotBarInitial2 <- ggplot(allTP_down, aes(x=avg_logFC, y=group, fill=direction))+
  theme_bw()+
  ggtitle("Downregulated")+
  geom_bar(position="dodge", stat="identity")+
  xlab("")+
  ylab("")+
  scale_fill_manual(name="DE genes",values=c("cornflowerblue"))+
  theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust = 0.5))+
  facet_wrap(vars(timePoint), scales="free", nrow=1)+
  tidytext::scale_y_reordered()+
  theme(legend.position="none")


####################
## Main Figure 7I ##
####################


### combine row1 and row2 panels
library(ggpubr)
combined_plot <- ggarrange(plotBarInitial2, 
                           plotBarInitial,
                           ncol = 1, heights = c(2,1), align = "hv")

# Combine the plots with the legend and add titles
final_plot <- ggarrange(combined_plot, legend, 
                        ncol = 2, widths = c(4, 1)) # Adjust width ratio as needed

# Annotate the final plot with axis titles and a main title

title <- "Genes contributing to GO enrichment in all TP (hRgl1)"
xaxis <- "Average logFC in hRgl1 (3D organoids) [Patient/Control]"
fig7I <- annotate_figure(final_plot,
			 top = text_grob(title, face = "bold", size = 14),
			 bottom = text_grob(xaxis, size = 12))

# To save the plot
ggsave(paste0(rootMain,"mainFigure7I.pdf"), fig7I, width = 10, height = 6)



#############
## hMidPre ##
#############

genes_intersection <- subset(allres, ctype=="hMidPre")
genes_intersection$evidence_codes <- NULL

##downregulated-Mito (day 120)
genes_intersection_down <- subset(genes_intersection, diffExpDir=="Downregulated" & timePoint=="120")
genes_intersection_down <- unique(unlist(strsplit(genes_intersection_down$intersection,",")))
res_downreg_g1 <- resDE[resDE$geneSymbol %in% genes_intersection_down,]
res_downreg_g1 <- subset(res_downreg_g1, cluster=="hMidPre" & direction=="down" & timePoint=="120")
res_downreg_g1$timePoint <- as.character(res_downreg_g1$timePoint)
res_downreg_g1$timePoint <- factor(res_downreg_g1$timePoint, levels=c("40","70","120"))

allTP <- res_downreg_g1
allTP$direction <- gsub("down","Downregulated",allTP$direction)
allTP$direction <- gsub("up","Upregulated", allTP$direction)

allTP$direction <- factor(allTP$direction, levels=c("Upregulated","Downregulated"))
allTP <- as.data.frame(allTP %>% group_by(timePoint) %>% arrange(-avg_logFC))

allTP$geneSymbol <- factor(allTP$geneSymbol, levels = rev(unique(allTP$geneSymbol)))


##############################
## Supplementary Figure 21D ##
##############################

figSupp21D <- ggplot(allTP, aes(x=avg_logFC, y=geneSymbol, fill=direction))+
  theme_bw()+
  ggtitle("Genes contributing to 'Mitochondrial' enrichment (hMidPre)")+
  geom_bar(position="dodge", stat="identity")+
  facet_grid(direction~timePoint, scales="free_y")+
  xlab("Average logFC in hMidPre (3D organoids) [Patient/Control]")+
  ylab("")+
  scale_fill_manual(name="DE genes",values=c("Upregulated"="brown1",
                                             "Downregulated"="cornflowerblue"))+
  theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust = 0.5),
        legend.position="top")

pdf(file=paste0(rootSupp,"suppFigure21D.pdf"), width=6, height = 10)
plot(figSupp21D)
dev.off()


##downregulated and upregulated - extracellular (day 40)
genes_intersection_up_40 <- subset(genes_intersection, diffExpDir=="Upregulated"  & timePoint=="40")
genes_intersection_up_40 <- unique(unlist(strsplit(genes_intersection_up_40$intersection,",")))
res_upreg_g1 <- resDE[resDE$geneSymbol %in% genes_intersection_up_40,]
res_upreg_g1 <- subset(res_upreg_g1, cluster=="hMidPre" & direction=="up" & timePoint=="40")
res_upreg_g1$timePoint <- as.character(res_upreg_g1$timePoint)
res_upreg_g1$timePoint <- factor(res_upreg_g1$timePoint, levels=c("40","70","120"))

genes_intersection_down_40 <- subset(genes_intersection, diffExpDir=="Downregulated"  & timePoint=="40")
genes_intersection_down_40 <- unique(unlist(strsplit(genes_intersection_down_40$intersection,",")))
res_downreg_g2 <- resDE[resDE$geneSymbol %in% genes_intersection_down_40,]
res_downreg_g2 <- subset(res_downreg_g2, cluster=="hMidPre" & direction=="down" & timePoint=="40")
res_downreg_g2$timePoint <- as.character(res_downreg_g2$timePoint)
res_downreg_g2$timePoint <- factor(res_downreg_g2$timePoint, levels=c("40","70","120"))


allTP <- rbind(res_upreg_g1, res_downreg_g2)
allTP$direction <- gsub("down","Downregulated",allTP$direction)
allTP$direction <- gsub("up","Upregulated", allTP$direction)

allTP$direction <- factor(allTP$direction, levels=c("Upregulated","Downregulated"))
allTP <- as.data.frame(allTP %>% group_by(timePoint) %>% arrange(-avg_logFC))

allTP$geneSymbol <- factor(allTP$geneSymbol, levels = rev(unique(allTP$geneSymbol)))


plotBarInitial <- ggplot(allTP, aes(x=avg_logFC, y=geneSymbol, fill=direction))+
  theme_bw()+
  ggtitle("Genes contributing to 'Extracellular' enrichment (hMidPre)")+
  geom_bar(position="dodge", stat="identity")+
  facet_grid(direction~timePoint, scales="free_y")+
  xlab("Average logFC in hMidPre (3D organoids) [Patient/Control]")+
  ylab("")+
  scale_fill_manual(name="DE genes",values=c("brown1","cornflowerblue"))+
  theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust = 0.5),
        axis.text.y=element_text(size=6),
        legend.position="top")

pdf(file=paste0(rootDir,"hMidPre_genesDE_enrich.pdf"), width=6, height = 10)
plot(plotBarInitial)
dev.off()


#########################################################################
#########################################################################
#########################################################################

library(ggupset)

geneUniverse <- read.table("saved/seurat/DEgeneUniverse.txt")
resDE <- read.table("saved/suppTables/TableS11.txt", header=TRUE)
stopifnot(all(resDE$geneSymbol %in% geneUniverse$V1))


neuroLineage <- c("hDA1a", "hDA1b", "hDA2", "hNbDA", "hPreDA")
resDE_neuroL <- resDE[resDE$cluster %in% neuroLineage,]

bothDirections <- sapply(unique(resDE_neuroL$direction), function(x){
  
  ## individual plots for either upregulated or downregulated
  tmp <- subset(resDE_neuroL, direction==x)
  
  ## facet for timepoint
  df_perTP <- sapply(unique(tmp$timePoint), function(y){
    
    tmp2 <- subset(tmp, timePoint==y)
    
    newdf <- data.frame(geneSymbol=names(split(tmp2$cluster, tmp2$geneSymbol)),
                        timePoint=rep(y, length(names(split(tmp2$cluster, tmp2$geneSymbol)))),
                        direction=x)
    
    combinations=unname(sapply(split(tmp2$cluster, tmp2$geneSymbol), function(x) sort(x)))
    newdf$combinations <- combinations
    
    return(newdf)

  }, simplify=F)
  
  df_perTP <- do.call("rbind", df_perTP)
  rownames(df_perTP) <- NULL
  
  return(df_perTP)
  
  
}, simplify=F)


bothDirections <- do.call("rbind", bothDirections)
rownames(bothDirections) <- NULL


upRegulated <- subset(bothDirections, direction=="up")
allUpset <- ggplot(upRegulated, aes(x=combinations)) +
  theme_bw()+
  geom_bar() +
  scale_x_upset()+
  facet_wrap(~timePoint, scales="free")+
  ggtitle("Intersection of upregulated genes in the neuronal lineage (Patients VS controls)")+
  xlab("")+
  ylab("Number of DE genes")+
  theme(plot.title=element_text(face="bold"))+
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1)+
  scale_y_continuous(limits=c(0,60), breaks=seq(0,60,10))

pdf(file=paste0(rootDir,"upsetAll_neuronLike_upreg.pdf"), height=4, width=12)
plot(allUpset)
dev.off()


downRegulated <- subset(bothDirections, direction=="down")
allUpset2 <- ggplot(downRegulated, aes(x=combinations)) +
  theme_bw()+
  geom_bar() +
  scale_x_upset()+
  facet_wrap(~timePoint, scales="free")+
  ggtitle("Intersection of downregulated genes in the neuronal lineage (Patients VS controls)")+
  xlab("")+
  ylab("Number of DE genes")+
  theme(plot.title=element_text(face="bold"))+
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1)+
  scale_y_continuous(limits=c(0,60), breaks=seq(0,60,10))

pdf(file=paste0(rootDir,"upsetAll_neuronLike_downreg.pdf"), height=4, width=12)
plot(allUpset2)
dev.off()



## GO enrichment test ##     

listGOtests <- list(c(40,"hDA1a","up"),c(40,"hDA1b","up"),c(40,"hDA2","up"),c(40,"hNbDA","up"),c(40,"hPreDA","up"),
                    c(40,"hDA1a","down"),c(40,"hDA1b","down"),c(40,"hDA2","down"),c(40,"hNbDA","down"),c(40,"hPreDA","down"),
                    c(70,"hDA1a","up"),c(70,"hDA1b","up"),c(70,"hDA2","up"),c(70,"hNbDA","up"),c(70,"hPreDA","up"),
                    c(70,"hDA1a","down"),c(70,"hDA1b","down"),c(70,"hDA2","down"),c(70,"hNbDA","down"),c(70,"hPreDA","down"),
                    c(120,"hDA1a","up"),c(120,"hDA1b","up"),c(120,"hDA2","up"),c(120,"hNbDA","up"),c(120,"hPreDA","up"),
                    c(120,"hDA1a","down"),c(120,"hDA1b","down"),c(120,"hDA2","down"),c(120,"hNbDA","down"),c(120,"hPreDA","down"))


names(listGOtests) <- sapply(listGOtests, function(x) paste(x, collapse="/"))

resultsGO <- sapply(1:length(listGOtests), function(x){
  
  geneList_DE <- subset(resDE, timePoint==listGOtests[[x]][1] & cluster==listGOtests[[x]][2] & direction==listGOtests[[x]][3])$geneSymbol
  
  gostres <- gost(query = geneList_DE, 
                  organism = "hsapiens", ordered_query = FALSE, 
                  multi_query = FALSE, significant = TRUE, exclude_iea = TRUE, 
                  measure_underrepresentation = FALSE, evcodes = TRUE, 
                  user_threshold = 0.05, correction_method = "g_SCS", 
                  domain_scope = "custom_annotated", custom_bg = geneUniverse$V1, 
                  numeric_ns = "", sources = c("GO:MF,GO:BP","GO:CC","KEGG"), as_short_link = FALSE, highlight = TRUE)
  
  dfResults <- gostres$result
  dfResults$test <- names(listGOtests)[x]
  return(dfResults)
  
}, simplify=F)


resultsGO_clean <- resultsGO[sapply(resultsGO, function(x) is.data.frame(x))]
resultsGO_clean <- sapply(resultsGO_clean, function(x) {x[x$source=="GO:CC",] }, simplify=F )
resultsGO_clean <- resultsGO_clean[sapply(resultsGO_clean, function(x) dim(x)[1]>0)]


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

allres$ctype <- sapply(strsplit(allres$test,"/"), function(x) x[2])
allres$timePoint <- as.numeric(sapply(strsplit(allres$test,"/"), function(x) x[1]))
allres$diffExpDir <- sapply(strsplit(allres$test,"/"), function(x) x[3])
allres$diffExpDir <- gsub("down","Downregulated",allres$diffExpDir)
allres$diffExpDir <- gsub("up","Upregulated", allres$diffExpDir)

allres$timePoint <- as.character(allres$timePoint)
allres$timePoint <- factor(allres$timePoint, levels=c("40","70","120"))



####################
## Main Figure 7F ##
####################

fig7F <- ggplot(data=allres, aes(x=timePoint, y=term_name, fill=OddsRatio))+
  theme_bw()+
  geom_tile()+
  facet_grid(diffExpDir~ctype, scales="free_y")+
  xlab("Timepoint [in days]")+
  ylab("GO enrichment: Cellular component terms")+
  ggtitle("Neuronal lineage (DE genes: patients vs ctrl)")+
  scale_fill_viridis_c()+
  theme(plot.title=element_text(hjust=0.5, face="bold"))


pdf(file=paste0(rootMain,"mainFigure7F.pdf"), width=8, height=8)
plot(fig7F)
dev.off()


