library(ggsignif)
library(ggpubr)
library(Seurat)
library(monocle3)
library(RColorBrewer)
library(ggthemes)
library(ggrastr)
library(base64enc)
library(ggplot2) 
library(readxl)
library(Biobase)
library(ggbeeswarm)
library(cowplot)
library(stringr)
library(ggridges)
library(tidyverse)
library(stringr)
library(ggsignif)
library(ggpubr)

rootMain <- "figures/main/"
rootSupp <- "figures/supp/"
rootDir <- "otherPlots/monocle3/"


setLast <- function(ctypesVec, lastCtype="Unk"){

  idLast <- match(lastCtype, ctypesVec)
  ordVec <- c(sort(ctypesVec[-idLast]), ctypesVec[idLast])
  return(ordVec)
}

mn.obj <- readRDS("saved/monocle3/monocle_integratedlogNormCountsAllGenes.RDS")
gene_module_df <- readRDS("saved/monocle3/monocle_geneModules_integratedlogNormCountsAllGenes.RDS")


module6 <- subset(gene_module_df, module==6)$id
module12 <- subset(gene_module_df, module==12)$id
module18 <- subset(gene_module_df, module==18)$id
module21 <- subset(gene_module_df, module==21)$id
module46 <- subset(gene_module_df, module==46)$id



modules_list <- list(module6, module12, module18, module21, module46)
names(modules_list) <- c("Module 6","Module 12","Module 18","Module 21", "Module 46")

geneUniveirse <- read.table("saved/seurat/DEgeneUniverse.txt")

modules_list <- sapply(modules_list, function(y){
  y <- y[y %in% geneUniverse$V1]
  all(y %in% geneUniverse$V1)
  return(y)
}, simplify=T)

resultsGO <- sapply(1:length(modules_list), function(x){

  geneList_DE <- modules_list[[x]]

  gostres <- gost(query = geneList_DE,
                  organism = "hsapiens", ordered_query = FALSE,
                  multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
                  measure_underrepresentation = FALSE, evcodes = TRUE,
                  user_threshold = 0.05, correction_method = "g_SCS",
                  domain_scope = "custom_annotated", custom_bg = geneUniverse$V1,
                  numeric_ns = "", sources = c("GO:MF,GO:BP","GO:CC","KEGG"), as_short_link = FALSE, highlight = TRUE)

  dfResults <- gostres$result
  dfResults$test <- names(modules_list)[x]
  return(dfResults)

}, simplify=F)



resultsGO_clean <- resultsGO[sapply(resultsGO, function(x) is.data.frame(x))]
resultsGO_clean <- sapply(resultsGO_clean, function(x) {x[x$highlighted==TRUE,] }, simplify=F )


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


#############################
## Supplemental Figure 11C ##
#############################


allres$short <- allres$term_name
mask_long<- nchar(allres$term_name)>30
allres$short[mask_long] <- sapply(allres[mask_long,]$term_name, function(x){

  mask <- unlist(sapply(strsplit(x, " "), function(y) cumsum(nchar(y))>30, simplify=F))
  #sapply(strsplit(x, " "), function(y) paste0(paste0(y[!mask], collapse=" "),"\n", y[mask]))
  sapply(strsplit(x, " "), function(y) paste0(paste0(y[!mask], collapse=" "),"\n",paste0(y[mask], collapse=" ")))

})
allres <- allres[order(allres$OddsRatio, decreasing=F),]

allres$short <- factor(allres$short, levels=unique(allres$short))

#allres$short <- factor(allres$short, levels=allres$short)

allres$test <- factor(allres$test, levels=c("Module 6","Module 12","Module 18"))


library(tidytext)
allres <- allres %>%
  mutate(group = tidytext::reorder_within(short, OddsRatio, within=test))


suppFig11C <- ggplot(allres, aes(x=OddsRatio, y=group, size=OddsRatio, col=p_value))+
  geom_point()+
  theme_bw()+
  scale_colour_gradient(low = "red", high = "blue", na.value = NA)+
  xlab("Odds Ratio")+
  ylab("GO enrichment [Cellular component]")+
  theme(axis.text.y=element_text(size=10),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(hjust=0.5, size=14, face="bold"))+
  tidytext::scale_y_reordered() +
  facet_wrap(vars(test), scales = "free_y", ncol=1)


pdf(file=paste0(rootSupp,"suppFigure11C.pdf"), width=12, height = 8)
plot(suppFig11C)
dev.off()




library(facefuns)
library(ggpubr)
library(dplyr)
library(coin)

## Compute module score for each cell

agg_matrix_perCell <- aggregate_gene_expression(mn.obj, gene_module_df)
row.names(agg_matrix_perCell) <- stringr::str_c("Module ", row.names(agg_matrix_perCell))

agg_matrix_perCell_sub <- agg_matrix_perCell[row.names(agg_matrix_perCell) %in% c("Module 6","Module 12", "Module 18", "Module 21","Module 46"),]
agg_matrix_perCell_sub <- as.data.frame(agg_matrix_perCell_sub)
agg_matrix_perCell_sub$Module <- row.names(agg_matrix_perCell_sub)
row.names(agg_matrix_perCell_sub) <- NULL


#####
#####


library(facefuns)
library(ggpubr)
library(dplyr)
library(coin)

## Compute module score for each cell

agg_matrix_perCell <- aggregate_gene_expression(mn.obj, gene_module_df)
row.names(agg_matrix_perCell) <- stringr::str_c("Module ", row.names(agg_matrix_perCell))

agg_matrix_perCell_sub <- agg_matrix_perCell[row.names(agg_matrix_perCell) %in% c("Module 6","Module 12", "Module 18", "Module 21","Module 46"),]
agg_matrix_perCell_sub <- as.data.frame(agg_matrix_perCell_sub)
agg_matrix_perCell_sub$Module <- row.names(agg_matrix_perCell_sub)
row.names(agg_matrix_perCell_sub) <- NULL


### Compute modules score for each cell in the corresponding groups

###################################
## foetal-ours vs foetal-LaManno ##
###################################

maskFoetalOursLaManno<- (colData(mn.obj)$Dataset2=="This work (2D, Organoids, Foetal)" & colData(mn.obj)$system=="Foetal") |
  colData(mn.obj)$Dataset2=="La Manno et al. 2016 (Foetal)"
mn.obj_sub <- mn.obj[,maskFoetalOursLaManno]

agg_matrix_perCell_long <- as.data.frame(agg_matrix_perCell_sub %>% pivot_longer(-c("Module"), names_to="cellId", values_to="ActivationScores"))
agg_matrix_perCell_long <- agg_matrix_perCell_long[agg_matrix_perCell_long$cellId %in% colnames(mn.obj_sub),]
agg_matrix_perCell_long <- cbind(agg_matrix_perCell_long, colData(mn.obj_sub)[match(agg_matrix_perCell_long$cellId, rownames(colData(mn.obj_sub))),])

agg_matrix_perCell_long$Dataset2 <- gsub("2D, Organoids, ","", agg_matrix_perCell_long$Dataset2)
agg_matrix_perCell_long$Dataset2 <- as.character(agg_matrix_perCell_long$Dataset2)
agg_matrix_perCell_long$Dataset2 <- factor(agg_matrix_perCell_long$Dataset2, levels=c("This work (Foetal)",
                                                                                      "La Manno et al. 2016 (Foetal)"))

agg_matrix_perCell_long_sub <- subset(agg_matrix_perCell_long, annotation_unified=="Neurons")

colVec_foetal_work_LaManno <- colVec
names(colVec_foetal_work_LaManno) <- gsub("2D, 3D, ","", names(colVec_foetal_work_LaManno))
colVec_foetal_work_LaManno <- colVec_foetal_work_LaManno[names(colVec_foetal_work_LaManno) %in% agg_matrix_perCell_long$Dataset2]

agg_matrix_perCell_long_sub$Module <- gsub("Module ","Mod",agg_matrix_perCell_long_sub$Module)
agg_matrix_perCell_long_sub$Module <- factor(agg_matrix_perCell_long_sub$Module,
                                             levels=unique(agg_matrix_perCell_long_sub$Module)[order(as.numeric(gsub("Mod","",unique(agg_matrix_perCell_long_sub$Module))))])


# Function to perform permutation test (only, no re-sampling due to differences on sample size between the two groups of comparison)

set.seed(123)
permutation_test <- function(data, group_var, response_var) {
  test_result <- oneway_test(as.formula(paste(response_var, "~", group_var)), data = data,
                             distribution = approximate(B = 10000))
  p_value <- pvalue(test_result)
  return(p_value)
}

# Function to perform permutation test with resampling (from the smaller group)
permutation_test_resample <- function(data, group_var, response_var, n_resamples = 100) {
  group_levels <- unique(data[[group_var]])
  if(length(group_levels) != 2) {
    stop("The group variable must have exactly 2 levels.")
  }

  group1 <- data[data[[group_var]] == group_levels[1], ]
  group2 <- data[data[[group_var]] == group_levels[2], ]

  n1 <- nrow(group1)
  n2 <- nrow(group2)

  resampled_p_values <- numeric(n_resamples)

  for (i in 1:n_resamples) {

    #print(i)

    if (n1 > n2) {
      group1_sample <- group1[sample(n1, n2, replace = FALSE), ]
      resampled_data <- rbind(group1_sample, group2)
    } else {
      group2_sample <- group2[sample(n2, n1, replace = FALSE), ]
      resampled_data <- rbind(group1, group2_sample)
    }
    test_result <- oneway_test(as.formula(paste(response_var, "~", group_var)),
                               data = resampled_data,
                               distribution = approximate(nresample = 10000))
    resampled_p_values[i] <- pvalue(test_result)
  }

  p_value <- mean(resampled_p_values)
  return(p_value)
}


# Perform permutation tests for each Module
set.seed(123)
perm_test_results <- agg_matrix_perCell_long_sub %>%
  group_by(Module) %>%
  summarise(p_value = permutation_test_resample(cur_data(), "Dataset2", "ActivationScores"))

# Adjust p-values using Bonferroni correction
perm_test_results <- perm_test_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Add significance labels
perm_test_results <- perm_test_results %>%
  mutate(p_signif = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01 ~ "**",
    p_adj < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

print(perm_test_results)

saveRDS(perm_test_results, "saved/permutationsDistr/perm_test_results_neuronalModule_distr_ours_LaManno.RDS")


#############################
## Supplemental Figure 11D ##
#############################


suppFig11D <- ggplot(agg_matrix_perCell_long_sub, aes(x = Module, y = ActivationScores, fill = Dataset2)) +
  theme_bw()+
  geom_split_violin(alpha = .4, trim = FALSE) +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape=NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F,
               position = position_dodge(.175)) +
  xlab("")+
  #scale_x_discrete(name = "", labels = c("Module 6","Module 12", "Module 18","","")) +
  scale_y_continuous(name = "Activation Score",
                     breaks = seq(-3, 4, 1),
                     limits = c(-3, 4)) +
  scale_fill_manual(values=colVec_foetal_work_LaManno, name = "") +
  theme(legend.position="top",
        axis.title=element_text(size=14),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=11))

for(i in 1:nrow(perm_test_results)) {
  suppFig11D <- suppFig11D + annotate("text",
                        x = perm_test_results$Module[i],
                        y = 4, # Adjust the y position based on your data
                        label = perm_test_results$p_signif[i],
                        size = 5,
                        color = "black")
}

pdf(file=paste0(rootSupp,"suppFigure11D.pdf"), height = 4.5, width=4.5)
plot(suppFig11D)
dev.off()



###################################
## foetal-ours vs foetal-Birtele ##
###################################

maskFoetalOursBirtele<- (colData(mn.obj)$Dataset2=="This work (2D, Organoids, Foetal)" & colData(mn.obj)$system=="Foetal") |
  colData(mn.obj)$Dataset2=="Birtele et al. 2022 (Foetal)"
mn.obj_sub <- mn.obj[,maskFoetalOursBirtele]


agg_matrix_perCell_long <- as.data.frame(agg_matrix_perCell_sub %>% pivot_longer(-c("Module"), names_to="cellId", values_to="ActivationScores"))
agg_matrix_perCell_long <- agg_matrix_perCell_long[agg_matrix_perCell_long$cellId %in% colnames(mn.obj_sub),]
agg_matrix_perCell_long <- cbind(agg_matrix_perCell_long, colData(mn.obj_sub)[match(agg_matrix_perCell_long$cellId, rownames(colData(mn.obj_sub))),])

agg_matrix_perCell_long$Dataset2 <- gsub("2D, Organoids, ","", agg_matrix_perCell_long$Dataset2)
agg_matrix_perCell_long$Dataset2 <- as.character(agg_matrix_perCell_long$Dataset2)
agg_matrix_perCell_long$Dataset2 <- factor(agg_matrix_perCell_long$Dataset2, levels=c("This work (Foetal)",
                                                                                      "Birtele et al. 2022 (Foetal)"))

agg_matrix_perCell_long_sub <- subset(agg_matrix_perCell_long, annotation_unified=="Neurons")

colVec_foetal_work_Birtele <- colVec
names(colVec_foetal_work_Birtele) <- gsub("2D, 3D, ","", names(colVec_foetal_work_Birtele))
colVec_foetal_work_Birtele <- colVec_foetal_work_Birtele[names(colVec_foetal_work_Birtele) %in% agg_matrix_perCell_long$Dataset2]

agg_matrix_perCell_long_sub$Module <- gsub("Module ","Mod",agg_matrix_perCell_long_sub$Module)
agg_matrix_perCell_long_sub$Module <- factor(agg_matrix_perCell_long_sub$Module,
                                             levels=unique(agg_matrix_perCell_long_sub$Module)[order(as.numeric(gsub("Mod","",unique(agg_matrix_perCell_long_sub$Module))))])

# Custom function for summarising with print statement
custom_summarise <- function(data, current_module, group_var, response_var) {
  # Print the current module
  print(paste("Processing module:", current_module))

  # Perform the permutation test with resampling
  p_value <- permutation_test_resample(data, group_var, response_var)

  return(p_value)
}

# Perform permutation tests with custom summarise function
set.seed(123)
perm_test_results <- agg_matrix_perCell_long_sub %>%
  group_by(Module) %>%
  summarise(p_value = custom_summarise(cur_data(), unique(Module), "Dataset2", "ActivationScores"))

# Adjust p-values using Bonferroni correction
perm_test_results <- perm_test_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Add significance labels
perm_test_results <- perm_test_results %>%
  mutate(p_signif = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01 ~ "**",
    p_adj < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

print(perm_test_results)

saveRDS(perm_test_results, "saved/permutationsDistr/perm_test_results_neuronalModule_distr_ours_Birtele.RDS")


####################
## Main Figure 4I ##
####################

fig4I <- ggplot(agg_matrix_perCell_long_sub, aes(x = Module, y = ActivationScores, fill = Dataset2)) +
  theme_bw()+
  geom_split_violin(alpha = .4, trim = FALSE) +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape=NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F,
               position = position_dodge(.175)) +
  xlab("")+
  #scale_x_discrete(name = "", labels = c("Module 6","Module 12", "Module 18","","")) +
  scale_y_continuous(name = "Activation Score",
                     breaks = seq(-3, 4, 1),
                     limits = c(-3, 4)) +
  scale_fill_manual(values=colVec_foetal_work_Birtele, name = "") +
  theme(legend.position="top",
        axis.title=element_text(size=14),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=11))

for(i in 1:nrow(perm_test_results)) {
  fig4I <- fig4I + annotate("text",
                        x = perm_test_results$Module[i],
                        y = 4, # Adjust the y position based on your data
                        label = perm_test_results$p_signif[i],
                        size = 5,
                        color = "black")
}

pdf(file=paste0(rootMain,"mainFigure4I.pdf"), height = 4.5, width=4.5)
plot(fig4I)
dev.off()



#################################
## foetal-ours vs foetal-Braun ##
#################################

maskFoetalOursBraun<- (colData(mn.obj)$Dataset2=="This work (2D, Organoids, Foetal)" & colData(mn.obj)$system=="Foetal") |
  colData(mn.obj)$Dataset2=="Braun et al. 2023 (Foetal)"
mn.obj_sub <- mn.obj[,maskFoetalOursBraun]

agg_matrix_perCell_long <- as.data.frame(agg_matrix_perCell_sub %>% pivot_longer(-c("Module"), names_to="cellId", values_to="ActivationScores"))
agg_matrix_perCell_long <- agg_matrix_perCell_long[agg_matrix_perCell_long$cellId %in% colnames(mn.obj_sub),]
agg_matrix_perCell_long <- cbind(agg_matrix_perCell_long, colData(mn.obj_sub)[match(agg_matrix_perCell_long$cellId, rownames(colData(mn.obj_sub))),])

agg_matrix_perCell_long$Dataset2 <- gsub("2D, Organoids, ","", agg_matrix_perCell_long$Dataset2)
agg_matrix_perCell_long$Dataset2 <- as.character(agg_matrix_perCell_long$Dataset2)
agg_matrix_perCell_long$Dataset2 <- factor(agg_matrix_perCell_long$Dataset2, levels=c("This work (Foetal)",
                                                                                      "Braun et al. 2023 (Foetal)"))

agg_matrix_perCell_long_sub <- subset(agg_matrix_perCell_long, annotation_unified=="Neurons")

colVec_foetal_work_Braun <- colVec
names(colVec_foetal_work_Braun) <- gsub("2D, 3D, ","", names(colVec_foetal_work_Braun))
colVec_foetal_work_Braun <- colVec_foetal_work_Braun[names(colVec_foetal_work_Braun) %in% agg_matrix_perCell_long$Dataset2]

agg_matrix_perCell_long_sub$Module <- gsub("Module ","Mod",agg_matrix_perCell_long_sub$Module)
agg_matrix_perCell_long_sub$Module <- factor(agg_matrix_perCell_long_sub$Module,
                                             levels=unique(agg_matrix_perCell_long_sub$Module)[order(as.numeric(gsub("Mod","",unique(agg_matrix_perCell_long_sub$Module))))])

# Custom function for summarising with print statement
custom_summarise <- function(data, current_module, group_var, response_var) {
  # Print the current module
  print(paste("Processing module:", current_module))

  # Perform the permutation test with resampling
  p_value <- permutation_test_resample(data, group_var, response_var)

  return(p_value)
}

# Perform permutation tests with custom summarise function
set.seed(123)

perm_test_results <- agg_matrix_perCell_long_sub %>%
  group_by(Module) %>%
  summarise(p_value = custom_summarise(cur_data(), unique(Module), "Dataset2", "ActivationScores"))

# Adjust p-values using Bonferroni correction
perm_test_results <- perm_test_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Add significance labels
perm_test_results <- perm_test_results %>%
  mutate(p_signif = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01 ~ "**",
    p_adj < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

print(perm_test_results)
saveRDS(perm_test_results, "saved/permutationsDistr/perm_test_results_neuronalModule_distr_ours_Braun.RDS")

####################
## Main Figure 4J ##
####################

fig4J <- ggplot(agg_matrix_perCell_long_sub, aes(x = Module, y = ActivationScores, fill = Dataset2)) +
  theme_bw()+
  geom_split_violin(alpha = .4, trim = FALSE) +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape=NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F,
               position = position_dodge(.175)) +
  xlab("")+
  #scale_x_discrete(name = "", labels = c("Module 6","Module 12", "Module 18","","")) +
  scale_y_continuous(name = "Activation Score",
                     breaks = seq(-3, 4, 1),
                     limits = c(-3, 4)) +
  scale_fill_manual(values=colVec_foetal_work_Braun, name = "") +
  theme(legend.position="top",
        axis.title=element_text(size=14),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=11))

for(i in 1:nrow(perm_test_results)) {
  fig4J <- fig4J + annotate("text",
                        x = perm_test_results$Module[i],
                        y = 4, # Adjust the y position based on your data
                        label = perm_test_results$p_signif[i],
                        size = 5,
                        color = "black")
}

pdf(file=paste0(rootMain,"mainFigure4J.pdf"), height = 4.5, width=4.5)
plot(fig4J)
dev.off()





###############################
## 2d/3d-ours vs foetal-ours ##
###############################

maskInVitro_Foetal<- (colData(mn.obj)$Dataset=="MLO" & colData(mn.obj)$Dataset!="Foetal") | (colData(mn.obj)$Dataset=="MLO" & colData(mn.obj)$system=="Foetal")
mn.obj_sub <- mn.obj[,maskInVitro_Foetal]
mn.obj_sub$Dataset3 <- "This work (in vitro)"
maskFoetal <- colData(mn.obj_sub)$system=="Foetal"
mn.obj_sub$Dataset3[maskFoetal] <- "This work (foetal)"

agg_matrix_perCell_long <- as.data.frame(agg_matrix_perCell_sub %>% pivot_longer(-c("Module"), names_to="cellId", values_to="ActivationScores"))
agg_matrix_perCell_long <- agg_matrix_perCell_long[agg_matrix_perCell_long$cellId %in% colnames(mn.obj_sub),]
agg_matrix_perCell_long <- cbind(agg_matrix_perCell_long, colData(mn.obj_sub)[match(agg_matrix_perCell_long$cellId, rownames(colData(mn.obj_sub))),])

agg_matrix_perCell_long$Dataset3 <- as.character(agg_matrix_perCell_long$Dataset3)
agg_matrix_perCell_long$Dataset3 <- factor(agg_matrix_perCell_long$Dataset3, levels=c("This work (in vitro)",
                                                                                      "This work (foetal)"))

agg_matrix_perCell_long_sub <- subset(agg_matrix_perCell_long, annotation_unified=="Neurons")


agg_matrix_perCell_long_sub$Module <- gsub("Module ","Mod",agg_matrix_perCell_long_sub$Module)
agg_matrix_perCell_long_sub$Module <- factor(agg_matrix_perCell_long_sub$Module,
                                             levels=unique(agg_matrix_perCell_long_sub$Module)[order(as.numeric(gsub("Mod","",unique(agg_matrix_perCell_long_sub$Module))))])

# Custom function for summarising with print statement
custom_summarise <- function(data, current_module, group_var, response_var) {
  # Print the current module
  print(paste("Processing module:", current_module))

  # Perform the permutation test with resampling
  p_value <- permutation_test_resample(data, group_var, response_var)

  return(p_value)
}

# Perform permutation tests with custom summarise function
set.seed(123)
perm_test_results <- agg_matrix_perCell_long_sub %>%
  group_by(Module) %>%
  summarise(p_value = custom_summarise(cur_data(), unique(Module), "Dataset3", "ActivationScores"))


# Adjust p-values using Bonferroni correction
perm_test_results <- perm_test_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))



# Add significance labels
perm_test_results <- perm_test_results %>%
  mutate(p_signif = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01 ~ "**",
    p_adj < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

print(perm_test_results)
saveRDS(perm_test_results, "saved/permutationsDistr/perm_test_results_oursInVitro_oursFoetal.RDS")



#############################
## Supplemental Figure 11E ##
#############################

library(viridis)

suppFig11E <- ggplot(agg_matrix_perCell_long_sub, aes(x = Module, y = ActivationScores, fill = Dataset3)) +
  theme_bw()+
  geom_split_violin(alpha = .4, trim = FALSE) +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape=NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F,
               position = position_dodge(.175)) +
  xlab("")+
  #scale_x_discrete(name = "", labels = c("Module 6","Module 12", "Module 18","","")) +
  scale_y_continuous(name = "Activation Score",
                     breaks = seq(-3, 4, 1),
                     limits = c(-3, 4)) +
  scale_fill_viridis_d()+
  theme(legend.position="top",
        axis.title=element_text(size=14),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=11),
        legend.title=element_blank())

for(i in 1:nrow(perm_test_results)) {
  suppFig11E <- suppFig11E + annotate("text",
                        x = perm_test_results$Module[i],
                        y = 4, # Adjust the y position based on your data
                        label = perm_test_results$p_signif[i],
                        size = 5,
                        color = "black")
}

pdf(file=paste0(rootSupp,"suppFigure11E.pdf"), height = 4.5, width=4.5)
plot(suppFig11E)
dev.off()



#################################
## 2d/3d-ours vs 3d-fiorenzano ##
#################################

maskInVitroOurs_Fiorenzano<- (colData(mn.obj)$Dataset=="MLO" & colData(mn.obj)$Dataset!="Foetal") | colData(mn.obj)$Dataset2=="Fiorenzano et al. 2021 (Organoids)"
mn.obj_sub <- mn.obj[,maskInVitroOurs_Fiorenzano]

agg_matrix_perCell_long <- as.data.frame(agg_matrix_perCell_sub %>% pivot_longer(-c("Module"), names_to="cellId", values_to="ActivationScores"))
agg_matrix_perCell_long <- agg_matrix_perCell_long[agg_matrix_perCell_long$cellId %in% colnames(mn.obj_sub),]
agg_matrix_perCell_long <- cbind(agg_matrix_perCell_long, colData(mn.obj_sub)[match(agg_matrix_perCell_long$cellId, rownames(colData(mn.obj_sub))),])

agg_matrix_perCell_long$Dataset2 <- gsub("2D, Organoids, Foetal","in vitro", agg_matrix_perCell_long$Dataset2)
agg_matrix_perCell_long$Dataset2 <- gsub("Organoids","3D", agg_matrix_perCell_long$Dataset2)

agg_matrix_perCell_long$Dataset2 <- as.character(agg_matrix_perCell_long$Dataset2)
agg_matrix_perCell_long$Dataset2 <- factor(agg_matrix_perCell_long$Dataset2, levels=c("This work (in vitro)",
                                                                                      "Fiorenzano et al. 2021 (3D)"))

agg_matrix_perCell_long_sub <- subset(agg_matrix_perCell_long, annotation_unified=="Neurons")

agg_matrix_perCell_long_sub$Module <- gsub("Module ","Mod",agg_matrix_perCell_long_sub$Module)
agg_matrix_perCell_long_sub$Module <- factor(agg_matrix_perCell_long_sub$Module,
                                             levels=unique(agg_matrix_perCell_long_sub$Module)[order(as.numeric(gsub("Mod","",unique(agg_matrix_perCell_long_sub$Module))))])

# Custom function for summarising with print statement
custom_summarise <- function(data, current_module, group_var, response_var) {
  # Print the current module
  print(paste("Processing module:", current_module))

  # Perform the permutation test with resampling
  p_value <- permutation_test_resample(data, group_var, response_var)

  return(p_value)
}

# Perform permutation tests with custom summarise function
set.seed(123)
perm_test_results <- agg_matrix_perCell_long_sub %>%
  group_by(Module) %>%
  summarise(p_value = custom_summarise(cur_data(), unique(Module), "Dataset2", "ActivationScores"))

# Adjust p-values using Bonferroni correction
perm_test_results <- perm_test_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Add significance labels
perm_test_results <- perm_test_results %>%
  mutate(p_signif = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01 ~ "**",
    p_adj < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

print(perm_test_results)
saveRDS(perm_test_results, "saved/permutationsDistr/perm_test_results_oursInVitro_Fiorenzano3D.RDS")



#############################
## Supplemental Figure 11F ##
#############################

library(viridis)

suppFig11F <- ggplot(agg_matrix_perCell_long_sub, aes(x = Module, y = ActivationScores, fill = Dataset2)) +
  theme_bw()+
  geom_split_violin(alpha = .4, trim = FALSE) +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape=NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F,
               position = position_dodge(.175)) +
  xlab("")+
  #scale_x_discrete(name = "", labels = c("Module 6","Module 12", "Module 18","","")) +
  scale_y_continuous(name = "Activation Score",
                     breaks = seq(-3, 4, 1),
                     limits = c(-3, 4)) +
  scale_fill_grey()+
  theme(legend.position="top",
        axis.title=element_text(size=14),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=11),
        legend.title=element_blank())

for(i in 1:nrow(perm_test_results)) {
  suppFig11F <- suppFig11F + annotate("text",
                        x = perm_test_results$Module[i],
                        y = 4, # Adjust the y position based on your data
                        label = perm_test_results$p_signif[i],
                        size = 5,
                        color = "black")
}

pdf(file=paste0(rootSupp,"suppFigure11F.pdf"), height = 4.5, width=4.5)
plot(suppFig11F)
dev.off()



##################################
## module specificity per ctype ##
##################################


#############################
## Supplemental Figure 11B ##
#############################



agg_matrix_perCell_sub2 <- agg_matrix_perCell[row.names(agg_matrix_perCell) %in% c("Module 1","Module 6","Module 12", "Module 18", "Module 21","Module 46"),]
agg_matrix_perCell_sub2 <- as.data.frame(agg_matrix_perCell_sub2)
agg_matrix_perCell_sub2$Module <- row.names(agg_matrix_perCell_sub2)
row.names(agg_matrix_perCell_sub2) <- NULL

agg_matrix_perCell_long2 <- as.data.frame(agg_matrix_perCell_sub2 %>% pivot_longer(-c("Module"), names_to="cellId", values_to="ActivationScores"))
agg_matrix_perCell_long2 <- agg_matrix_perCell_long2[agg_matrix_perCell_long2$cellId %in% colnames(mn.obj),]
agg_matrix_perCell_long2 <- cbind(agg_matrix_perCell_long2, colData(mn.obj)[match(agg_matrix_perCell_long2$cellId, rownames(colData(mn.obj))),])


colVec_ctypes <- readRDS(file="saved/colVec_ctypesColours.RDS")
colVec_ctypes[!names(colVec_ctypes) %in% c("Neurons","Neuroblasts","Progenitors","Precursors","IPC")] <- "grey"

agg_matrix_perCell_long2$Module <- factor(agg_matrix_perCell_long2$Module, levels=unique(agg_matrix_perCell_long2$Module)[order(as.numeric(gsub("Module ", "",unique(agg_matrix_perCell_long2$Module))))])
agg_matrix_perCell_long2$annotation_unified <- factor(agg_matrix_perCell_long2$annotation_unified, levels=rev(names(colVec_ctypes)))

suppFigure11B <- ggplot(agg_matrix_perCell_long2, aes(x = ActivationScores, y = annotation_unified, fill = annotation_unified)) +
  theme_bw()+
  ggridges::geom_density_ridges(scale = 2, show.legend = FALSE, alpha=0.9) +
  facet_wrap(~Module)+
  scale_x_continuous(name = "Activation Scores",
                      limits = c(-2, 4)) +
  scale_fill_manual(values=colVec_ctypes)+
  ylab("")


pdf(file=paste0(rootSupp,"suppFigure11B.pdf"), height = 6, width=6)
plot(suppFigure11B)
dev.off()


#########
#########

## now just showing all modules for neurons only

agg_matrix_perCell_neurons <- as.data.frame(agg_matrix_perCell)
agg_matrix_perCell_neurons$Module <- rownames(agg_matrix_perCell_neurons)
rownames(agg_matrix_perCell_neurons) <- NULL

ids_neurons <- rownames(colData(mn.obj))[colData(mn.obj)$annotation_unified=="Neurons"]
to_retain <- colnames(agg_matrix_perCell_neurons) %in% ids_neurons
to_retain[length(to_retain)] <- TRUE

agg_matrix_perCell_neurons <- agg_matrix_perCell_neurons[,to_retain]

agg_matrix_perCell_neurons_long <- as.data.frame(agg_matrix_perCell_neurons %>% pivot_longer(-c("Module"), names_to="cellId", values_to="ActivationScores"))
agg_matrix_perCell_neurons_long <- agg_matrix_perCell_neurons_long[agg_matrix_perCell_neurons_long$cellId %in% colnames(mn.obj),]
agg_matrix_perCell_neurons_long <- cbind(agg_matrix_perCell_neurons_long, colData(mn.obj)[match(agg_matrix_perCell_neurons_long$cellId, rownames(colData(mn.obj))),])

colVec_Modules <- rep("grey", length(unique(agg_matrix_perCell_neurons_long$Module)))
names(colVec_Modules) <- paste0("Module ", 1:50)

colVec_Modules[names(colVec_Modules) %in% c("Module 1","Module 6","Module 12", "Module 18", "Module 21","Module 46")] <- "red"

agg_matrix_perCell_neurons_long$Module <- factor(agg_matrix_perCell_neurons_long$Module , levels=rev(paste0("Module ", 1:50)))

ctypeActScoresNeuronsOnlyPlot <- ggplot(agg_matrix_perCell_neurons_long, aes(x = ActivationScores, y = Module, fill = Module)) +
  theme_bw()+
  ggridges::geom_density_ridges(scale = 2, show.legend = FALSE, alpha=0.9) +
  scale_x_continuous(name = "Activation Scores",
                      limits = c(-3, 3.5)) +
  scale_fill_manual(values=colVec_Modules)+
  ylab("")


pdf(file=paste0(rootDir,"ctypeActScoresDistrOnlyNeurons.pdf"), height = 6, width=4)
plot(ctypeActScoresNeuronsOnlyPlot)
dev.off()


## compute average and sd

standard_dev <- function(x) {
  return(sd(x))
}

# Compute mean and standard dev for each module
summary_stats <- agg_matrix_perCell_neurons_long %>%
  group_by(Module) %>%
  summarise(
    mean_activation = mean(ActivationScores, na.rm = TRUE),
    sd_activation = standard_dev(ActivationScores)
  )

summary_stats <- as.data.frame(summary_stats)


## pick comparisons of interest

p <-ggplot(summary_stats, aes(x=mean_activation, y=Module, col=Module)) +
  theme_bw()+
  geom_point()+
  geom_errorbar(aes(xmin=mean_activation-sd_activation, xmax=mean_activation+sd_activation))+
  scale_color_manual(values=colVec_Modules)+
  theme(legend.position="none")+
  labs(y="",
       x=expression(paste("Mean Activation Score ", "\u00B1", " SD")))+

pdf(file=paste0(rootDir,"ctypeActScoresDistrOnlyNeurons_errorbar.pdf"), height = 6, width=4)
plot(p)
dev.off()


###########
###########


## produce plot per Module 

agg_matrix_perCell_neurons_long_sub <- agg_matrix_perCell_neurons_long[agg_matrix_perCell_neurons_long$Module %in% c("Module 6","Module 12", "Module 18", "Module 21","Module 46"),]
agg_matrix_perCell_neurons_long_sub$Module <- droplevels(agg_matrix_perCell_neurons_long_sub$Module)


## Module 46 ##

#############################
## Supplemental Figure 11G ##
#############################

modName <- "Module 46"
colVec_dataset <- readRDS(file="saved/colVec_datasetColours.RDS")

tmp <- subset(agg_matrix_perCell_neurons_long_sub, Module==tt)
tmp$samplesToPseudobulk2 <- droplevels(tmp$samplesToPseudobulk2)


summary_stats <- tmp %>%
    group_by(samplesToPseudobulk2) %>%
    summarise(
      mean_activation = mean(ActivationScores, na.rm = TRUE),
      sd_activation = standard_dev(ActivationScores))

summary_stats <- as.data.frame(summary_stats)
summary_stats$color <- gsub("^ ","",gsub(".+-","", summary_stats$samplesToPseudobulk2))

summary_stats[grepl("This work", summary_stats$color),]$color <- "This work (2D, 3D, Foetal)"
stopifnot(all(summary_stats$color %in% names(colVec_dataset)))

summary_stats$samplesToPseudobulk2 <- factor(summary_stats$samplesToPseudobulk2, levels=rev(levels(summary_stats$samplesToPseudobulk2)))

suppFig11G <-ggplot(summary_stats, aes(x=mean_activation, y=samplesToPseudobulk2, col=color)) +
    theme_bw()+
    geom_point()+
    geom_errorbar(aes(xmin=mean_activation-sd_activation, xmax=mean_activation+sd_activation))+
    scale_color_manual(values=colVec_dataset)+
    theme(legend.position="none")+
    labs(y="",
         x=expression(paste("Mean Activation Score ", "\u00B1", " SD")))+
    ggtitle(tt)

pdf(file=paste0(rootSupp,"suppFigure11G.pdf"), height = 6, width=5)
plot(suppFig11G)
dev.off()



## Module 12 ##

####################
## Main Figure 4K ##
####################


modName <- "Module 12"
colVec_dataset <- readRDS(file="saved/colVec_datasetColours.RDS")

tmp <- subset(agg_matrix_perCell_neurons_long_sub, Module==tt)
tmp$samplesToPseudobulk2 <- droplevels(tmp$samplesToPseudobulk2)

summary_stats <- tmp %>%
    group_by(samplesToPseudobulk2) %>%
    summarise(
      mean_activation = mean(ActivationScores, na.rm = TRUE),
      sd_activation = standard_dev(ActivationScores))

summary_stats <- as.data.frame(summary_stats)
summary_stats$color <- gsub("^ ","",gsub(".+-","", summary_stats$samplesToPseudobulk2))

summary_stats[grepl("This work", summary_stats$color),]$color <- "This work (2D, 3D, Foetal)"
stopifnot(all(summary_stats$color %in% names(colVec_dataset)))

summary_stats$samplesToPseudobulk2 <- factor(summary_stats$samplesToPseudobulk2, levels=rev(levels(summary_stats$samplesToPseudobulk2)))

fig4K <- ggplot(summary_stats, aes(x=mean_activation, y=samplesToPseudobulk2, col=color)) +
    theme_bw()+
    geom_point()+
    geom_errorbar(aes(xmin=mean_activation-sd_activation, xmax=mean_activation+sd_activation))+
    scale_color_manual(values=colVec_dataset)+
    theme(legend.position="none")+
    labs(y="",
         x=expression(paste("Mean Activation Score ", "\u00B1", " SD")))+
    ggtitle(tt)

pdf(file=paste0(rootMain,"mainFigure4K.pdf"), height = 6, width=5)
plot(fig4K)
dev.off()





















