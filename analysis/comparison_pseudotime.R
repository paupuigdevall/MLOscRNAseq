library(tidyverse)
library(ggrastr)


rootSupp <- "figures/supp/"


## correlate cell types

pseudotimePerCell <- readRDS("saved/pseudotime/pseudotimePerCellSlingshot.RDS")
pseudotimePerCell$samplesToPseudobulk <- str_replace_all(pseudotimePerCell$samplesToPseudobulk, "Organoids", "3D")
pseudotimePerCell$ranked <- rank(pseudotimePerCell$pseudotimeMaxWeight)


pseudotimePerCellDestiny <- readRDS("saved/pseudotime/pseudotimePerCellDestiny.RDS")
pseudotimePerCellMonocle3 <- readRDS("saved/pseudotime/pseudotimePerCellMonocle3.RDS")


## set to NA those cells belonging not to the main partition
pseudotimePerCellMonocle3[is.infinite(pseudotimePerCellMonocle3$pseudotime_df),]$pseudotime_df <- NA
pseudotimePerCellMonocle3$pseudotime_ranked <- rank(pseudotimePerCellMonocle3$pseudotime_df, na.last="keep")


## Also, only in slinghot, cells belonging to rare cell types (nCells<50) were previously discarded, that explains the different dimensions with destiny and monocle3.
##QC
stopifnot(all(rownames(pseudotimePerCellDestiny)==rownames(pseudotimePerCellMonocle3)))

pseudotimeRaw <- data.frame(cellid=rownames(pseudotimePerCellDestiny),
                            slingshot=pseudotimePerCell[match(rownames(pseudotimePerCellDestiny), rownames(pseudotimePerCell)),]$pseudotimeMaxWeight,
                            destiny=pseudotimePerCellDestiny$pseudotime_MNN1_raw,
                            monocle3=pseudotimePerCellMonocle3[match(rownames(pseudotimePerCellDestiny), rownames(pseudotimePerCellMonocle3)),]$pseudotime_df)


pseudotimeOrderCells <- data.frame(cellid=rownames(pseudotimePerCellDestiny),
                                   slingshot=pseudotimePerCell[match(rownames(pseudotimePerCellDestiny), rownames(pseudotimePerCell)),]$ranked,
                                   destiny=pseudotimePerCellDestiny$pseudotime_MNN1_ordCells,
                                   monocle3=pseudotimePerCellMonocle3[match(rownames(pseudotimePerCellDestiny), rownames(pseudotimePerCellMonocle3)),]$pseudotime_ranked)

## Destiny with monocle3 ##
cor(pseudotimeRaw$destiny, pseudotimeRaw$monocle3, use="complete.obs", method="pearson")
#[1] 0.4844381
cor(pseudotimeOrderCells$destiny, pseudotimeOrderCells$monocle3, use="complete.obs", method="pearson")
# [1] 0.8586005


## Destiny with slingshot ##
cor(pseudotimeRaw$destiny, pseudotimeRaw$slingshot, use="complete.obs", method="pearson")
#[1] 0.3069006
cor(pseudotimeOrderCells$destiny, pseudotimeOrderCells$slingshot, use="complete.obs", method="pearson")
# [1] 0.2544053


tmp <- pseudotimeOrderCells
colnames(tmp)[-1] <- paste0(colnames(pseudotimeOrderCells)[-1],"_","CellsRanked")
tmp2 <- pseudotimeRaw
colnames(tmp2)[-1] <- c("slinghot_ptMaxWeight","destiny_firstEigenVec","monocle3_rawpt")
suppTab7 <- merge(tmp2, tmp, by="cellid")
stopifnot(dim(suppTab7)[1]==dim(tmp)[1])

## Export Table S7 ##

write.table(suppTab7, "saved/suppTables/TableS7.txt",
            quote=F, col.names=T, row.names=F, sep="\t")


## Slingshot with monocle3##
cor(pseudotimeRaw$slingshot, pseudotimeRaw$monocle3, use="complete.obs", method="pearson")
# [1] 0.4
cor(pseudotimeOrderCells$slingshot, pseudotimeOrderCells$monocle3, use="complete.obs", method="pearson")
# [1] 0.3956378



####################
## Supp Figure 8B ##
####################

library(tidyr)
#df <- as.data.frame(pseudotimeOrderCells %>% pivot_longer(-c("cellid"), names_to="Method", values_to="PseudotimeOrdNumCells"))
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

combinations <- combn(colnames(pseudotimeOrderCells)[2:length(colnames(pseudotimeOrderCells))], 2)

list_comb <- sapply(1:ncol(combinations), function(x) {
  c("cellid",combinations[,x])
}, simplify=F)


df <- sapply(list_comb, function(y){
  #print(y)
  tmp <- pseudotimeOrderCells[,y]
  
  if (length(grep("slingshot", y))>0){
    tmp$slingshot <- as.integer(tmp$slingshot)
  } 
  colnames(tmp) <- c("cellid","rep.x","rep.y")
  tmp$correlation <- paste(firstup(y[2:3]), collapse="-")
  return(tmp)
  
}, simplify=F)

df <- do.call("rbind", df)


reproPlot <- ggplot(df, aes(x=rep.x, y=rep.y))+
  geom_point_rast(alpha=0.1, size=0.1)+
  #geom_point(alpha=0.5, size=2)+
  facet_wrap(~correlation)+
  xlab("Ordered number of cells - Method 1")+
  ylab("Ordered number of cells - Method 2")+
  #ggtitle("Reproducibility on line proportion")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, size=14, face="bold"),
        axis.title=element_text(size=14),
        axis.text=element_text(size=10),
        strip.text = element_text(size = 14))+
  geom_smooth(method="lm", col="grey40", alpha=0.25, size=0.5)


corr <- sapply(unique(df$correlation), function(x){
  
  tmp <- subset(df, correlation==x)
  fit11 <- lm(rep.y ~ rep.x, data = tmp)
  
  data.frame(correlation=x,
             adj.r.squared=signif(summary(fit11)$adj.r.squared, 3),
             p=signif(summary(fit11)$coef[2,4], 3))
  
}, simplify=F)

corr <- do.call("rbind", corr)
rownames(corr) <- NULL
#rownames(corr) <- corr$type

suppFig8B <- reproPlot +
  geom_text(data=corr,
            aes(label=paste0("italic(R) ^ 2 == ", adj.r.squared)), x=125000, y=20000, col="blue",size=4, parse= T)


pdf(file=paste0(rootSupp,"suppFigure8B.pdf"), width=8, height = 6)
plot(suppFig8B)
dev.off()






