##### Pearson computation among superclones and bulk samples
library(dplyr)
library(SingleCellExperiment)
library(ggplot2)
library(copykit)
### TNBCellines samples processing pipeline
TNSLoc <- "/volumes/seq/projects/CNA_projects/DT_CNA/cell_line/" 
BT20 <- "BT20/BT20_P1_P2_P3/final_result"
MDAMB231 <- "231/MDAMB231/MDAMB231/MDAMB231_P1_P2_P3/output/final_result"
MDAMB157 <- "MDAMB157/MB157_P1_P2_P3/final_result" 
MDAMB453 <- "MDAMB453/MB453_P1_P2_P3/final_result"
HCC70 <- "HCC70/HCC70_P2_P3_P4_combine/output/final_result" 
BT549 <- "BT549/output/final_result"
SUM159 <- "SUM159/SUM159_Parental_P2_P3_P4/output/final_result"
TNSpath <- c(BT20,MDAMB231,MDAMB157,MDAMB453,HCC70,BT549,SUM159)
TNSCases <- c("BT20","MDAMB231","MDAMB157","MDAMB453","HCC70","BT549","SUM159")
FilReso <- c(0.85, 0.85, 0.85, 0.9, 0.85, 0.9, 0.85)

DF6 <- NULL # A superclone-wised seg-ratio dataframe for TNBC cell lines, used for computing Pearson corr.
DF7 <- NULL # A bulk seg-ratio dataframe for TNBC cell lines, used for computing Pearson corr.
for(W in 1:length(TNSCases)){
  TNS <- copykit::readVarbinCNA(paste(TNSLoc, TNSpath[W],sep = ""), remove_Y = TRUE)
  if(W == 2){
    TNS <- TNS[ , 1:986]
  }
  SummarizedExperiment::rowRanges(TNS)
  TNS <- copykit::runMetrics(TNS)
  #plotHeatmap(TNS, label = "is_normal", row_split = "is_normal")
  TNS <- filterCells(TNS, resolution = FilReso[W])
  #plotHeatmap(TNS, label = "filtered", row_split = "filtered")
  TNS <- TNS[,colData(TNS)$filtered == "kept"]
  OriCellN <- ncol(TNS)
  if(ncol(TNS) < 1000){
    Addi <- sample(1:ncol(TNS), 1000-ncol(TNS), replace = F)
    MOCK <- TNS[ , Addi]
    rownames(MOCK@colData) <- paste("ADDon_Mock", 1:(1000-ncol(TNS)), sep = "")
    MOCK@colData$sample <- paste("ADDon_Mock", 1:(1000-ncol(TNS)), sep = "")
    TNS <- cbind(TNS, MOCK)
  }else{
    TNS <- TNS[ , sample(1:ncol(TNS), 1000, replace = F)]
  }
  head(SummarizedExperiment::colData(TNS))
  TNS <- copykit::runDistMat(TNS)
  TNS <- copykit::runUmap(TNS, min_dist = 0, n_neighbors = 40)
  SingleCellExperiment::reducedDim(TNS, 'umap', withDimnames = FALSE)
  opt_k <- findOptimalK(TNS)$k
  TNS <- copykit::findClusters(TNS, k_superclones = 35, k_subclones = opt_k)
  NofSup <- TNS@colData$superclones %>% unique() %>% length()
  NofSub <- TNS@colData$subclones %>% unique() %>% length()
  
  ### Convert factors to numeric and obtain the ordered colData matrix
  TNS@colData$superclones <- as.character(TNS@colData$superclones)
  TNS@colData$subclones <- as.character(TNS@colData$subclones)
  TES <- TNS@colData[order(TNS@colData[ ,6], TNS@colData[ ,7]), ]
  # Sort segment ratios according to ordered colData sample names
  seg_data <- t(copykit::segment_ratios(TNS))
  SEQ <- NULL
  for (i in 1:length(TES$sample)) { SEQ <- c(SEQ, which(rownames(seg_data)%in%rownames(TES)[i])) }
  SEQ
  seg_data <- seg_data[SEQ, ]
  NewsegMat <- cbind.data.frame(TES$superclones, TES$subclones, seg_data)
  colnames(NewsegMat)[1:2] <- c("Superclone","Subclone")
  
  ### CNAs class color band bottom annotation from subclones' point of view
  subcN <- length(table(TES$superclones))
  subc_lonli <- split(NewsegMat, NewsegMat$Superclone)
  subc_consli <- parallel::mclapply(subc_lonli, function(x) { apply(x[,-c(1:2)], 2, median) }, mc.cores = 30)
  subc_csdf <- as.data.frame(do.call(rbind, subc_consli))
  
  ### Subclone ratios matrix
  rownames(subc_csdf) <- paste(TNSCases[W], "_", rownames(subc_csdf), sep = "")
  DF6 <- rbind(subc_csdf, DF6)
  
  NewsegMat <- cbind(rep(TNSCases[W], nrow(NewsegMat)), NewsegMat)
  colnames(NewsegMat)[1] <- "Sample"
  All_lonli <- split(NewsegMat, NewsegMat$Sample)
  All_consli <- parallel::mclapply(All_lonli, function(x) { apply(x[,-c(1:3)], 2, median) }, mc.cores = 30)
  All_csdf <- as.data.frame(do.call(rbind, All_consli))
  ### Whole sample ratios matrix
  rownames(All_csdf) <- paste(TNSCases[W], sep = "")
  DF7 <- rbind(All_csdf, DF7)
}

### Tumors samples processing pipeline (TN titles)
TNSLoc <- "/volumes/seq/projects/CNA_projects/DT_CNA/snap_frozen/Breast/TNBC/"
TN4 <- "TN4/output/final_result"
TN7 <- "TN7/output/final_result"
TN17 <- "TN17/output/final_result"
TN20 <- "TN20/TN20_2020_06_19_combine_new_run_output/output/final_result"
TN21 <- "TN21/output/final_result"
TN26 <- "TN26/output/final_result"
TN27 <- "TN27/output/final_result"
TN28 <- "TN28/TN28_2020_06_19_combine_new_run_output/output/final_result"
TNSpath <- c(TN4,TN7,TN17,TN20,TN21,TN26,TN27,TN28)
TNSCases <- c("TN4","TN7","TN17","TN20","TN21","TN26","TN27","TN28")

df6 <- NULL # A superclone-wised seg-ratio dataframe for TNBC tumor cases, used for computing Pearson corr.
df7 <- NULL # A bulk seg-ratio dataframe for TNBC tumor cases, used for computing Pearson corr.
for(C in 1:length(TNSCases)){
  TNS <- copykit::readVarbinCNA(paste(TNSLoc, TNSpath[C],sep = ""), remove_Y = TRUE)
  SummarizedExperiment::rowRanges(TNS)
  TNS <- copykit::runMetrics(TNS)
  TNS <- findNormalCells(TNS)
  #plotHeatmap(TNS, label = "is_normal", row_split = "is_normal")
  TNS <- TNS[,colData(TNS)$is_normal == FALSE]
  TNS <- filterCells(TNS, resolution = 0.9)
  #plotHeatmap(TNS, label = "filtered", row_split = "filtered")
  TNS <- TNS[,colData(TNS)$filtered == "kept"]
  OriCellN <- ncol(TNS)
  if(ncol(TNS) < 1000){
    Addi <- sample(1:ncol(TNS), 1000-ncol(TNS), replace = F)
    MOCK <- TNS[ , Addi]
    rownames(MOCK@colData) <- paste("ADDon_Mock", 1:(1000-ncol(TNS)), sep = "")
    MOCK@colData$sample <- paste("ADDon_Mock", 1:(1000-ncol(TNS)), sep = "")
    TNS <- cbind(TNS, MOCK)
  }else{
    TNS <- TNS[ , sample(1:ncol(TNS), 1000, replace = F)]
  }
  head(SummarizedExperiment::colData(TNS))
  TNS <- copykit::runDistMat(TNS)
  TNS <- copykit::runUmap(TNS, min_dist = 0, n_neighbors = 40)
  SingleCellExperiment::reducedDim(TNS, 'umap', withDimnames = FALSE)
  opt_k <- findOptimalK(TNS)$k
  TNS <- copykit::findClusters(TNS, k_superclones = 35, k_subclones = opt_k)
  NofSup <- TNS@colData$superclones %>% unique() %>% length()
  NofSub <- TNS@colData$subclones %>% unique() %>% length()
  
  ### Convert factors to numeric and obtain the ordered colData matrix
  TNS@colData$superclones <- as.character(TNS@colData$superclones)
  TNS@colData$subclones <- as.character(TNS@colData$subclones)
  TES <- TNS@colData[order(TNS@colData[ ,8], TNS@colData[ ,9]), ]
  # Sort segment ratios according to ordered colData sample names
  seg_data <- t(copykit::segment_ratios(TNS))
  SEQ <- NULL
  for (i in 1:length(TES$sample)) { SEQ <- c(SEQ, which(rownames(seg_data)%in%rownames(TES)[i])) }
  SEQ
  seg_data <- seg_data[SEQ, ]
  NewsegMat <- cbind.data.frame(TES$superclones, TES$subclones, seg_data)
  colnames(NewsegMat)[1:2] <- c("Superclone","Subclone")
  
  ### CNAs class color band bottom annotation from subclones' point of view
  subcN <- length(table(TES$superclones))
  subc_lonli <- split(NewsegMat, NewsegMat$Superclone)
  subc_consli <- parallel::mclapply(subc_lonli, function(x) { apply(x[,-c(1:2)], 2, median) }, mc.cores = 30)
  subc_csdf <- as.data.frame(do.call(rbind, subc_consli))
  
  ### Subclone ratios matrix
  rownames(subc_csdf) <- paste(TNSCases[C], "_", rownames(subc_csdf), sep = "")
  df6 <- rbind(subc_csdf, df6)
  
  NewsegMat <- cbind(rep(TNSCases[C], nrow(NewsegMat)), NewsegMat)
  colnames(NewsegMat)[1] <- "Sample"
  All_lonli <- split(NewsegMat, NewsegMat$Sample)
  All_consli <- parallel::mclapply(All_lonli, function(x) { apply(x[,-c(1:3)], 2, median) }, mc.cores = 30)
  All_csdf <- as.data.frame(do.call(rbind, All_consli))
  ### Whole sample ratios matrix
  rownames(All_csdf) <- paste(TNSCases[C], sep = "")
  df7 <- rbind(All_csdf, df7)
}
### Make superclone correlation map using build-in functions in R
DFdf6 <- rbind(DF6, df6)
cormat <- round(cor(t(DFdf6)),2)
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
library(reshape2)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap
tiff("Superclone_PearsonCorr.tiff", width=2300, height=2300, compression="lzw", res=300)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 5, hjust = 1, face="bold"),
        axis.text.y = element_text(vjust = 1, size = 5, hjust = 1, face="bold"),
        plot.title = element_text(size = 20, hjust = 0.5, face="bold"))+
  coord_fixed()+ ggtitle("Correlation_among_superclones")
dev.off()

### Make corrplot to nicely summarize
DFdf7 <- rbind(DF7, df7)
cormat2 <- round(cor(t(DFdf7)),2)
library(corrplot)
tiff("scCNAsClonal_Correlation.tiff", width=2000, height=2000, compression="lzw", res=300)
corrplot.mixed(cormat2, order="hclust", tl.col="black", title = "scCNAs Correlation among Samples",
               lower.col=colorRampPalette(c("blue","white","red"))(200), mar=c(0, 0, 1, 0), tl.pos = "d",
               upper.col=colorRampPalette(c("blue","white","red"))(200), number.cex = 0.7, tl.cex = 0.7)
dev.off()

