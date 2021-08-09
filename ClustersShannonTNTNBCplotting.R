##### Compare all aspects between TNBC cell lines and TNBC tumors(TN titles).
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

DF1 <- NULL # A dataframe collecting all numbers of superclones and subclones in TNBC cell lines.
DF2 <- NULL # A dataframe collecting all Shannon diversity scores of subclones in TNBC cell lines.
DF3 <- NULL # A dataframe collecting the varied CNA-classes total counts in each TNBC cell lines.
DF4 <- NULL # A subclone-wised seg-ratio dataframe for TNBC cell lines, used for computing Pearson corr.
DF5 <- NULL # A dataframe collecting the gain-loss attributed summaries for each TNBC cell lines.
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
  SSMat <- cbind(rep(TNSCases[W], 2), rep(OriCellN, 2), c("Superclone", "Subclone"), c(NofSup, NofSub))
  colnames(SSMat) <- c("Sample", "Total_Cell", "Type_of_cluster","Num_of_cluster")
  DF1 <- rbind.data.frame(DF1, SSMat)
  colnames(DF1) <- c("Sample", "Total_Cell", "Type_of_cluster","Num_of_cluster")
  
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
  
  ### Subclone perspective ShannonIndex
  SubMT <- as.data.frame(table(TES$subclones))
  Pi <- SubMT$Freq/sum(SubMT$Freq)
  SubMT <- cbind(SubMT, Pi, log(Pi), -(Pi*log(Pi))) 
  ddf2 <- cbind(TNSCases[W], sum(SubMT$`-(Pi * log(Pi))`))
  colnames(ddf2) <- c("Sample","ShannonIndex")
  DF2 <- rbind(DF2, ddf2)
  colnames(DF2) <- c("Sample","ShannonIndex")
  
  ### CNAs class color band bottom annotation from subclones' point of view
  subcN <- length(table(TES$subclones))
  subc_lonli <- split(NewsegMat, NewsegMat$Subclone)
  subc_consli <- parallel::mclapply(subc_lonli, function(x) { apply(x[,-c(1:2)], 2, median) }, mc.cores = 30)
  subc_csdf <- as.data.frame(do.call(rbind, subc_consli))
  # Let's set up the difference criteria between column mean and each subclonal median seg_ratio.
  subcVec <- NULL
  for (i in 1:ncol(subc_csdf)) {
    if(all(log2(subc_csdf[,i])<=0.25)&all(log2(subc_csdf[,i])>=(-0.25))){
      subcVec[i] <- 0
      next
    }
    SR <- 1-sum(abs(log2(median(subc_csdf[,i]))-log2(subc_csdf[,i])) >= 0.25)/subcN
    #print(SR)
    if(SR==(1-(1/subcN))){
      subcVec[i] <- 3
      next
    }else if(SR==1){
      subcVec[i] <- 1
      next
    }else{
      subcVec[i] <- 2
    }
  }
  table(subcVec)
  ddf3 <- cbind(rep(TNSCases[W], 3), c("cCNA", "sCNA", "uCNA"), c("orange", "black", "purple"), table(subcVec)[names(table(subcVec))!=0])
  colnames(ddf3) <- c("Sample","CNA_Class","Color","Count")
  DF3 <- rbind(DF3, ddf3)
  colnames(DF3) <- c("Sample","CNA_Class","Color","Count")
  ### Subclone ratios matrix
  rownames(subc_csdf) <- paste(TNSCases[W], "_", rownames(subc_csdf), sep = "")
  DF4 <- rbind(subc_csdf, DF4)
  
  ### Gain and loss statistics
  NewsegMat <- cbind(rep(TNSCases[W], nrow(NewsegMat)), NewsegMat)
  colnames(NewsegMat)[1] <- "Sample"
  All_lonli <- split(NewsegMat, NewsegMat$Sample)
  All_consli <- parallel::mclapply(All_lonli, function(x) { apply(x[,-c(1:3)], 2, median) }, mc.cores = 30)
  All_csdf <- as.data.frame(do.call(rbind, All_consli))
  All_csdf <- log2(All_csdf)
  
  G1 <- length(All_csdf[All_csdf>=0.25])/ncol(All_csdf)*100
  N1 <- length(All_csdf[All_csdf<0.25&All_csdf>(-0.25)])/ncol(All_csdf)*100
  L1 <- length(All_csdf[All_csdf<=(-0.25)])/ncol(All_csdf)*100
  ddf5 <- cbind(rep(TNSCases[W], 3), c("Gain", "Neutral", "Loss"), c(G1,N1,L1))
  colnames(ddf5) <- c("Sample","Alteration","Percentage")
  DF5 <- rbind(DF5, ddf5)
  colnames(DF5) <- c("Sample","Alteration","Percentage")
}

### TN samples processing pipeline
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

df1 <- NULL # A dataframe collecting all numbers of superclones and subclones in tumor samples.
df2 <- NULL # A dataframe collecting all Shannon diversity scores of subclones in tumor samples.
df3 <- NULL # A dataframe collecting the varied CNA-classes total counts in each tumor samples.
df4 <- NULL # A subclone-wised seg-ratio dataframe for tumor samples, used for computing Pearson corr.
df5 <- NULL # A dataframe collecting the gain-loss attributed summaries for each tumor samples.

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
  SSMat <- cbind(rep(TNSCases[C], 2), rep(OriCellN, 2), c("Superclone", "Subclone"), c(NofSup, NofSub))
  colnames(SSMat) <- c("Sample", "Total_Cell", "Type_of_cluster","Num_of_cluster")
  df1 <- rbind.data.frame(df1, SSMat)
  colnames(df1) <- c("Sample", "Total_Cell", "Type_of_cluster","Num_of_cluster")
  
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
  
  ### Subclone perspective ShannonIndex
  SubMT <- as.data.frame(table(TES$subclones))
  Pi <- SubMT$Freq/sum(SubMT$Freq)
  SubMT <- cbind(SubMT, Pi, log(Pi), -(Pi*log(Pi))) 
  ddf2 <- cbind(TNSCases[C], sum(SubMT$`-(Pi * log(Pi))`))
  colnames(ddf2) <- c("Sample","ShannonIndex")
  df2 <- rbind(df2, ddf2)
  colnames(df2) <- c("Sample","ShannonIndex")
  
  ### CNAs class color band bottom annotation from subclones' point of view
  subcN <- length(table(TES$subclones))
  subc_lonli <- split(NewsegMat, NewsegMat$Subclone)
  subc_consli <- parallel::mclapply(subc_lonli, function(x) { apply(x[,-c(1:2)], 2, median) }, mc.cores = 30)
  subc_csdf <- as.data.frame(do.call(rbind, subc_consli))
  # Let's set up the difference criteria between column mean and each subclonal median seg_ratio.
  subcVec <- NULL
  for (i in 1:ncol(subc_csdf)) {
    if(all(log2(subc_csdf[,i])<=0.25)&all(log2(subc_csdf[,i])>=(-0.25))){
      subcVec[i] <- 0
      next
    }
    SR <- 1-sum(abs(log2(median(subc_csdf[,i]))-log2(subc_csdf[,i])) >= 0.25)/subcN
    #print(SR)
    if(SR==(1-(1/subcN))){
      subcVec[i] <- 3
      next
    }else if(SR==1){
      subcVec[i] <- 1
      next
    }else{
      subcVec[i] <- 2
    }
  }
  table(subcVec)
  ddf3 <- cbind(rep(TNSCases[C], 3), c("cCNA", "sCNA", "uCNA"), c("orange", "black", "purple"), table(subcVec)[names(table(subcVec))!=0])
  colnames(ddf3) <- c("Sample","CNA_Class","Color","Count")
  df3 <- rbind(df3, ddf3)
  colnames(df3) <- c("Sample","CNA_Class","Color","Count")
  ### Subclone ratios matrix
  rownames(subc_csdf) <- paste(TNSCases[C], "_", rownames(subc_csdf), sep = "")
  df4 <- rbind(subc_csdf, df4)
  
  ### Gain and loss statistics
  NewsegMat <- cbind(rep(TNSCases[C], nrow(NewsegMat)), NewsegMat)
  colnames(NewsegMat)[1] <- "Sample"
  All_lonli <- split(NewsegMat, NewsegMat$Sample)
  All_consli <- parallel::mclapply(All_lonli, function(x) { apply(x[,-c(1:3)], 2, median) }, mc.cores = 30)
  All_csdf <- as.data.frame(do.call(rbind, All_consli))
  All_csdf <- log2(All_csdf)
  
  G1 <- length(All_csdf[All_csdf>=0.25])/ncol(All_csdf)*100
  N1 <- length(All_csdf[All_csdf<0.25&All_csdf>(-0.25)])/ncol(All_csdf)*100
  L1 <- length(All_csdf[All_csdf<=(-0.25)])/ncol(All_csdf)*100
  ddf5 <- cbind(rep(TNSCases[C], 3), c("Gain", "Neutral", "Loss"), c(G1,N1,L1))
  colnames(ddf5) <- c("Sample","Alteration","Percentage")
  df5 <- rbind(df5, ddf5)
  colnames(df5) <- c("Sample","Alteration","Percentage")
}

### Make statistically summarized plots
# Colors Customization
library(cluster) #General color Idex
colors = c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabebe"
           , "#469990", "#e6beff", "#9A6324", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9", "#000000")

DFdf1 <- rbind(DF1, df1)
DFdf1$Num_of_cluster <- as.numeric(as.character(DFdf1$Num_of_cluster))
DFdf1 <- DFdf1[ order(DFdf1$Num_of_cluster, decreasing = F),]
tiff("Clonality_Statistic.tiff", width=2000, height=1000, compression="lzw", res=300)
ggplot(DFdf1, aes(x = Sample, y = Num_of_cluster, fill = Type_of_cluster)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_y_continuous(breaks=seq(0,31,2), limits = c(0, 31)) +
  theme(axis.text.x = element_text(size=6, angle = 45, vjust = 0.5), axis.text.y = element_text(size=12)) +
  scale_x_discrete()
dev.off()

DFdf2 <- rbind(DF2, df2) %>% as.data.frame()
DFdf2$ShannonIndex <- as.numeric(as.character(DFdf2$ShannonIndex))
DFdf2 <- DFdf2[ order(DFdf2$ShannonIndex, decreasing = F), ]
level_order <- DFdf2$Sample
tiff("ShannonIndexing.tiff", width=2000, height=1000, compression="lzw", res=300)
ggplot(DFdf2, aes(x = factor(Sample, level = level_order), y = ShannonIndex, fill = colors[1:nrow(DFdf2)])) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_y_continuous(breaks=seq(0,4,0.25), limits = c(0, 4)) +
  theme(axis.text.x = element_text(size=6, angle = 45, vjust = 0.5), axis.text.y = element_text(size=12), legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete() + ggtitle("Subclone_Diversity") +
  xlab("Sample") + ylab("Shannon_Diversity")
dev.off()

DFdf3 <- rbind(DF3, df3) %>% as.data.frame()
DFdf3$Count <- as.numeric(as.character(DFdf3$Count))
tiff("CNAs_Statistics.tiff", width=2000, height=1000, compression="lzw", res=300)
ggplot(DFdf3, aes(x = Sample, y = Count, fill = CNA_Class)) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  scale_y_continuous() + scale_fill_manual(breaks = c("cCNA", "sCNA", "uCNA"), values=c("orange", "black", "purple")) +
  theme(axis.text.x = element_text(size=6, angle = 45, vjust = 0.5), axis.text.y = element_text(size=12), plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete() + ggtitle("Total_CNAs_Attribution") +
  xlab("Sample") + ylab("CNAs_Counts")
dev.off()

DFdf5 <- rbind(DF5, df5) %>% as.data.frame()
DFdf5$Percentage <- as.numeric(as.character(DFdf5$Percentage))
tiff("GainANDLoss.tiff", width=2000, height=1000, compression="lzw", res=300)
ggplot(DFdf5, aes(x = Sample, y = Percentage, fill = Alteration)) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  scale_y_continuous(breaks=seq(0,101,10), limits = c(0, 101)) + scale_fill_manual(breaks = c("Gain", "Neutral", "Loss"), values=c("red", "gray", "blue")) +
  theme(axis.text.x = element_text(size=6, angle = 45, vjust = 0.5), axis.text.y = element_text(size=12), plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete() + ggtitle("Gain_and_Loss_Percentage") +
  xlab("Sample") + ylab("Percentage_across_Genome")
dev.off()

### Make subclone correlation map using build-in functions in R
DFdf4 <- rbind(DF4, df4)
cormat <- round(cor(t(DFdf4)),2)
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
tiff("Subclone_PearsonCorr.tiff", width=4000, height=4000, compression="lzw", res=300)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 5, hjust = 1, face="bold"),
        axis.text.y = element_text(vjust = 1, size = 5, hjust = 1, face="bold"),
        plot.title = element_text(size = 20, hjust = 0.5, face="bold"))+
  coord_fixed()+ ggtitle("Correlation_among_subclones")
dev.off()
