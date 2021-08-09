### Check bin counts of cells in messy profiling in CNV heatmaps
TNBCLoc <- "/volumes/seq/projects/CNA_projects/DT_CNA/cell_line/" 
BT20 <- "BT20/BT20_P1_P2_P3/final_result"
MDAMB231 <- "231/MDAMB231/MDAMB231/MDAMB231_P1_P2_P3/output/final_result"
MDAMB157 <- "MDAMB157/MB157_P1_P2_P3/final_result" 
MDAMB453 <- "MDAMB453/MB453_P1_P2_P3/final_result"
HCC70 <- "HCC70/HCC70_P2_P3_P4_combine/output/final_result" 
BT549 <- "BT549/output/final_result"
SUM159 <- "SUM159/SUM159_Parental_P2_P3_P4/output/final_result"
TNBCpath <- c(BT20,MDAMB231,MDAMB157,MDAMB453,HCC70,BT549,SUM159)
TNBCases <- c("BT20","MDAMB231","MDAMB157","MDAMB453","HCC70","BT549","SUM159")
FilReso <- c(0.85, 0.85, 0.85, 0.90, 0.85, 0.9, 0.85)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
library(copykit)
BeginT <- Sys.time()
FListVecs <- NULL
for (M in 1:length(TNBCases)) {
  Cas <- TNBCases[M]
  print(Cas)
  TNBCelline <- copykit::readVarbinCNA(paste(TNBCLoc, TNBCpath[M],sep = ""), remove_Y = TRUE)
  if(M == 2){
    TNBCelline <- TNBCelline[ , 1:986]
  }
  SummarizedExperiment::rowRanges(TNBCelline)
  TNBCelline <- copykit::runMetrics(TNBCelline)
  TNBCelline <- copykit::filterCells(TNBCelline,resolution = FilReso[M])
  #plotHeatmap(TNBCelline, label = "filtered", row_split = "filtered")
  #head(SummarizedExperiment::colData(TNBCelline))
  TNBCelline <- TNBCelline[ , which(SummarizedExperiment::colData(TNBCelline)$filtered == "kept")]
  TNBCelline <- copykit::runDistMat(TNBCelline)
  TNBCelline <- copykit::runUmap(TNBCelline, min_dist = 0, n_neighbors = 40)
  SingleCellExperiment::reducedDim(TNBCelline, 'umap', withDimnames = FALSE)
  
  CluVec <- NULL
  for (i in 1:ceiling(ncol(TNBCelline)/100)) {
    if(i == ceiling(ncol(TNBCelline)/100)){
      i <- ncol(TNBCelline)/100
    }
    SAMa <- TNBCelline[, sample.int(ncol(TNBCelline), i*100)]
    opt_k <- findOptimalK(SAMa)$k
    CluVec100 <- NULL
    for (j in 1:50) {
      SAMa <- TNBCelline[, sample.int(ncol(TNBCelline), i*100)]
      SAMa <- copykit::runDistMat(SAMa)
      SAMa <- copykit::runUmap(SAMa, min_dist = 0, n_neighbors = 40)
      SingleCellExperiment::reducedDim(SAMa, 'umap', withDimnames = FALSE)
      SAMa <- copykit::findClusters(SAMa, k_superclones = 35, k_subclones = opt_k)
      Nsub <- unique(SAMa@colData$subclones) %>% length()
      CluVec100 <- c(CluVec100, Nsub)
    }
    CluVec <- c(CluVec, round(sum(CluVec100)/50, 1))
  }
  FListVecs[[M]] <- CluVec
}
TimeSpent <- Sys.time() - BeginT

plotHeatmap(SAMa, label = c("superclones", "subclones"))
### Line plot with all 7 TNBCellines dissected clustering vectors
# Colors Customization
library(cluster) #General color Idex
colors = c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabebe"
           , "#469990", "#e6beff", "#9A6324", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9", "#000000")
plot(1:7, pch = 19, col=colors[1:7], cex = 5)

tiff("Clustering_of_downsizedcell_sampling.tiff", width=3000, height=2500, compression="lzw", res=300)
plot(FListVecs[[1]], type = "b", pch = 19, col = colors[1], xlab = "Hundred Cell", ylab = "Num_of_Subclone",
     main = "Clustering of downsized cell sampling", xlim = c(1, 14), ylim = c(1, 18))
axis(1, at = seq(1, 14, by = 1))
axis(2, at = seq(1, 18, by = 1))
lines(FListVecs[[2]], type = "b", pch = 19, col = colors[2])
lines(FListVecs[[3]], type = "b", pch = 19, col = colors[3])
lines(FListVecs[[4]], type = "b", pch = 19, col = colors[4])
lines(FListVecs[[5]], type = "b", pch = 19, col = colors[5])
lines(FListVecs[[6]], type = "b", pch = 19, col = colors[6])
lines(FListVecs[[7]], type = "b", pch = 19, col = colors[7])
legend("topleft", legend = c("BT20", "MDAMB231", "MDAMB157", "MDAMB453", "HCC70", "BT549", "SUM159"),
       col = colors[1:7],
       title="TNBCelline", text.font=4, bg="aliceblue", pch = 19, bty = "o", pt.cex = 2, 
       cex = 1, text.col = "black", horiz = F , inset = c(0.05, 0.05))
dev.off()