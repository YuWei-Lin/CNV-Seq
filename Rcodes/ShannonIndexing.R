### Calculate Shannon entropy in TNBC cell lines
TNBCLoc <- "/volumes/seq/projects/CNA_projects/DT_CNA/cell_line/" 
MDAMB157 <- "MDAMB157/MB157_P1_P2_P3/final_result"
BT20 <- "BT20/BT20_P1_P2_P3/final_result"
MDAMB453 <- "MDAMB453/MB453_P1_P2_P3/final_result" 
HCC70 <- "HCC70P1/output/final_result"
BT549 <- "BT549/BT549ParentalP4_PM2133/output/final_result"
HCC38 <- "HCC38/output_200k/final_result"
HCC1806 <- "HCC1806/output/final_result" 
Hs587T <- "HS587/HS587T_Parental_P1/output/final_result"

TNBCases <- c(MDAMB157,BT20,MDAMB453,HCC70,BT549,HCC38,HCC1806,Hs587T)
TNBCasess <- c("MDAMB157","BT20","MDAMB453","HCC70","BT549","HCC38","HCC1806","Hs587T")
FilReso <- c(0.75, 0.75, 0.95, 0.9, 0.6, 0.7,0.7,0.8)
library(dplyr)
library(ggplot2)

### Quantify Shannon entropy from the perspective of superclone and subclone
SuperC.DIV <- NULL
SubC.DIV <- NULL
for (i in 1:length(TNBCasess)) {
  TNBCelline <- copykit::readVarbinCNA(paste(TNBCLoc, TNBCases[i],sep = ""), remove_Y = TRUE)
  TNBCelline <- copykit::runMetrics(TNBCelline)
  TNBCelline <- copykit::filterCells(TNBCelline,resolution = FilReso[i], plot_heatmap = F) 
  TNBCelline <- TNBCelline[ ,which(SummarizedExperiment::colData(TNBCelline)$filtered == "kept")] #"kept" and "removed"
  TNBCelline <- copykit::runDistMat(TNBCelline)
  TNBCelline <- copykit::runUmap(TNBCelline,min_dist = 0.001,n_neighbors = 30)
  SingleCellExperiment::reducedDim(TNBCelline, 'umap', withDimnames = FALSE)
  TNBCelline <- copykit::findClusters(TNBCelline, k_major = 35, k_minor = 21)
  ### Convert factors to numeric and obtain the ordered colData matrix
  TNBCelline@colData$subclones <- as.numeric(as.character(TNBCelline@colData$subclones))
  TES <- TNBCelline@colData[order(TNBCelline@colData[ ,6], TNBCelline@colData[ ,7]), ]
  ### Superclone perspective
  SupMT <- as.data.frame(table(TES$superclones))
  Pi <- SupMT$Freq/sum(SupMT$Freq)
  SupMT <- cbind(SupMT, Pi, log(Pi), -(Pi*log(Pi)))
  OutputM <- cbind(rep(TNBCasess[i], nrow(SupMT)+2), c(rep("Hi", nrow(SupMT)), "Hmax", "Even"), c(as.character(SupMT$Var1), "X", "Y"), 
                   c(SupMT$`-(Pi * log(Pi))`, log(nrow(SupMT)), sum(SupMT$`-(Pi * log(Pi))`)/log(nrow(SupMT))))
  colnames(OutputM) <- c("Celline","Class","Superclone","ShannonIndex")
  write.csv(OutputM, paste("/volumes/lab/users/yuwei/", TNBCasess[i],"_Superclone_ShannonMatrix.csv", sep = ""), row.names = F)
  SuperC.DIV <- rbind.data.frame(SuperC.DIV, OutputM)
  ### Subclone perspective
  SubMT <- as.data.frame(table(TES$subclones))
  Pi <- SubMT$Freq/sum(SubMT$Freq)
  SubMT <- cbind(SubMT, Pi, log(Pi), -(Pi*log(Pi)))
  OutputMM <- cbind(rep(TNBCasess[i], nrow(SubMT)+2), c(rep("Hi", nrow(SubMT)), "Hmax", "Even"), c(as.character(SubMT$Var1), "X", "Y"), 
                   c(SubMT$`-(Pi * log(Pi))`, log(nrow(SubMT)), sum(SubMT$`-(Pi * log(Pi))`)/log(nrow(SubMT))))
  colnames(OutputMM) <- c("Celline","Class","Subclone","ShannonIndex")
  write.csv(OutputMM, paste("/volumes/lab/users/yuwei/", TNBCasess[i],"_Subclone_ShannonMatrix.csv", sep = ""), row.names = F)
  SubC.DIV <- rbind.data.frame(SubC.DIV, OutputMM)
}
SuperC.DIV$ShannonIndex <- round(as.numeric(as.character(SuperC.DIV$ShannonIndex)),2)
SubC.DIV$ShannonIndex <- round(as.numeric(as.character(SubC.DIV$ShannonIndex)),2)
positions <- c("Hi", "Hmax", "Even")

### Stacked barchart plotting by ggplot2
tiff("ShannonIndex_SupercloneBarChart.tiff", width=2000, height=1000, compression="lzw", res=300)
ggplot(SuperC.DIV, aes(x = Class, y = ShannonIndex, fill = Superclone)) + 
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ Celline) + 
  scale_y_continuous(breaks=seq(0,3,0.25), limits = c(0, 3)) + 
  geom_text(aes(label = ShannonIndex), position = position_stack(vjust = 0.5), size = 5) +
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12)) + 
  scale_x_discrete(limits = positions)
dev.off()

tiff("ShannonIndex_SubcloneBarChart.tiff", width=2000, height=1000, compression="lzw", res=300)
ggplot(SubC.DIV, aes(x = Class, y = ShannonIndex, fill = Subclone)) + 
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ Celline) + 
  scale_y_continuous(breaks=seq(0,3,0.25), limits = c(0, 3)) + 
  geom_text(aes(label = ShannonIndex), position = position_stack(vjust = 0.5), size = 5) +
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12)) + 
  scale_x_discrete(limits = positions)
dev.off()



'''
test  <- data.frame(person=c("A", "B", "C", "D", "E"), 
                    value1=c(100,150,120,80,150),     
                    value2=c(25,30,45,30,30) , 
                    value3=c(100,120,150,150,200)) 
library(reshape2) # for melt
melted <- melt(test, "person")
melted$cat <- ''
melted[melted$variable == 'value1',]$cat <- "first"
melted[melted$variable != 'value1',]$cat <- "second"
ggplot(melted, aes(x = cat, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ person)
'''

