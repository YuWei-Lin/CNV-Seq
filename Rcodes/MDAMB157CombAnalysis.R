### Libraries 
library(dplyr)
library(SingleCellExperiment)
library(ggplot2)
library(copykit)

### Merge study of MDAMB157 single cells and bulks of its 84 subclones
# Import MDAMB157 scCNV:
Loc <- "/volumes/seq/projects/CNA_projects/DT_CNA/cell_line/MDAMB157/MB157_P1_P2_P3/final_result"
MDAMB157 <- copykit::readVarbinCNA(Loc, remove_Y = TRUE)

# Import MDAMB157Subs bulkCNV:
SEC <- "/volumes/lab/users/yuwei/MDAMB/output/final_result"
Sub157 <- copykit::readVarbinCNA(SEC, remove_Y = TRUE)

# Combine scCNV objects, RMSE filtering and hdbscan clustering
CaseName <- "Merged_MDAMB157"
scCNVcomb <- cbind(MDAMB157, Sub157)
SummarizedExperiment::rowRanges(scCNVcomb)
scCNVcomb <- copykit::runMetrics(scCNVcomb)
scCNVcomb <- copykit::filterCells(scCNVcomb, resolution = 0.85)
tiff("MDAMB157_Merged_scCNV_Fil(0.85).tiff", width=2000, height=2000, compression="lzw", res=300)
plotHeatmap(scCNVcomb, label = "filtered", row_split = "filtered")
dev.off()
scCNVcomb <- scCNVcomb[ ,which(SummarizedExperiment::colData(scCNVcomb)$filtered == "kept")]
scCNVcomb <- copykit::runDistMat(scCNVcomb)
scCNVcomb <- copykit::runUmap(scCNVcomb, min_dist = 0, n_neighbors = 40)
SingleCellExperiment::reducedDim(scCNVcomb, 'umap', withDimnames = FALSE)
opt_k <- findOptimalK(scCNVcomb)$k
scCNVcomb <- copykit::findClusters(scCNVcomb, k_superclones = 35, k_subclones = opt_k)

### Umap and Heapmap plotting
# Adding new columns of umap coorfinates, identities and source of data points
head(scCNVcomb@colData)
GetMiddleName <- function(x) {
  y <- strsplit(x, "_")[[1]][2]
  y <- gsub('[0-9]+', '', toupper(y))
  return(y)
}
scCNVcomb@colData$"Umap_X" <- SingleCellExperiment::reducedDim(scCNVcomb, 'umap', withDimnames = FALSE)[, 1]
scCNVcomb@colData$"Umap_Y" <- SingleCellExperiment::reducedDim(scCNVcomb, 'umap', withDimnames = FALSE)[, 2]
scCNVcomb@colData$"Ident" <- unlist(lapply(rownames(scCNVcomb@colData), GetMiddleName))
scCNVcomb@colData$"Source" <- c(rep("Single-Cell", nrow(scCNVcomb@colData)-84), rep("Subclones_Bulk", 84))
# Convert factors to numeric and obtain the ordered colData matrix
scCNVcomb@colData$superclones <- as.character(scCNVcomb@colData$superclones)
scCNVcomb@colData$subclones <- as.character(scCNVcomb@colData$subclones)
TES <- scCNVcomb@colData[order(scCNVcomb@colData[ ,6], scCNVcomb@colData[ ,7], scCNVcomb@colData[ ,11]), ]


#### Do you want to include CNAs events? Need do it manually like before. Can't use copykit!!!
# Sort segment ratios according to ordered colData sample names
seg_data <- t(copykit::segment_ratios(scCNVcomb))
SEQ <- NULL
for (i in 1:length(TES$sample)) { SEQ <- c(SEQ, which(rownames(seg_data)%in%rownames(TES)[i])) }
SEQ
seg_data <- seg_data[SEQ, ]

### chromosome bar aesthetic
chr_ranges <- as.data.frame(SummarizedExperiment::rowRanges(scCNVcomb))
chr_lengths <- rle(as.numeric(chr_ranges$seqnames))$lengths
if (any(chr_ranges$seqnames == "24") ||
    any(chr_ranges$seqnames == "Y") ||
    any(chr_ranges$seqnames == "chrY")) {
  chr_binary <- rep(c(2, 1), length(chr_lengths) / 2)
} else {
  chr_binary <- c(rep(c(2, 1), (length(chr_lengths) / 2)), 2)
}
chr <- data.frame(chr = rep.int(x = chr_binary, times = chr_lengths))
chr_rl_c <- c(1, cumsum(chr_lengths)) # getting lengths for chr numbers annotation
# creating a data frame to calculate rowMeans
chr_df <- data.frame(a = chr_rl_c[1:length(chr_rl_c) - 1], b = chr_rl_c[2:length(chr_rl_c)])
chr_l_means <- round(rowMeans(chr_df))
chrom.names <- c(1:22, "X", "Y")
# creating the vector for chr number annotations
v <- vector(length = sum(chr_lengths), mode = "character")
suppressWarnings(v[chr_l_means] <- chrom.names)
v[is.na(v)] <- ""
# chr bar with the chr names
chr_bar <-
  ComplexHeatmap::HeatmapAnnotation(
    chr_text = ComplexHeatmap::anno_text(v[1:ncol(seg_data)], gp = grid::gpar(fontsize = 14)),
    df = as.character(chr[1:nrow(chr),]),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    which = "column",
    col = list(df = c("1" = "grey88", "2" = "black"))
  )
# Cluster Annotations
metadata <- SummarizedExperiment::colData(scCNVcomb) %>% as.data.frame()
metadata <- metadata[rownames(seg_data), ]
metadata_anno_df <- metadata %>% dplyr::select(superclones,subclones,Source) 
library(cluster) # Colors Customization and color Idex
colors = c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#42d4f4", 
           "#f032e6", "#bfef45", "#fabebe", "#469990", "#e6beff", "#9A6324", "#800000", 
           "#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9", "#000000")
BR <- c("cornsilk4", "deeppink2")
names(BR) <- c("Single-Cell", "Subclones_Bulk")
cluster_anno <-
  ComplexHeatmap::rowAnnotation(
    df = metadata_anno_df,
    col = list(superclones = copykit:::superclones_pal(),
               subclones = copykit:::subclones_pal(), Source = BR),
    show_annotation_name = FALSE
  )
### Extract coordinates of desired genes in SYMPOL format
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#BCDri <- c("ESR1","TP53","NF1","AKT1","KMT2C","PTEN","PIK3CA","ARID1A","CDKN1B","ERBB2","CBFB","CDH1","FOXA1","GATA3","GPS2","MAP2K4","MAP3K1","NCOR1","RB1","RUNX1","TBX3")
BCDri <- c("ATM", "BARD1", "BRCA1", "BRCA2", "CDH1", "CHEK2", "NF1", "PALB2", "PTEN", "RAD51C", "RAD51D", "TP53")
myGeneSymbols <- select(org.Hs.eg.db, keys = BCDri, columns = c("SYMBOL","ENTREZID"), keytype = "SYMBOL")
myGeneSymbolsTx <- select(TxDb.Hsapiens.UCSC.hg19.knownGene,
                          keys = myGeneSymbols$ENTREZID,
                          columns = c("GENEID", "TXID", "TXCHROM", "TXSTART", "TXEND"),
                          keytype = "GENEID")
res <- merge(myGeneSymbols, myGeneSymbolsTx, by.x = "ENTREZID", by.y = "GENEID")
# Overlap gene ranges to find positions
BCDriGR <- makeGRangesFromDataFrame(res)
BCDriGR$symbol <- res$SYMBOL
VarbinGR <- makeGRangesFromDataFrame(scCNVcomb@rowRanges)
olaps <- findOverlaps(BCDriGR, VarbinGR)
mk_df <- tibble(gene = BCDriGR$symbol[queryHits(olaps)], pos = subjectHits(olaps)) %>% dplyr::distinct(gene, .keep_all = TRUE) 
NewsegMat <- cbind.data.frame(TES$superclones, TES$subclones, seg_data)
colnames(NewsegMat)[1:2] <- c("Superclone","Subclone")

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
mk_df$"CNA class" <- subcVec[mk_df$pos]
for (c in 1:nrow(mk_df)) {
  if(mk_df$`CNA class`[c]==0){
    mk_df$"color"[c] <- "grey88" 
    next
  }else if(mk_df$`CNA class`[c]==1){
    mk_df$"color"[c] <- "orange" 
    next
  }else if(mk_df$`CNA class`[c]==2){mk_df$"color"[c] <- "black" 
  next
  }else{mk_df$"color"[c] <- "purple"}
}
Bottom_bar <-
  ComplexHeatmap::columnAnnotation(
    df = as.character(subcVec),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    col = list(df = c("0" = "grey88", "1" = "black", "2" = "orange", "3" = "purple")),
    #txtCol <- as.character(mk_df$`CNA class`),
    foo = ComplexHeatmap::anno_mark(
      at = mk_df$pos,
      labels = mk_df$gene,
      side = "bottom",
      link_height = grid::unit(8, "mm"),
      labels_gp = grid::gpar(fontsize = 8, col= mk_df$color)
    )
  )

### plotting
tiff(paste(CaseName, "_scCNV_Heatmap.tiff", sep = ""), width=3000, height=2500, compression="lzw", res=300)
ComplexHeatmap::Heatmap(
  log2(seg_data + 1e-3),
  use_raster = TRUE,
  left_annotation = cluster_anno,
  #column_title = "Genomic coordinates",
  column_title_gp = grid::gpar(fontsize = 18),
  #column_title_side = "bottom",
  row_title = paste(nrow(seg_data), " cells", sep = ""),
  row_title_gp = grid::gpar(fontsize = 18),
  heatmap_legend_param = list(title = "log2(segratio)"),
  top_annotation = chr_bar,
  bottom_annotation = Bottom_bar,
  col = circlize::colorRamp2(breaks = c(-2,0.1,2),c("dodgerblue3", "white", "firebrick3")),
  cluster_rows = F,
  border = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = F,
  show_heatmap_legend = TRUE,
  width = NULL,
)
dev.off()

### Percentage DataFrame
PerC1 <- table(TES[TES$Source=="Single-Cell", 6])
PerC2 <- table(TES[TES$Source=="Subclones_Bulk", 6])
PerC3 <- table(TES[TES$Source=="Single-Cell", 7])
PerC4 <- table(TES[TES$Source=="Subclones_Bulk", 7])
PerD.Super <- cbind(PerC1, round(PerC1/1335*100, 2), PerC2, round(PerC2/84*100, 2))
PerD.Sub <-  cbind(PerC3, round(PerC3/1335*100, 2), PerC4, round(PerC4/84*100, 2))
colnames(PerD.Super) <- c("NofParental", "P.Perct", "NofSubc", "S.Perct")
colnames(PerD.Sub) <- c("NofParental", "P.Perct", "NofSubc", "S.Perct")
PerD.Sub <- rbind(PerD.Sub[c(1, 8:15), ], PerD.Sub[2:7, ])
### Plot a barplot showing percentages of cell counts from both super- and sub- clusters perspectives.
# Creat ggplot_df: (Superclones) 
NewDSper <- rbind(PerD.Super[,1:2], PerD.Super[, 3:4])
NewDSper <- cbind(c(rep("Parental", 6),rep("SubC", 6)), rownames(rbind(PerD.Super[,1:2], PerD.Super[, 3:4])), as.numeric(NewDSper[,2]))
NewDSper <- as.data.frame(NewDSper)
colnames(NewDSper) <- c("Source", "Cluster", "Percent")
NewDSper$Percent <- NewDSper$Percent %>% as.character() %>% as.numeric()
ggplot(data=NewDSper, aes(x=Cluster, y=Percent, fill=Source)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_y_continuous(breaks=seq(0,60,5), limits = c(0, 60)) +
  geom_text(aes(label = paste(Percent, "%", sep = "")), vjust=1.6, position = position_dodge(0.9), size = 3.5)
# Creat ggplot_df: (Subclones) 
NewDSub <- rbind(PerD.Sub[,1:2], PerD.Sub[, 3:4])
NewDSub <- cbind(c(rep("Parental", 15),rep("SubC", 15)), rownames(rbind(PerD.Sub[,1:2], PerD.Sub[, 3:4])), as.numeric(NewDSub[,2]))
NewDSub<- as.data.frame(NewDSub)
colnames(NewDSub) <- c("Source", "Cluster", "Percent")
NewDSub$Percent <- NewDSub$Percent %>% as.character() %>% as.numeric()
ggplot(data=NewDSub, aes(x=Cluster, y=Percent, fill=Source)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_y_continuous(breaks=seq(0,30,5), limits = c(0, 30)) +
  geom_text(aes(label = paste(Percent, "%", sep = "")), vjust=1.6, position = position_dodge(0.9), size = 3.5)
# Create custom Umap using ggplot2
head(scCNVcomb@colData)
ggUmapdf_SC <- as.data.frame(scCNVcomb@colData[1:1335,c(1,6,7,8,9,11)])
ggUmapdf_bulk <- as.data.frame(scCNVcomb@colData[1336:1419,c(1,6,7,8,9,11)])
write.csv(ggUmapdf_bulk, "MDAMB157_Subclones_ClusterList.csv", row.names = FALSE)
MatchSuperColor <- function(x) {
  Colpal <- as.data.frame(superclones_pal())
  return(Colpal$`superclones_pal()`[which(x == rownames(Colpal))])
} # element only. So please use lapply.
MatchSubColor <- function(x) {
  Colpal <- as.data.frame(subclones_pal())
  return(Colpal$`subclones_pal()`[which(x == rownames(Colpal))])
} # element only. So please use lapply.
ggUmapdf_bulk$"SuperColor" <- lapply(ggUmapdf_bulk$superclones, MatchSuperColor) %>% unlist() 
ggUmapdf_bulk$"SubcColor" <- lapply(ggUmapdf_bulk$subclones, MatchSubColor) %>% unlist() 
ggUmapdf_SC$"SuperColor" <- lapply(ggUmapdf_SC$superclones, MatchSuperColor) %>% unlist() 
ggUmapdf_SC$"SubcColor" <- lapply(ggUmapdf_SC$subclones, MatchSubColor) %>% unlist()

ggUmapdf <- rbind(ggUmapdf_SC, ggUmapdf_bulk) %>% as.data.frame()

q <- ggplot(data=ggUmapdf_SC, aes(x=Umap_X, y=Umap_Y)) +
      scale_shape_identity() +
      geom_point(aes(fill = subclones, color = superclones, shape = 21, stroke = 1)) + 
      scale_fill_manual(values=subclones_pal()) +
      scale_color_manual(values=superclones_pal()) +
      geom_point(data=ggUmapdf_bulk, aes(x=Umap_X, y=Umap_Y, fill = subclones, color = superclones, size = 1, shape=23, stroke = 1)) +
      guides(fill = guide_legend(override.aes = list(shape = 21))) 

Base <- ggplot(ggUmapdf, aes(x=Umap_X, y=Umap_Y))+
         geom_point(aes(color = superclones, size = 3, shape = Source)) + 
         scale_color_manual(values=superclones_pal()) +
         scale_shape_manual(values = c(19,18))
Base2 <- Base + ggnewscale::new_scale_colour() + 
                geom_point(aes(color = subclones, shape = Source), fill = NA) + 
                scale_color_manual(values = subclones_pal())
##### You may consider using Rplot to do your Plotting. ggplot2 is GGing twice!!!

tiff("UMAP of MDA-MB-157.tiff", width=3000, height=2500, compression="lzw", res=300)
plot(ggUmapdf_SC$Umap_X, ggUmapdf_SC$Umap_Y, pch = 19, cex = 3, col = adjustcolor(ggUmapdf_SC$SuperColor, alpha.f = 0.5), 
     main="UMAP of MDA-MB-157", xlab = "Umap1", ylab = "Umap2")
points(ggUmapdf_SC$Umap_X, ggUmapdf_SC$Umap_Y, pch = 16, cex = 1, col = as.character(ggUmapdf_SC$SubcColor))
legend("bottomleft", legend = paste("s", 1:6, sep = ""), 
       col = adjustcolor(superclones_pal()[1:6], alpha.f = 0.8),
       title="Superclone", text.font=4, bg="aliceblue", pch = 19, bty = "o", pt.cex = 2, 
       cex = 1.2, text.col = "black", horiz = F , inset = c(0.05, 0.05))
legend("topright", legend = paste("c", 1:15, sep = ""), 
       col = subclones_pal()[1:15],
       title="Subclone", text.font=4, bg="aliceblue", pch = 19, bty = "o", pt.cex = 2, 
       cex = 1.2, text.col = "black", horiz = F , inset = c(0.05, 0.05))
text(c(9,2.5,-11,-7,-4,-14),c(-8,9,5,4.5,-8,2), labels = paste("s", 1:6, sep = ""), cex = 2, lwd = 2)
dev.off()

tiff("UMAP of MDA-MB-157 Expansions.tiff", width=3000, height=2500, compression="lzw", res=300)
plot(ggUmapdf_SC$Umap_X, ggUmapdf_SC$Umap_Y, pch = 19, cex = 3, col = adjustcolor(ggUmapdf_SC$SuperColor, alpha.f = 0.05), 
     main="UMAP of MDA-MB-157 Expansions", xlab = "Umap1", ylab = "Umap2")
points(ggUmapdf_SC$Umap_X, ggUmapdf_SC$Umap_Y, pch = 16, cex = 1, col = adjustcolor(as.character(ggUmapdf_SC$SubcColor), alpha.f = 0.01))

#points(ggUmapdf_bulk$Umap_X, ggUmapdf_bulk$Umap_Y, pch = 19, cex = 3, col = adjustcolor(ggUmapdf_bulk$SuperColor, alpha.f = 0.3))               
points(ggUmapdf_bulk$Umap_X, ggUmapdf_bulk$Umap_Y, pch = 8, cex = 1.75, col = as.character(ggUmapdf_bulk$SubcColor), lwd = 1.35)

legend("bottomleft", legend = paste("s", 1:6, sep = ""), 
       col = adjustcolor(superclones_pal()[1:6], alpha.f = 0.7),
       title="Superclone", text.font=4, bg="aliceblue", pch = 19, bty = "o", pt.cex = 2, 
       cex = 1.2, text.col = "black", horiz = F , inset = c(0.05, 0.05))
legend("topright", legend = paste("c", c(1:3, 5:13, 15), sep = ""), 
       col = subclones_pal()[c(1:3, 5:13, 15)],
       title="Subclone", text.font=4, bg="aliceblue", pch = 8, bty = "o", pt.cex = 2, 
       cex = 1.2, text.col = "black", horiz = F , inset = c(0.05, 0.05))
legend("topleft", legend = c("Parental_SC", "Expansion_Bulk"),
       title="Source", text.font=4, bg="aliceblue", pch = c(19,8), bty = "o", pt.cex = 2, 
       cex = 1, text.col = "black", horiz = F , inset = c(0.05, 0.05))
text(c(9,2.5,-11,-7,-4,-14),c(-8,9,5,4.5,-8,2), labels = paste("s", 1:6, sep = ""), cex = 2, lwd = 2)
dev.off()

plot(1:15, pch = 19, cex=5, col=subclones_pal()[1:15])

### S2C4 and S5C6 sampling pValue
#H0: Missing in S2C4 and S5C6 single clones is due to sampling error.
#H1: Missing in S2C4 and S5C6 single clones is NOT due to sampling error.

Popu <- TES[TES$Source == "Single-Cell", c(6,7,11)]
hitsC4 <- 0
for (w in 1:100000) {
  Trial <- sample(Popu$subclones, 84, replace = T)
  if("c4"%in%Trial == F){
    hitsC4 + 1
  }
}
pValC4 <- (hitsC4/100000)*100

hitsC6 <- 0
for (w in 1:100000) {
  Trial <- sample(Popu$subclones, 84, replace = T)
  if("c6"%in%Trial == F){
    hitsC6 + 1
  }
}
pValC6 <- (hitsC6/100000)*100
# The result rejects H0 and supports H1.