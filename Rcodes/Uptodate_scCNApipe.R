### Downstrem pipeline for individual TNBC cell line
### Will generate UMAP, clustered heatmap with CNAs annotations, and consensus NJ tree.
DataLoc <- "/volumes/seq/projects/CNA_projects/DT_CNA/cell_line/HCC70/HCC70_P2_P3_P4_combine/output/final_result"
CaseName <- "HCC70"
library(dplyr)
library(SingleCellExperiment)
library(ggplot2)
library(copykit)
library(ggtree)
library(ape)
library(phangorn)
TNBCelline <- copykit::readVarbinCNA(DataLoc, remove_Y = TRUE)
#TNBCelline <- TNBCelline[ , 1:986]
SummarizedExperiment::rowRanges(TNBCelline)
TNBCelline <- copykit::runMetrics(TNBCelline)
TNBCelline <- copykit::filterCells(TNBCelline, resolution = 0.85) # Adjust the resolution.
tiff(paste(CaseName, "_scCNV_filtered.tiff", sep = ""), width=3000, height=2500, compression="lzw", res=300)
plotHeatmap(TNBCelline, label = "filtered", row_split = "filtered")
dev.off()
head(SummarizedExperiment::colData(TNBCelline))
table(TNBCelline@colData$filtered)
TNBCelline <- TNBCelline[ ,which(SummarizedExperiment::colData(TNBCelline)$filtered == "kept")]
TNBCelline <- copykit::runDistMat(TNBCelline)
TNBCelline <- copykit::runUmap(TNBCelline, min_dist = 0, n_neighbors = 40)
SingleCellExperiment::reducedDim(TNBCelline, 'umap', withDimnames = FALSE)
opt_k <- findOptimalK(TNBCelline)$k
TNBCelline <- copykit::findClusters(TNBCelline, k_superclones = 35, k_subclones = opt_k)
TNBCelline@colData$"Umap_X" <- SingleCellExperiment::reducedDim(TNBCelline, 'umap', withDimnames = FALSE)[, 1]
TNBCelline@colData$"Umap_Y" <- SingleCellExperiment::reducedDim(TNBCelline, 'umap', withDimnames = FALSE)[, 2]

# Convert factors to numeric and count numbers of superclones and subclones
TNBCelline@colData$superclones <- as.character(TNBCelline@colData$superclones)
TNBCelline@colData$subclones <- as.character(TNBCelline@colData$subclones)
N.Superc <- length(unique(TNBCelline@colData$superclones)) 
N.Subc <- length(unique(TNBCelline@colData$subclones)) 

# Functions to assign colors
MatchSuperColor <- function(x) {
  Colpal <- as.data.frame(superclones_pal())
  return(Colpal$`superclones_pal()`[which(x == rownames(Colpal))])
} # element only. So please use lapply.
MatchSubColor <- function(x) {
  Colpal <- as.data.frame(subclones_pal())
  return(Colpal$`subclones_pal()`[which(x == rownames(Colpal))])
} # element only. So please use lapply.
TNBCelline@colData$"SuperColor" <- lapply(TNBCelline@colData$superclones, MatchSuperColor) %>% unlist() 
TNBCelline@colData$"SubcColor" <- lapply(TNBCelline@colData$subclones, MatchSubColor) %>% unlist()

### Generate customised UMAP. Adjust legends position!
tiff(paste(CaseName, "_scCNV_HDBSCAN_UMAP.tiff", sep = ""), width=3000, height=2500, compression="lzw", res=300)
plot(TNBCelline@colData$Umap_X, TNBCelline@colData$Umap_Y, pch = 19, cex = 3, col = adjustcolor(TNBCelline@colData$SuperColor, alpha.f = 0.5), 
     main=paste(CaseName, "_UMAP", sep = ""), xlab = "Umap1", ylab = "Umap2")
points(TNBCelline@colData$Umap_X, TNBCelline@colData$Umap_Y, pch = 16, cex = 1, col = as.character(TNBCelline@colData$SubcColor))
legend("topleft", legend = paste("s", 1:(N.Superc), sep = ""), 
       col = adjustcolor(superclones_pal()[1:(N.Superc)], alpha.f = 0.8),
       title="Superclone", text.font=4, bg="aliceblue", pch = 19, bty = "o", pt.cex = 2, 
       cex = 1.2, text.col = "black", horiz = F , inset = c(0.05, 0.05))
legend(-9, 3, legend = paste("c", 1:(N.Subc), sep = ""), 
       col = subclones_pal()[1:(N.Subc)],
       title="Subclone", text.font=4, bg="aliceblue", pch = 19, bty = "o", pt.cex = 2, 
       cex = 1.2, text.col = "black", horiz = F , inset = c(0.05, 0.05))
#text(c(9,2.5,-11,-7,-4,-14),c(-8,9,5,4.5,-8,2), labels = paste("s", 1:6, sep = ""), cex = 2, lwd = 2)
dev.off()

# Quick copykit functions to take a look.
#copykit::plotUmap(TNBCelline)
#plotRatio(TNBCelline)

# Get ordered colData matrix
head(SummarizedExperiment::colData(TNBCelline))
TES <- TNBCelline@colData[order(TNBCelline@colData[ ,6], TNBCelline@colData[ ,7]), ]
TES[1:10, ]

# Sort segment ratios according to ordered colData sample names
seg_data <- t(copykit::segment_ratios(TNBCelline))
SEQ <- NULL
for (i in 1:length(TES$sample)) { SEQ <- c(SEQ, which(rownames(seg_data)%in%rownames(TES)[i])) }
SEQ
seg_data <- seg_data[SEQ, ]

### chromosome bar aesthetic
chr_ranges <- as.data.frame(SummarizedExperiment::rowRanges(TNBCelline))
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
metadata <- SummarizedExperiment::colData(TNBCelline) %>% as.data.frame()
metadata <- metadata[rownames(seg_data), ]
metadata_anno_df <- metadata %>% dplyr::select(superclones,subclones) 
library(cluster) # Colors Customization and color Idex
colors = c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#42d4f4", 
           "#f032e6", "#bfef45", "#fabebe", "#469990", "#e6beff", "#9A6324", "#800000", 
           "#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9", "#000000")
cluster_anno <-
  ComplexHeatmap::rowAnnotation(
    df = metadata_anno_df,
    col = list(superclones = copykit:::superclones_pal(),
               subclones = copykit:::subclones_pal()),
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
VarbinGR <- makeGRangesFromDataFrame(TNBCelline@rowRanges)
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
    col = list(df = c("0" = "grey88", "1" = "orange", "2" = "black", "3" = "purple")),
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

### More Legends for Figures
tiff("CNA_class_legend.tiff", width=3000, height=2500, compression="lzw", res=300)
plot(1:10)
legend("bottomright", legend = c("nsCNA","cCNA","sCNA","uCNA"), title = "CNA Classes",
       col = c("grey88","orange","black","purple"), text.font=4, bg="white", pch = 15, bty = "o", pt.cex = 2, 
       cex = 1.2, text.col = "black", horiz = F , inset = c(0.05, 0.05))
dev.off()
tiff("HorizonClusterLegend.tiff", width=3500, height=2500, compression="lzw", res=300)
plot(1:10)
legend("top", legend = paste("s", 1:9, sep = ""), title = "Superclone",
       col = copykit:::superclones_pal()[1:9], text.font=4, bg="white", pch = 15, bty = "o", pt.cex = 2, 
       cex = 1.2, text.col = "black", horiz = T , inset = c(0.05, 0.05))
legend("bottom", legend = paste("c", 1:14, sep = ""), title = "Subclone",
       col = copykit:::subclones_pal()[1:14], text.font=4, bg="white", pch = 15, bty = "o", pt.cex = 2, 
       cex = 1.2, text.col = "black", horiz = T , inset = c(0.05, 0.05))
dev.off()

### Integer Copy Number Heatmap 
library(RColorBrewer) 
heatVal <- round(seg_data*4.52, digits = 0)
table(heatVal)
grcol <- colorRampPalette(c("#2166AC","#F0F0F0","#B2182B"))(9) # Modify when needed
heatVal[heatVal>=8] <- 8
#heatVal[heatVal==0] <- 1
table(heatVal)
tiff(paste(CaseName, "_scCNV_IntegCN.tiff", sep = ""), width=3000, height=2500, compression="lzw", res=300)
ComplexHeatmap::Heatmap(
  heatVal,
  use_raster = TRUE,
  left_annotation = cluster_anno,
  #column_title = "Genomic coordinates",
  column_title_gp = grid::gpar(fontsize = 18),
  #column_title_side = "bottom",
  row_title = paste(nrow(seg_data), " cells", sep = ""),
  row_title_gp = grid::gpar(fontsize = 18),
  heatmap_legend_param = list(title = "Integer Copy Number"),
  top_annotation = chr_bar,
  bottom_annotation = Bottom_bar,
  col = grcol,
  cluster_rows = F,
  border = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = F,
  show_heatmap_legend = TRUE,
  width = NULL,
)
dev.off()

### CopyKit ConsensusLines
head(SummarizedExperiment::colData(TNBCelline))
TNBCelline <- calcConsensus(TNBCelline)
plotConsensusLine(TNBCelline)

### Concensus HeatMap
long_list <- split(NewsegMat, NewsegMat$Subclone) 
consensus_list <- parallel::mclapply(long_list, function(x) { apply(x[,-c(1:2)], 2, median) }, mc.cores = 30)
cs_df <- do.call(rbind, consensus_list)
# Covert cs_df to IntegCN values
cs_df <- round(cs_df*4.52, digits = 0)
table(cs_df)
grcol2 <- colorRampPalette(c("#2166AC","#F0F0F0","#B2182B"))(7) # Get the right colors range.
cs_df[cs_df>=7] <- 7
#cs_df[cs_df==0] <- 1
table(cs_df)
df2 <- metadata_anno_df$subclones %>% unique() %>% as.data.frame() # Change column for getting Subclones
colnames(df2) <- "subclones"
df3 <- metadata_anno_df %>% distinct(superclones, subclones)
ConSenAnno <- ComplexHeatmap::rowAnnotation(
  df = df3,
  col = list(superclones = copykit:::superclones_pal(), subclones = copykit:::subclones_pal()),
  foo = ComplexHeatmap::anno_text(df3$superclones),
  foo2 = ComplexHeatmap::anno_text(df3$subclones),
  show_annotation_name = FALSE)
cs_df <- cs_df[match(df3$subclones, rownames(cs_df)), ] 

tiff(paste(CaseName, "_scCNV_ConsensusHeatmap.tiff", sep = ""), width=3000, height=2500, compression="lzw", res=300)
ComplexHeatmap::Heatmap(
  cs_df,
  use_raster = T,
  left_annotation = ConSenAnno,
  #column_title = "Genomic coordinates",
  column_title_gp = grid::gpar(fontsize = 18),
  #column_title_side = "bottom",
  row_title = "Subclone",
  row_title_gp = grid::gpar(fontsize = 18),
  heatmap_legend_param = list(title = "Integer Copy Number"),
  top_annotation = chr_bar,
  bottom_annotation = Bottom_bar,
  col = grcol2,
  cluster_rows = F,
  border = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = F,
  show_heatmap_legend = TRUE,
  row_split = 1:nrow(df2),
  width = NULL
)
dev.off()

### Phylogenetic trees
library(phangorn)
SUB_DIST <- dist(cs_df)
#SUB_UPGMA <- upgma(SUB_DIST)
SUB_NJ  <- NJ(SUB_DIST)
SUB_NJ$tip.label <- paste(df3$superclones, SUB_NJ$tip.label, sep = "")
#SUB_UPGMA$tip.label <- paste(df3$superclones, SUB_UPGMA$tip.label, sep = "")
tiff(paste(CaseName, "_scCNV_NJ_Tree.tiff", sep = ""), width=2500, height=2500, compression="lzw", res=300)
plot(SUB_NJ, main=paste(CaseName, "_NJ_Tree", sep = ""), cex=2.5)
dev.off()