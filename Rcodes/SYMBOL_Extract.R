TNBCLoc <- "/volumes/seq/projects/CNA_projects/DT_CNA/cell_line/"
HCC70 <- "HCC70P1/output/final_result"
library(dplyr)
TNBCelline <- copykit::readVarbinCNA(paste(TNBCLoc, HCC70,sep = ""), remove_Y = TRUE)
SummarizedExperiment::rowRanges(TNBCelline)
TNBCelline <- copykit::runMetrics(TNBCelline)
#copykit::plotMetrics(TNBCelline, label = "rmse") *This is not working.
TNBCelline <- copykit::filterCells(TNBCelline,resolution = 0.9)
tiff("HCC70_scCNV_Fil(0.9).tiff", width=3000, height=2500, compression="lzw", res=300)
copykit::filterCells(TNBCelline, resolution = 0.9)
dev.off()
head(SummarizedExperiment::colData(TNBCelline))
TNBCelline <- TNBCelline[ ,which(SummarizedExperiment::colData(TNBCelline)$filtered == "kept")]
TNBCelline <- copykit::runDistMat(TNBCelline)
TNBCelline <- copykit::runUmap(TNBCelline,min_dist = 0.001,n_neighbors = 30)
SingleCellExperiment::reducedDim(TNBCelline, 'umap', withDimnames = FALSE)
TNBCelline <- copykit::findClusters(TNBCelline, k_major = 35, k_minor = 21)
tiff("HCC70_scCNV_Fil(0.9)_Umap.tiff", width=2000, height=2000, compression="lzw", res=300)
copykit::plotUmap(TNBCelline)
dev.off()
spatial_info <- as.data.frame(SummarizedExperiment::colData(TNBCelline))
spatial_info$spatial_location <- stringr::str_extract(spatial_info$sample,"(s[0-9]){1}")
SummarizedExperiment::colData(TNBCelline)$spatial_location <- spatial_info$spatial_location
#copykit::plotUmap(TNBCelline,label = "spatial_location")

# Convert factors to numeric and obtain the ordered colData matrix
TNBCelline@colData$subclones <- as.numeric(as.character(TNBCelline@colData$subclones))
TES <- TNBCelline@colData[order(TNBCelline@colData[ ,7], TNBCelline@colData[ ,8]), ]
TES$sample[1:5]

#obtaining data and sort cell rows using ordered position
seg_data <- t(copykit::segment_ratios(TNBCelline))
SEQ <- NULL
for (i in 1:length(TES$sample)) { SEQ <- c(SEQ, which(rownames(seg_data)%in%rownames(TES)[i])) }
SEQ
seg_data <- seg_data[SEQ, ]
#chromosome bar aesthetic
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
# getting lengths for chr numbers annotation
chr_rl_c <- c(1, cumsum(chr_lengths))
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
#hc <- fastcluster::hclust(copykit::distMat(TNBCelline), method = "ward.D2")
#seg_data <- seg_data[hc$order, ]
TNBCelline@colData$subclones <- as.character(TNBCelline@colData$subclones)
metadata <- SummarizedExperiment::colData(TNBCelline) %>% as.data.frame()
metadata <- metadata[rownames(seg_data), ]
metadata_anno_df <- metadata %>% dplyr::select(superclones,subclones) 
# Colors Customization
library(cluster) #General color Idex
colors = c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#42d4f4", 
           "#f032e6", "#bfef45", "#fabebe", "#469990", "#e6beff", "#9A6324", "#800000", 
           "#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9", "#000000")
cluster_anno <-
  ComplexHeatmap::rowAnnotation(
    df = metadata_anno_df,
    col = list(superclones = copykit:::major_palette,
               subclones = copykit:::minor_palette),
    show_annotation_name = FALSE
  )
### Extract coordinates of desired genes in SYMPOL format
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
BCDri <- c("ESR1","TP53","NF1","AKT1","KMT2C","PTEN","PIK3CA","ARID1A","CDKN1B",
           "ERBB2","CBFB","CDH1","FOXA1","GATA3","GPS2","MAP2K4","MAP3K1","NCOR1","RB1","RUNX1","TBX3")
myGeneSymbols <- select(org.Hs.eg.db, keys = BCDri, columns = c("SYMBOL","ENTREZID"), keytype = "SYMBOL")
myGeneSymbolsTx <- select(TxDb.Hsapiens.UCSC.hg19.knownGene,
                          keys = myGeneSymbols$ENTREZID,
                          columns = c("GENEID", "TXID", "TXCHROM", "TXSTART", "TXEND"),
                          keytype = "GENEID")
res <- merge(myGeneSymbols, myGeneSymbolsTx, by.x = "ENTREZID", by.y = "GENEID")
### Overlap gene ranges to find positions
BCDriGR <- makeGRangesFromDataFrame(res)
BCDriGR$symbol <- res$SYMBOL
VarbinGR <- makeGRangesFromDataFrame(TNBCelline@rowRanges)
olaps <- findOverlaps(BCDriGR, VarbinGR)
mk_df <- tibble(gene = BCDriGR$symbol[queryHits(olaps)], pos = subjectHits(olaps)) %>% dplyr::distinct(gene, .keep_all = TRUE) 
NewsegMat <- cbind.data.frame(TES$superclones, TES$subclones, seg_data)
colnames(NewsegMat)[1:2] <- c("Superclone","Subclone")
### CNAs class color band bottom annotation from subclones' point of view
subcN <- max(TES$subclones)
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
      at = mk_df$pos[c(1:5, 10:11, 13,16:17,21)],
      labels = mk_df$gene[c(1:5, 10:11, 13,16:17,21)],
      side = "bottom",
      link_height = grid::unit(8, "mm"),
      labels_gp = grid::gpar(fontsize = 8, col= mk_df$color[c(1:5, 10:11, 13,16:17,21)])
    )
  )
Bottom_bar_Super <-
  ComplexHeatmap::columnAnnotation(
    df = as.character(subcVec),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    col = list(df = c("0" = "grey88", "1" = "grey88", "2" = "grey88", "3" = "grey88")),
    foo = ComplexHeatmap::anno_mark(
      at = mk_df$pos[c(1:5, 10:11, 13,16:17,21)],
      labels = mk_df$gene[c(1:5, 10:11, 13,16:17,21)],
      side = "bottom",
      link_height = grid::unit(8, "mm"),
      labels_gp = grid::gpar(fontsize = 8)
    )
  )
#plotting
tiff("HCC70_scCNV_dbscan_Fil(0.9).tiff", width=3000, height=2500, compression="lzw", res=300)
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

### Superclonal Concensus Map
long_list <- split(NewsegMat, NewsegMat$Superclone) 
consensus_list <- parallel::mclapply(long_list, function(x) { apply(x[,-c(1:2)], 2, median) }, mc.cores = 30)
cs_df <- as.data.frame(do.call(rbind, consensus_list))
df2 <- metadata_anno_df$superclones %>% unique() %>% as.data.frame() 
colnames(df2) <- "superclones"
ConSenAnno <- ComplexHeatmap::rowAnnotation(
  df = df2,
  col = list(superclones = copykit:::major_palette),
  show_annotation_name = FALSE)
tiff("HCC70_scCNV_SuperConsensus.tiff", width=3000, height=2500, compression="lzw", res=300)
ComplexHeatmap::Heatmap(
  log2(as.matrix(cs_df) + 1e-3),
  #use_raster = TRUE,
  left_annotation = ConSenAnno,
  column_title = "Genomic coordinates",
  column_title_gp = grid::gpar(fontsize = 18),
  column_title_side = "bottom",
  row_title = "Superclone",
  row_title_gp = grid::gpar(fontsize = 18),
  heatmap_legend_param = list(title = "log2(segratio)"),
  top_annotation = chr_bar,
  bottom_annotation = Bottom_bar_Super,
  col = circlize::colorRamp2(breaks = c(-2,0.1,2),c("dodgerblue3", "white", "firebrick3")),
  cluster_rows = T,
  border = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = T,
  show_heatmap_legend = TRUE,
  row_split = rownames(cs_df))
dev.off()
### Subclonal Concensus Map
long_list <- split(NewsegMat, NewsegMat$Subclone) 
consensus_list <- parallel::mclapply(long_list, function(x) { apply(x[,-c(1:2)], 2, median) }, mc.cores = 30)
cs_df <- as.data.frame(do.call(rbind, consensus_list))
df2 <- metadata_anno_df$subclones %>% unique() %>% as.data.frame() # Change column for getting Subclones
colnames(df2) <- "subclones"
ConSenAnno <- ComplexHeatmap::rowAnnotation(
  df = df2,
  col = list(#superclones = copykit:::major_palette,
    subclones = copykit:::minor_palette),
  show_annotation_name = FALSE)
tiff("HCC70_scCNV_SubConsensus.tiff", width=3000, height=2500, compression="lzw", res=300)
ComplexHeatmap::Heatmap(
  log2(as.matrix(cs_df) + 1e-3),
  #use_raster = TRUE,
  left_annotation = ConSenAnno,
  column_title = "Genomic coordinates",
  column_title_gp = grid::gpar(fontsize = 18),
  column_title_side = "bottom",
  row_title = "Subclone",
  row_title_gp = grid::gpar(fontsize = 18),
  heatmap_legend_param = list(title = "log2(segratio)"),
  top_annotation = chr_bar,
  bottom_annotation = Bottom_bar,
  col = circlize::colorRamp2(breaks = c(-2,0.1,2),c("dodgerblue3", "white", "firebrick3")),
  cluster_rows = T,
  border = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = T,
  show_heatmap_legend = TRUE,
  row_split = rownames(cs_df)
)
dev.off()