# ==========================================
# All annotation is under GRCH38 build
# ==========================================
library(GenomicRanges)
library(IRanges)
library(AnnotationHub)
library(ensembldb)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) # Provides gene range information


#' Add annotation under GRCH38 build 
#' 
#' @param final_zscore_matrix the input Z-score matrix
#' @returns Dataframe of genes + variant + ENTREZID
annotation_df <- function(final_zscore_matrix){
  ah <- AnnotationHub(ask = FALSE) # ah <- AnnotationHub()
  ensdb <- query(ah, "EnsDb.Hsapiens.v101")[[1]]# ensdb <- query(ah, "EnsDb.Hsapiens.v105")[[1]]  # Adjust index based on query results # # Find Ensembl database for GRCh38 # metadata(ensdb) # genome_build:GRCh38
  
  geneRanges <- function(db, column="symbol") {
    g <- genes(db, columns = column)
    col <- mcols(g)[[column]]
    genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
    mcols(genes)[[column]] <- as.character(unlist(col))
    genes
  }
  
  # Extract gene ranges with symbol column
  gns <- geneRanges(ensdb, column="symbol") 
  
  # Convert gene ranges to a data frame and extract relevant columns
  df_gns <- data.frame(gns) %>%
    dplyr::rename(SYMBOL = symbol)
  gtf.gene <- df_gns %>%
    mutate(chr = gsub("chr", "", seqnames)) %>%
    dplyr::select(chr, start, end, SYMBOL)
  
  # Extract gene ranges with entrezid column
  df_entrez <- geneRanges(ensdb, column="entrezid") %>%
    data.frame() %>%
    mutate(chr = gsub("chr", "", seqnames)) %>%
    mutate(ChrPos = paste(chr, start, sep=":")) %>% 
    dplyr::rename(ENTREZID = entrezid) %>% 
    dplyr::select(ChrPos, ENTREZID)
  
  #' Convert from string to range
  #' 
  #' @param pos A vector of strings ex. chr1 2938302 2938329
  #' @param delim Delimiter for string splitting
  #' @param region Boolean of whether region or just one position
  #'
  #' @returns Dataframe of ranges
  #' 
  string2range <- function(pos, delim=' ', region=TRUE) {
    posp <- as.data.frame(do.call(rbind, strsplit(pos, delim)))
    posp[,1] <- posp[,1]
    posp[,2] <- as.numeric(as.character(posp[,2]))
    if(region) {
      posp[,3] <- as.numeric(as.character(posp[,3]))
    } else {
      posp[,3] <- posp[,2]
    }
    return(posp)
  }
  
  #' Convert from ranges to GRanges
  #' 
  #' @param df Dataframe with columns as sequence name, start, and end
  #' 
  #' @returns GRanges version 
  range2GRanges <- function(df) {
    require(GenomicRanges)
    require(IRanges)
    gr <- GenomicRanges::GRanges(
      seqnames = df[,1],
      ranges=IRanges(start = df[,2], end = df[,3])
    )
    return(gr)
  }
  
  # convert SNPs to GRanges
  snps <- rownames(final_zscore_matrix)
  snps.ranges <- string2range(snps, delim=":", region=FALSE)
  snps.granges <- range2GRanges(snps.ranges)
  names(snps.granges) <- snps
  
  # convert genes to GRanges
  gtf.granges <- range2GRanges(gtf.gene)
  names(gtf.granges) <-  gtf.gene$SYMBOL  #gene.names
  
  hits <- GenomicRanges::nearest(snps.granges, gtf.granges)
  # make vector of SNPs to gene
  
  df_gns <- df_gns %>% 
    mutate(chr = gsub("chr","",seqnames)) %>%
    mutate(ChrPos = paste(chr,start,sep=":")) %>%
    dplyr::select(chr, start, end, ChrPos, SYMBOL) %>%
    merge(df_entrez, by="ChrPos")
  
  df_hits <- data.frame(gene=names(gtf.granges)[hits]) %>%
    mutate(variant = names(snps.granges)) %>%
    merge(df_gns[,c("SYMBOL","ENTREZID")], by.x="gene", by.y="SYMBOL") %>%
    dplyr::filter(!duplicated(variant))
  
  return(df_hits)
}



#' Exclude the genes with annotation starting with xxxx
#' 
#' @param with_text A string specifying the text that, if found in annotations, will result in exclusion.
#' @param final_zscore_matrix the input Z-score matrix
#'
#' @returns excluded version of final_zscore_matrix
exclude_genes_with_annotation <- function(final_zscore_matrix,with_text){
  annotation_data <- annotation_df(final_zscore_matrix) # Generate the annotation data frame
  genes_to_exclude <- annotation_data$gene[grepl(with_text, annotation_data$gene)] # Identify rows where the `gene` column contains the specified text
  variants_to_exclude <- annotation_data$variant[annotation_data$gene %in% genes_to_exclude] # Identify corresponding `variant` rows to exclude in the Z-score matrix
  rows_to_keep <- !(rownames(final_zscore_matrix) %in% variants_to_exclude) # Keep only rows in the Z-score matrix that are not in the `variants_to_exclude`
  filtered_matrix <- final_zscore_matrix[rows_to_keep, , drop = FALSE] # Filter the Z-score matrix
  return(filtered_matrix)
}




