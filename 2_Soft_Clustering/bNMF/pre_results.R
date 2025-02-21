library(LDlinkR)
library(readxl)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(AnnotationHub)
library(ensembldb)
library(ComplexHeatmap)
library(cowplot)
library(gridExtra)
library(stringr)
library(forcats)
library(pheatmap)
library(RColorBrewer)

load('./test_results/pipeline_data_wo_corr.RData')
source("./prep_bNMF.R")  # fetch_summary_stats & prep_z_matrix
source("./bnmf-clustering/scripts/run_bNMF.R")  # run_bNMF & summarize_bNMF
source('./choose_variants.R') # ld_pruning, count_traits_per_variant, fina_variants_needing_proxies, & choose_potential_proxies
source("./exclude_genes.R")
project_final_dir_wo_chr6_20_40 <- "./final_test_results_wo_corrImmCell_wo_chr6_20_40"

main_dir <- project_final_dir
weight_cutoff <- 0.87848 # project_final_dir
# main_dir <- project_final_dir_wo_chr6_20_40 

# weight_cutoff <- 1.1119 # project_final_dir_wo_chr6_20_40
align_dir <- project_dir
k <- k
loci_file<-"query"
GTEx<-F
my_traits<-gwas_traits
eval_GTEx <- GTEx

df_run_summary <- fread(file.path(main_dir,"run_summary.txt"),
                        stringsAsFactors = F, data.table = F)

k_counts <- df_run_summary %>% dplyr::count(K) %>%
  mutate(perc = n/nrow(df_run_summary))
# knitr::kable(caption="Frequency of each K converged to...") %>%
# kable_styling(full_width = F)

if (is.null(k)){
  k <- k_counts$K[which.max(k_counts$n)]
}

w <- fread(sprintf("%s/L2EU.W.mat.%i.txt",main_dir, k),
           stringsAsFactors = FALSE, data.table = F)
n_variants <- nrow(w)

h <- fread(sprintf("%s/L2EU.H.mat.%i.txt",main_dir,k),
           stringsAsFactors = FALSE, data.table = F) %>%
  rename_at(
    vars(starts_with("X2hr")), function(x) {gsub("X","",x) })

df_gwas <- fread(sprintf("%s/alignment_GWAS_summStats.csv",align_dir),
                 stringsAsFactors = F, data.table=F)

orig_traits <- colnames(h)
traits_short <- str_before_last(names(h), "[_]")
df_traits <- data.frame(Trait=unique(traits_short))


##########################################
file_path <- file.path(main_dir,"sorted_cluster_weights_K6_rev.xlsx")
data <- read_excel(file_path, sheet = "cluster1")
data_part1 <- data[, 1:which(names(data) == "...10") - 1]
data_part2 <- data[, which(names(data) == "trait"):ncol(data)]
names(data_part1) <- sub("\\.\\.\\..*", "", names(data_part1))
names(data_part2) <- sub("\\.\\.\\..*", "", names(data_part2))


for (i in 1:6) {
  column_name <- paste0("X", i)
  if (column_name %in% names(data_part1)) {
    # Filter rows and select columns
    data_part1_filtered <- data_part1[data_part1[[column_name]] > weight_cutoff, c("gene", column_name)]
    data_part1_filtered <- data_part1_filtered %>% arrange(gene)
    data_part2_filtered <- data_part2[data_part2[[column_name]] > weight_cutoff, c("trait", column_name)]
    data_part2_filtered <- na.omit(data_part2_filtered)
    data_part2_filtered <- data_part2_filtered %>%
      mutate(positive_trait = case_when(
        str_detect(trait, "_pos$") ~ 1,  # trait ends with "_pos" : 1
        str_detect(trait, "_neg$") ~ 0,  # trait ends with "_neg" : 0
        TRUE ~ NA_real_  # otherwise: NA
      )) %>%
      mutate(trait = str_replace(trait, "_(pos|neg)$", ""))%>% # delete _(pos|neg)
      mutate(trait_sort_key = gsub(".*_panel", "", trait)) %>%  
      arrange(trait_sort_key)
    
    p1 <- ggplot(data_part1_filtered, aes(x = forcats::fct_inorder(gene), y = get(column_name))) +
      geom_bar(stat = "identity", fill = "skyblue") +
      geom_text(aes(label = round(get(column_name), 2)),
                vjust = -0.5, size = 5) +
      theme_minimal() +
      labs(x = "Gene", y = "Weights", title = paste("Cluster", i, "- Gene Weights")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
    
    p2 <- ggplot(data_part2_filtered, 
                 aes(x = forcats::fct_inorder(str_trunc(trait, 55)), # truncate the length of the trait character 
                     y = get(column_name), 
                     fill = factor(positive_trait))) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = round(get(column_name), 2)),
                vjust = -0.5, size = 5) +
      theme_minimal() +
      labs(x = "Trait", y = "Weights", title = paste("Cluster", i, "- Trait Weights"), fill = "Trait Type") +
      scale_fill_manual(values = c("1" = "gold", "0" = "steelblue"),
                        labels = c("1" = "Traits w/ increasing values", "0" = "Traits w/ decreasing values")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
            legend.position = c(0.90, 0.85),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 12),
            # legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
            legend.margin = margin(t = 5, r = 5, b = 5, l = 5))
  }
  file_path <- file.path(main_dir, paste0("cluster", i, "_barplot.png"))
  ggsave(file_path, grid.arrange(p1, p2, ncol = 1), width = 16, height = 16, dpi=300)
  # file_path_gene <- file.path(main_dir, paste0("cluster", i, "_barplot_gene.png"))
  # ggsave(file_path_gene, plot = p1, width = 10, height = 6, dpi = 300)
  # file_path_trait <- file.path(main_dir, paste0("cluster", i, "_barplot_trait.png"))
  # ggsave(file_path_trait, plot = p2, width = 12, height = 8, dpi = 300)
}


################################
################################
################################
main_dir <- project_final_dir_wo_chr6_20_40 
weight_cutoff <- 1.1119

align_dir <- project_dir
k <- k
loci_file<-"query"
GTEx<-F
my_traits<-gwas_traits
eval_GTEx <- GTEx

df_run_summary <- fread(file.path(main_dir,"run_summary.txt"),
                        stringsAsFactors = F, data.table = F)

k_counts <- df_run_summary %>% dplyr::count(K) %>%
  mutate(perc = n/nrow(df_run_summary))
# knitr::kable(caption="Frequency of each K converged to...") %>%
# kable_styling(full_width = F)

if (is.null(k)){
  k <- k_counts$K[which.max(k_counts$n)]
}

w <- fread(sprintf("%s/L2EU.W.mat.%i.txt",main_dir, k),
           stringsAsFactors = FALSE, data.table = F)
n_variants <- nrow(w)

h <- fread(sprintf("%s/L2EU.H.mat.%i.txt",main_dir,k),
           stringsAsFactors = FALSE, data.table = F) %>%
  rename_at(
    vars(starts_with("X2hr")), function(x) {gsub("X","",x) })

df_gwas <- fread(sprintf("%s/alignment_GWAS_summStats.csv",align_dir),
                 stringsAsFactors = F, data.table=F)

orig_traits <- colnames(h)
traits_short <- str_before_last(names(h), "[_]")
df_traits <- data.frame(Trait=unique(traits_short))


##########################################
file_path <- file.path(main_dir,"sorted_cluster_weights_K4_rev.xlsx")
data <- read_excel(file_path, sheet = "cluster1")
data_part1 <- data[, 1:which(names(data) == "...8") - 1]
data_part2 <- data[, which(names(data) == "trait"):ncol(data)]
names(data_part1) <- sub("\\.\\.\\..*", "", names(data_part1))
names(data_part2) <- sub("\\.\\.\\..*", "", names(data_part2))
for (i in 1:4) {
  column_name <- paste0("X", i)
  if (column_name %in% names(data_part1)) {
    # Filter rows and select columns
    data_part1_filtered <- data_part1[data_part1[[column_name]] > weight_cutoff, c("gene", column_name)]
    data_part1_filtered <- data_part1_filtered %>% arrange(gene)
    data_part2_filtered <- data_part2[data_part2[[column_name]] > weight_cutoff, c("trait", column_name)]
    data_part2_filtered <- na.omit(data_part2_filtered)
    data_part2_filtered <- data_part2_filtered %>%
      mutate(positive_trait = case_when(
        str_detect(trait, "_pos$") ~ 1,  # trait ends with "_pos" : 1
        str_detect(trait, "_neg$") ~ 0,  # trait ends with "_neg" : 0
        TRUE ~ NA_real_  # otherwise: NA
      )) %>%
      mutate(trait = str_replace(trait, "_(pos|neg)$", ""))%>% # delete _(pos|neg)
      mutate(trait_sort_key = gsub(".*_panel", "", trait)) %>%  
      arrange(trait_sort_key)
    
    p1 <- ggplot(data_part1_filtered, aes(x = forcats::fct_inorder(gene), y = get(column_name))) +
      geom_bar(stat = "identity", fill = "skyblue") +
      geom_text(aes(label = round(get(column_name), 2)),
                vjust = -0.5, size = 5) +
      theme_minimal() +
      labs(x = "Gene", y = "Weights", title = paste("Cluster", i, "- Gene Weights")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
    
    p2 <- ggplot(data_part2_filtered, 
                 aes(x = forcats::fct_inorder(str_trunc(trait, 55)), # truncate the length of the trait character 
                     y = get(column_name), 
                     fill = factor(positive_trait))) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = round(get(column_name), 2)),
                vjust = -0.5, size = 5) +
      theme_minimal() +
      labs(x = "Trait", y = "Weights", title = paste("Cluster", i, "- Trait Weights"), fill = "Trait Type") +
      scale_fill_manual(values = c("1" = "gold", "0" = "steelblue"),
                        labels = c("1" = "Traits w/ increasing values", "0" = "Traits w/ decreasing values")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
            legend.position = c(0.90, 0.85),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 12),
            # legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
            legend.margin = margin(t = 5, r = 5, b = 5, l = 5))
  }
  file_path <- file.path(main_dir, paste0("cluster", i, "_barplot.png"))
  ggsave(file_path, grid.arrange(p1, p2, ncol = 1), width = 16, height = 16, dpi=300)
  # file_path_gene <- file.path(main_dir, paste0("cluster", i, "_barplot_gene.png"))
  # ggsave(file_path_gene, plot = p1, width = 10, height = 6, dpi = 300)
  # file_path_trait <- file.path(main_dir, paste0("cluster", i, "_barplot_trait.png"))
  # ggsave(file_path_trait, plot = p2, width = 12, height = 8, dpi = 300)
}


######
# Assuming `z_mat` is already loaded as a data frame or matrix
# Convert to matrix if it's a data frame with row names
z_mat <- as.matrix(z_mat)
# Set a color palette: blue for negative, white for zero, red for positive values
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
# Create a color palette for NA values (missingness)
na_color <- "grey"  # Grey color for missing values
z_mat_for_clustering <- z_mat
z_mat_for_clustering[is.na(z_mat_for_clustering)] <- 0
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
# Define breaks so that 0 is always white
breaks <- seq(min(z_mat_for_clustering, na.rm = TRUE), max(z_mat_for_clustering, na.rm = TRUE), length.out = 101)
zero_index <- which.min(abs(breaks - 0))  # Find where zero is in the breaks
breaks[zero_index] <- 0
colnames(z_mat_for_clustering) <- substr(colnames(z_mat_for_clustering), 1, 30)
pheatmap(
  z_mat_for_clustering,
  color = color_palette,        # Custom color palette
  breaks = breaks,              # Custom breaks to ensure 0 = white
  cluster_rows = TRUE,          # Cluster loci (rows)
  cluster_cols = TRUE,          # Cluster traits (columns)
  show_rownames = FALSE,        # Hide row names
  show_colnames = T,        # Hide column names
  fontsize_col = 10,            # Column name font size
  na_col = na_color,            # Grey for missing values
  scale = "none",               # No scaling (Z-scores are standardized)
  main = "Heatmap of Z-Scores (Zero = White, Non-zero = Colored)"
)

missingness_matrix <- ifelse(is.na(z_mat), 1, 0)
missing_colors <- c("white", "grey")
colnames(missingness_matrix) <- substr(colnames(missingness_matrix), 1, 30)
pheatmap(
  missingness_matrix,
  color = missing_colors,
  cluster_rows = TRUE,          # Cluster rows to show missingness patterns
  cluster_cols = TRUE,          # Cluster columns for trait patterns
  show_rownames = FALSE,        # Hide row names
  show_colnames = T,         # Show column names (traits)
  fontsize_col = 10,            # Font size for columns
  legend_breaks = c(0, 1),      # Legend showing 0 and 1
  legend_labels = c("Present", "Missing"),  # Legend labels
  main = "Heatmap of Missingness (Grey = Missing, White = Present)"
)