
#----

start=Sys.time()

# load requires packages
install.packages("devtools")
install.packages("pacman")
install.packages("cowplot")
pacman::p_load(tidyverse, data.table, readxl, magrittr, dplyr, strex,
               rstudioapi, DT, kableExtra, GenomicRanges)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
BiocManager::install("Homo.sapiens")
BiocManager::install("ComplexHeatmap")
BiocManager::install('ensembldb')
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

library(devtools)
devtools::install_github("CBIIT/LDlinkR")
library(LDlinkR)
library(readxl)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(AnnotationHub)
library(ensembldb)
library(ComplexHeatmap)
library(cowplot)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")
# load('./test_results/pipeline_data.RData') # to load the saved environment variables (for debugging)
load('./test_results/pipeline_data_wo_corr.RData')

# load project scripts containing bNMF functions
source("./prep_bNMF.R")  # fetch_summary_stats & prep_z_matrix
source("./bnmf-clustering/scripts/run_bNMF.R")  # run_bNMF & summarize_bNMF
source('./choose_variants.R') # ld_pruning, count_traits_per_variant, fina_variants_needing_proxies, & choose_potential_proxies
source("./exclude_genes.R")

# USER INPUTS
project_dir = './test_results' # path to where you want results saved
# project_final_dir = './final_test_results'
project_final_dir = './final_test_results_wo_corrImmCell'
user_token = 'ac846b46f561' # 'cb5457b210a6' # token for LDlinkR api
rsID_map_file = '../Data/GWAS summary stats/rsid_map_fromMainGWAS.txt' # '/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/rsid_map_fromMainGWAS.txt'

# create project folder 
dir.create(project_dir)

#------------------------------------------------

# SECTION 1: PULL IN GWAS INFORMATION

data_dir =  '../Data/GWAS summary stats' # '/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats'
# no need for rsID_map_file since GWAS data contains the rsIDs

# GWAS for main trait
gwas <- read_excel(file.path(data_dir, "gwas_traits.xlsx"),
                   sheet="main_gwas") %>% data.frame()

# GWAS for clustering traits
gwas_traits <- read_excel(file.path(data_dir, "gwas_traits.xlsx"),
                          sheet="trait_gwas") %>% data.frame()

# GWAS to be used for final allele alignment
main_ss_filepath <- gwas %>% pull(full_path)
gwas_ss_files <- setNames(gwas$full_path, gwas$study)
trait_ss_files <- setNames(gwas_traits$full_path, gwas_traits$trait)
trait_ss_size <- setNames(gwas_traits$N, gwas_traits$trait)
metabolic_traits <- c("diamante_T2D-European","diamante_T2Dbmiadj-European","childhood-bmi_7years",
            "childhood-bmi_3years", "MAGIC_HbA1c-EUR","MAGIC_FI-EUR(negative control trait)",
            "MAGIC_FG-EUR","cardio-ukbb_CAD(negative control trait)")
metabolic_trait_ss_files <- trait_ss_files[metabolic_traits]

#------------------------------------------------

# SECTION 2: PULL SIGNIFICANT VARIANTS FROM MAIN TRAIT GWAS

# P-value threshold for variants in main trait
PVCUTOFF = 5e-8
n_gwas <- length(gwas_ss_files)
vars_sig = data.frame(VAR_ID = as.character(),
                      P_VALUE = as.numeric(),
                      Risk_Allele=as.character(),
                      GWAS=as.character())
print(sprintf("Pulling significant SNPs w/ pval<%.1e from %i T1D GWAS...", PVCUTOFF, n_gwas))

for(i in 1:n_gwas) {
  print(paste0("...Reading ", names(gwas_ss_files)[i], "..."))
  
  vars <- fread(gwas_ss_files[i], data.table = F, stringsAsFactors=F)
  
  if (!"hm_beta" %in% colnames(vars)){
    print("Converting Odds Ratio to Log Odds Ratio...")
    vars <- vars %>%
      mutate(hm_beta = log(as.numeric(hm_odds_ratio)))
  }
  vars <- vars %>%
    filter(as.numeric(p_value) <= PVCUTOFF) %>%
    subset(grepl("^[0-9]+_[0-9]+_[ACGT]+_[ACGT]", hm_variant_id)) %>%
    separate(hm_variant_id, into=c("CHR", "POS", "REF", "ALT"), sep="_", remove = F) %>%  # effect allele is ALT
    mutate(Risk_Allele = ifelse(hm_beta>=0, ALT, REF)) %>%
    mutate(GWAS = gwas$study[i]) %>%
    select(hm_variant_id, p_value, Risk_Allele, GWAS)
  
  print(nrow(vars))
  vars_sig = rbind(vars_sig, vars)
}
print(paste("No. total SNPs below pval cutoff:",nrow(vars_sig)))

# remove duplicates
vars_sig_uniq <- vars_sig %>%
  arrange(hm_variant_id, p_value) %>%
  filter(!duplicated(hm_variant_id)) %>% # so we remove duplicates with the higher pvalue
  rename(PVALUE = p_value)
print(paste("No. unique SNPs:",nrow(vars_sig_uniq)))

# remove indels
vars_sig_noIndels <- vars_sig_uniq %>%
  separate(hm_variant_id, into=c("CHR","POS","REF","ALT"),sep="_",remove = F) %>%
  mutate(alleles = paste0(REF,ALT)) %>%
  subset(nchar(alleles)==2 | (nchar(alleles)<=4 & grepl(",",alleles))) %>%
  rename(VAR_ID = hm_variant_id) %>%
  select(VAR_ID, PVALUE, Risk_Allele, GWAS)
print(paste("No. SNPs excluding indels:",nrow(vars_sig_noIndels)))

# save.image(file = file.path(project_dir, "pipeline_data.RData"))
save.image(file = file.path(project_dir, "pipeline_data_wo_corr.RData"))

#------------------------------------------------

# SECTION 3: VARIANT PRUNING (LD-BASED) for MAIN TRAIT GWAS

# LD pruning
print("LD-pruning using EUR panel in LDlinkR::SNPclip...")
ld_prune(df_snps = vars_sig_noIndels,
         pop = "EUR",
         output_dir = project_dir,
         r2 = 0.05,
         maf=0.001,
         my_token = user_token,
         genome_assembly = 'grch38',
         chr = c(1:22))

# combine LD-pruning results
print("Combining SNP.clip results...")
ld_files <- list.files(path = project_dir,
                       pattern = "^snpClip_results",
                       full.names = T)

df_clipped_res = data.frame("RS_Number"=as.character(),
                            "Position"=as.character(),
                            "Alleles"= as.character(),
                            "Details"=as.character())

rename_cols_clipped <- c(RS_Number="RS Number",
                         Position_grch37="Position")

for (ld_file in ld_files){
  df <- fread(ld_file, stringsAsFactors = F, data.table = F) %>%
    dplyr::rename(any_of(rename_cols_clipped))
  df_clipped_res <- rbind(df_clipped_res, df)
}

df_clipped_kept <- df_clipped_res %>%
  filter(Details=="Variant kept.")

pruned_vars <- vars_sig_noIndels %>%
  separate(VAR_ID, into=c("CHR","POS","REF","ALT"), sep = "_",remove = F) %>%
  mutate(ChrPos = paste0("chr", CHR, ":", POS)) %>%
  filter(ChrPos %in% df_clipped_kept$Position)
print(sprintf("T1D SNPs pruned from %i to %i...", nrow(vars_sig_noIndels), nrow(pruned_vars)))

# save.image(file = file.path(project_dir, "pipeline_data.RData"))
save.image(file = file.path(project_dir, "pipeline_data_wo_corr.RData"))

#------------------------------------------------

# SECTION 4: VARIANT MISSINGNESS: connect Main GWAS and Trait GWASs - time consuming
# The missing ratio represents the fraction of selected SNPs from the main GWAS that are present across all trait GWAS for each SNP.
# Calculate the proportion of non-missing values in each row.

sample_size_traitGWAS <- gwas_traits$N # trait_ss_size

print("Searching for variants in trait GWAS...")
gwas_variants <- pruned_vars$VAR_ID
# Use parallel programming # modified Ns for the eight metabolic traits!!!!!! 
df_Ns <- count_traits_per_variant(gwas_variants,
                                  ss_files = trait_ss_files,
                                  sample_size = sample_size_traitGWAS, # trait_ss_size
                                  savepath_varid="./test_results/all_snps_varids.tmp",
                                  savepath_varid_inverse="./test_results/all_snps_varids_inverse.tmp")
# i=1-737: Autoimmune_Rheumatoid_Arthritis & immune cell # w_corr
# i=738-745: Metabolic # w_corr
# using a different LD pruning strategy on the main T1D GWAS, as only 371 SNPs remain in the main GWAS right now. (?)
# save.image(file = file.path(project_dir, "pipeline_data.RData"))
save.image(file = file.path(project_dir, "pipeline_data_wo_corr.RData"))

# fix column names
df_Ns_rev <- df_Ns %>%
  pivot_wider(names_from="trait", values_from="N") %>%
  data.frame() %>%
  column_to_rownames("hm_variant_id") %>%
  set_colnames(names(trait_ss_files))

print("Calculating variant missingess in traits...")
variant_counts_df <- data.frame(VAR_ID=rownames(df_Ns_rev),
                                frac=rowSums(!is.na(df_Ns_rev[,names(trait_ss_files)]))/length(trait_ss_files))
var_nonmissingness <- ifelse(
  gwas_variants %in% variant_counts_df$VAR_ID,
  # if in counts data frame, take the non-missing fraction:
  variant_counts_df$frac[match(gwas_variants, variant_counts_df$VAR_ID)],
  # else not in data frame, so non-missing fraction is 0:
  0
)
var_nonmissingness <- setNames(var_nonmissingness, gwas_variants)

# save.image(file = file.path(project_dir, "pipeline_data.RData"))
save.image(file = file.path(project_dir, "pipeline_data_wo_corr.RData"))

#----

# SECTION 5: DETERMINE VARIANTS NEEDING PROXIES

print("Identifying variants needing proxies...")
proxies_needed_df <- find_variants_needing_proxies(gwas_variant_df=pruned_vars,
                                                   var_nonmissingness=var_nonmissingness,
                                                   rsID_map_file = rsID_map_file,
                                                   missing_cutoff = 0.3) # 0.8 or 0.5 or 0.3?
# [1] "Choosing variants in need of proxies..."
# [1] "...63 strand-ambiguous variants"
# [1] "...0 multi-allelic variants"
# [1] "...81 variants with excessive missingness" # 371-81=290---I guess this might be caused by immune cell GWAS(highly similar SNPs)
# [1] "...129 unique variants in total in need of proxies"
# [1] "...129 of these are mapped to rsIDs"
# save.image(file = file.path(project_dir, "pipeline_data.RData"))
save.image(file = file.path(project_dir, "pipeline_data_wo_corr.RData"))


#----
# SECTION 6: PROXY SEARCH - time consuming

print("Searching for proxies with TopLD API...")
proxy_search_results <- choose_proxies(need_proxies = proxies_needed_df,
                                       method="LDlink",
                                       LDlink_token = user_token,
                                       rsID_map_file = rsID_map_file,
                                       trait_ss_files = trait_ss_files,
                                       pruned_variants = pruned_vars,
                                       population="EUR")
# proxy_df saved at "./combined_query_snp_list_grch38.txt"
write.csv(as.data.frame(proxy_search_results),  file.path(project_dir, "proxy_search_results.csv"), row.names = FALSE)
# proxy_search_results <- read.csv(file.path(project_dir, "proxy_search_results.csv"), stringsAsFactors = FALSE)

df_proxies <- proxy_search_results %>%
  dplyr::select(hm_variant_id, proxy_VAR_ID) %>%
  dplyr::inner_join(pruned_vars[,c("VAR_ID","GWAS")], by=c("hm_variant_id"="VAR_ID")) %>%
  mutate(Risk_Allele=NA, PVALUE=NA) # add two new NA columns: Risk_Allele and PVALUE

# save.image(file = file.path(project_dir, "pipeline_data.RData"))
save.image(file = file.path(project_dir, "pipeline_data_wo_corr.RData"))


#----

# SECTION 7: Fetch summary statistics for SNPs in trait GWAS - time consuming - fread

# proxy_search_results <- read.csv(file.path(project_dir, "proxy_search_results.csv"), stringsAsFactors = FALSE)
print("Prepping input for fetch_summary_stats...")
df_orig_snps <- pruned_vars %>% # 371 SNPs in pruned_vars
  filter(!VAR_ID %in% proxies_needed_df$hm_variant_id) %>%
  rename(hm_variant_id = VAR_ID)

df_input_snps <- rbind(df_orig_snps[,c("hm_variant_id","PVALUE", "Risk_Allele", "GWAS")],
                       df_proxies[,c("hm_variant_id","PVALUE", "Risk_Allele", "GWAS")]) %>%
  arrange(PVALUE) %>%
  filter(!duplicated(hm_variant_id))%>%
  rename(VAR_ID = hm_variant_id)

df_orig_snps <- df_orig_snps %>%
  rename(VAR_ID = hm_variant_id)

cat(sprintf("\n%i original SNPs...\n", nrow(df_orig_snps)))
cat(sprintf("\n%i proxy SNPs...\n", nrow(df_proxies)))
cat(sprintf("\n%i total unique SNPs!\n", nrow(df_input_snps)))
# 242 original SNPs...
# 41 proxy SNPs...
# 283 total unique SNPs!

initial_zscore_matrices <- fetch_summary_stats(
  df_variants=df_input_snps,
  gwas_ss_file=main_ss_filepath,
  trait_ss_files=trait_ss_files,
  trait_ss_size=trait_ss_size,
  pval_cutoff=0.05
)
saveRDS(initial_zscore_matrices, file = file.path(project_dir,"initial_zscore_matrices.rds"))

# save.image(file = file.path(project_dir, "pipeline_data.RData"))
save.image(file = file.path(project_dir, "pipeline_data_wo_corr.RData"))
# system(sprintf("mv alignment_GWAS_summStats.csv '%s'", project_dir))

#----

# Section 8: get rsIDs for final variant set

print("Getting rsIDs for final snps and saving to results...")
z_mat <- initial_zscore_matrices$z_mat
N_mat <- initial_zscore_matrices$N_mat

df_var_ids <- df_input_snps %>%
  separate(VAR_ID, into=c("Chr","Pos","Ref","Alt"),sep="_",remove = F) %>%
  mutate(ChrPos=paste(Chr,Pos,sep = ":")) %>%
  subset(ChrPos %in% rownames(z_mat))
write(df_var_ids$VAR_ID,file.path(project_dir, 'my_snps.tmp'))

system(sprintf("grep -wFf '%s' '%s' > '%s'",file.path(project_dir, 'my_snps.tmp'),
               rsID_map_file, file.path(project_dir, "rsID_map.txt")))

df_rsIDs <- fread(cmd=sprintf("grep -wFf '%s' '%s'",file.path(project_dir, 'my_snps.tmp'),rsID_map_file),
                  header = F,
                  col.names = c("VAR_ID","rsID"))
print(sprintf("rsIDs found for %i of %i SNPs...", nrow(df_rsIDs), nrow(df_var_ids)))

df_rsIDs_final <- df_rsIDs %>%
  filter(VAR_ID %in% df_var_ids$VAR_ID)
write_delim(x=df_rsIDs_final,
            file = file.path(project_dir, "rsID_map.txt"),
            col_names = T)

# save.image(file = file.path(project_dir, "pipeline_data.RData"))
save.image(file = file.path(project_dir, "pipeline_data_wo_corr.RData"))


#----

# Section 9: Fill missing data in z-score and N matrices

df_snps <- df_input_snps %>%
  inner_join(df_rsIDs_final, by="VAR_ID") %>%
  data.frame()

print("Searching for cover proxies for missing z-scores...")
initial_zscore_matrices_final <- fill_missing_zscores(initial_zscore_matrices,
                                                      df_snps,
                                                      trait_ss_files,
                                                      trait_ss_size,
                                                      main_ss_filepath,
                                                      rsID_map_file,
                                                      method_fill="median",
                                                      population="EUR")
# save.image(file = file.path(project_dir, "pipeline_data.RData"))
save.image(file = file.path(project_dir, "pipeline_data_wo_corr.RData"))


#----

# Section 10.) Generate non-negative z-score matrix

# This section is to filter traits w/ no pvalues below cutoff, prune traits by correlation (remove traits with Pearson |r| > 0.80), perform sample size adjustment

section10_output_file_path <- file.path(project_dir, "section10_prep_z_printoutput.txt") # Set file path to save the output
sink(section10_output_file_path) # Start redirecting console output to the file
prep_z_output <- prep_z_matrix(z_mat = initial_zscore_matrices_final$z_mat,
                               N_mat = initial_zscore_matrices_final$N_mat,
                               corr_cutoff = 0.8)
sink() # Stop redirecting output to the file and restore to the console
cat("All output in section 10 has been saved to:", section10_output_file_path, "\n")

# prep_z_output has two outputs:

#   1.) The scaled, non-negative z-score matrix
final_zscore_matrix <- prep_z_output$final_z_mat

#   2.) Results from the trait filtering
df_traits_filtered <- prep_z_output$df_traits
write_csv(x = df_traits_filtered,
          file = file.path(project_dir,"df_traits.csv"))

# prep_z_matrix also save trait correlation matrix to working dir, so move to project dir
system(sprintf("mv trait_cor_mat.txt '%s'", project_dir))

print(sprintf("Final matrix: %i SNPs x %i traits",
              nrow(final_zscore_matrix),
              ncol(final_zscore_matrix)/2))

# save.image(file = file.path(project_dir, "pipeline_data.RData"))
save.image(file = file.path(project_dir, "pipeline_data_wo_corr.RData"))

#----

# Section 11.) Run bNMF 
bnmf_reps <- run_bNMF(final_zscore_matrix,
                      n_reps=25,
                      tolerance = 1e-6)
summarize_bNMF(bnmf_reps, dir_save=project_final_dir)
# K: Kth cluster
# W: feature-to-cluster matrices.
# H: cluster-to-gene matrices.
# Use the internal function make_run_summary() to calculate the number of clusters (K) and the corresponding likelihood for each run, generating a summary table.

# files
# run_summary.txt: Contains a summary table of each run, listing the selected K value and the corresponding likelihood for each iteration.
# L2EU.W.mat.K.txt and L2EU.H.mat.K.txt: The cleaned versions of the W and H matrices, representing the contributions of features to clusters and clusters to traits.
# W_plot_K.pdf and H_plot_K.pdf: Heatmaps for each K, visualizing the associations in the W and H matrices.

# The summarized result  identify the most meaningful number of clusters (K), and visualize the W and H matrices. 
# save.image(file = file.path(project_dir, "pipeline_data.RData"))
save.image(file = file.path(project_dir, "pipeline_data_wo_corr.RData"))

# end=Sys.time()
# print("Total pipeline runtime:")
# print(end-start)

#----

# format results
k <- NULL
if (is.null(k)){
  html_filename <- file.path(project_final_dir, "results_for_maxK.html")
} else {
  html_filename <- sprintf(file.path(project_final_dir, "results_for_K_%i.html"), k)
}

rmarkdown::render(
  './format_bNMF_results.Rmd',
  output_file = html_filename,
  params = list(main_dir = project_final_dir,
                align_dir = project_dir,
                k = k,
                loci_file="query",
                GTEx=F,
                my_traits=gwas_traits)
)


#----
# Section 12.) Run bNMF using final_zscore_matrix excluding genes w/ text: "HLA"
# project_final_dir_wo_HLA <- "./final_test_results_wo_HLA"
project_final_dir_wo_HLA <- "./final_test_results_wo_corrImmCell_wo_HLA"
final_zscore_matrix_wo_HLA <- exclude_genes_with_annotation(final_zscore_matrix,with_text="HLA")
bnmf_reps_wo_HLA <- run_bNMF(final_zscore_matrix_wo_HLA,
                      n_reps=25,
                      tolerance = 1e-6)
summarize_bNMF(bnmf_reps_wo_HLA, dir_save=project_final_dir_wo_HLA)

# format results
k <- NULL
if (is.null(k)){
  html_filename <- file.path(project_final_dir_wo_HLA, "results_for_maxK.html")
} else {
  html_filename <- sprintf(file.path(project_final_dir_wo_HLA, "results_for_K_%i.html"), k)
}

rmarkdown::render(
  './format_bNMF_results.Rmd',
  output_file = html_filename,
  params = list(main_dir = project_final_dir_wo_HLA,
                align_dir = project_dir,
                k = k,
                loci_file="query",
                GTEx=F,
                my_traits=gwas_traits)
)


#----
# Section 13.) Run bNMF using final_zscore_matrix excluding genes within chr6:20Mb-40Mb
project_final_dir_wo_chr6_20_40 <- "./final_test_results_wo_corrImmCell_wo_chr6_20_40"

rownames_df <- data.frame(
  rownames = rownames(final_zscore_matrix),
  stringsAsFactors = FALSE
)

rownames_df <- transform(
  rownames_df,
  chr = sub(":.*", "", rownames),
  pos = as.numeric(sub(".*:", "", rownames))
)

filtered_rows <- with(
  rownames_df,
  !(chr == "6" & pos >= 20000000 & pos <= 40000000)
)

final_zscore_matrix_wo_chr6_20_40 <- final_zscore_matrix[filtered_rows, ]
bnmf_reps_wo_chr6_20_40 <- run_bNMF(final_zscore_matrix_wo_chr6_20_40,
                             n_reps=25,
                             tolerance = 1e-6)
summarize_bNMF(bnmf_reps_wo_chr6_20_40, dir_save=project_final_dir_wo_chr6_20_40)

# format results
k <- NULL
if (is.null(k)){
  html_filename <- file.path(project_final_dir_wo_chr6_20_40, "results_for_maxK_wo_chr6_20_40.html")
} else {
  html_filename <- sprintf(file.path(project_final_dir_wo_chr6_20_40, "results_for_K_%i_wo_chr6_20_40.html"), k)
}

rmarkdown::render(
  './format_bNMF_results.Rmd',
  output_file = html_filename,
  params = list(main_dir = project_final_dir_wo_chr6_20_40,
                align_dir = project_dir,
                k = k,
                loci_file="query",
                GTEx=F,
                my_traits=gwas_traits)
)





