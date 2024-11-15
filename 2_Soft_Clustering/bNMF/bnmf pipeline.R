
#----

start=Sys.time()

# load requires packages
install.packages("devtools")
install.packages("pacman")
pacman::p_load(tidyverse, data.table, readxl, magrittr, dplyr, strex,
               rstudioapi, DT, kableExtra, GenomicRanges)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
BiocManager::install("Homo.sapiens")
library(devtools)
devtools::install_github("CBIIT/LDlinkR")
library(LDlinkR)
library(readxl)
# load('/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/bNMF/test_results/pipeline_data.RData') # to load the saved environment variables (for debugging)

# load project scripts containing bNMF functions
source('/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/bNMF/choose_variants.R') # ld_pruning, count_traits_per_variant, fina_variants_needing_proxies, & choose_potential_proxies
source("./bnmf-clustering/scripts/prep_bNMF.R")  # fetch_summary_stats & prep_z_matrix
source("./bnmf-clustering/scripts/run_bNMF.R")  # run_bNMF & summarize_bNMF

# USER INPUTS
project_dir = '/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/bNMF/test_results' # path to where you want results saved
user_token = 'cb5457b210a6' # 'ac846b46f561' # token for LDlinkR api

# create project folder 
dir.create(project_dir)

#------------------------------------------------

# SECTION 1: PULL IN GWAS INFORMATION

data_dir =  '/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats'
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

save.image(file = file.path(project_dir, "pipeline_data.RData"))

#------------------------------------------------

# SECTION 3: VARIANT PRUNING (LD-BASED)

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

save.image(file = file.path(project_dir, "pipeline_data.RData"))


#------------------------------------------------

# SECTION 4: VARIANT MISSINGNESS

print("Searching for variants in trait GWAS...")
gwas_variants <- pruned_vars$VAR_ID
df_Ns <- count_traits_per_variant(gwas_variants,
                                  ss_files = trait_ss_files,
                                  savepath_varid="/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/bNMF/test_results/all_snps_varids.tmp")
# fix column names
df_Ns_rev <- df_Ns %>%
  column_to_rownames("VAR_ID") %>%
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

save.image(file = file.path(project_dir, "pipeline_data.RData"))


