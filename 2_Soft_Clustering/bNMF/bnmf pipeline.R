
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

# load project scripts containing bNMF functions
source("./bnmf-clustering/scripts/choose_variants.R")  # ld_pruning, count_traits_per_variant, fina_variants_needing_proxies, & choose_potential_proxies
source("./bnmf-clustering/scripts/prep_bNMF.R")  # fetch_summary_stats & prep_z_matrix
source("./bnmf-clustering/scripts/run_bNMF.R")  # run_bNMF & summarize_bNMF

# USER INPUTS
project_dir = '/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/bNMF/test_results' # path to where you want results saved
user_token = 'ac846b46f561' # token for LDlinkR api

# create project folder 
dir.create(project_dir)

#----

# SECTION 1: PULL IN GWAS INFORMATION
data_dir =  '/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data'
# no need for rsID_map_file since GWAS data contains the rsIDs

# GWAS for main trait
gwas <- read_excel(file.path(data_dir, "clustering_data_sources_example.xlsx"),
                   sheet="main_gwas") %>%
  data.frame()

# GWAS for clustering traits
gwas_traits <- read_excel(file.path(data_dir, "clustering_data_sources_example.xlsx"),
                          sheet="trait_gwas")


