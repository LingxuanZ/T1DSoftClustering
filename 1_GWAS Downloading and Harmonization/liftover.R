# install liftOver from BiocManager
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
BiocManager::install("liftOver")
  # Do you want to install from sources the package which needs compilation? (Yes/no/cancel) no
  # Update all/some/none? [a/s/n]: a
  # browseVignettes("liftOver") # documentation
suppressPackageStartupMessages({
library('liftOver')
})

library(gwascat)
library(liftOver)
library(data.table) # install.packages("data.table")
library(R.utils)
library(rtracklayer)

input_directory <- '/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data'
file_name <- 'Immune_Cell_100_32929287-GCST90001490-EFO_0007937.h.tsv.gz'
data <- fread(file.path(input_directory, file_name)) # install.packages("R.utils")

# Import the chain file
# chain_file <- "/scratch/scjp_root/scjp0/ccrober/reference/hg19ToHg38.over.chain.gz" # cannot access the folder and the file 
chain_file <- "/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Harmonized Data/Program/hg19ToHg38.over.chain.gz"
chain <- import.chain(chain_file)

seqlevelsStyle(data) = "UCSC"  # necessary

