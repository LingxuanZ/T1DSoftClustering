# 1 Downloading GWAS data

## 1.1 GWAS data source

[Genotype Imputation Information](https://docs.google.com/spreadsheets/d/1l0u3dfIYQImFPsy1rGyCjebZeVbtcU6L0zyITQKowDk/edit?gid=0#gid=0) is extracted from all the articles below.

- Immune cell phenotypes
  - [Orru et al., 2021](https://www.nature.com/articles/s41588-020-0684-4) (multiple GWAS catalog accessions: from GCST0001391 to GCST0002121)
    - GWAS Catalog with accession numbers from GCST0001391 (https://www.ebi.ac.uk/gwas/studies/GCST0001391) to GCST0002121 (https://www.ebi.ac.uk/gwas/studies/GCST0002121)
    - The accession number for each trait is reported in Supplementary Table [1B](https://www.nature.com/articles/s41588-020-0684-4#MOESM3). ([Supplementary Tables](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-020-0684-4/MediaObjects/41588_2020_684_MOESM3_ESM.xlsx))
- T1D
  - [Chiou et al., 2021](https://www.nature.com/articles/s41586-021-03552-w) (GCST90014023)
    - Full summary statistics for the T1D GWAS have been deposited into the NHGRI-EBI GWAS catalogue with accession number GCST90014023 and can be downloaded from http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90014001-GCST90015000/GCST90014023/. 
  - [McGrail et al., 2024](https://www.medrxiv.org/content/10.1101/2024.07.31.24311310v1) (stats not available yet) 
- Autoimmune traits
  - Multiple sclerosis - [Patsopoulus et al., 2019](https://www.science.org/doi/full/10.1126/science.aav7188) (requested summary stats through the IMSGC website https://imsgc.net/)
  - Rheumatoid arthritis - [Ishigaki et al., 2022](https://www.nature.com/articles/s41588-022-01213-w) (deposited six sets of summary statistics [multi-ancestry, EUR-GWAS and EAS-GWAS for all RA and seropositive RA] to the GWAS Catalog under accession IDs GCST90132222, GCST90132223, GCST90132224, GCST90132225, GCST90132226 and GCST90132227)
  - Ankylosing spondylitis
  - Crohn’s disease
  - Ulcerative colitis
  - Autoimmune thyroid disease
  - Celiac
  - Juvenile idiopathic arthritis (GWAS underpowered, small sample size)
- Metabolic traits (already hosted on Great Lakes: /scratch/scjp_root/scjp1/xiaoouw/Muscle_snATAC_nucQTL/data/GWAS)
  - Some traits to start with
    - diamante_T2D-European
    - diamante_T2Dbmiadj-European
    - childhood-bmi_7years
    - childhood-bmi_3years
    - MAGIC_HbA1c-EUR
    - MAGIC_FI-EUR (negative control trait)
    - MAGIC_FG-EUR
    - cardio-ukbb_CAD (negative control trait)

## 1.2 GWAS Original data Downloading Codes

The code is written to automatically download GWAS files (ending in .h.tsv.gz) from the [GWAS Catalog](https://www.ebi.ac.uk/gwas/home) after obtaining the accession numbers provided in the papers.

accession_numbers_dict = {

​    "T1D": "GCST90014023",  # T1D

​    "Autoimmune_Rheumatoid_Arthritis_1": "GCST90132222", 

​    "Autoimmune_Rheumatoid_Arthritis_2": "GCST90132223",

​    "Autoimmune_Rheumatoid_Arthritis_3": "GCST90132224",

​    "Autoimmune_Rheumatoid_Arthritis_4": "GCST90132225",

​    "Autoimmune_Rheumatoid_Arthritis_5": "GCST90132226",

​    "Autoimmune_Rheumatoid_Arthritis_6": "GCST90132227",

} and  Immune_Cell_1 to Immune_Cell_731: GCST0001391 to GCST0002121

728 .h.tsv.gz files in total are saved in the directory: **'/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data'**

See see details in **download_GWAS.py file**

# 2 Harmonize GWAS

## 2.1  Liftover

Before merging, ensure that allele directions are consistent across all datasets.

Ensure Consistent Reference Panels.

### 2.1.1 Install liftover

```sh
# {sh}

# install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh # default path = $HOME/miniconda3
# run Miniconda
export PATH="/home/zhulx/miniconda3/bin:$PATH"
conda --version
# Install liftover using conda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c bioconda ucsc-liftover
```

### 2.1.2 use liftover

I used the LifeOver software under [LiftOver userguidde][https://gatk.broadinstitute.org/hc/en-us/articles/360037060932-LiftoverVcf-Picard]. 

first: convert https://genome.sph.umich.edu/wiki/LiftOver#Binary_liftOver_tool =====convert  extract the chromosome.

```sh
# {sh}
source activate
conda activate /home/zhulx/miniconda3
liftOver # check if liftOver is usable

# 

java -jar picard.jar LiftoverVcf \\ 
     I=input.vcf \\ # input.vcf
     O=lifted_over.vcf \\ # lifted_over.vcf
     CHAIN=b37tohg38.chain \\ # b37tohg38.chain
     REJECT=rejected_variants.vcf \\ # rejected_variants.vcf
     R=reference_sequence.fasta # reference_sequence.fasta

```



## 2.2 Variant allele alighment (rel-alt and alt-ref)

