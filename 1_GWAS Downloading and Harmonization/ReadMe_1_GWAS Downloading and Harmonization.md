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

} and  Immune_Cell_1 to Immune_Cell_731: GCST90001391 to GCST90002121

738 .h.tsv.gz files in total are saved in the directory: **'/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data'**

See see details in **download_GWAS.py file**

# 2 Harmonise GWAS ---- no need to harmonize, since the data from GWAS Catalog harmonised folder has already been hormonised.

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

#### For vcf

The chain file is downloaded from USCS:  [hg19ToHg38.over.chain.gz][https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/], and saved at "**/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Harmonized Data/Program/hg19ToHg38.over.chain.gz**".



I used the LifeOver software under [LiftOver userguide][https://gatk.broadinstitute.org/hc/en-us/articles/360037060932-LiftoverVcf-Picard]. 

1. Install the LiftOver software

```sh
# sh
java -version
wget https://github.com/broadinstitute/picard/releases/download/2.26.10/picard.jar
alias picard='java -jar /home/zhulx/software/picard.jar'
```

2. Use the LiftOver software

```sh
# {sh}
cd "/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data"
echo "##fileformat=VCFv4.2" > /home/zhulx/header.vcf
echo "##source=GWASCatalog" >> /home/zhulx/header.vcf
echo "##INFO=<ID=HM_VARIANT_ID,Number=1,Type=String,Description=\"Variant ID in Harmonization\">" >> /home/zhulx/header.vcf
echo "##INFO=<ID=HM_RSID,Number=1,Type=String,Description=\"RS ID\">" >> /home/zhulx/header.vcf
echo "##INFO=<ID=HM_EFFECT_ALLELE_FREQ,Number=1,Type=Float,Description=\"Effect Allele Frequency\">" >> /home/zhulx/header.vcf
echo "##INFO=<ID=BETA,Number=1,Type=Float,Description=\"Beta\">" >> /home/zhulx/header.vcf
echo "##INFO=<ID=ODDS_RATIO,Number=1,Type=Float,Description=\"Odds Ratio\">" >> /home/zhulx/header.vcf
echo "##INFO=<ID=P_VALUE,Number=1,Type=Float,Description=\"P-value\">" >> /home/zhulx/header.vcf
echo "##INFO=<ID=CI_LOWER,Number=1,Type=Float,Description=\"Confidence Interval Lower Bound\">" >> /home/zhulx/header.vcf
echo "##INFO=<ID=CI_UPPER,Number=1,Type=Float,Description=\"Confidence Interval Upper Bound\">" >> /home/zhulx/header.vcf
echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> /home/zhulx/header.vcf
zcat 	Immune_Cell_100_32929287-GCST90001490-EFO_0007937.h.tsv.gz | awk -v OFS='\t' 'NR > 1 {print $1, $2, ".", $3, $4, ".", "PASS", "AF=" $5}' >> variants.vcf
cat header.vcf variants.vcf > final.vcf
java -jar picard.jar LiftoverVcf \\ 
     I=final.vcf \\ # input.vcf
     O=harmonized_example.vcf \\ # lifted_over.vcf
     CHAIN=b37tohg38.chain \\ # b37tohg38.chain
     REJECT=rejected_variants.vcf \\ # rejected_variants.vcf
     R=reference_sequence.fasta # reference_sequence.fasta
rm header.vcf variants.vcf final.vcf

```

#### For tsv files

The chain file is downloaded from USCS:  [hg19ToHg38.over.chain.gz][https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/], and saved at "**/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Harmonized Data/Program/hg19ToHg38.over.chain.gz**".

I used the LifeOver software under [Bioconfuctor/liftOver][https://github.com/freeseek/score?tab=readme-ov-file#liftover-vcfs].

- Install the segment_liftover

#### For tsv files

Install liftover from [BCFtools/liftover][https://github.com/freeseek/score?tab=readme-ov-file]

```sh
# {sh}
module load plink/2.4a4
module load Bioinformatics
mkdir bcftools-1.21
cd bcftools-1.21
wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2
tar -xjf bcftools-1.12.tar.bz2
cd bcftools-1.21
./configure
make
make install
export PATH=/home/zhulx/software/bcftools-1.21/bin:$PATH
module load bcftools-1.21
export BCFTOOLS_PLUGINS=/sw/spack/bio/pkgs/gcc-10.3.0/bcftools/1.12-g4b275ez/libexec/bcftools # found at "which bcftools /sw/spack/bio/pkgs/gcc-10.3.0/bcftools/1.12-g4b275ez/bin/bcftools" 

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz{,.tbi}
bcftools +liftover --no-version -Ou ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz -- \
  -s $HOME/GRCh37/human_g1k_v37.fasta \
  -f $HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  -c $HOME/GRCh38/hg19ToHg38.over.chain.gz \
  --reject ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.reject.bcf \
  --reject-type b \
  --write-src | \
bcftools sort -o ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.hg38.bcf -Ob --write-index
```





## 2.2 Variant allele alighment (rel-alt and alt-ref)----no need, been harmonised already



