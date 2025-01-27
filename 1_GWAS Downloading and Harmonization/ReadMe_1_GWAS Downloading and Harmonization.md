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
  
  **==Every file is under GRCh38 - Confirmed==**
  
  - Some traits to start with  (8 files in total)
    - diamante_T2D-European: diamante_T2D-European.bed.gz
    - diamante_T2Dbmiadj-European: diamante_T2Dbmiadj-European.bed.gz
    - childhood-bmi_7years: childhood-bmi_7years.bed.gz
    - childhood-bmi_3years: childhood-bmi_3years.bed.gz
    - MAGIC_HbA1c-EUR: MAGIC_HbA1c-EUR.bed.gz 
    - MAGIC_FI-EUR (negative control trait): MAGIC_FI-EUR.bed.gz
    - MAGIC_FG-EUR: MAGIC_FG-EUR.bed.gz
    - cardio-ukbb_CAD (negative control trait): cardio-ukbb_CAD.bed.gz

Effect_allele = ALT = EA

other_allele = REF = NEA

## 1.2 GWAS Original data Downloading Codes

The code is written to automatically download GWAS files (ending in .h.tsv.gz) from the [GWAS Catalog](https://www.ebi.ac.uk/gwas/home) after obtaining the accession numbers provided in the papers.

accession_numbers_dict = {

​    "T1D": "GCST90014023",  # T1D: hm_variant_id and hm_rsid are harmonized columns

​    "Autoimmune_Rheumatoid_Arthritis_1": "GCST90132222",  # Autoimmune_Rheumatoid_Arthritis: rsid is the only one harmonized column; the first 4 columns  (chromosome,base_pair_location,effect_allele,other_allele) are in grch38: chr_loc_other_effect

​    "Autoimmune_Rheumatoid_Arthritis_2": "GCST90132223",

​    "Autoimmune_Rheumatoid_Arthritis_3": "GCST90132224",

​    "Autoimmune_Rheumatoid_Arthritis_4": "GCST90132225",

​    "Autoimmune_Rheumatoid_Arthritis_5": "GCST90132226",

​    "Autoimmune_Rheumatoid_Arthritis_6": "GCST90132227",

} and  Immune_Cell_1 to Immune_Cell_731: GCST90001391 to GCST90002121

738 .h.tsv.gz files in total are saved in the directory: **'/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data'**

See see details in **download_GWAS.py file** 

==**For Metabolic traits, remember to modify the `N` values, as each SNP has a different sample size.**==



# 2 Harmonise GWAS ---- no need to harmonize, since the data from GWAS Catalog harmonised folder has already been hormonised.

## 2.0 Add a new column: hm_variant_id (required)

The datasets of  Immune_Cell and T1D already contain the columns of hm_variant_id and hm_rsid

We need to add a new column of hm_variant_id into Autoimmune_Rheumatoid_Arthritis data sets and rename rsid into hm_rsid.

```sh
# {sh}
#!/bin/bash

# Define file paths
source_dir="/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data"
old_dir="${source_dir}/old"

# Create a directory to store old files
mkdir -p "$old_dir"

# Find all .tsv.gz files starting with "Autoimmune_Rheumatoid_Arthritis" and process them
for file in "$source_dir"/Autoimmune_Rheumatoid_Arthritis*.tsv.gz; do
    # Get the filename (without path)
    filename=$(basename "$file")
    # Move the file to the old directory
    mv "$file" "$old_dir"
done

# Find all .tsv.gz files starting with "Autoimmune_Rheumatoid_Arthritis" and process them
for file in "$old_dir"/Autoimmune_Rheumatoid_Arthritis*.tsv.gz; do
    # Get the filename (without path)
    filename=$(basename "$file")
    # Unzip, add the new column, and save the file in the original directory with the same filename
    zcat "$old_dir/$filename" | \
    tr -d '\r' | \
    awk 'BEGIN {OFS="\t"} 
        NR==1 {gsub("rsid", "hm_rsid", $0); print $0, "hm_variant_id"} 
        NR>1 {print $0, $1"_"$2"_"$4"_"$3}' \
    | gzip > "$source_dir/$filename"
done
```

Metabolic Traits

```sh
# {sh}
#!/bin/bash

source_dir="/scratch/scjp_root/scjp1/xiaoouw/Muscle_snATAC_nucQTL/data/GWAS"
save_dir="/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data"

files=(
    "diamante_T2D-European.bed.gz"
    "diamante_T2Dbmiadj-European.bed.gz"
    "childhood-bmi_7years.bed.gz"
    "childhood-bmi_3years.bed.gz"
    "MAGIC_HbA1c-EUR.bed.gz"
    "MAGIC_FI-EUR.bed.gz"
    "MAGIC_FG-EUR.bed.gz"
)

mkdir -p "$save_dir"

# Iterate over each file for copying and modification
for file in "${files[@]}"; do
  echo "Processing $file..."
  
  # Copy the file from the source directory to the target directory
  cp "$source_dir/$file" "$save_dir/"

  # Decompress the file for modification
  gunzip -c "$save_dir/$file" > "$save_dir/${file%.gz}"

  # Modify the file to add and rename columns
  awk -v OFS="\t" '
  BEGIN {header_processed = 0}
  {
    # Process the header row
    if (header_processed == 0) {
      for (i=1; i<=NF; i++) {
        if ($i == "SNP") {
          $i = "hm_rsid" # Rename "SNP" to "hm_rsid"
        }
      }
      print $0, "hm_variant_id" # Add the new column "hm_variant_id"
      header_processed = 1
    } else {
      chrom_numeric = $1
      gsub(/^chr/, "", chrom_numeric) # Remove "chr" from the start of snp_chrom
      print $0, chrom_numeric"_"$3"_"$11"_"$10
    }
  }' "$save_dir/${file%.gz}" > "$save_dir/${file%.gz}_temp"

  # Overwrite the original file with the modified content
  mv "$save_dir/${file%.gz}_temp" "$save_dir/${file%.gz}"

  # Recompress the modified file
  gzip -f "$save_dir/${file%.gz}"
done

echo "All files processed and saved to $save_dir"


# Since the EA values in the "cardio-ukbb_CAD.bed.gz" have the form: 1:569406_G_A, not A
# I need to process cardio-ukbb_CAD.bed.gz to get the correct form of hm_variant_id.
file = "cardio-ukbb_CAD.bed.gz"
echo "Processing $file..."
cp "$source_dir/$file" "$save_dir/" # Copy the file from the source directory to the target directory
gunzip -c "$save_dir/$file" > "$save_dir/${file%.gz}" # Decompress the file for modification
awk -v OFS="\t" ' # Modify the file to add and rename columns
BEGIN {header_processed = 0}
{
  # Process the header row
  if (header_processed == 0) {
    for (i=1; i<=NF; i++) {
      if ($i == "SNP") {
        $i = "hm_rsid" # Rename "SNP" to "hm_rsid"
      }
    }
    print $0, "hm_variant_id" # Add the new column "hm_variant_id"
    header_processed = 1
  } else {
    chrom_numeric = $1
    gsub(/^chr/, "", chrom_numeric) # Remove "chr" from the start of snp_chrom
    allele = $10
    sub(/^[^_]*_/, "", allele) 
    print $0, chrom_numeric"_"$3"_"allele
  }
}' "$save_dir/${file%.gz}" > "$save_dir/${file%.gz}_temp"

mv "$save_dir/${file%.gz}_temp" "$save_dir/${file%.gz}" # Overwrite the original file with the modified content
gzip -f "$save_dir/${file%.gz}"    # Recompress the modified file
```

Autoimmune_Inflammatory_bowel_disease: Rename variant_id and rsid into hm_variant_id and hm_rsid (grch38 confirmed).

```sh
#!/bin/bash

# Define file paths
source_dir="/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data"
old_dir="${source_dir}/old"

# Create a directory to store old files
mkdir -p "$old_dir"

# Find all .tsv.gz files starting with "Autoimmune_Rheumatoid_Arthritis" and process them
for file in "$source_dir"/Autoimmune_Inflammatory_bowel_disease*.tsv.gz; do
    # Get the filename (without path)
    filename=$(basename "$file")
    # Move the file to the old directory
    mv "$file" "$old_dir"
done


# Find all .tsv.gz files starting with "Autoimmune_Rheumatoid_Arthritis" and process them
for file in "$old_dir"/Autoimmune_Inflammatory_bowel_disease*.tsv.gz; do
    # Get the filename (without path)
    filename=$(basename "$file")
    # Unzip, rename the column, and save the file in the original directory with the same filename
    zcat "$old_dir/$filename" | \
    awk 'BEGIN {OFS="\t"}
        NR==1 {
            for (i=1; i<=NF; i++) {
                if ($i == "variant_id") $i = "hm_variant_id";
                else if ($i == "rsid") $i = "hm_rsid";
            }
        }
        {print $0}' | \
    gzip > "$source_dir/$filename"
done
```

## 2.1 Delete duplicated SNPs

the row for Autoimmune_Ankylosing_spondylitis with 16_11285918_T_C is duplicated.---remain the first row and delete the second row.

```sh
source_dir="/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data"
old_dir="${source_dir}/old"
file="${source_dir}/Autoimmune_Ankylosing_spondylitis_23749187-GCST005529-EFO_0003898.h.tsv.gz"
mkdir -p "$old_dir"
mv "$file" "$old_dir"

input_file="/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data/old/Autoimmune_Ankylosing_spondylitis_23749187-GCST005529-EFO_0003898.h.tsv.gz"
output_file="/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data/Autoimmune_Ankylosing_spondylitis_23749187-GCST005529-EFO_0003898.h.tsv.gz"

zcat "$input_file" | awk -F'\t' '
BEGIN { OFS = FS }
NR == 1 { print; next }  # print the first row (colnames)
$1 == "16_11285918_T_C" {
    if (!seen[$1]++) print
    next
}
{ print }
' | gzip > "$output_file"
```



## 2.2  Liftover (optional)

Before merging, ensure that allele directions are consistent across all datasets.

Ensure Consistent Reference Panels.

### 2.2.1 Install liftover

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

### 2.2.2 use liftover

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





## 2.3 Variant allele alighment (rel-alt and alt-ref)----no need, been harmonised already (optional)



