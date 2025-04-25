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



## 2.2  Liftover 

See GWAS summary stats/generate_pathsource_gwas_traits.R: ######## Add Beta Cell Function trait GWASs ### Rename file and harmonization: 

```bash
# try
cd "/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data"
input_tsv="meta1xinsG30.filthetmafn.rsid.selectedcolumns"
output_bed="Beta_cell_function_xinsG30.bed.gz"
rejected_bed="Beta_cell_function_xinsG30_unmapped.hg19.bed.gz"
```

Slurs script

```bash
#!/bin/bash
#SBATCH --time=00:25:00
#SBATCH --mem=4G
#SBATCH --account=scjp1
#SBATCH --output='../BetaCellFunc_Harmonization.log'
#SBATCH --mail-user=zhulx@umich.edu
#SBATCH --mail-type=END,FAIL

module load Bioinformatics UCSC-utilities/4.0.0
module load htslib

# === Input variables ===
input_tsv=$1          # .selectedcolumns (hg19) file
output_bed=$2         # output .bed.gz file (hg38)
rejected_bed=$3       # rejected (unmapped) .bed.gz file
prefix="${output_bed%%.*}"

# === Tools and files ===
chain="/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/hg19ToHg38.over.chain.gz"
liftover_tool=$(which liftOver)

# === Temporary files ===
input_tsv_renames=$(mktemp)
bed_tmp=$(mktemp)
bed4_only_bed=$(mktemp)
mapped_bed=$(mktemp)
rejected_tmp=$(mktemp)
mapped_full_with_id=$(mktemp)
mapped_full=$(mktemp)
rejected_full=$(mktemp)

# === Saved Temporar files ===
lookup_file="${prefix}.lookup.tsv"
lookup_sorted_file="${prefix}.lookup.sorted.tsv"
joined_file="${prefix}.joined.tsv"
lookup_rejected_file="${prefix}.lookup_rejected.tsv"
lookup_rejected_sorted_file="${prefix}.lookup_rejected.sorted.tsv"
bed_tmp_sorted_file="${prefix}.bed_tmp.sorted.tsv"
rejected_joined_file="${prefix}.rejected_joined.tsv"


echo "Step 1: Renaming header columns..."
awk 'NR==1 {
  for (i=1; i<=NF; i++) {
    if ($i=="P-value") $i="p_value"
    else if ($i=="Effect") $i="hm_beta"
    else if ($i=="StdErr") $i="standard_error"
    else if ($i=="Allele1") $i="ALT"
    else if ($i=="Allele2") $i="REF"
    else if ($i=="MarkerName") $i="CHRPOS_grch37"
    else if ($i=="TotalSampleSize") $i="n_complete_samples"
  }
}
{ print }' OFS='\t' "$input_tsv" > "$input_tsv_renames"

# ================================================================================

echo "Step 2: Converting to BED format (with all columns)..."
awk 'BEGIN{OFS="\t"} NR==1 {
  print "chrom", "start", "end", $0     # head：add colnames in BED
  next
}
{
  split($1, a, ":");
  print "chr"a[1], a[2]-1, a[2], $0          # output all columns
}' "$input_tsv_renames" > "$bed_tmp"
tail -n +2 "$bed_tmp" | cut -f1-4 > "$bed4_only_bed"

# ================================================================================

echo "Step 3: Running liftOver..."
$liftover_tool "$bed4_only_bed" "$chain" "$mapped_bed" "$rejected_tmp"

# ================================================================================

echo "Step 4: Compressing output..."
LC_ALL=C sort -k4,4 "$bed_tmp" > "$bed_tmp_sorted_file"

cut -f1,3,4 "$mapped_bed" | awk -F'\t' 'BEGIN{OFS="\t"} {gsub("chr", "", $1); print $3, $1, $2}' > "$lookup_file"
LC_ALL=C sort -k1,1 "$lookup_file" > "$lookup_sorted_file"
join -t $'\t' -1 1 -2 4 "$lookup_sorted_file" "$bed_tmp_sorted_file" > "$joined_file"
awk -F'\t' 'BEGIN{OFS=FS}{
    $7 = toupper($7)  # 将第7列（ALT）转大写
    $8 = toupper($8)  # 将第8列（REF）转大写
    print
}' "$joined_file" > "${joined_file}.tmp" && mv "${joined_file}.tmp" "$joined_file"
header=$(echo -e "hm_variant_id\tchrom\thm_position\tCHRPOS_grch37\tstart_hg19\tend_hg19\tALT\tREF\tFreq1\tFreqSE\tMinFreq\tMaxFreq\thm_beta\tstandard_error\tp_value\tn_complete_samples")
awk -F'\t' 'BEGIN{OFS="\t"}
{
  chr=$2; pos=$3; ref=toupper($8); alt=toupper($7)
  id = chr"_"pos"_"ref"_"alt
  $1 = id  # replace CHRPOS_grch37 with hm_variant_id
  print
}' "$joined_file" > "$mapped_full_with_id"
echo -e "$header" | cat - "$mapped_full_with_id" | gzip -c > "$output_bed" # echo -e "$header" | cat - "$mapped_full_with_id" | bgzip -c > "$output_bed"

# rej
cut -f1,3,4 "$rejected_tmp" | \
awk -F'\t' 'BEGIN{OFS="\t"}
{
  gsub("chr", "", $1)
  if ($1 != "" && $2 != "" && $3 != "") {
    print $3, $1, $2
  }
}' > "$lookup_rejected_file"
LC_ALL=C sort -k1,1 "$lookup_rejected_file" > "$lookup_rejected_sorted_file"
join -t $'\t' -1 1 -2 4 "$lookup_rejected_sorted_file" "$bed_tmp_sorted_file" | \
awk -F'\t' 'BEGIN{OFS="\t"} {
  $7 = toupper($7);  # ALT
  $8 = toupper($8);  # REF
  print
}' > "$rejected_joined_file"
header_rej=$(echo -e "CHRPOS_grch37\tchrom\tposition_19\tchrom\tstart_hg19\tend_hg19\tALT\tREF\tFreq1\tFreqSE\tMinFreq\tMaxFreq\thm_beta\tstandard_error\tp_value\tn_complete_samples")
echo -e "$header_rej" | cat - "$rejected_joined_file" | gzip -c > "$rejected_bed"

rm -f "$lookup_file" \
      "$lookup_sorted_file" \
      "$joined_file" \
      "$lookup_rejected_file" \
      "$lookup_rejected_sorted_file" \
      "$bed_tmp_sorted_file" \
      "$rejected_joined_file"
echo "Done!"
```

```bash
cd "/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data"

sbatch "../BetaCellFunc_Harmonization.sh" \
    meta1bigair.filthetmafn.rsid.selectedcolumns \
    Beta_cell_function_BIGTT-AIR.bed.gz \
    Beta_cell_function_BIGTT-AIR_Rej.hg19.unmapped.bed.gz # Submitted batch job 24598278

sbatch "../BetaCellFunc_Harmonization.sh" \
    meta1cir.filthetmafn.rsid.selectedcolumns \
    Beta_cell_function_CIR.bed.gz \
    Beta_cell_function_CIR_Rej.hg19.unmapped.bed.gz # Submitted batch job 24598281

sbatch "../BetaCellFunc_Harmonization.sh" \
    meta1di.filthetmafn.rsid.selectedcolumns \
    Beta_cell_function_DI.bed.gz \
    Beta_cell_function_DI_Rej.hg19.unmapped.bed.gz # Submitted batch job 24598282

sbatch "../BetaCellFunc_Harmonization.sh" \
    meta1dibig.filthetmafn.rsid.selectedcolumns \
    Beta_cell_function_DIBIG.bed.gz \
    Beta_cell_function_DIBIG_Rej.hg19.unmapped.bed.gz # Submitted batch job 24598283

sbatch "../BetaCellFunc_Harmonization.sh" \
    meta1homab.filthetmafn.rsid.selectedcolumns \
    Beta_cell_function_HOMA-beta.bed.gz \
    Beta_cell_function_HOMA-beta_Rej.hg19.unmapped.bed.gz # Submitted batch job 24598284

sbatch "../BetaCellFunc_Harmonization.sh" \
    meta1stumvoll.filthetmafn.rsid.selectedcolumns \
    Beta_cell_function_Stumvoll.bed.gz \
    Beta_cell_function_Stumvoll_Rej.hg19.unmapped.bed.gz # Submitted batch job 24598287

sbatch "../BetaCellFunc_Harmonization.sh" \
    meta1xinsdG30.filthetmafn.rsid.selectedcolumns \
    Beta_cell_function_xinsdG30.bed.gz \
    Beta_cell_function_xinsdG30_Rej.hg19.unmapped.bed.gz # Submitted batch job 24598289

sbatch "../BetaCellFunc_Harmonization.sh" \
    meta1xinsG30.filthetmafn.rsid.selectedcolumns \
    Beta_cell_function_xinsG30.bed.gz \
    Beta_cell_function_xinsG30_Rej.hg19.unmapped.bed.gz # Submitted batch job 24598294
```

