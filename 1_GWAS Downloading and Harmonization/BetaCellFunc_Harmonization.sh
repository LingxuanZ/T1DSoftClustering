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