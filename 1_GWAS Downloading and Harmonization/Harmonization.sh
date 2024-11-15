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