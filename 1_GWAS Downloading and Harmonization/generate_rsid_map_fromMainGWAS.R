# ---------------------Illustration---------------------
# This R script is designed to extract the two columns (hm_variant_id and hm_rsid), and save the results in text file.
# input: main GWAS summary statistics path
# output: the path of the varid_grch38 - rsid map file.
# Since only the SNPs from the main GWAS are used, the IDs of other SNPs are unnecessary.
# ------------------------------------------------------

library(data.table)

# Input file path
input_file <- "/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data/T1D_34012112-GCST90014023-EFO_0001359.h.tsv.gz"

# Output file path
output_file <- "/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/rsid_map_fromMainGWAS.txt"

# Read data
cat("Reading input file...\n")
data <- fread(input_file)

# Check if required columns exist
if (!all(c("hm_variant_id", "hm_rsid") %in% colnames(data))) {
  stop("Error: The input file does not contain 'hm_variant_id' or 'hm_rsid' columns.")
}

# Extract required columns
cat("Extracting hm_variant_id and hm_rsid columns...\n")
subset_data <- data[, .(hm_variant_id, hm_rsid)]

# Check for empty rows (rows where all values are NA or empty strings)
empty_rows <- rowSums(is.na(subset_data) | subset_data == "") == ncol(subset_data)
print(paste("Number of completely empty rows:", sum(empty_rows)))
# Check for rows containing any NA values
rows_with_na <- rowSums(is.na(subset_data)) > 0
print(paste("Number of rows with at least one NA:", sum(rows_with_na)))
# delete rows with NA
subset_data <- subset_data[!rows_with_na, ]

# Save to output file
cat("Writing to output file...\n")
fwrite(subset_data, file = output_file, sep = "\t", quote = FALSE)
cat("Process completed successfully. Data saved to:", output_file, "\n")
