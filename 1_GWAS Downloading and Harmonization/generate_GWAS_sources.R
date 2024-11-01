install.packages("xlsx")
# Load necessary libraries
library(httr)
library(jsonlite)
library(xlsx)
library(dplyr)

# Set the directory path
directory_path <- "/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data"

# List files in the directory
files <- list.files(directory_path, pattern = "GCST[0-9]+", full.names = FALSE)

# Extract Accession Numbers and Types from file names
accession_numbers <- gsub(".*(GCST[0-9]+).*", "\\1", files)
types <- gsub("^(\\D+).*", "\\1", files) # Extract the characters before _(the first number) in the filename
types <- ifelse(types == "T", "T1D_", types)
types <- substr(types, 1, nchar(types) - 1) # Remove the last character from each entry in `types`

# Function to get Trait(s) from GWAS Catalog using accession number
get_trait_from_gwas <- function(accession) {
  url <- paste0("https://www.ebi.ac.uk/gwas/rest/api/studies/", accession,"/efoTraits")
  response <- GET(url)
  
  if (status_code(response) == 200) {
    content <- content(response, as = "text", encoding = "UTF-8")
    data <- fromJSON(content)
    return(data$`_embedded`$efoTraits$trait)
  } else {
    return(NA)
  }
}
trait <- as.vector(sapply(accession_numbers, get_trait_from_gwas)) # some accession numbers have multiple traits
trait <- lapply(trait, function(x) paste(x, collapse = "; "))

# Create a data frame to store Accession Numbers and Traits
gwas_data <- cbind(
  Type = types,
  Accession = accession_numbers,
  Trait = as.vector(trait),
  full_path = files
)

# Print or save the results
print(gwas_data)
# To save as a CSV file
write.csv(gwas_data, "/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/gwas_traits.csv", row.names = FALSE)
