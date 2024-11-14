# ---------------------Illustration---------------------
# This .R file generates trait information from the supplementary tables provided by papers and the corresponding GWAS paths.
# The trait information is saved as "gwas_traits.xlsx" under "/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats".

# Autoimmune_Rheumatoid_Arthritis: The six traits for the six summary statistics are multi-ancestry, EUR-, and EAS-GWAS for all RA and seropositive RA. Here I simply saved the traits as 'rheumatoid arthritis' and 'rheumatoid arthritis; anti-citrullinated protein antibody seropositivity; rheumatoid factor seropositivity measurement'.
# Immune_Cell: The traits are described in more detail in the supplementary table. Thus, I used the information provided in the table directly.
# ------------------------------------------------------


# Load the necessary library
library(openxlsx) # install.packages('openxlsx')

## Construct two sheets of paths for both main GWAS and other trait GWAS
# Main GWAS
directory_path <- '/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data'   # Define the directory path where the GWAS summary statistics files are stored
main_gwas <- data.frame(study = "T1D_GWAS") # Create a data frame for the main GWAS study, starting with a column 'study' containing "T1D_GWAS"
main_gwas$root_path <- directory_path # Add a column 'root_path' to the data frame, setting it to the directory path defined above
files <- list.files(directory_path, pattern = "T1D", full.names = FALSE) # List all files in the specified directory that contain "T1D" in their names; only filenames are returned
main_gwas$file <- files # Add a new column 'file' for the file names
main_gwas$full_path <- file.path(main_gwas$root_path,main_gwas$file) # Add a new columnn of 'full_path' for full path

# Other trait GWAS 
types <- c('Autoimmune_Rheumatoid_Arthritis','Immune_Cell') # Define the trait types
files <- c()
n_types <- c()
traits <- c()
for(type in types){
  a <- list.files(directory_path,pattern=type, full.names = FALSE)
  n_a <- length(a)
  files <- c(files,a) 
  n_types <- c(n_types,n_a)
  trait_dir <- '/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/traits_information'
  if(type=='Immune_Cell'){
    data_trait <- read_excel(file.path(trait_dir,"Immune cell phenotypes Orru et al. Supplementary Table 1B.xlsx"), sheet = "Supplementary Table 1B",skip = 2)
    traits <- c(traits,data_trait$Trait)
  }else if(type=='Autoimmune_Rheumatoid_Arthritis'){
    traits <- c(traits,rep(c('all rheumatoid arthritis', 'rheumatoid arthritis; anti-citrullinated protein antibody seropositivity; rheumatoid factor seropositivity measurement'), each = 3))
  } # according to the paper, the six traits for the six summary statistics are multi-ancestry, EUR-, and EAS-GWAS for all RA and seropositive RA.
} # 731 files for Immune cell; 6 files for Autoimmune_Rheumatoid_Arthritis;  
trait_gwas <- cbind(
  type = rep(types, times = n_types),
  trait = traits,
  root_path = directory_path,
  file = files,
  full_path = file.path(directory_path,files)
)


## Save the generated excel book for the paths for all GWAS traits
wb <- createWorkbook() # Create a new workbook

# Add sheets and write data frames to respective sheets
addWorksheet(wb, "main_gwas")
writeData(wb, "main_gwas", main_gwas)

addWorksheet(wb, "trait_gwas")
writeData(wb, "trait_gwas", trait_gwas)

# Save the workbook
save_dir <- '/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats'
saveWorkbook(wb, file.path(save_dir,"gwas_traits.xlsx"), overwrite = TRUE)
