# ---------------------Illustration---------------------
# This .R file generates trait information from the supplementary tables provided by papers and the corresponding GWAS paths.
# This .R file also automatically download the sample size N from Catalog urls and save N in the trait information.
# The trait information is saved as "gwas_traits.xlsx" under "/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats".

# Autoimmune_Rheumatoid_Arthritis: The six traits for the six summary statistics are multi-ancestry, EUR-, and EAS-GWAS for all RA and seropositive RA. 
# Immune_Cell: The traits are described in more detail in the supplementary table. Thus, I used the information provided in the table directly.
# ------------------------------------------------------


# Load the necessary library
library(openxlsx) # install.packages('openxlsx')
library(httr)
library(jsonlite)
library(dplyr)
library(readxl)

# Define the function to get the sample size N (number_of_individuals)
get_number_of_individuals <- function(accession) {
  # Construct the URL
  url <- paste0("https://www.ebi.ac.uk/gwas/rest/api/studies/", accession)
  
  # Make a GET request to the API
  response <- GET(url)
  
  # Check if the request was successful
  if (status_code(response) != 200) {
    stop("Failed to retrieve data. Check the accession code or network connection.")
  }
  
  # Parse the JSON content
  content <- content(response, as = "text", encoding = "UTF-8")
  data <- fromJSON(content)
  
  # Extract the sum of numberOfIndividuals value: "Using sum helps to avoid issues with GWAS studies on populations with different ancestries, as data$ancestries$numberOfIndividuals is a vector."
  if (!is.null(data$ancestries$numberOfIndividuals)) {
    return(sum(data$ancestries$numberOfIndividuals))
  } else {
    stop("numberOfIndividuals not found in the response.")
  }
}


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
size_N <- c()
for(type in types){
  a <- list.files(directory_path,pattern=type, full.names = FALSE)
  accession_numbers <- gsub(".*(GCST[0-9]+).*", "\\1", a)
  n_a <- length(a)
  files <- c(files,a) 
  n_types <- c(n_types,n_a)
  size_N <- c(size_N,as.vector(sapply(accession_numbers, get_number_of_individuals)))
  trait_dir <- '/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/traits_information'
  if(type=='Immune_Cell'){
    data_trait <- read_excel(file.path(trait_dir,"Immune cell phenotypes Orru et al. Supplementary Table 1B.xlsx"), sheet = "Supplementary Table 1B",skip = 2)
    data_trait <- data_trait %>%
      filter(`GWAS Catalog Accession Number` %in% accession_numbers) %>%
      arrange(match(`GWAS Catalog Accession Number`, accession_numbers))
    ordered_traits <- data_trait$Trait
    ordered_traits <- paste0(ordered_traits, "_panel_", data_trait$Panel) # Add panel info to all traits
    
    # # only for duplicated ordered_traits: Modify ordered_traits by appending _panel_{Panel} to duplicated elements
    # duplicated_elements <- ordered_traits[duplicated(ordered_traits)] # Identify duplicated elements
    # ordered_traits <- ifelse( # Modify ordered_traits by appending _panel_{Panel} to duplicated elements
    #   ordered_traits %in% duplicated_elements, 
    #   paste0(ordered_traits, "_panel_", data_trait$Panel),
    #   ordered_traits
    # )
    
    traits <- c(traits,ordered_traits)
  }else if(type=='Autoimmune_Rheumatoid_Arthritis'){
    traits <- c(traits,'all rheumatoid arthritis multi-ancestry','all rheumatoid arthritis EUR','all rheumatoid arthritis EAS',
                'rheumatoid arthritis; anti-citrullinated protein antibody seropositivity; rheumatoid factor seropositivity measurement multi-ancestry',
                'rheumatoid arthritis; anti-citrullinated protein antibody seropositivity; rheumatoid factor seropositivity measurement EUR',
                'rheumatoid arthritis; anti-citrullinated protein antibody seropositivity; rheumatoid factor seropositivity measurement EAS')
    # traits <- c(traits,rep(c('all rheumatoid arthritis', 'rheumatoid arthritis; anti-citrullinated protein antibody seropositivity; rheumatoid factor seropositivity measurement'), each = 3))
  } # according to the paper, the six traits for the six summary statistics are multi-ancestry, EUR-, and EAS-GWAS for all RA and seropositive RA.
} # 731 files for Immune cell; 6 files for Autoimmune_Rheumatoid_Arthritis;  
trait_gwas <- data.frame(
  type = rep(types, times = n_types),
  trait = traits,
  root_path = directory_path,
  file = files,
  full_path = file.path(directory_path,files),
  N = size_N
) # dim = 737 * 6

unique_traits <- unique(traits) # length(unique_traits)  is 735


## Save the generated excel book for the paths for all GWAS traits
wb <- createWorkbook() # Create a new workbook

# Add sheets and write data frames to respective sheets
addWorksheet(wb, "main_gwas")
writeData(wb, "main_gwas", main_gwas)

addWorksheet(wb, "trait_gwas_w_correlatedImmuneCe")
writeData(wb, "trait_gwas_w_correlatedImmuneCe", trait_gwas)

# Save the workbook
save_dir <- '/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats'
saveWorkbook(wb, file.path(save_dir,"gwas_traits.xlsx"), overwrite = TRUE)


######## Add Metabolic trait GWASs
save_dir <- '/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats'
file_path <- file.path(save_dir, "gwas_traits.xlsx")
traits <- c("diamante_T2D-European","diamante_T2Dbmiadj-European","childhood-bmi_7years",
            "childhood-bmi_3years", "MAGIC_HbA1c-EUR","MAGIC_FI-EUR(negative control trait)",
            "MAGIC_FG-EUR","cardio-ukbb_CAD(negative control trait)")
file <- c("diamante_T2D-European.bed.gz","diamante_T2Dbmiadj-European.bed.gz","childhood-bmi_7years.bed.gz",
          "childhood-bmi_3years.bed.gz","MAGIC_HbA1c-EUR.bed.gz","MAGIC_FI-EUR.bed.gz",
          "MAGIC_FG-EUR.bed.gz","cardio-ukbb_CAD.bed.gz")
file <- paste0("Metabolic_", file)
root_path <- rep("/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data", length(traits))
N <- rep(0, length(traits))
type <- rep("Metabolic", length(traits))
full_path = file.path(root_path,file)
new_data <- data.frame(
  type = type,
  trait = traits,
  root_path = root_path,
  file = file,
  full_path = full_path,
  N = N,
  stringsAsFactors = FALSE
)
wb <- loadWorkbook(file_path) # Load the existing workbook
existing_data <- read.xlsx(wb, sheet = "trait_gwas_w_correlatedImmuneCe") # Write the new data to the "trait_gwas" sheet (append to existing data)
combined_data <- rbind(existing_data, new_data)
writeData(wb, sheet = "trait_gwas_w_correlatedImmuneCe", x = combined_data) # Overwrite the "trait_gwas" sheet with the updated data
saveWorkbook(wb, file_path, overwrite = TRUE) # Save the workbook
cat("Data successfully added to the 'trait_gwas' sheet.\n")
# remember to modify N: because each SNP has different sample size!!!!


######## Add Autoimmune traits: Ankylosing spondylitis, Inflammatory bowel disease, thyroid Grave’s disease, thyroid Hashimoto’s disease, Celiac
directory_path <- '/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data'   # Define the directory path where the GWAS summary statistics files are stored
types <- c('Autoimmune_Ankylosing_spondylitis','Autoimmune_Inflammatory_bowel_disease',
           'Autoimmune_thyroid_Graves_disease','Autoimmune_thyroid_Hashimotos_disease',
           'Autoimmune_Celiac') # Define the trait types
files <- c()
n_types <- c()
traits <- types
size_N <- c()
for(type in types){
  a <- list.files(directory_path,pattern=type, full.names = FALSE)
  accession_numbers <- gsub(".*(GCST[0-9]+).*", "\\1", a)
  n_a <- length(a)
  files <- c(files,a) 
  n_types <- c(n_types,n_a)
  size_N <- c(size_N,as.vector(sapply(accession_numbers, get_number_of_individuals)))
} 
new_data <- data.frame(
  type = types,
  trait = traits,
  root_path = directory_path,
  file = files,
  full_path = file.path(directory_path,files),
  N = size_N
) # dim = 5 * 6
unique_traits <- unique(traits) # length(unique_traits)  is 5
# save data
wb <- loadWorkbook(file_path) # Load the existing workbook
existing_data <- read.xlsx(wb, sheet = "trait_gwas_w_correlatedImmuneCe") # Write the new data to the "trait_gwas" sheet (append to existing data)
combined_data <- rbind(existing_data, new_data)
writeData(wb, sheet = "trait_gwas_w_correlatedImmuneCe", x = combined_data) # Overwrite the "trait_gwas" sheet with the updated data
saveWorkbook(wb, file_path, overwrite = TRUE) # Save the workbook
cat("Data successfully added to the 'trait_gwas_w_correlatedImmuneCe' sheet.\n")
# remember to modify N: because each SNP has different sample size!!!!


######## Add sheet excluding gentic correlated GWASs: trait_gwas
wb <- loadWorkbook(file_path)
existing_data <- read.xlsx(wb, sheet = "trait_gwas_w_correlatedImmuneCe") 
selected_immunecellTraits <- c(
  "CD20 on memory B cell", "CD40 on monocytes","CD127 on CD4+","CD25 on CD24+ CD27+",
  "CD24 on memory B cell","HLA DR on  HLA DR+ T cell","CD45RA on naive CD4+","CD38 on IgD+ CD38dim",
  "CD19 on memory B cell","CD25 on CD4 Treg","CD11c on granulocyte","HLA DR on monocyte",
  "CD64 on CD14+ CD16- monocyte","CCR7 on naive CD4+","CD45 on lymphocyte_panel_TBNK","BAFF-R on memory B cell",
  "CD14 on CD14+ CD16+  monocyte","IgD on IgD+ CD38dim","CD86 on granulocyte","CD4 on activated Treg",
  "CD3 on CD8br","CD8 on  CD8br","CD39 on CD39+ CD4+","BAFF-R on CD20- CD38-",
  "CD28 on CD4 Treg","CD123 on plasmacytoid DC","CD3 on CD4 Treg","CD25 on CD39+ resting Treg",
  "CD27 on sw mem","CCR2 on CD62L+ myeloid DC","HLA DR on  HLA DR+ CD4+","CX3CR1 on monocyte",
  "HVEM on CD4+","CD19 on CD20-" 
) # selected 34 representative immune cell traits using gentic correlation matrix 
# HLA DR on  HLA DR+ CD4+ # not HLA DR on HLA DR+ CD4+
# CD14 on CD14+ CD16+  monocyte # not CD14 on CD14+ CD16+ monocyte
# CD8 on  CD8br # not CD8 on CD8br
# HLA DR on  HLA DR+ T cell # not HLA DR on HLA DR+ T cell
# CD45 on lymphocyte_panel_TBNK # not CD45 on lymphocyte
immuneCellTraits <- existing_data[existing_data$type=='Immune_Cell',]
otherTraits <- existing_data[existing_data$type!='Immune_Cell',]
selected_immunecellTraits_lower <- tolower(selected_immunecellTraits)
filtered_immuneCellTraits <- immuneCellTraits %>%
  filter(sapply(trait, function(x) {
    trait_lower <- tolower(x)
    prefix <- tolower(sub("_.*", "", x))
    prefix %in% selected_immunecellTraits_lower || trait_lower %in% selected_immunecellTraits_lower
  }))
# check whether all selected_immunecellTraits are selected
trait_prefixes <- sub("_.*", "", filtered_immuneCellTraits$trait)
non_matching_traits <- selected_immunecellTraits[(!selected_immunecellTraits %in% trait_prefixes) & (!selected_immunecellTraits %in%filtered_immuneCellTraits$trait)]
non_matching_traits
# check duplicated traits started with selected_immunecellTraits
filtered_immuneCellTraits <- filtered_immuneCellTraits %>%
  mutate(prefix = tolower(sub("_.*", "", trait)))
duplicated_prefixes <- filtered_immuneCellTraits$prefix[duplicated(filtered_immuneCellTraits$prefix)]
duplicate_rows <- filtered_immuneCellTraits %>%
  filter(prefix %in% duplicated_prefixes)
print(duplicate_rows)
filtered_immuneCellTraits <- filtered_immuneCellTraits %>%
  select(-prefix)
# save new sheet: trait_gwas
combined_data <- rbind(filtered_immuneCellTraits, otherTraits)
addWorksheet(wb, "trait_gwas")
writeData(wb, sheet = "trait_gwas", x = combined_data) # Overwrite the "trait_gwas" sheet with the updated data
saveWorkbook(wb, file_path, overwrite = TRUE) # Save the workbook
cat("Data successfully added to the 'trait_gwas' sheet.\n")

