# ---------------------Illustration---------------------
# This .py file automatically retrieves harmonized GWAS data from the GWAS Catalog using the accession numbers provided in the papers.
# The harmonized data is downloaded from the URL: f"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST{accession_range_start}-GCST{accession_range_end}/{accession_number}/harmonised/"
# This URL contains the folder for harmonized GWAS data using the genome assembly GRCh38.
# ------------------------------------------------------

import os
import requests
import re
from bs4 import BeautifulSoup

############### accession numbers input ############### 
accession_numbers_dict = {
    "T1D": "GCST90014023",  # T1D
    "Autoimmune_Rheumatoid_Arthritis_1": "GCST90132222", 
    "Autoimmune_Rheumatoid_Arthritis_2": "GCST90132223",
    "Autoimmune_Rheumatoid_Arthritis_3": "GCST90132224",
    "Autoimmune_Rheumatoid_Arthritis_4": "GCST90132225",
    "Autoimmune_Rheumatoid_Arthritis_5": "GCST90132226",
    "Autoimmune_Rheumatoid_Arthritis_6": "GCST90132227",
}
# Add a range of accession numbers for Immune_cell
immune_cell_start = 90001391
immune_cell_end = 90002121  # Based on the description (731 numbers)
# Add all the immune cell accession numbers to the dictionary
for i in range(immune_cell_start, immune_cell_end + 1):
    accession_number = f"GCST{i:08d}"  # Ensure the accession number format with leading zeros
    key = f"Immune_Cell_{i - immune_cell_start + 1}"  # Sequential key for each immune cell
    accession_numbers_dict[key] = accession_number
########################################################



##### define the functions to automatically downloading GWAS data ##### 
def file_exists(url):
    """
    Check if a file exists at the given URL by sending a HEAD request.

    :param url: The URL of the file to check
    :return: True if the file exists, False otherwise
    """
    response = requests.head(url)
    return response.status_code == 200

def find_matching_file(url, suffix=".h.tsv.gz"):
    """
    Find the first file in the directory at the given URL that ends with the specified suffix.
    Note: .h denotes harmonized files.

    :param url: The URL of the directory to search
    :param suffix: The suffix that the file should end with (default is ".h.tsv.gz")
    :return: The full URL of the matching file, or None if no such file is found
    """
    response = requests.get(url)
    
    # If the request is successful, parse the HTML
    if response.status_code == 200:
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Find all <a> tags (which represent links in the HTML)
        for link in soup.find_all('a'):
            file_name = link.get('href')
            
            # Check if the file name ends with the specified suffix
            if file_name and file_name.endswith(suffix):
                return os.path.join(url, file_name)
    
    # Return None if no matching file is found
    return None

def download_file(url, local_filename):
    """
    Downloads a file from a given URL and saves it locally.

    :param url: URL to the file
    :param local_filename: local path where the file will be saved
    """
    # Send a GET request to the URL
    with requests.get(url, stream=True) as r:
        r.raise_for_status()  # Raises an HTTPError for bad responses
        # Open the local file for writing in binary mode
        with open(local_filename, 'wb') as f:
            # Write the contents of the response to the file
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    print(f"File downloaded: {local_filename}")
################################################################### 



################## Path to save the downloaded files ################## 
save_directory = '/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data' # '/home/zhulx/Parker Lab/T1D Soft Clustering/Data/GWAS summary stats/Original Data'
# Create the directory if it does not exist
if not os.path.exists(save_directory):
    os.makedirs(save_directory)
#######################################################################



####################### Download the GWAS data #######################  
# Iterate over the accession numbers and find the corresponding files
for key, accession_number in accession_numbers_dict.items():
    # Define the directory URL based on the accession number range
    accession_range_start = int(accession_number[4:]) // 1000 * 1000 + 1
    if accession_range_start > int(accession_number[4:]):
        accession_range_start = accession_range_start-1000
    accession_range_end = accession_range_start + 999
    directory_url = f"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST{accession_range_start}-GCST{accession_range_end}/{accession_number}/harmonised/"
    
    if file_exists(directory_url):
        # Find the file in the directory that ends with .h.tsv.gz
        matching_file_url = find_matching_file(directory_url, suffix=".h.tsv.gz")
        
        if matching_file_url:
            # Extract the file name from the URL
            file_name_from_url = os.path.basename(matching_file_url)
            
            # Define the local filename as {key}_{file_name_from_url}
            filename = os.path.join(save_directory, f"{key}_{file_name_from_url}")
            
            # Download the file
            download_file(matching_file_url, filename)
        else:
            # Print out the directory URL if no matching file is found
            print(f"No matching file found in: {directory_url}")
    else:
        print(f"The directort: {directory_url} doesn't exist")
######################################################################



####################### Download the additional Autoimmune GWAS data #######################  
accession_numbers_dict_new1 = {
    "Autoimmune_Ankylosing_spondylitis": "GCST005529", 
    "Autoimmune_Inflammatory_bowel_disease": "GCST90446792",
    "Autoimmune_thyroid_Graves_disease": "GCST90018847",
    "Autoimmune_thyroid_Hashimotos_disease": "GCST90018855",
    "Autoimmune_Celiac": "GCST005523"
}
save_directory = '/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data' # '/home/zhulx/Parker Lab/T1D Soft Clustering/Data/GWAS summary stats/Original Data'
# save_directory = '/nfs/turbo/umms-scjp/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data'
# Iterate over the accession numbers and find the corresponding files
for key, accession_number in accession_numbers_dict_new1.items():
    if accession_number[4:][0] != "0": 
        # Define the directory URL based on the accession number range
        accession_range_start = int(accession_number[4:]) // 1000 * 1000 + 1
        if accession_range_start > int(accession_number[4:]):
            accession_range_start = accession_range_start-1000
        accession_range_end = accession_range_start + 999
        directory_url = f"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST{accession_range_start}-GCST{accession_range_end}/{accession_number}/harmonised/"
    else:
        character = accession_number[4:]
        match = re.match(r"0*", character)
        leading_zeros = match.group() if match else ""
        accession_range_start = int(accession_number[4:]) // 1000 * 1000 + 1
        if accession_range_start > int(accession_number[4:]):
            accession_range_start = accession_range_start-1000
        accession_range_end = accession_range_start + 999
        directory_url = f"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST{leading_zeros}{accession_range_start}-GCST{leading_zeros}{accession_range_end}/{accession_number}/harmonised/"

    if file_exists(directory_url):
        # Find the file in the directory that ends with .h.tsv.gz
        matching_file_url = find_matching_file(directory_url, suffix=".h.tsv.gz")
        
        if matching_file_url:
            # Extract the file name from the URL
            file_name_from_url = os.path.basename(matching_file_url)
            
            # Define the local filename as {key}_{file_name_from_url}
            filename = os.path.join(save_directory, f"{key}_{file_name_from_url}")
            
            # Download the file
            download_file(matching_file_url, filename)
        else:
            # Print out the directory URL if no matching file is found
            print(f"No matching file found in: {directory_url}")
    else:
        print(f"The directort: {directory_url} doesn't exist")
######################################################################  

