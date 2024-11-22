library(parallel) # Load the parallel package
library(dplyr)
library(tidyr)

# The following ld_prund function is a modified version of the ld_prune function from the package available at https://github.com/gwas-partitioning/bnmf-clustering.
ld_prune <- function(df_snps,
                     pop,
                     my_token,
                     r2=0.1,
                     maf=0.01,
                     chr=1:22,
                     genome_assembly = 'grch38',
                     output_dir="./") {
  
  snp_clip_input <- df_snps %>% 
    separate(VAR_ID, into=c("CHR","POS","REF","ALT"),sep="_",remove = F) %>%
    mutate(ChrPos = paste0("chr",CHR,":",POS)) %>%
    arrange(PVALUE)
  
  df_clipped <- data.frame(ChrPos=as.character(),
                           rsID=as.character())
  for (i in chr){
    start = Sys.time()
    cur_chr <- snp_clip_input %>%
      filter(CHR==i)
    print(sprintf("Chr %i (%i SNPs)",i, nrow(cur_chr)))
    
    
    if (nrow(cur_chr) == 0) {
      next
    }
    else if (nrow(cur_chr) == 1) {
      
      cur_snps <- cur_chr %>%
        pull(ChrPos)
      # if only one SNP, use LDhap to get the variant info
      clipped_res <- LDlinkR::LDhap(snps = cur_snps, 
                                    pop = pop, 
                                    token = my_token,
                                    genome_build = genome_assembly,
                                    table_type = "variant"
      ) %>%
        rename(Alleles=Allele_Frequency) %>%
        mutate(Details="Variant kept.")
    } else { # >1 SNP
      if (between(nrow(cur_chr), 1, 5000)) { 
        cur_snps <- cur_chr %>%
          pull(ChrPos)
        
      } else {
        print("Chromosome has >5000 SNPs; breaking into sections...")
        var_df_list <- split(cur_chr, (seq(nrow(cur_chr))-1) %/% 5000)
        cur_snps <- c()
        
        for (j in 1:length(var_df_list)){
          
          print(sprintf("Pruning subset %i for chromosome %i...", j, i))
          var_df <- var_df_list[[j]]
          cur_snps_j <- var_df %>%
            pull(ChrPos)
          
          clipped_res_split <- LDlinkR::SNPclip(
            cur_snps_j,
            pop = pop,
            r2_threshold = r2,
            maf_threshold = maf,
            token = my_token,
            file = FALSE,
            genome_build = genome_assembly)
          
          clipped_snps <- clipped_res_split %>%
            filter(Details=="Variant kept.") %>%
            pull(RS_Number)
          print(sprintf("Subset %i pruned to %i SNPs...", j, length(clipped_snps)))
          cur_snps <- c(cur_snps, clipped_snps)
        }
      }
      
      print(sprintf("Performing final chromosomal pruning for %i SNPs...", length(cur_snps)))
      clipped_res <- LDlinkR::SNPclip(
        cur_snps,
        pop = pop,
        r2_threshold = r2,
        maf_threshold = maf,
        token = my_token,
        file = FALSE,
        genome_build = genome_assembly)
    }
    
    fwrite(x=clipped_res,
           file = file.path(output_dir, sprintf("snpClip_results_%s_chr%i.txt", pop, i)),
           quote = F,
           sep = "\t")
    
    df_clipped_final <- clipped_res %>%
      filter(Details=="Variant kept.")
    print(sprintf("Chr%i pruned from %i to %i SNPs...",i, nrow(cur_chr), nrow(df_clipped_final)))
    
    end = Sys.time()
    print(end-start)
    
  }
  print("Done!")
  
}


count_traits_per_variant <- function(gwas_variants, ss_files,sample_size,savepath_varid) {
  
  # Given a vector of variants, a vector of sample_size, the savepath of selected var_id_grch38 for main GWAS 
  # and a named vector of summary statistics files for traits to be clustered, output a vector of non-missing 
  # trait fractions per variant
  
  # Since the number of GWAS files is very large, we opted not to use the original 
  # sequential algorithm from https://github.com/gwas-partitioning/bnmf-clustering.
  # Instead, we implemented a parallel algorithm to accelerate the processing speed.
  # Additionally, I optimized the code and data structure, reducing some redundant steps,
  # which more than doubled the speed and reduced memory usage.
  
  
  print("Assessing variant missingness across traits...")
  write(gwas_variants, savepath_varid)
  
  numCores <- detectCores() - 1 # Define the number of cores to use, leaving one core for the system
  variant_df_list <- mclapply(1:length(ss_files), function(i) { # Use mclapply for parallel processing; mclapply is the parallel version of lapply.
    # print(sprintf("...Reading %s...", names(ss_files)[i])) 
    library(dplyr)
    library(data.table)
    start_time_i <- Sys.time()  # Start time
    headers <- as.character(fread(ss_files[i], nrows=1,
                                  data.table=F, stringsAsFactors=F, header=F))
    hm_variant_id_index <- which(headers == "hm_variant_id")
    if (length(hm_variant_id_index) == 0) {
      stop("Error: 'hm_variant_id' column not found in headers.")
    }
    
    if (endsWith(ss_files[i],".gz")) {
      df <- fread(
        cmd = sprintf("gzip -cd '%s' | cut -f%d | fgrep -wf '%s'", ss_files[i], hm_variant_id_index, savepath_varid),
        header = F,
        col.names = "hm_variant_id",
        data.table = F,
        stringsAsFactors = F
      )
    } else {
      df <- fread(cmd=sprintf("cut -f%d |fgrep -wf %s %s", hm_variant_id_index, savepath_varid, ss_files[i]),
                  header=F,
                  col.names="hm_variant_id",
                  data.table=F,
                  stringsAsFactors=F)
    }
    end_time_i <- Sys.time()  # End time
    time_taken <- as.numeric(end_time_i - start_time_i, units = "mins")  # Calculate duration
    print(sprintf("...Reading %s...: %d rows, Time taken: %s min", names(ss_files)[i], nrow(df), time_taken))
    # print(nrow(df))
    return(df) # A df contains the hm_variant_id column only.
  }, mc.cores = numCores)  # mc.cores parameter specifies the number of parallel processes to run.
  
  # make dataframe of Ns
  df_N <- lapply(seq_along(variant_df_list), function(i) {
    variant_df_list[[i]]$N <- sample_size[i]  # Add a new column 'N' with the value sample_size[i]
    variant_df_list[[i]]  # Return the modified data frame
  })
  df_N <- df_N %>%
    setNames(names(ss_files)) %>%
    bind_rows(.id="trait") # %>%
    # pivot_wider(names_from="trait", values_from="N") %>%
    # data.frame()
}


find_variants_needing_proxies <- function(gwas_variant_df, var_nonmissingness,
                                          rsID_map_file, missing_cutoff=0.8) {
  
  # Given a data frame containing GWAS variants and alleles as well as a vector
  # of trait missingness fractions per variant (from count_traits_per_variant),
  # output a vector of variants that need proxies
  # Criteria (any of the following):
  #   Strand-ambiguous (AT or GC)
  #   Multi-allelic
  #   Low-count (available in < 80% of traits)
  # rsID_map_file should point to a whitespace-delimited file with columns
  # corresponding to VAR_ID and rsID
  
  print("Choosing variants in need of proxies...")
  
  gwas_variant_df <- gwas_variant_df %>%
    rename(hm_variant_id = VAR_ID) %>%
    separate(hm_variant_id, into=c("CHR", "POS", "REF", "ALT"),
             sep="_", remove=F)
  
  need_proxies_varid <- with(gwas_variant_df, {
    strand_ambig <- hm_variant_id[paste0(REF, ALT) %in% c("AT", "TA", "CG", "GC")]
    print(paste0("...", length(strand_ambig), " strand-ambiguous variants"))
    
    multi_allelic <- grep("^[0-9]+_[0-9]+_[ACGT]+_[ACGT]+,[ACGT]+$", hm_variant_id, value=T)  # i.e. ALT allele has a comma
    print(paste0("...", length(multi_allelic), " multi-allelic variants"))
    
    low_cnt <- hm_variant_id[!(hm_variant_id %in% names(var_nonmissingness)) |
                        var_nonmissingness[hm_variant_id] < missing_cutoff]
    print(paste0("...", length(low_cnt), " variants with excessive missingness"))
    
    unique(c(strand_ambig, multi_allelic, low_cnt)) 
  })
  print(paste0("...", length(need_proxies_varid), " unique variants in total"))
  
  if (length(need_proxies_varid) == 0) return(tibble(hm_variant_id=c(), hm_rsid=c()))
  
  write(need_proxies_varid, "need_proxies_varid.tmp")
  varid_rsid_map <- fread(cmd= paste0("grep -wFf need_proxies_varid.tmp \"", rsID_map_file, "\""),#paste0("grep -wFf need_proxies_varid.tmp ",rsID_map_file),
                          header=F, col.names=c("hm_variant_id", "hm_rsid"),
                          data.table=F, stringsAsFactors=F)
  need_proxies_rsid <- varid_rsid_map$hm_rsid[match(need_proxies_varid, 
                                                 varid_rsid_map$hm_variant_id)]
  print(paste0("...", length(unique(varid_rsid_map$hm_rsid)), 
               " of these are mapped to rsIDs"))
  system("rm need_proxies_varid.tmp")
  
  tibble(hm_variant_id=need_proxies_varid) %>%
    left_join(varid_rsid_map, by="hm_variant_id") %>%
    left_join(gwas_variant_df[,c("hm_variant_id","PVALUE")], by="hm_variant_id")
}


choose_proxies <- function(need_proxies, 
                           rsID_map_file,
                           trait_ss_files,
                           pruned_variants,
                           method="LDlink",
                           LDlink_token=NULL,
                           population="EUR",
                           frac_nonmissing_num=0.8,
                           sample_size_traitGWAS = sample_size_traitGWAS,
                           r2_num=0.8) {
  # Modification: I have abandoned the topLD method; here we choose to use only the LDlinkR method.
  
  # Given a vector of variants (rsIDs) needing proxies
  # (from find_variants_needing_proxies) and an LD reference file,
  # output a data frame linking each variant to a data frame containing possible
  # proxies (variant ID + r^2 + alleles)
  # Criteria for eligibility:
  #   Not strand-ambiguous
  #   Trait fraction >= 80%
  #   r^2 >= 0.8 with the index variant
  # Choose based on first trait count, then r^2
  
  # First, run "/path/to/tabix /path/to/LDfile rsID_1 rsID_2 ...
  print(paste("Num rows need_proxies:",nrow(need_proxies)))
  if (method %in% c("LDlink","LDlinkR","LDproxy")) {
    
    print("Using LDlinkR:LDproxy_batch to find proxies...")
    need_proxies <- need_proxies %>%
      separate(hm_variant_id, into=c("CHR","POS","REF","ALT"),sep = "_",remove = F) %>%
      mutate(query_snp = paste0("chr", CHR, ":", POS)) %>%
      select(-c(CHR, POS, REF, ALT))
    need_proxies_snps <- need_proxies$query_snp
    
    LDlinkR::LDproxy_batch(need_proxies_snps,
                           pop = population,
                           r2d = "r2",
                           token = LDlink_token,
                           append = T,
                           genome_build = "grch38")
    proxy_df <- read.table("./combined_query_snp_list_grch38.txt",sep = "\t",row.names = NULL) %>%
      filter(R2>r2_num) %>%
      filter(!Coord %in% need_proxies_snps) %>%
      inner_join(need_proxies, by = "query_snp") %>%
      arrange(PVALUE) %>%
      filter(!duplicated(RS_Number)) %>%
      dplyr::select(hm_rsid, proxy_rsID=RS_Number, r2=R2) 
    need_proxies <- need_proxies %>%
      select(-c(query_snp))
    
  } else if (method=="TopLD") { # use TopLD
    # print(sprintf("Using TopLD to find proxies for %s!", population))
    # if (nrow(need_proxies)<100) {
    #   write(need_proxies$rsID, "need_proxies_rsIDs.tmp")
    #   system(sprintf("%s -thres %.1f -pop %s -maf 0.01 -inFile need_proxies_rsIDs.tmp -outputLD outputLD.txt -outputInfo outputInfo.txt", topLD_path, r2_num, population))
    # } else { # need to split up
    #   print("Splitting proxy df into subsets (more than 100 SNPs)...")
    #   proxy_df_list <- split(need_proxies, (seq(nrow(need_proxies))-1) %/% 100)
    #   
    #   system("touch outputLD.txt")
    #   print("Running TopLD for proxy df segments...")
    #   for (j in 1:length(proxy_df_list)){
    #     print(sprintf("Querying LD subset %i/%i",j, length(proxy_df_list)))
    #     df <- proxy_df_list[[j]]
    #     write(df$rsID, "need_proxies_rsIDs.tmp")
    #     system(sprintf("%s -thres %.1f -pop %s -maf 0.01 -inFile need_proxies_rsIDs.tmp -outputLD outputLD_temp.txt -outputInfo outputInfo.txt", topLD_path, r2_num, population))
    #     system("cat outputLD_temp.txt >> outputLD.txt")
    #   }
    # }
    # proxy_df <- fread("outputLD.txt", stringsAsFactors = F, data.table = F) %>%
    #   select(rsID=rsID1, proxy_rsID=rsID2, r2=R2) %>%
    #   subset(proxy_rsID %like% "rs")
  } 
  else {
    stop("Enter appropriate proxy search method: LDlinkR !") # Using stop function
    
  }
  print(paste("No. possible proxies found:",nrow(proxy_df))) # proxy_df should have columns (rsID, proxy_rsID, r2)
  write(proxy_df$proxy_rsID, "./test_results/potential_proxies_rsid.tmp")
  
  if (nrow(proxy_df)>0) {
    print("Creating proxy rsID map...")
    potential_proxies_map <- fread(cmd=paste0("grep -wFf ./test_results/potential_proxies_rsid.tmp \"", rsID_map_file, "\""),
                                   header=F, col.names=c("proxy_VAR_ID", "proxy_rsID"),
                                   data.table=F, stringsAsFactors=F) %>% distinct()
    print(head(potential_proxies_map))
    
    proxy_variants <- potential_proxies_map$proxy_VAR_ID
    
    proxy_missingness <- count_traits_per_variant(
      proxy_variants,
      trait_ss_files,
      sample_size = sample_size_traitGWAS,
      savepath_varid="./test_results/all_snps_varids_proxies.tmp"
    )

    # get proxy missingness
    df_Ns_rev <- proxy_missingness %>%
      pivot_wider(names_from="trait", values_from="N") %>%
      data.frame() %>%
      column_to_rownames("hm_variant_id") %>%
      set_colnames(names(trait_ss_files)) # proxy_missingness %>% column_to_rownames("hm_variant_id")
    df_Ns_rev[df_Ns_rev == 'NULL'] <- NA
    
    # get variant counts
    variant_counts_df <- data.frame(VAR_ID=rownames(df_Ns_rev),
                                    frac=rowSums(!is.na(df_Ns_rev))/length(trait_ss_files)) # the matrix for missingness ratio
    
    proxy_missingness_new <- ifelse(
      proxy_variants %in% variant_counts_df$VAR_ID,
      variant_counts_df$frac[match(proxy_variants, variant_counts_df$VAR_ID)],  # If in counts data frame, take the non-missing fraction
      0  # If not in data frame, then the non-missing fraction is 0
    )
    proxy_missingness_new <- setNames(proxy_missingness_new, proxy_variants)
    
    proxy_missingness_df <- tibble(
      proxy_VAR_ID=names(proxy_missingness_new),
      frac_nonmissing=proxy_missingness_new
    )
    
    final_proxy_df <- proxy_df %>%
      inner_join(potential_proxies_map, by="proxy_rsID") %>%
      separate(proxy_VAR_ID, into=c("CHR", "POS", "REF", "ALT"), 
               sep="_", remove=F) %>%
      inner_join(proxy_missingness_df, by="proxy_VAR_ID") %>%
      filter(
        !(paste0(REF, ALT) %in% c("AT", "TA", "CG", "GC")),  # Not strand-ambiguous
        !grepl("^[0-9]+_[0-9]+_[ACGT]+_[ACGT]+,[ACGT]+$", proxy_VAR_ID),  # Not multi-allelic
        frac_nonmissing >= frac_nonmissing_num,  # Sufficient fraction of traits non-missing
        r2 >= r2_num  # Sufficient LD with the proxied variant
      ) %>%
      group_by(hm_rsid) %>%
      arrange(desc(frac_nonmissing),
              desc(r2),
              CHR) %>%  # Arbitrary sort for reproducibility in case of missingness + r2 ties
      dplyr::slice(1) %>%
      ungroup() %>%
      inner_join(need_proxies, by="hm_rsid") %>%  # added to include orig VAR_ID in output 
      data.frame()
  } else {
    final_proxy_df <- NULL
  }
  proxies_found <- final_proxy_df$hm_rsid
  
  no_proxies_found <- setdiff(need_proxies$hm_rsid, proxies_found)
  print(paste0("No proxies needed for ", 
               length(setdiff(pruned_variants$VAR_ID, need_proxies$hm_variant_id)), 
               " variants."))
  print(paste0("Proxies found for ", length(proxies_found), " variants."))
  print(paste0("No adequate proxies found for ", length(no_proxies_found), 
               " variants."))
  
  return(final_proxy_df)
  
}

