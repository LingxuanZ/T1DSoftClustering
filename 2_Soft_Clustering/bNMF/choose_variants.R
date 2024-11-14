
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


count_traits_per_variant <- function(gwas_variants, ss_files) {
  
  # Given a vector of variants and a named vector of summary statistics files
  # for traits to be clustered, output a vector of non-missing trait fractions
  # per variant
  
  print("Assessing variant missingness across traits...")
  write(gwas_variants, "all_snps_varids.tmp")
  
  rename_cols <- c(N_PH="N")
  
  variant_df_list <- lapply(1:length(ss_files), function(i) {
    print(sprintf("...Reading %s...", names(ss_files)[i]))
    
    headers <- as.character(fread(ss_files[i], nrows=1,
                                  data.table=F, stringsAsFactors=F, header=F))
    
    if (endsWith(ss_files[i],".gz")) {
      df <- fread(cmd=sprintf("gzip -cd %s | fgrep -wf all_snps_varids.tmp ",ss_files[i]),
                  header=F,
                  col.names=headers,
                  data.table=F,
                  stringsAsFactors=F) %>%
        rename(any_of(rename_cols))
      
    } else {
      df <- fread(cmd=sprintf("fgrep -wf all_snps_varids.tmp %s ",ss_files[i]),
                  header=F,
                  col.names=headers,
                  data.table=F,
                  stringsAsFactors=F) %>%
        rename(any_of(rename_cols)) 
      
    }
    print(nrow(df))
    return(df)
  })
  
  # make dataframe of Ns
  df_N <- variant_df_list %>%
    setNames(names(ss_files)) %>%
    bind_rows(.id="trait") %>%
    select(trait, VAR_ID, N_PH) %>%
    pivot_wider(names_from="trait", values_from="N_PH") %>%
    data.frame()
}


