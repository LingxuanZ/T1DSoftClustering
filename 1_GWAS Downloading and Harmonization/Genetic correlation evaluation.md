This file is intended to quantify and evaluate the genetic correlation between different trait GWASs, as the number of immune cell GWASs is too large and the descriptions of the traits appear to be highly correlated.

Based on methods summarized in the paper "[Comparison of methods for estimating genetic correlation between complex traits using GWAS summary statistics](https://pmc.ncbi.nlm.nih.gov/articles/PMC8425307)", I chose to use the methods LDSC, GNOVA and HDL to quantify the genetic correlation between immune cell GWASs and between Autoimmune_Rheumatoid_Arthritis GWASs.

- LDSC is the most widely used method in terms of genome-wide genetic correlation estimation
- GNOVA estimator tends to be more accurate than LDSC estimator.
- High-definition likelihood (HDL), a recently proposed method, also outperformed LDSC in simulations. 

# 1 LDSC method

ReproGen Consortium, ReproGen Consortium, Psychiatric Genomics Consortium, et al. An atlas of genetic correlations across human diseases and traits. Nat Genet 2015;47(11):1236.[[Google Scholar](https://scholar.google.com/scholar_lookup?journal=Nat Genet&title=An atlas of genetic correlations across human diseases and traits&volume=47&issue=11&publication_year=2015&pages=1236&pmid=26414676&doi=10.1038/ng.3406&)]

**Github**: https://github.com/bulik/ldsc # Not applicable because the original software on Github was written in Python 2, while Python 3 is currently used. 

**PyPi**: https://pypi.org/project/ldsc/#files # The latest LDSC software can be used in Python 3.

**Limitation**:  When there is genetic heterogeneity between the actual sample and reference data from which LD scores are estimated, the accuracy of LDSC decreases further.

```sh
### {sh}
cd "/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Genetic Corr/software"
cd ldsc

module load python
conda install pip
pip install scipy

# python3 -m venv tutorial_env
# source tutorial_env/bin/activate
# python3 -m venv ldsc
# source ldsc/bin/activate


python3 -m pip install ldsc
python3 -m ldsc --help 
python3 -m munge_sumstats --help 
```

```sh
### {sh}
python3 -m munge_sumstats \
    --sumstats /scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data/Immune_Cell_1_32929287-GCST90001391-EFO_0007937.h.tsv.gz \
    --out formatted_sumstats_1 \
    --merge-alleles w_hm3.snplist #####
    
python3 -m ldsc \
    --rg "/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data/Immune_Cell_1_32929287-GCST90001391-EFO_0007937.h.tsv.gz","/scratch/scjp_root/scjp0/zhulx/T1D Soft Clustering/Data/GWAS summary stats/Original Data/Immune_Cell_10_32929287-GCST90001400-EFO_0007937.h.tsv.gz" \
    --ref-ld-chr eur_w_ld_chr/ \ ####
    --w-ld-chr eur_w_ld_chr/ \ ######
    --out genetic_corr_1
```

still dealing with the software

# 2 GNOVA