## Load plink before starting R
# module load plink/1.90b6.6
# R

## Now run R code
library('data.table')
library('SummarizedExperiment')
library("here")
library('sessioninfo')

## To avoid issues with running this code on qsub
data.table::setDTthreads(threads = 1)

## Find the samples for this project
load("/dcl01/lieber/ajaffe/lab/Nicotine/NAc/RNAseq/paired_end_n239/count_data/NAc_Nicotine_hg38_rseGene_rawCounts_allSamples_n239.rda", verbose = TRUE)
stopifnot(length(unique(rse_gene$BrNum)) == ncol(rse_gene))

## For converting BrNum's into numbers
brnumerical <- function(x) {
    as.integer(gsub("Br|_.*", "", x))
}

libd_bfile <- "/dcl01/lieber/ajaffe/Brain/Imputation/Subj_Cleaned/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38"

## Read the LIBD fam data
libd_fam <- fread(
    paste0(libd_bfile, ".fam"),
    col.names = c('famid', 'w_famid', 'w_famid_fa', 'w_famid_mo', 'sex_code', 'phenotype'))
libd_fam$brnumerical <- brnumerical(libd_fam$famid)
setkey(libd_fam, 'brnumerical')

## Filter the LIBD data to the one specific to this project
# region <- "NAc_Nicotine"
message(paste(Sys.time(), 'processing', "NAc_Nicotine"))
samp_file <- paste0('samples_to_extract_', "NAc_Nicotine", '.txt')

## Which NAc samples have genotype data?
filter_m <- match(brnumerical(rse_gene$BrNum), libd_fam$brnumerical)
message(Sys.time(), " NAc BrNums that have genotype data")
table(!is.na(filter_m))

fwrite(
    libd_fam[filter_m[!is.na(filter_m)], 1:2], ## can be more involved
    file = samp_file,
    sep = '\t', col.names = FALSE
)
newbfile_root <- 'LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38_filtered_NAc_Nicotine'

dir.create("duplicate_snps_bim", showWarnings = FALSE)
newbfile <- here::here("twas", "filter_data", "duplicate_snps_bim", paste0(
    newbfile_root,
    "_duplicateSNPs"
))

## Extract
message(paste(Sys.time(), 'running bfile extract for', newbfile))
system(paste("plink --bfile", libd_bfile,
	"--keep", samp_file, "--make-bed --out",
	newbfile, " --memory 100000 --biallelic-only"))

## Check that we have the right data
newbfile_fam <- fread(paste0(newbfile, ".fam"),
    col.names = c("famid", "w_famid", "w_famid_fa", "w_famid_mo", "sex_code", "phenotype")
)
check_m <- match(brnumerical(newbfile_fam$famid), brnumerical(colData(rse_gene)$BrNum))
stopifnot(all(!is.na(check_m)))


## Re-run but now make the SNV names unique
dir.create("unique_snps_bim", showWarnings = FALSE)
newbfile_unique <- here::here("twas", "filter_data", "unique_snps_bim", paste0(
    newbfile_root,
    "_uniqueSNPs"
))

## Extract again (could also copy and rename, but it's fast this way)
message(paste(Sys.time(), 'running bfile extract for', newbfile_unique))
system(paste("plink --bfile", libd_bfile,
             "--keep", samp_file, "--make-bed --out",
             newbfile_unique, " --memory 100000 --biallelic-only"))


message(paste(Sys.time(), "reading the bim file", newbfile_unique))
bim <- fread(
    paste0(newbfile_unique, ".bim"),
    col.names = c("chr", "snp", "position", "basepair", "allele1", "allele2")
)

table(duplicated(bim$snp))
#    FALSE     TRUE
# 10943065    44114

## Make names unique
message(Sys.time(), " making the variant names unique")
bim$snp <- make.names(bim$snp, unique = TRUE)
stopifnot(all(!duplicated(bim$snp)))

## Ovewrite the PLINK bim file
fwrite(bim, file = paste0(newbfile_unique, ".bim"), sep = " ", col.names = FALSE)


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value                                      
# version  R version 4.0.2 Patched (2020-06-24 r78746)
# os       CentOS Linux 7 (Core)                      
# system   x86_64, linux-gnu                          
# ui       X11                                        
# language (EN)                                       
# collate  en_US.UTF-8                                
# ctype    en_US.UTF-8                                
# tz       US/Eastern                                 
# date     2020-09-04                                 
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date       lib source                            
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.0.0)                    
# backports              1.1.6    2020-04-05 [1] CRAN (R 4.0.0)                    
# Biobase              * 2.48.0   2020-04-27 [1] Bioconductor                      
# BiocGenerics         * 0.34.0   2020-04-27 [1] Bioconductor                      
# bitops                 1.0-6    2013-08-17 [2] CRAN (R 4.0.0)                    
# cli                    2.0.2    2020-02-28 [1] CRAN (R 4.0.0)                    
# colorout             * 1.2-2    2020-05-08 [1] Github (jalvesaq/colorout@726d681)
# colorspace             1.4-1    2019-03-18 [2] CRAN (R 4.0.0)                    
# crayon                 1.3.4    2017-09-16 [1] CRAN (R 4.0.0)                    
# data.table           * 1.12.8   2019-12-09 [1] CRAN (R 4.0.0)                    
# DelayedArray         * 0.14.0   2020-04-27 [1] Bioconductor                      
# digest                 0.6.25   2020-02-23 [1] CRAN (R 4.0.0)                    
# dplyr                  0.8.5    2020-03-07 [1] CRAN (R 4.0.0)                    
# ellipsis               0.3.0    2019-09-20 [1] CRAN (R 4.0.0)                    
# fansi                  0.4.1    2020-01-08 [1] CRAN (R 4.0.0)                    
# GenomeInfoDb         * 1.24.0   2020-04-27 [1] Bioconductor                      
# GenomeInfoDbData       1.2.3    2020-05-18 [2] Bioconductor                      
# GenomicRanges        * 1.40.0   2020-04-27 [1] Bioconductor                      
# ggplot2                3.3.0    2020-03-05 [1] CRAN (R 4.0.0)                    
# glue                   1.4.0    2020-04-03 [1] CRAN (R 4.0.0)                    
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.0.0)                    
# here                 * 0.1      2017-05-28 [1] CRAN (R 4.0.0)                    
# htmltools              0.4.0    2019-10-04 [1] CRAN (R 4.0.0)                    
# htmlwidgets            1.5.1    2019-10-08 [1] CRAN (R 4.0.0)                    
# httpuv                 1.5.2    2019-09-11 [1] CRAN (R 4.0.0)                    
# IRanges              * 2.22.1   2020-04-28 [1] Bioconductor                      
# jsonlite               1.6.1    2020-02-02 [2] CRAN (R 4.0.0)                    
# later                  1.0.0    2019-10-04 [1] CRAN (R 4.0.0)                    
# lattice                0.20-41  2020-04-02 [3] CRAN (R 4.0.2)                    
# lifecycle              0.2.0    2020-03-06 [1] CRAN (R 4.0.0)                    
# magrittr               1.5      2014-11-22 [1] CRAN (R 4.0.0)                    
# Matrix                 1.2-18   2019-11-27 [3] CRAN (R 4.0.2)                    
# matrixStats          * 0.56.0   2020-03-13 [1] CRAN (R 4.0.0)                    
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.0.0)                    
# pillar                 1.4.4    2020-05-05 [1] CRAN (R 4.0.0)                    
# pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 4.0.0)                    
# png                    0.1-7    2013-12-03 [2] CRAN (R 4.0.0)                    
# promises               1.1.0    2019-10-04 [1] CRAN (R 4.0.0)                    
# purrr                  0.3.4    2020-04-17 [1] CRAN (R 4.0.0)                    
# R6                     2.4.1    2019-11-12 [2] CRAN (R 4.0.0)                    
# Rcpp                   1.0.4.6  2020-04-09 [1] CRAN (R 4.0.0)                    
# RCurl                  1.98-1.2 2020-04-18 [2] CRAN (R 4.0.0)                    
# rlang                  0.4.6    2020-05-02 [1] CRAN (R 4.0.0)                    
# rmote                * 0.3.4    2020-05-08 [1] Github (cloudyr/rmote@fbce611)    
# rprojroot              1.3-2    2018-01-03 [2] CRAN (R 4.0.0)                    
# S4Vectors            * 0.26.0   2020-04-27 [1] Bioconductor                      
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.0.0)                    
# servr                  0.16     2020-03-02 [1] CRAN (R 4.0.0)                    
# sessioninfo          * 1.1.1    2018-11-05 [1] CRAN (R 4.0.0)                    
# SummarizedExperiment * 1.18.1   2020-04-30 [1] Bioconductor                      
# tibble                 3.0.1    2020-04-20 [1] CRAN (R 4.0.0)                    
# tidyselect             1.1.0    2020-05-11 [2] CRAN (R 4.0.0)                    
# vctrs                  0.2.4    2020-03-10 [1] CRAN (R 4.0.0)                    
# withr                  2.2.0    2020-04-20 [1] CRAN (R 4.0.0)                    
# xfun                   0.13     2020-04-13 [1] CRAN (R 4.0.0)                    
# XVector                0.28.0   2020-04-27 [1] Bioconductor                      
# zlibbioc               1.34.0   2020-04-27 [1] Bioconductor                      
# 
# [1] /users/lcollado/R/4.0
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0/R/4.0/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0/R/4.0/lib64/R/library
system('plink --version')
# PLINK v1.90b6.6 64-bit (10 Oct 2018)
