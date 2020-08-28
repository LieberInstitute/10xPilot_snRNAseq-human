## Load plink before starting R
# module load plink/1.90b6.6
# R

## Now run R code
library('data.table')
library('SummarizedExperiment')
library('sessioninfo')
library("tidyverse")

# system("wget https://www.dropbox.com/s/f60w61ifky3wkwh/20116_0.gwas.imputed_v3.both_sexes.tsv.bgz?dl=1")
# system("mv 20116_0.gwas.imputed_v3.both_sexes.tsv.bgz?dl=1 20116_0.gwas.imputed_v3.both_sexes.tsv.bgz")
# system("gunzip -c 20116_0.gwas.imputed_v3.both_sexes.tsv.bgz > 20116_0.gwas.imputed_v3.both_sexes.tsv")

gwas_test <- fread("/users/aseyedia/NAc_TWAS/20116_0.gwas.imputed_v3.both_sexes.tsv", sep="\t")

## get the snpMap that contains the hg38 positions
message(paste(Sys.time(), 'loading LIBD genotype data'))
load('/dcl01/lieber/ajaffe/Brain/Imputation/Subj_Cleaned/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38.Rdata', verbose = TRUE)

colnames(snpMap) <- tolower(colnames(snpMap))
colnames(snpMap)[colnames(snpMap) == 'pos'] <- 'basepair'
snpMap <- data.table(snpMap)
setkey(snpMap, 'chr', 'basepair')

table(gwas_test$variant %in% snpMap$snp)
# > table(gwas_test$variant %in% snpMap$snp)
#
#    FALSE     TRUE
# 13352549   438918

gwas_test <- tidyr::separate(gwas_test, variant, c("chr", "basepair", "main", "alt"), sep = ":")

## Change 23 to X; Y is not there
as.character(1:24) %in% unique(gwas_test$chr)
gwas_test$chr[gwas_test$chr == '23'] <- 'X'
stopifnot(all(c(as.character(1:22), 'X') %in% unique(gwas_test$chr)))

ids <- with(gwas_test, paste(chr, basepair, sep = ':'))

# this command returns false
identical(length(unique(ids)), length(ids))
# > table(duplicated(ids))
#
#    FALSE     TRUE
# 13773232    18235

ids_snpmap <- with(snpMap, paste(chr, basepair, sep = ':'))
length(unique(ids_snpmap))
length(ids_snpmap)

# > length(unique(ids_snpmap))
# [1] 10933669
# > length(ids_snpmap)
# [1] 10988411


table(ids %in% ids_snpmap)
# > table(ids %in% ids_snpmap)
#
#   FALSE    TRUE
# 5627759 8163708

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# - Session info -------------------------------------------------------------------------------------------------------
#  setting  value
#  version  R version 4.0.2 Patched (2020-06-24 r78746)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  C
#  ctype    C
#  tz       US/Eastern
#  date     2020-07-20
#
# - Packages -----------------------------------------------------------------------------------------------------------
#  package              * version  date       lib source
#  assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.0.0)
#  backports              1.1.8    2020-06-17 [2] CRAN (R 4.0.2)
#  Biobase              * 2.49.0   2020-04-27 [2] Bioconductor
#  BiocGenerics         * 0.35.4   2020-06-04 [2] Bioconductor
#  bitops                 1.0-6    2013-08-17 [2] CRAN (R 4.0.0)
#  blob                   1.2.1    2020-01-20 [2] CRAN (R 4.0.0)
#  broom                  0.7.0    2020-07-09 [2] CRAN (R 4.0.2)
#  cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.0.0)
#  cli                    2.0.2    2020-02-28 [2] CRAN (R 4.0.0)
#  colorspace             1.4-1    2019-03-18 [2] CRAN (R 4.0.0)
#  crayon                 1.3.4    2017-09-16 [2] CRAN (R 4.0.0)
#  data.table           * 1.12.8   2019-12-09 [2] CRAN (R 4.0.0)
#  DBI                    1.1.0    2019-12-15 [2] CRAN (R 4.0.0)
#  dbplyr                 1.4.4    2020-05-27 [2] CRAN (R 4.0.2)
#  DelayedArray         * 0.15.7   2020-07-14 [2] Bioconductor
#  dplyr                * 1.0.0    2020-05-29 [2] CRAN (R 4.0.2)
#  ellipsis               0.3.1    2020-05-15 [2] CRAN (R 4.0.0)
#  fansi                  0.4.1    2020-01-08 [2] CRAN (R 4.0.0)
#  forcats              * 0.5.0    2020-03-01 [2] CRAN (R 4.0.0)
#  fs                     1.4.2    2020-06-30 [2] CRAN (R 4.0.2)
#  generics               0.0.2    2018-11-29 [2] CRAN (R 4.0.0)
#  GenomeInfoDb         * 1.25.8   2020-07-03 [2] Bioconductor
#  GenomeInfoDbData       1.2.3    2020-05-13 [2] Bioconductor
#  GenomicRanges        * 1.41.5   2020-06-09 [2] Bioconductor
#  ggplot2              * 3.3.2    2020-06-19 [2] CRAN (R 4.0.2)
#  glue                   1.4.1    2020-05-13 [2] CRAN (R 4.0.0)
#  gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.0.0)
#  haven                  2.3.1    2020-06-01 [2] CRAN (R 4.0.2)
#  hms                    0.5.3    2020-01-08 [2] CRAN (R 4.0.0)
#  httr                   1.4.1    2019-08-05 [2] CRAN (R 4.0.0)
#  IRanges              * 2.23.10  2020-06-13 [2] Bioconductor
#  jsonlite               1.7.0    2020-06-25 [2] CRAN (R 4.0.2)
#  lattice                0.20-41  2020-04-02 [3] CRAN (R 4.0.2)
#  lifecycle              0.2.0    2020-03-06 [2] CRAN (R 4.0.0)
#  lubridate              1.7.9    2020-06-08 [1] CRAN (R 4.0.2)
#  magrittr               1.5      2014-11-22 [2] CRAN (R 4.0.0)
#  Matrix               * 1.2-18   2019-11-27 [3] CRAN (R 4.0.2)
#  matrixStats          * 0.56.0   2020-03-13 [2] CRAN (R 4.0.0)
#  modelr                 0.1.8    2020-05-19 [1] CRAN (R 4.0.2)
#  munsell                0.5.0    2018-06-12 [2] CRAN (R 4.0.0)
#  pillar                 1.4.6    2020-07-10 [2] CRAN (R 4.0.2)
#  pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.0.0)
#  purrr                * 0.3.4    2020-04-17 [2] CRAN (R 4.0.0)
#  R6                     2.4.1    2019-11-12 [2] CRAN (R 4.0.0)
#  Rcpp                   1.0.5    2020-07-06 [2] CRAN (R 4.0.2)
#  RCurl                  1.98-1.2 2020-04-18 [2] CRAN (R 4.0.0)
#  readr                * 1.3.1    2018-12-21 [2] CRAN (R 4.0.0)
#  readxl                 1.3.1    2019-03-13 [2] CRAN (R 4.0.0)
#  reprex                 0.3.0    2019-05-16 [1] CRAN (R 4.0.2)
#  rlang                  0.4.7    2020-07-09 [2] CRAN (R 4.0.2)
#  rstudioapi             0.11     2020-02-07 [2] CRAN (R 4.0.0)
#  rvest                  0.3.5    2019-11-08 [2] CRAN (R 4.0.0)
#  S4Vectors            * 0.27.12  2020-06-09 [2] Bioconductor
#  scales                 1.1.1    2020-05-11 [2] CRAN (R 4.0.0)
#  sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 4.0.0)
#  stringi                1.4.6    2020-02-17 [2] CRAN (R 4.0.0)
#  stringr              * 1.4.0    2019-02-10 [2] CRAN (R 4.0.0)
#  SummarizedExperiment * 1.19.6   2020-07-09 [2] Bioconductor
#  tibble               * 3.0.3    2020-07-10 [2] CRAN (R 4.0.2)
#  tidyr                * 1.1.0    2020-05-20 [2] CRAN (R 4.0.2)
#  tidyselect             1.1.0    2020-05-11 [2] CRAN (R 4.0.0)
#  tidyverse            * 1.3.0    2019-11-21 [1] CRAN (R 4.0.2)
#  vctrs                  0.3.2    2020-07-15 [2] CRAN (R 4.0.2)
#  withr                  2.2.0    2020-04-20 [2] CRAN (R 4.0.0)
#  xml2                   1.3.2    2020-04-23 [2] CRAN (R 4.0.0)
#  XVector                0.29.3   2020-06-25 [2] Bioconductor
#  zlibbioc               1.35.0   2020-04-27 [2] Bioconductor
#
# [1] /users/aseyedia/R/4.0.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/library


