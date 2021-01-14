library(data.table)
library(dplyr)
library(devtools)

setDTthreads(1)

# Load snpMap
load(
    "/dcl01/lieber/ajaffe/Brain/Imputation/Subj_Cleaned/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38.Rdata"
)

snpMap <- as.data.table(snpMap)

# Read the GWAS files
hg19_gwas_si <-
    fread("NAc_GWAS/SmokingInitiation.txt")

hg19_gwas_sc <-
    fread("NAc_GWAS/SmokingCessation.txt")

hg19_gwas_dpw <-
    fread("NAc_GWAS/DrinksPerWeek.txt")

hg19_gwas_cpd <-
    fread("NAc_GWAS/CigarettesPerDay.txt")

hg19_gwas_aoi <-
    fread("NAc_GWAS/AgeofInitiation.txt")

# Make key
snpMap$hg19_key <- paste0(snpMap$CHR, "_", snpMap$POS)

hg19_gwas_si$hg19_key <- paste0(hg19_gwas_si$CHROM, "_", hg19_gwas_si$POS)

hg19_gwas_sc$hg19_key <- paste0(hg19_gwas_sc$CHROM, "_", hg19_gwas_sc$POS)

hg19_gwas_dpw$hg19_key <- paste0(hg19_gwas_dpw$CHROM, "_", hg19_gwas_dpw$POS)

hg19_gwas_cpd$hg19_key <- paste0(hg19_gwas_cpd$CHROM, "_", hg19_gwas_cpd$POS)

hg19_gwas_aoi$hg19_key <- paste0(hg19_gwas_aoi$CHROM, "_", hg19_gwas_aoi$POS)

# Merge in the hg38 coordinates based on hg19 coordinates
hg38_gwas_si <-
    merge(hg19_gwas_si, snpMap[, c("chr_hg38", "pos_hg38", "hg19_key")], by = "hg19_key")

hg38_gwas_sc <-
    merge(hg19_gwas_sc, snpMap[, c("chr_hg38", "pos_hg38", "hg19_key")], by = "hg19_key")

hg38_gwas_dpw <-
    merge(hg19_gwas_dpw, snpMap[, c("chr_hg38", "pos_hg38", "hg19_key")], by = "hg19_key")

hg38_gwas_cpd <-
    merge(hg19_gwas_cpd, snpMap[, c("chr_hg38", "pos_hg38", "hg19_key")], by = "hg19_key")

hg38_gwas_aoi <-
    merge(hg19_gwas_aoi, snpMap[, c("chr_hg38", "pos_hg38", "hg19_key")], by = "hg19_key")

# Drop hg19 coords from gwas
hg38_gwas_si <- hg38_gwas_si[, -c("hg19_key", "CHROM", "POS")]
hg38_gwas_sc <- hg38_gwas_sc[, -c("hg19_key", "CHROM", "POS")]
hg38_gwas_dpw <- hg38_gwas_dpw[, -c("hg19_key", "CHROM", "POS")]
hg38_gwas_cpd <- hg38_gwas_cpd[, -c("hg19_key", "CHROM", "POS")]
hg38_gwas_aoi <- hg38_gwas_aoi[, -c("hg19_key", "CHROM", "POS")]

# Calculate Odds Ratio
hg38_gwas_si[, OR:= exp(BETA)]
hg38_gwas_sc[, OR:= exp(BETA)]
hg38_gwas_dpw[, OR := exp(BETA)]
hg38_gwas_cpd[, OR := exp(BETA)]
hg38_gwas_aoi[, OR := exp(BETA)]

# reorder columns
col_order <- c(
    "chr_hg38",
    "RSID",
    "pos_hg38",
    "REF",
    "ALT",
    "AF",
    "STAT",
    "PVALUE",
    "BETA",
    "SE",
    "N",
    "EFFECTIVE_N",
    "Number_of_Studies",
    "ANNO",
    "ANNOFULL"
)

hg38_gwas_si <- hg38_gwas_si[, ..col_order]
hg38_gwas_sc <- hg38_gwas_sc[, ..col_order]
hg38_gwas_dpw <- hg38_gwas_dpw[, ..col_order]
hg38_gwas_cpd <- hg38_gwas_cpd[, ..col_order]
hg38_gwas_aoi <- hg38_gwas_aoi[, ..col_order]

names(hg38_gwas_si)[1] <- "CHR"
names(hg38_gwas_sc)[1] <- "CHR"
names(hg38_gwas_dpw)[1] <- "CHR"
names(hg38_gwas_cpd)[1] <- "CHR"
names(hg38_gwas_aoi)[1] <- "CHR"

names(hg38_gwas_si)[3] <- "BP"
names(hg38_gwas_sc)[3] <- "BP"
names(hg38_gwas_dpw)[3] <- "BP"
names(hg38_gwas_cpd)[3] <- "BP"
names(hg38_gwas_aoi)[3] <- "BP"

# read bim file with unique rsIDs
uniq_bim <-
    fread(
        "/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/twas/filter_data/unique_snps_bim/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38_filtered_NAc_Nicotine_uniqueSNPs.bim"
    )

colnames(uniq_bim) <- c("CHR", "SNP", "dummy", "BP", "A1", "A2")

# check for overlap between gwas and unique bim file
uniq_bim$new_test <- paste0(uniq_bim$CHR, "_", uniq_bim$BP)

hg38_gwas_si$new_test <-
    paste0(gsub("chr", "" , hg38_gwas_si$CHR), "_", hg38_gwas_si$BP)

hg38_gwas_sc$new_test <-
    paste0(gsub("chr", "" , hg38_gwas_sc$CHR), "_", hg38_gwas_sc$BP)

hg38_gwas_dpw$new_test <-
    paste0(gsub("chr", "" , hg38_gwas_dpw$CHR), "_", hg38_gwas_dpw$BP)

hg38_gwas_cpd$new_test <-
    paste0(gsub("chr", "" , hg38_gwas_cpd$CHR), "_", hg38_gwas_cpd$BP)

hg38_gwas_aoi$new_test <-
    paste0(gsub("chr", "" , hg38_gwas_aoi$CHR), "_", hg38_gwas_aoi$BP)

# table(uniq_bim$new_test %in% hg38_gwas_si$new_test)
#
#   FALSE    TRUE
# 4154352 6832827

# It does not seem that doing this (lifting over each GWAS file to hg38) has
# done much for the overlap between the unique bim file and the GWAS. In fact,
# now I am doubting if it was hg19 to begin with.

# table(uniq_bim$new_test %in% hg38_gwas_sc$new_test)
#
#   FALSE    TRUE
# 4154352 6832827

# table(uniq_bim$new_test %in% hg38_gwas_sc$new_test)
#
#   FALSE    TRUE
# 3925017 7062162

# table(uniq_bim$new_test %in% hg38_gwas_dpw$new_test)
#
#   FALSE    TRUE
# 4095180 6891999

# table(uniq_bim$new_test %in% hg38_gwas_cpd$new_test)
#
#   FALSE    TRUE
# 4025959 6961220

# table(uniq_bim$new_test %in% hg38_gwas_aoi$new_test)
#
#   FALSE    TRUE
# 4033185 6953994

# merge in unique rsIDs
hg38_gwas_si <-
    merge(hg38_gwas_si[, -"RSID"], uniq_bim[, c("SNP", "new_test")], by = "new_test")

hg38_gwas_sc <-
    merge(hg38_gwas_sc[, -"RSID"], uniq_bim[, c("SNP", "new_test")], by = "new_test")

hg38_gwas_dpw <-
    merge(hg38_gwas_dpw[, -"RSID"], uniq_bim[, c("SNP", "new_test")], by = "new_test")

hg38_gwas_cpd <-
    merge(hg38_gwas_cpd[, -"RSID"], uniq_bim[, c("SNP", "new_test")], by = "new_test")

hg38_gwas_aoi <-
    merge(hg38_gwas_aoi[, -"RSID"], uniq_bim[, c("SNP", "new_test")], by = "new_test")

hg38_gwas_si$new_test <- NULL

hg38_gwas_sc$new_test <- NULL

hg38_gwas_dpw$new_test <- NULL

hg38_gwas_cpd$new_test <- NULL

hg38_gwas_aoi$new_test <- NULL

col_order <- c(
    "CHR",
    "SNP",
    "BP",
    "REF",
    "ALT",
    "AF",
    "STAT",
    "PVALUE",
    "BETA",
    "SE",
    "N",
    "EFFECTIVE_N",
    "Number_of_Studies",
    "ANNO",
    "ANNOFULL"
)

hg38_gwas_si <- hg38_gwas_si[, ..col_order]
hg38_gwas_sc <- hg38_gwas_sc[, ..col_order]
hg38_gwas_dpw <- hg38_gwas_dpw[, ..col_order]
hg38_gwas_cpd <- hg38_gwas_cpd[, ..col_order]
hg38_gwas_aoi <- hg38_gwas_aoi[, ..col_order]


# hg38_gwas_si$effect <- log(hg38_gwas_si$OR)
# hg38_gwas_sc$effect <- log(hg38_gwas_sc$OR)
# hg38_gwas_dpw$effect <- log(hg38_gwas_dpw$OR)
# hg38_gwas_cpd$effect <- log(hg38_gwas_cpd$OR)
# hg38_gwas_aoi$effect <- log(hg38_gwas_aoi$OR)
#
# hg38_gwas$Z <- hg38_gwas$effect / hg38_gwas$SE

# pdf(file = "PGC_BIP_hg38_hist.png.pdf", useDingbats = FALSE)
#
# hist(hg38_gwas$effect, color = "gold")
# hist(hg38_gwas$Z, color = "darkorange")
#
# dev.off()

hg38_gwas_si <-
    hg38_gwas_si[, c("SNP", "REF", "ALT", "EFFECTIVE_N", "PVALUE")]
colnames(hg38_gwas_si) <- c("SNP", "A1", "A2", "Neff", "PVALUE")

hg38_gwas_sc <-
    hg38_gwas_sc[, c("SNP", "REF", "ALT", "EFFECTIVE_N", "PVALUE")]
colnames(hg38_gwas_sc) <- c("SNP", "A1", "A2", "Neff", "PVALUE")

hg38_gwas_dpw <-
    hg38_gwas_dpw[, c("SNP", "REF", "ALT", "EFFECTIVE_N", "PVALUE")]
colnames(hg38_gwas_dpw) <- c("SNP", "A1", "A2", "Neff", "PVALUE")

hg38_gwas_cpd <-
    hg38_gwas_cpd[, c("SNP", "REF", "ALT", "EFFECTIVE_N", "PVALUE")]
colnames(hg38_gwas_cpd) <- c("SNP", "A1", "A2", "Neff", "PVALUE")

hg38_gwas_aoi <-
    hg38_gwas_aoi[, c("SNP", "REF", "ALT", "EFFECTIVE_N", "PVALUE")]
colnames(hg38_gwas_aoi) <- c("SNP", "A1", "A2", "Neff", "PVALUE")

dir.create("clean_gwas/", showWarnings = FALSE)

write.table(
    hg38_gwas_si,
    file = "clean_gwas/SmokingInitiation_Clean.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)

write.table(
    hg38_gwas_sc,
    file = "clean_gwas/SmokingCessation_Clean.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)

write.table(
    hg38_gwas_dpw,
    file = "clean_gwas/DrinksPerWeek_Clean.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)

write.table(
    hg38_gwas_cpd,
    file = "clean_gwas/CigarettesPerDay_Clean.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)

write.table(
    hg38_gwas_aoi,
    file = "clean_gwas/AgeofInitiation_Clean.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2021-01-04 17:01:29 EST"
# > proc.time()
#     user   system  elapsed
#  919.864   82.686 4641.836
# > options(width = 120)
# > session_info()
# ??? Session info ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#  setting  value
#  version  R version 4.0.3 Patched (2020-11-29 r79529)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2021-01-04
#
# ??? Packages ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#  package     * version date       lib source
#  assertthat    0.2.1   2019-03-21 [2] CRAN (R 4.0.3)
#  callr         3.5.1   2020-10-13 [2] CRAN (R 4.0.3)
#  cli           2.2.0   2020-11-20 [2] CRAN (R 4.0.3)
#  crayon        1.3.4   2017-09-16 [2] CRAN (R 4.0.3)
#  data.table  * 1.13.6  2020-12-30 [2] CRAN (R 4.0.3)
#  desc          1.2.0   2018-05-01 [2] CRAN (R 4.0.3)
#  devtools    * 2.3.2   2020-09-18 [2] CRAN (R 4.0.3)
#  digest        0.6.27  2020-10-24 [2] CRAN (R 4.0.3)
#  dplyr       * 1.0.2   2020-08-18 [2] CRAN (R 4.0.3)
#  ellipsis      0.3.1   2020-05-15 [2] CRAN (R 4.0.3)
#  fansi         0.4.1   2020-01-08 [2] CRAN (R 4.0.3)
#  fs            1.5.0   2020-07-31 [2] CRAN (R 4.0.3)
#  generics      0.1.0   2020-10-31 [2] CRAN (R 4.0.3)
#  glue          1.4.2   2020-08-27 [2] CRAN (R 4.0.3)
#  lifecycle     0.2.0   2020-03-06 [2] CRAN (R 4.0.3)
#  magrittr      2.0.1   2020-11-17 [2] CRAN (R 4.0.3)
#  memoise       1.1.0   2017-04-21 [2] CRAN (R 4.0.3)
#  pillar        1.4.7   2020-11-20 [2] CRAN (R 4.0.3)
#  pkgbuild      1.2.0   2020-12-15 [2] CRAN (R 4.0.3)
#  pkgconfig     2.0.3   2019-09-22 [2] CRAN (R 4.0.3)
#  pkgload       1.1.0   2020-05-29 [2] CRAN (R 4.0.3)
#  prettyunits   1.1.1   2020-01-24 [2] CRAN (R 4.0.3)
#  processx      3.4.5   2020-11-30 [2] CRAN (R 4.0.3)
#  ps            1.5.0   2020-12-05 [2] CRAN (R 4.0.3)
#  purrr         0.3.4   2020-04-17 [2] CRAN (R 4.0.3)
#  R6            2.5.0   2020-10-28 [2] CRAN (R 4.0.3)
#  remotes       2.2.0   2020-07-21 [2] CRAN (R 4.0.3)
#  rlang         0.4.10  2020-12-30 [2] CRAN (R 4.0.3)
#  rprojroot     2.0.2   2020-11-15 [2] CRAN (R 4.0.3)
#  sessioninfo   1.1.1   2018-11-05 [2] CRAN (R 4.0.3)
#  testthat      3.0.1   2020-12-17 [2] CRAN (R 4.0.3)
#  tibble        3.0.4   2020-10-12 [2] CRAN (R 4.0.3)
#  tidyselect    1.1.0   2020-05-11 [2] CRAN (R 4.0.3)
#  usethis     * 2.0.0   2020-12-10 [2] CRAN (R 4.0.3)
#  vctrs         0.3.6   2020-12-17 [2] CRAN (R 4.0.3)
#  withr         2.3.0   2020-09-22 [2] CRAN (R 4.0.3)
#
# [1] /users/aseyedia/R/4.0.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/library
