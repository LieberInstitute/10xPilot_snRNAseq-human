## Load plink before starting R
# module load plink/1.90b6.6
# R

## Now run R code
library("data.table")
library("SummarizedExperiment")
library("here")
library("recount")
library("sva")
library("sessioninfo")

dir.create("rda", showWarnings = FALSE)

## To avoid issues with running this code on qsub
data.table::setDTthreads(threads = 1)

## Find the samples for this project
load("/dcl01/lieber/ajaffe/lab/Nicotine/NAc/RNAseq/paired_end_n239/count_data/NAc_Nicotine_hg38_rseGene_rawCounts_allSamples_n239.rda", verbose = TRUE)
stopifnot(length(unique(rse_gene$BrNum)) == ncol(rse_gene))

## Note that some rse_gene$BrNums have 0s that the rownames(mds)
## don't have
load("/dcl01/lieber/ajaffe/lab/Nicotine/NAc/RNAseq/paired_end_n239/genotype_data/Nicotine_NAc_Genotypes_n206_mds.rda", verbose = TRUE)

## For converting BrNum's into numbers
brnumerical <- function(x) {
    as.integer(gsub("Br|_.*", "", x))
}

libd_bfile <- "/dcl01/lieber/ajaffe/Brain/Imputation/Subj_Cleaned/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38"

## Read the LIBD fam data
libd_fam <- fread(
    paste0(libd_bfile, ".fam"),
    col.names = c("famid", "w_famid", "w_famid_fa", "w_famid_mo", "sex_code", "phenotype")
)
libd_fam$brnumerical <- brnumerical(libd_fam$famid)
setkey(libd_fam, "brnumerical")

## Filter the LIBD data to the one specific to this project
# region <- "NAc_Nicotine"
message(paste(Sys.time(), "processing", "NAc_Nicotine"))
samp_file <- paste0("samples_to_extract_", "NAc_Nicotine", ".txt")

## Which NAc samples have genotype data and MDS data?
samples_in_all <- intersect(
    intersect(brnumerical(rse_gene$BrNum), libd_fam$brnumerical),
    brnumerical(rownames(mds))
)

################################
## Subset and save all key files
################################
rse_gene <- rse_gene[, brnumerical(rse_gene$BrNum) %in% samples_in_all]

## Match rse_gene to mds
m_to_mds <- match(brnumerical(rse_gene$BrNum), brnumerical(rownames(mds)))
stopifnot(all(!is.na(m_to_mds)))
mds <- mds[m_to_mds, ]

## Append mds to colData
colData(rse_gene) <- cbind(colData(rse_gene), mds)

## Compute RPKM
assays(rse_gene)$RPKM <- getRPKM(rse_gene, "Length")

## Compute gene PCs
message(Sys.time(), " computing gene PCs on log2(RPKM + 1)")
pcaGene <- prcomp(t(log2(assays(rse_gene)$RPKM + 1)))
save(pcaGene, file = "rda/pcaGene.Rdata")

message(Sys.time(), " determine how many gene PCs to adjust for")
mod <- model.matrix(~ Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, data = colData(rse_gene))
kGene <- num.sv(log2(assays(rse_gene)$RPKM + 1), mod)
stopifnot(kGene > 0)
genePCs <- pcaGene$x[, seq_len(kGene)]
save(genePCs, file = "rda/genePCs.Rdata")

## Add gene PCs to rse_gene
colData(rse_gene) <- cbind(colData(rse_gene), genePCs)

## Save for later
save(rse_gene, file = paste0("rda/NAc_Nicotine_hg38_rseGene_rawCounts_allSamples_n", ncol(rse_gene), ".Rdata"))

## Now extract the genotype data too
filter_m <- match(brnumerical(rse_gene$BrNum), libd_fam$brnumerical)
stopifnot(all(!is.na(filter_m)))
fwrite(
    libd_fam[filter_m, 1:2], ## can be more involved
    file = samp_file,
    sep = "\t", col.names = FALSE
)
newbfile_root <- "LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38_filtered_NAc_Nicotine"

dir.create("duplicate_snps_bim", showWarnings = FALSE)
newbfile <- here::here("twas", "filter_data", "duplicate_snps_bim", paste0(
    newbfile_root,
    "_duplicateSNPs"
))

## Extract
message(paste(Sys.time(), "running bfile extract for", newbfile))
system(paste(
    "plink --bfile", libd_bfile,
    "--keep", samp_file, "--make-bed --out",
    newbfile, " --memory 100000 --biallelic-only"
))

# PLINK v1.90b6.6 64-bit (10 Oct 2018)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to /dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/twas/filter_data/duplicate_snps_bim/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38_filtered_NAc_Nicotine_duplicateSNPs.log.
# Options in effect:
#   --bfile /dcl01/lieber/ajaffe/Brain/Imputation/Subj_Cleaned/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38
#   --biallelic-only
#   --keep samples_to_extract_NAc_Nicotine.txt
#   --make-bed
#   --memory 100000
#   --out /dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/twas/filter_data/duplicate_snps_bim/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38_filtered_NAc_Nicotine_duplicateSNPs
#
# 257850 MB RAM detected; reserving 100000 MB for main workspace.
# 10987179 variants loaded from .bim file.
# 2054 people (1291 males, 763 females) loaded from .fam.
# --keep: 205 people remaining.
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 205 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 45 het. haploid genotypes present (see
# /dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/twas/filter_data/duplicate_snps_bim/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38_filtered_NAc_Nicotine_duplicateSNPs.hh
# ); many commands treat these as missing.
# Total genotyping rate in remaining samples is 0.978247.
# 10987179 variants and 205 people pass filters and QC.
# Note: No phenotypes present.
# --make-bed to
# /dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/twas/filter_data/duplicate_snps_bim/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38_filtered_NAc_Nicotine_duplicateSNPs.bed
# +
# /dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/twas/filter_data/duplicate_snps_bim/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38_filtered_NAc_Nicotine_duplicateSNPs.bim
# +
# /dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/twas/filter_data/duplicate_snps_bim/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38_filtered_NAc_Nicotine_duplicateSNPs.fam
# ... done.


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
message(paste(Sys.time(), "running bfile extract for", newbfile_unique))
system(paste(
    "plink --bfile", libd_bfile,
    "--keep", samp_file, "--make-bed --out",
    newbfile_unique, " --memory 100000 --biallelic-only"
))


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
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.0.2 Patched (2020-06-24 r78746)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2020-09-04
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version  date       lib source
#  acepack                1.4.1    2016-10-29 [2] CRAN (R 4.0.0)
#  annotate               1.66.0   2020-04-27 [1] Bioconductor
#  AnnotationDbi          1.50.0   2020-04-27 [1] Bioconductor
#  askpass                1.1      2019-01-13 [1] CRAN (R 4.0.0)
#  assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.0.0)
#  backports              1.1.6    2020-04-05 [1] CRAN (R 4.0.0)
#  base64enc              0.1-3    2015-07-28 [2] CRAN (R 4.0.0)
#  Biobase              * 2.48.0   2020-04-27 [1] Bioconductor
#  BiocFileCache          1.12.0   2020-04-27 [1] Bioconductor
#  BiocGenerics         * 0.34.0   2020-04-27 [1] Bioconductor
#  BiocParallel         * 1.22.0   2020-04-27 [1] Bioconductor
#  biomaRt                2.44.0   2020-04-27 [1] Bioconductor
#  Biostrings             2.56.0   2020-04-27 [1] Bioconductor
#  bit                    1.1-15.2 2020-02-10 [2] CRAN (R 4.0.0)
#  bit64                  0.9-7    2017-05-08 [2] CRAN (R 4.0.0)
#  bitops                 1.0-6    2013-08-17 [2] CRAN (R 4.0.0)
#  blob                   1.2.1    2020-01-20 [2] CRAN (R 4.0.0)
#  BSgenome               1.56.0   2020-04-27 [1] Bioconductor
#  bumphunter             1.30.0   2020-04-27 [1] Bioconductor
#  checkmate              2.0.0    2020-02-06 [1] CRAN (R 4.0.0)
#  cli                    2.0.2    2020-02-28 [1] CRAN (R 4.0.0)
#  cluster                2.1.0    2019-06-19 [3] CRAN (R 4.0.2)
#  codetools              0.2-16   2018-12-24 [3] CRAN (R 4.0.2)
#  colorout             * 1.2-2    2020-05-08 [1] Github (jalvesaq/colorout@726d681)
#  colorspace             1.4-1    2019-03-18 [2] CRAN (R 4.0.0)
#  crayon                 1.3.4    2017-09-16 [1] CRAN (R 4.0.0)
#  curl                   4.3      2019-12-02 [1] CRAN (R 4.0.0)
#  data.table           * 1.12.8   2019-12-09 [1] CRAN (R 4.0.0)
#  DBI                    1.1.0    2019-12-15 [2] CRAN (R 4.0.0)
#  dbplyr                 1.4.3    2020-04-19 [1] CRAN (R 4.0.0)
#  DelayedArray         * 0.14.0   2020-04-27 [1] Bioconductor
#  derfinder              1.22.0   2020-04-27 [1] Bioconductor
#  derfinderHelper        1.22.0   2020-04-27 [1] Bioconductor
#  digest                 0.6.25   2020-02-23 [1] CRAN (R 4.0.0)
#  doRNG                  1.8.2    2020-01-27 [2] CRAN (R 4.0.0)
#  downloader             0.4      2015-07-09 [1] CRAN (R 4.0.0)
#  dplyr                  0.8.5    2020-03-07 [1] CRAN (R 4.0.0)
#  edgeR                  3.30.0   2020-04-27 [1] Bioconductor
#  ellipsis               0.3.0    2019-09-20 [1] CRAN (R 4.0.0)
#  fansi                  0.4.1    2020-01-08 [1] CRAN (R 4.0.0)
#  foreach                1.5.0    2020-03-30 [2] CRAN (R 4.0.0)
#  foreign                0.8-80   2020-05-24 [3] CRAN (R 4.0.2)
#  Formula                1.2-3    2018-05-03 [2] CRAN (R 4.0.0)
#  genefilter           * 1.70.0   2020-04-27 [1] Bioconductor
#  GenomeInfoDb         * 1.24.0   2020-04-27 [1] Bioconductor
#  GenomeInfoDbData       1.2.3    2020-05-18 [2] Bioconductor
#  GenomicAlignments      1.24.0   2020-04-27 [1] Bioconductor
#  GenomicFeatures        1.40.0   2020-04-27 [1] Bioconductor
#  GenomicFiles           1.24.0   2020-04-27 [1] Bioconductor
#  GenomicRanges        * 1.40.0   2020-04-27 [1] Bioconductor
#  GEOquery               2.56.0   2020-04-27 [1] Bioconductor
#  ggplot2                3.3.0    2020-03-05 [1] CRAN (R 4.0.0)
#  glue                   1.4.0    2020-04-03 [1] CRAN (R 4.0.0)
#  gridExtra              2.3      2017-09-09 [2] CRAN (R 4.0.0)
#  gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.0.0)
#  here                 * 0.1      2017-05-28 [1] CRAN (R 4.0.0)
#  Hmisc                  4.4-0    2020-03-23 [1] CRAN (R 4.0.0)
#  hms                    0.5.3    2020-01-08 [2] CRAN (R 4.0.0)
#  htmlTable              1.13.3   2019-12-04 [1] CRAN (R 4.0.0)
#  htmltools              0.4.0    2019-10-04 [1] CRAN (R 4.0.0)
#  htmlwidgets            1.5.1    2019-10-08 [1] CRAN (R 4.0.0)
#  httpuv                 1.5.2    2019-09-11 [1] CRAN (R 4.0.0)
#  httr                   1.4.1    2019-08-05 [1] CRAN (R 4.0.0)
#  IRanges              * 2.22.1   2020-04-28 [1] Bioconductor
#  iterators              1.0.12   2019-07-26 [2] CRAN (R 4.0.0)
#  jpeg                   0.1-8.1  2019-10-24 [2] CRAN (R 4.0.0)
#  jsonlite               1.6.1    2020-02-02 [2] CRAN (R 4.0.0)
#  knitr                  1.28     2020-02-06 [1] CRAN (R 4.0.0)
#  later                  1.0.0    2019-10-04 [1] CRAN (R 4.0.0)
#  lattice                0.20-41  2020-04-02 [3] CRAN (R 4.0.2)
#  latticeExtra           0.6-29   2019-12-19 [2] CRAN (R 4.0.0)
#  lifecycle              0.2.0    2020-03-06 [1] CRAN (R 4.0.0)
#  limma                  3.44.1   2020-04-28 [1] Bioconductor
#  locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.0.0)
#  magrittr               1.5      2014-11-22 [1] CRAN (R 4.0.0)
#  Matrix                 1.2-18   2019-11-27 [3] CRAN (R 4.0.2)
#  matrixStats          * 0.56.0   2020-03-13 [1] CRAN (R 4.0.0)
#  memoise                1.1.0    2017-04-21 [2] CRAN (R 4.0.0)
#  mgcv                 * 1.8-31   2019-11-09 [3] CRAN (R 4.0.2)
#  munsell                0.5.0    2018-06-12 [2] CRAN (R 4.0.0)
#  nlme                 * 3.1-148  2020-05-24 [3] CRAN (R 4.0.2)
#  nnet                   7.3-14   2020-04-26 [3] CRAN (R 4.0.2)
#  openssl                1.4.1    2019-07-18 [1] CRAN (R 4.0.0)
#  pillar                 1.4.4    2020-05-05 [1] CRAN (R 4.0.0)
#  pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 4.0.0)
#  plyr                   1.8.6    2020-03-03 [2] CRAN (R 4.0.0)
#  png                    0.1-7    2013-12-03 [2] CRAN (R 4.0.0)
#  prettyunits            1.1.1    2020-01-24 [1] CRAN (R 4.0.0)
#  progress               1.2.2    2019-05-16 [1] CRAN (R 4.0.0)
#  promises               1.1.0    2019-10-04 [1] CRAN (R 4.0.0)
#  purrr                  0.3.4    2020-04-17 [1] CRAN (R 4.0.0)
#  qvalue                 2.20.0   2020-04-27 [1] Bioconductor
#  R6                     2.4.1    2019-11-12 [2] CRAN (R 4.0.0)
#  rappdirs               0.3.1    2016-03-28 [1] CRAN (R 4.0.0)
#  RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.0.0)
#  Rcpp                   1.0.4.6  2020-04-09 [1] CRAN (R 4.0.0)
#  RCurl                  1.98-1.2 2020-04-18 [2] CRAN (R 4.0.0)
#  readr                  1.3.1    2018-12-21 [1] CRAN (R 4.0.0)
#  recount              * 1.14.0   2020-04-27 [1] Bioconductor
#  rentrez                1.2.2    2019-05-02 [1] CRAN (R 4.0.0)
#  reshape2               1.4.4    2020-04-09 [2] CRAN (R 4.0.0)
#  rlang                  0.4.6    2020-05-02 [1] CRAN (R 4.0.0)
#  rmote                * 0.3.4    2020-05-08 [1] Github (cloudyr/rmote@fbce611)
#  rngtools               1.5      2020-01-23 [2] CRAN (R 4.0.0)
#  rpart                  4.1-15   2019-04-12 [3] CRAN (R 4.0.2)
#  rprojroot              1.3-2    2018-01-03 [2] CRAN (R 4.0.0)
#  Rsamtools              2.4.0    2020-04-27 [1] Bioconductor
#  RSQLite                2.2.0    2020-01-07 [2] CRAN (R 4.0.0)
#  rstudioapi             0.11     2020-02-07 [2] CRAN (R 4.0.0)
#  rtracklayer            1.48.0   2020-04-27 [1] Bioconductor
#  S4Vectors            * 0.26.0   2020-04-27 [1] Bioconductor
#  scales                 1.1.1    2020-05-11 [2] CRAN (R 4.0.0)
#  servr                  0.16     2020-03-02 [1] CRAN (R 4.0.0)
#  sessioninfo          * 1.1.1    2018-11-05 [1] CRAN (R 4.0.0)
#  stringi                1.4.6    2020-02-17 [2] CRAN (R 4.0.0)
#  stringr                1.4.0    2019-02-10 [1] CRAN (R 4.0.0)
#  SummarizedExperiment * 1.18.1   2020-04-30 [1] Bioconductor
#  survival               3.1-12   2020-04-10 [1] CRAN (R 4.0.0)
#  sva                  * 3.36.0   2020-04-27 [2] Bioconductor
#  tibble                 3.0.1    2020-04-20 [1] CRAN (R 4.0.0)
#  tidyr                  1.0.3    2020-05-07 [1] CRAN (R 4.0.0)
#  tidyselect             1.1.0    2020-05-11 [2] CRAN (R 4.0.0)
#  VariantAnnotation      1.34.0   2020-04-27 [1] Bioconductor
#  vctrs                  0.2.4    2020-03-10 [1] CRAN (R 4.0.0)
#  withr                  2.2.0    2020-04-20 [1] CRAN (R 4.0.0)
#  xfun                   0.13     2020-04-13 [1] CRAN (R 4.0.0)
#  XML                    3.99-0.3 2020-01-20 [2] CRAN (R 4.0.0)
#  xml2                   1.3.2    2020-04-23 [1] CRAN (R 4.0.0)
#  xtable                 1.8-4    2019-04-21 [2] CRAN (R 4.0.0)
#  XVector                0.28.0   2020-04-27 [1] Bioconductor
#  zlibbioc               1.34.0   2020-04-27 [1] Bioconductor
#
# [1] /users/lcollado/R/4.0
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0/R/4.0/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0/R/4.0/lib64/R/library
system("plink --version")
# PLINK v1.90b6.6 64-bit (10 Oct 2018)
