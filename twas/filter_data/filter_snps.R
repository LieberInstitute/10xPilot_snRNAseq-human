## Load plink before starting R
# module load plink/1.90b6.6
# R

## Now run R code
library('data.table')
library('SummarizedExperiment')
library('sessioninfo')

## get the snpMap that contains the hg38 positions
message(paste(Sys.time(), 'loading BSP2 genotype data'))
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda', verbose = TRUE)

colnames(snpMap) <- tolower(colnames(snpMap))
colnames(snpMap)[colnames(snpMap) == 'pos'] <- 'basepair'
snpMap <- data.table(snpMap)
setkey(snpMap, 'chr', 'basepair')


## Get their LD reference data
## Following the instructions at http://gusevlab.org/projects/fusion/#computing-your-own-functional-weights
# message(paste(Sys.time(), 'Downloading LD reference data'))
# system('wget https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2')
# # http://how-to.wikia.com/wiki/How_to_untar_a_tar_file_or_gzip-bz2_tar_file
# system('tar xvjf LDREF.tar.bz2')
#
# ## There's 22 bim files
# dir('LDREF', '.*bim$')
# stopifnot(length(dir('LDREF', '.*bim$')) == 22)
#
# ## Looks like it's one per chromosome
# system('head LDREF/1000G.EUR.1.bim')
# system('head LDREF/1000G.EUR.2.bim')
#
# their_bims <- dir('LDREF', '.*bim$', full.names = TRUE)
# names(their_bims) <- dir('LDREF', '.*bim$')

## Use their LD reference already ported to hg38 (well, hg19 coordinates but subsetted to those present in hg38)
their_bims_hg38 <- dir('/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38', '.*bim.original$', full.names = TRUE)
names(their_bims_hg38) <- dir('/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38', '.*bim.original$')
their_bims <- their_bims_hg38

ldref_bim <- do.call(rbind, lapply(their_bims, function(input_bim) {
    message(paste(Sys.time(), 'reading file', input_bim))
    fread(input_bim, col.names = c('chr', 'snp', 'position', 'basepair', 'allele1', 'allele2'))
}))


## Load our bim file
message(paste(Sys.time(), 'Loading BSP2 bim file'))
bfile <- '/dcl01/lieber/ajaffe/Brain/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2.bim'
bsp2_bim <- fread(
    bfile,
    col.names = c('chr', 'snp', 'position', 'basepair', 'allele1', 'allele2')
)

## Filter ours to their
message(paste(Sys.time(), 'Preparing to filter our bim file'))
setkey(ldref_bim, 'chr', 'basepair')
setkey(bsp2_bim, 'chr', 'basepair')

message(paste(Sys.time(), 'Filtering our bim file'))
bsp2_bim_filt <- bsp2_bim[.(ldref_bim$chr, ldref_bim$basepair)]

table(is.na(bsp2_bim_filt$snp))
#
#   FALSE
# 1022566

# This is where we realized that we need to drop the biallelic SNPs
paste('Out of the original', nrow(ldref_bim), 'SNPs,', nrow(bsp2_bim_filt), 'are present in BSP2. That is,', round(nrow(bsp2_bim_filt) / nrow(ldref_bim) * 100, 2), 'percent.')

pairs <- with(bsp2_bim_filt, paste(chr, basepair, sep = '_'))
table(duplicated(pairs))
#   FALSE    TRUE
# 1022547      20
options(width = 120)
bsp2_bim_filt[pairs %in% pairs[which(duplicated(pairs))], ]
dim(bsp2_bim_filt[pairs %in% pairs[which(duplicated(pairs))], ])
# [1] 39  6 ### Originally (when using hg19 and their snp names): # [1] 54  6

dim(bsp2_bim_filt)
# [1] 1022567       6 ### Originally: # [1] 1190349       6

## Drop these duplicated ones
bsp2_bim_filt <- bsp2_bim_filt[!pairs %in% pairs[which(duplicated(pairs))], ]
dim(bsp2_bim_filt)
# [1] 1022528       6 ### Originally: # [1] 1190295       6

## Are they in our snpMap? Do they have hg38 positions?
table(is.na(snpMap$pos_hg38))
#   FALSE    TRUE
# 7023286     574
snpMap <- snpMap[!is.na(snpMap$pos_hg38), ]

m_hg38 <- match(with(bsp2_bim_filt, paste(chr, basepair, sep = '-')), with(snpMap, paste(chr, basepair, sep = '-')))
table(is.na(m_hg38))
#   FALSE
# 1022527
bsp2_bim_filt <- bsp2_bim_filt[!is.na(m_hg38), ]
dim(bsp2_bim_filt)
# [1] 1022527       6

message(paste(Sys.time(), 'Write new filtered bim & filtered SNP info files'))
fwrite(bsp2_bim_filt,
    file = 'LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_LDfiltered.bim',
    sep = '\t', col.names = FALSE
)
filt_snps <- 'LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_LDfiltered_snps.txt'
fwrite(
    bsp2_bim_filt[, 'snp'],
    file = filt_snps,
    sep = '\t', col.names = FALSE
)

## Read in BSP2 fam data
bsp2_fam <- fread(
    '/dcl01/lieber/ajaffe/Brain/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2.fam', 
    col.names = c('famid', 'w_famid', 'w_famid_fa', 'w_famid_mo', 'sex_code', 'phenotype'))
setkey(bsp2_fam, 'famid')

## Load in BSP2 data, subset to those samples we used for the eQTL analyses
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
keepInd = which(colData(rse_gene)$Age > 13)
rse_gene = rse_gene[,keepInd]

## Write a new fam file per brain region
for(region in unique(colData(rse_gene)$Region)) {
    # region <- 'HIPPO'
    message(paste(Sys.time(), 'processing', region))
    samp_file <- paste0('samples_to_extract_', region, '.txt')
    fwrite(
        bsp2_fam[.(colData(rse_gene)$BrNum[colData(rse_gene)$Region == region]), 1:2],
        file = samp_file,
        sep = '\t', col.names = FALSE
    )
    newbfile <- paste0(
        '/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/filter_data/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_LDfiltered_',
        region
    )
    
    ## Extract
    message(paste(Sys.time(), 'running bfile extract'))
    ## --biallelic-only doesn't really do anything at this point
    ## since we already dropped the duplicated snps
    system(paste("plink --bfile", gsub('.bim', '', bfile), '--extract', filt_snps, 
    	"--keep", samp_file, "--make-bed --out", 
    	newbfile, " --memory 225000 --biallelic-only"))
    
    ## Change to hg38 positions
    
    ## Read the new file
    message(paste(Sys.time(), 'read the new bim file'))
    final_bim <- fread(paste0(newbfile, '.bim'),
        col.names = c('chr', 'snp', 'position', 'basepair', 'allele1', 'allele2'),
        colClasses = c('character', 'character', 'numeric', 'integer', 'character', 'character')
    )
    ## Keep a copy of the original one
    system(paste0('mv ', newbfile, '.bim ', newbfile, '.bim.original'))
    
    ## Complete matching code later
    m_final <- match(with(final_bim, paste(chr, basepair, sep = '-')), with(snpMap, paste(chr, basepair, sep = '-')))
    stopifnot(!any(is.na(m_final)))
    final_bim$basepair <- snpMap$pos_hg38[m_final]
    ## Also use our snp names
    final_bim$snp <- snpMap$snp[m_final]
    
    message(paste(Sys.time(), 'Write new filtered bim file for region', region))
    fwrite(
        final_bim,
        file = paste0(newbfile, '.bim'),
        sep = '\t', col.names = FALSE
    )
}


## Checks that lead us to drop the duplicated snps from the filt_snps
## file since --biallelic-only was taking second priority to the
## --extract list
hippo_bim <- fread(
    paste0(newbfile, '.bim.original'),
    col.names = c('chr', 'snp', 'position', 'basepair', 'allele1', 'allele2')
)
setkey(hippo_bim, 'chr', 'basepair')

dim(hippo_bim)
dim(ldref_bim)
dim(hippo_bim[.(ldref_bim$chr, ldref_bim$basepair)])
dim(ldref_bim[.(hippo_bim$chr, hippo_bim$basepair)])
pairs_hippo <- with(hippo_bim, paste(chr, basepair, sep = '-'))
pairs_ldref <- with(ldref_bim, paste(chr, basepair, sep = '-'))
length(unique(pairs_ldref))
length(unique(pairs_hippo))

m <- match(pairs_ldref, pairs_hippo)
table(is.na(m))
#   FALSE    TRUE
# 1022527      19
m2 <- match(pairs_hippo, pairs_ldref)
table(is.na(m2))
#   FALSE
# 1022527
gr_ldref <- GRanges(seqnames = ldref_bim$chr, IRanges(start = ldref_bim$basepair, width = 1))
gr_hippo <- GRanges(seqnames = hippo_bim$chr, IRanges(start = hippo_bim$basepair, width = 1))
table(countOverlaps(gr_hippo, gr_ldref))
## Originally:
#       1
# 1182028
## Now after dropping duplicated snps:
#       1
# 1181974
## After moving to hg38:
#       1
# 1022527
table(countOverlaps(gr_ldref, gr_hippo))
## Originally:
#    0       1       2       3
# 8321 1181974      24       2
## Now after dropping duplicated snps:
#    0       1
# 8347 1181974
## After moving to hg38:
#  0       1
# 19 1022527


hippo_bim[subjectHits(findOverlaps(gr_ldref[countOverlaps(gr_ldref, gr_hippo) > 1], gr_hippo)), ]
## Originally:
#     chr                                   snp position  basepair        allele1          allele2
#  1:   1             1:1627805:A:<CN0>:1673809        0   1627805          <CN0>                A
#  2:   1             1:1627805:A:<CN2>:1673809        0   1627805          <CN2>                A
#  3:   1                  rs601141:1627805:A:C        0   1627805              A                C
#  4:   1         1:108656210:C:<CN0>:108657767        0 108656210          <CN0>                C
#  5:   1                rs614456:108656210:C:T        0 108656210              T                C
#  6:   1         1:189785833:G:<CN0>:189787567        0 189785833          <CN0>                G
#  7:   1                             rs3931586        0 189785833              G                A
#  8:   1            rs36022107:236026148:C:CAT        0 236026148              C              CAT
#  9:   1               rs7513729:236026148:C:T        0 236026148              T                C
# 10:   2        2:202023575:CATTATTATTATTATT:C        0 202023575              C CATTATTATTATTATT
# 11:   2               rs7582581:202023575:C:T        0 202023575              T                C
# 12:   3                  3:74565355:C:<CN0>:0        0  74565355          <CN0>                C
# 13:   3                             rs4676973        0  74565355              C                A
# 14:   4    4:45884781:C:<INS:ME:ALU>:45885041        0  45884781   <INS:ME:ALU>                C
# 15:   4                rs4345219:45884781:C:T        0  45884781              T                C
# 16:   4                 4:155982149:A:<CN0>:0        0 155982149          <CN0>                A
# 17:   4              rs41333846:155982149:A:G        0 155982149              G                A
# 18:   5                  5:52639700:G:<CN0>:0        0  52639700          <CN0>                G
# 19:   5                            rs10491443        0  52639700              A                G
# 20:   5                       5:138647829:G:T        0 138647829              T                G
# 21:   5              rs12186596:138647829:G:A        0 138647829              A                G
# 22:   6                rs9374661:98025180:A:C        0  98025180              C                A
# 23:   6                rs9374661:98025180:A:G        0  98025180              G                A
# 24:   6  6:116143075:T:<INS:ME:ALU>:116143354        0 116143075   <INS:ME:ALU>                T
# 25:   6               rs6908592:116143075:T:C        0 116143075              C                T
# 26:   8           8:91144724:T:<CN0>:91185759        0  91144724          <CN0>                T
# 27:   8           8:91144724:T:<CN2>:91185759        0  91144724          <CN2>                T
# 28:   8                rs7831184:91144724:T:C        0  91144724              T                C
# 29:  10                 10:30323389:G:<CN0>:0        0  30323389          <CN0>                G
# 30:  10               rs11595639:30323389:G:A        0  30323389              A                G
# 31:  10          10:50282772:T:<CN0>:50283171        0  50282772          <CN0>                T
# 32:  10               rs11597087:50282772:T:G        0  50282772              G                T
# 33:  11                     11:103123999:G:GT        0 103123999             GT                G
# 34:  11                rs648387:103123999:G:T        0 103123999              T                G
# 35:  12                   12:77900397:AAAAC:A        0  77900397              A            AAAAC
# 36:  12                             rs1491050        0  77900397              A                C
# 37:  12                      12:107063617:G:T        0 107063617              T                G
# 38:  12                    12:107063617:GGT:G        0 107063617              G              GGT
# 39:  13                   13:38045165:C:CACAA        0  38045165          CACAA                C
# 40:  13                             rs7324856        0  38045165              A                C
# 41:  15          15:93836167:G:<CN0>:93838577        0  93836167          <CN0>                G
# 42:  15                             rs1483324        0  93836167              G                A
# 43:  17              17:789237:A:<CN2>:806860        0    789237          <CN2>                A
# 44:  17                             rs6598829        0    789237              A                G
# 45:  19            19:9863483:G:<CN0>:9872020        0   9863483          <CN0>                G
# 46:  19                rs10415132:9863483:G:A        0   9863483              A                G
# 47:  19          19:51373190:T:<CN2>:51384645        0  51373190          <CN2>                T
# 48:  19                rs1997563:51373190:T:C        0  51373190              C                T
# 49:  20                  20:7951997:A:<CN0>:0        0   7951997          <CN0>                A
# 50:  20                 rs6118035:7951997:A:C        0   7951997              C                A
# 51:  21 21:41740153:A:<INS:ME:LINE1>:41746105        0  41740153 <INS:ME:LINE1>                A
# 52:  21                             rs2837584        0  41740153              A                G
# 53:  22                             rs6007805        0  48529399              A                G
# 54:  22                rs6007805:48529399:G:T        0  48529399              T                G
#     chr                                   snp position  basepair        allele1          allele2

## Now after dropping duplicated snps:
# Empty data.table (0 rows) of 6 cols: chr,snp,position,basepair,allele1,allele2
    
## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.5.1 Patched (2018-10-29 r75535)
#  os       Red Hat Enterprise Linux Server release 6.9 (Santiago)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2019-01-28
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  assertthat             0.2.0     2017-04-11 [2] CRAN (R 3.5.0)
#  bindr                  0.1.1     2018-03-13 [1] CRAN (R 3.5.0)
#  bindrcpp               0.2.2     2018-03-29 [1] CRAN (R 3.5.0)
#  Biobase              * 2.42.0    2018-10-30 [2] Bioconductor
#  BiocGenerics         * 0.28.0    2018-10-30 [1] Bioconductor
#  BiocParallel         * 1.16.5    2019-01-04 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  cli                    1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.4-0     2019-01-13 [2] CRAN (R 3.5.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  data.table           * 1.12.0    2019-01-13 [1] CRAN (R 3.5.1)
#  DelayedArray         * 0.8.0     2018-10-30 [2] Bioconductor
#  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr                  0.7.8     2018-11-10 [1] CRAN (R 3.5.1)
#  GenomeInfoDb         * 1.18.1    2018-11-12 [1] Bioconductor
#  GenomeInfoDbData       1.2.0     2018-11-02 [2] Bioconductor
#  GenomicRanges        * 1.34.0    2018-10-30 [1] Bioconductor
#  ggplot2                3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  glue                   1.3.0     2018-07-17 [1] CRAN (R 3.5.1)
#  gtable                 0.2.0     2016-02-26 [2] CRAN (R 3.5.0)
#  htmltools              0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv                 1.4.5.1   2018-12-18 [2] CRAN (R 3.5.1)
#  IRanges              * 2.16.0    2018-10-30 [1] Bioconductor
#  later                  0.7.5     2018-09-18 [2] CRAN (R 3.5.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval               0.2.1     2017-10-29 [2] CRAN (R 3.5.0)
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  Matrix                 1.2-15    2018-11-01 [3] CRAN (R 3.5.1)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.0)
#  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr                  0.2.5     2018-05-29 [2] CRAN (R 3.5.0)
#  R6                     2.3.0     2018-10-04 [2] CRAN (R 3.5.1)
#  Rcpp                   1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl                  1.95-4.11 2018-07-15 [2] CRAN (R 3.5.1)
#  rlang                  0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  servr                  0.11      2018-10-23 [1] CRAN (R 3.5.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  SummarizedExperiment * 1.12.0    2018-10-30 [1] Bioconductor
#  tibble                 2.0.1     2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun                   0.4       2018-10-23 [1] CRAN (R 3.5.1)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
#
system('plink --version')
# PLINK v1.90b6.6 64-bit (10 Oct 2018)
