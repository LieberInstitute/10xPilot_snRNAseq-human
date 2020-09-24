library("SingleCellExperiment")
library("iSEE")
library("pryr")
library("here")
library("sessioninfo")

load(here("rdas", "regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda"), verbose = TRUE)

dim(sce.nac.all)
# [1] 33538 13241
length(unique(colnames(sce.nac.all)))
# [1] 13224

col_tab <- table(colnames(sce.nac.all))
col_tab[col_tab > 1]
# AGAAATGCACTTCCTG-1 AGGACTTTCTAGCAAC-1 ATGAGTCTCGACATCA-1 CACGGGTGTGAGCAGT-1 CATTCATGTCCTACGG-1
# 2                  2                  2                  2                  2
# CCGTTCAGTCCATAGT-1 CTATAGGTCACTGCTC-1 CTCTCGACAGTTAAAG-1 GACCGTGGTTGTCTAG-1 GACTGATGTCCGATCG-1
# 2                  2                  2                  2                  2
# GATTTCTAGTGGCGAT-1 GCAGTTAAGTAAGACT-1 GGGTTATGTGCGAGTA-1 TAAGCACAGTCATGAA-1 TCGCAGGTCCGTTTCG-1
# 2                  2                  2                  2                  2
# TGATCAGGTAGGATAT-1 TGGGATTGTGATGAAT-1
# 2                  2

length(unique(rownames(sce.nac.all)))
# [1] 33538

sce_nac_small <- sce.nac.all
assays(sce_nac_small) <- assays(sce_nac_small)["logcounts"]
rowData(sce_nac_small)
# sce_nac_small$Sample <- factor(gsub(".*premRNA/|/outs/.*", "", sce_nac_small$Sample))
sce_nac_small$sample <- factor(sce_nac_small$sample)
sce_nac_small$region <- factor(sce_nac_small$region)
sce_nac_small$donor <- factor(sce_nac_small$donor)
sce_nac_small$processDate <- factor(sce_nac_small$processDate)
sce_nac_small$protocol <- factor(sce_nac_small$protocol)
colData(sce_nac_small) <- colData(sce_nac_small)[, !colnames(colData(sce_nac_small)) %in% c("Sample", "prelimCluster", "collapsedCluster", "lessCollapsed", "cellType", "cellType.split", "prelimCluster.split", "cellType.moreSplit", "cellType.final")]
sce_nac_small$cell_type <- factor(sce.nac.all$cellType.final)
colData(sce_nac_small)

## Make the rows more browsable
colnames(sce_nac_small) <- paste0(sce_nac_small$sample, '_', sce_nac_small$Barcode)
sce_nac_small$Barcode <- make.names(sce_nac_small$Barcode, unique = TRUE)
metadata(sce_nac_small) <- list()
sce_nac_small

## It's all "Gene Expression", so we can remove it
rowData(sce_nac_small)$Type <- NULL

## Test and get the "initial" code
# iSEE(sce_nac_small)


object_size(sce.nac.all)
# 2.04 GB
object_size(sce_nac_small)
# 1.29 GB

dir.create(here("shiny_apps", "tran2020_NAc"), showWarnings = FALSE)
saveRDS(sce_nac_small, file = here("shiny_apps", "tran2020_NAc", "sce_nac_small.rds"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# - Session info -------------------------------------------------------------------------------------------------------
#     setting  value
# version  R version 4.0.2 (2020-06-22)
# os       Windows 10 x64
# system   x86_64, mingw32
# ui       RStudio
# language (EN)
# collate  English_United States.1252
# ctype    English_United States.1252
# tz       America/New_York
# date     2020-09-23
#
# - Packages -----------------------------------------------------------------------------------------------------------
#     package              * version  date       lib source
# assertthat             0.2.1    2019-03-21 [1] CRAN (R 4.0.2)
# Biobase              * 2.48.0   2020-04-27 [1] Bioconductor
# BiocGenerics         * 0.34.0   2020-04-27 [1] Bioconductor
# BiocManager            1.30.10  2019-11-16 [1] CRAN (R 4.0.2)
# bitops                 1.0-6    2013-08-17 [1] CRAN (R 4.0.0)
# circlize               0.4.10   2020-06-15 [1] CRAN (R 4.0.2)
# cli                    2.0.2    2020-02-28 [1] CRAN (R 4.0.2)
# clue                   0.3-57   2019-02-25 [1] CRAN (R 4.0.2)
# cluster                2.1.0    2019-06-19 [1] CRAN (R 4.0.2)
# codetools              0.2-16   2018-12-24 [1] CRAN (R 4.0.2)
# colorspace             1.4-1    2019-03-18 [1] CRAN (R 4.0.2)
# colourpicker           1.1.0    2020-09-14 [1] CRAN (R 4.0.2)
# ComplexHeatmap         2.4.3    2020-07-25 [1] Bioconductor
# crayon                 1.3.4    2017-09-16 [1] CRAN (R 4.0.2)
# DelayedArray         * 0.14.1   2020-07-15 [1] Bioconductor
# digest                 0.6.25   2020-02-23 [1] CRAN (R 4.0.2)
# dplyr                  1.0.2    2020-08-18 [1] CRAN (R 4.0.2)
# DT                     0.15     2020-08-05 [1] CRAN (R 4.0.2)
# ellipsis               0.3.1    2020-05-15 [1] CRAN (R 4.0.2)
# fansi                  0.4.1    2020-01-08 [1] CRAN (R 4.0.2)
# fastmap                1.0.1    2019-10-08 [1] CRAN (R 4.0.2)
# generics               0.0.2    2018-11-29 [1] CRAN (R 4.0.2)
# GenomeInfoDb         * 1.24.2   2020-06-15 [1] Bioconductor
# GenomeInfoDbData       1.2.3    2020-06-30 [1] Bioconductor
# GenomicRanges        * 1.40.0   2020-04-27 [1] Bioconductor
# GetoptLong             1.0.2    2020-07-06 [1] CRAN (R 4.0.2)
# ggplot2                3.3.2    2020-06-19 [1] CRAN (R 4.0.2)
# GlobalOptions          0.1.2    2020-06-10 [1] CRAN (R 4.0.2)
# glue                   1.4.2    2020-08-27 [1] CRAN (R 4.0.2)
# gtable                 0.3.0    2019-03-25 [1] CRAN (R 4.0.2)
# htmltools              0.5.0    2020-06-16 [1] CRAN (R 4.0.2)
# htmlwidgets            1.5.1    2019-10-08 [1] CRAN (R 4.0.2)
# httpuv                 1.5.4    2020-06-06 [1] CRAN (R 4.0.2)
# igraph                 1.2.5    2020-03-19 [1] CRAN (R 4.0.2)
# IRanges              * 2.22.2   2020-05-21 [1] Bioconductor
# iSEE                 * 2.0.0    2020-04-27 [1] Bioconductor
# jsonlite               1.7.1    2020-09-07 [1] CRAN (R 4.0.2)
# later                  1.1.0.1  2020-06-05 [1] CRAN (R 4.0.2)
# lattice                0.20-41  2020-04-02 [1] CRAN (R 4.0.2)
# lifecycle              0.2.0    2020-03-06 [1] CRAN (R 4.0.2)
# magrittr               1.5      2014-11-22 [1] CRAN (R 4.0.2)
# Matrix                 1.2-18   2019-11-27 [1] CRAN (R 4.0.2)
# matrixStats          * 0.56.0   2020-03-13 [1] CRAN (R 4.0.2)
# mgcv                   1.8-33   2020-08-27 [1] CRAN (R 4.0.2)
# mime                   0.9      2020-02-04 [1] CRAN (R 4.0.0)
# miniUI                 0.1.1.1  2018-05-18 [1] CRAN (R 4.0.2)
# munsell                0.5.0    2018-06-12 [1] CRAN (R 4.0.2)
# nlme                   3.1-149  2020-08-23 [1] CRAN (R 4.0.2)
# pillar                 1.4.6    2020-07-10 [1] CRAN (R 4.0.2)
# pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 4.0.2)
# png                    0.1-7    2013-12-03 [1] CRAN (R 4.0.0)
# promises               1.1.1    2020-06-09 [1] CRAN (R 4.0.2)
# pryr                 * 0.1.4    2018-02-18 [1] CRAN (R 4.0.2)
# purrr                  0.3.4    2020-04-17 [1] CRAN (R 4.0.2)
# R6                     2.4.1    2019-11-12 [1] CRAN (R 4.0.2)
# RColorBrewer           1.1-2    2014-12-07 [1] CRAN (R 4.0.0)
# Rcpp                   1.0.5    2020-07-06 [1] CRAN (R 4.0.2)
# RCurl                  1.98-1.2 2020-04-18 [1] CRAN (R 4.0.0)
# rintrojs               0.2.2    2019-05-29 [1] CRAN (R 4.0.2)
# rjson                  0.2.20   2018-06-08 [1] CRAN (R 4.0.0)
# rlang                  0.4.7    2020-07-09 [1] CRAN (R 4.0.2)
# rstudioapi             0.11     2020-02-07 [1] CRAN (R 4.0.2)
# S4Vectors            * 0.26.1   2020-05-16 [1] Bioconductor
# scales                 1.1.1    2020-05-11 [1] CRAN (R 4.0.2)
# sessioninfo          * 1.1.1    2018-11-05 [1] CRAN (R 4.0.2)
# shape                  1.4.5    2020-09-13 [1] CRAN (R 4.0.2)
# shiny                  1.5.0    2020-06-23 [1] CRAN (R 4.0.2)
# shinyAce               0.4.1    2019-09-24 [1] CRAN (R 4.0.2)
# shinydashboard         0.7.1    2018-10-17 [1] CRAN (R 4.0.2)
# shinyjs                2.0.0    2020-09-09 [1] CRAN (R 4.0.2)
# shinyWidgets           0.5.3    2020-06-01 [1] CRAN (R 4.0.2)
# SingleCellExperiment * 1.10.1   2020-04-28 [1] Bioconductor
# stringi                1.5.3    2020-09-09 [1] CRAN (R 4.0.2)
# stringr                1.4.0    2019-02-10 [1] CRAN (R 4.0.2)
# SummarizedExperiment * 1.18.2   2020-07-09 [1] Bioconductor
# tibble                 3.0.3    2020-07-10 [1] CRAN (R 4.0.2)
# tidyselect             1.1.0    2020-05-11 [1] CRAN (R 4.0.2)
# vctrs                  0.3.4    2020-08-29 [1] CRAN (R 4.0.2)
# vipor                  0.4.5    2017-03-22 [1] CRAN (R 4.0.2)
# viridisLite            0.3.0    2018-02-01 [1] CRAN (R 4.0.2)
# withr                  2.3.0    2020-09-22 [1] CRAN (R 4.0.2)
# xtable                 1.8-4    2019-04-21 [1] CRAN (R 4.0.2)
# XVector                0.28.0   2020-04-27 [1] Bioconductor
# zlibbioc               1.34.0   2020-04-27 [1] Bioconductor
#
# [1] D:/R/R-4.0.2/library
