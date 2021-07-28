library("SingleCellExperiment")
library("rafalib")
library("iSEE")
library("lobstr")
library("here")
library("whisker")
library("usethis")
library("withr")
library("rsconnect")
library("sessioninfo")

load(here("rdas", "revision", "regionSpecific_Amyg-n5_cleaned-combined_SCE_MNT2021.rda"), verbose = TRUE)

source(here("shiny_apps", "00_clean_functions.R"))

explore_sce_original(sce.amy)
# [1] "Dimensions:"
# [1] 33538 15177
# [1] "Number of unique cell names:"
# [1] 15138
# [1] "Repeated cell names:"
#
# AAATGGATCCCGTAAA-1 AACCAACGTATCGCGC-1 AACTTCTTCGCTCTAC-1 AAGTACCTCTTCTTCC-1
#                  2                  2                  2                  2
# ACTATTCCAGCCTATA-1 AGACCCGGTGGCCCAT-1 AGCTTCCGTTCAATCG-1 ATCGCCTGTGGTACAG-1
#                  2                  2                  2                  2
# ATGAAAGAGACACACG-1 ATTGTTCTCATGACAC-1 CAACAGTCAGTCCCGA-1 CAACCAACATCAGCGC-1
#                  2                  2                  2                  2
# CAACCTCGTTCTCACC-1 CAACCTCTCTGATTCT-1 CACAGGCCATTAAAGG-1 CAGTTCCTCTATTTCG-1
#                  2                  2                  2                  2
# CCCTCTCGTATCCCAA-1 CGATCGGTCGGCATAT-1 CGCCATTGTGCTGCAC-1 CGGACACCACGTACTA-1
#                  2                  2                  2                  2
# CGTAATGGTCACTCTC-1 CGTGAATAGACCGTTT-1 GAAGCGATCGGTAGAG-1 GAGGCAAGTCTCAGGC-1
#                  2                  2                  2                  2
# GATGAGGTCCCATACC-1 GCACGTGCAAGCTGTT-1 GCTACAACATGTAACC-1 GGACGTCTCTTCGGTC-1
#                  2                  2                  2                  2
# GGTAATCAGGATATAC-1 GTCACTCTCTCCCATG-1 TCAAGACCAAATGATG-1 TCACATTTCTCACTCG-1
#                  2                  2                  2                  2
# TCATCCGTCGATTGAC-1 TCCTAATTCGATAACC-1 TCCTTCTAGGTATTGA-1 TCGAAGTTCTTGCAAG-1
#                  2                  2                  2                  2
# TGGTTAGCAGTCAACT-1 TTCTGTACATCTATCT-1 TTCTTCCAGTTGCTCA-1
#                  2                  2                  2
# [1] "Number of unique genes names:"
# [1] 33538

sce_small <- create_small_sce_2021(sce.amy)
# * 1.1673461 GB
# * 0.6082883 GB
dim(sce_small)
# [1] 33538 14039
colData(sce_small)
rowData(sce_small)

save_sce_small(sce_small, "AMY", prefix = "tran2021_")
save_cell_colors(cell_colors.amy, "AMY")

create_app(sce_small, "AMY", cellmarkers = cellmarkers_fig_s7_2021, prefix = "tran2021_")

withr::with_dir(here("shiny_apps", "tran2021_AMY"), source("deploy.R"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

## For the portion done on the cluster. For deploying, check 06_clean_NAc_2021.R
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.1.0 Patched (2021-05-18 r80330)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2021-07-28
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version  date       lib source
#  assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
#  Biobase              * 2.52.0   2021-05-19 [2] Bioconductor
#  BiocGenerics         * 0.38.0   2021-05-19 [2] Bioconductor
#  bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
#  bslib                  0.2.5.1  2021-05-18 [2] CRAN (R 4.1.0)
#  Cairo                  1.5-12.2 2020-07-07 [1] CRAN (R 4.1.0)
#  circlize               0.4.13   2021-06-09 [1] CRAN (R 4.1.0)
#  cli                    3.0.1    2021-07-17 [2] CRAN (R 4.1.0)
#  clue                   0.3-59   2021-04-16 [2] CRAN (R 4.1.0)
#  cluster                2.1.2    2021-04-17 [3] CRAN (R 4.1.0)
#  codetools              0.2-18   2020-11-04 [2] CRAN (R 4.1.0)
#  colorout               1.2-2    2021-05-25 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace             2.0-2    2021-06-24 [2] CRAN (R 4.1.0)
#  colourpicker           1.1.0    2020-09-14 [2] CRAN (R 4.1.0)
#  ComplexHeatmap         2.8.0    2021-05-19 [1] Bioconductor
#  crayon                 1.4.1    2021-02-08 [2] CRAN (R 4.1.0)
#  DBI                    1.1.1    2021-01-15 [2] CRAN (R 4.1.0)
#  DelayedArray           0.18.0   2021-05-19 [2] Bioconductor
#  digest                 0.6.27   2020-10-24 [2] CRAN (R 4.1.0)
#  doParallel             1.0.16   2020-10-16 [2] CRAN (R 4.1.0)
#  dplyr                  1.0.7    2021-06-18 [2] CRAN (R 4.1.0)
#  DT                     0.18     2021-04-14 [2] CRAN (R 4.1.0)
#  ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
#  fansi                  0.5.0    2021-05-25 [2] CRAN (R 4.1.0)
#  fastmap                1.1.0    2021-01-25 [2] CRAN (R 4.1.0)
#  foreach                1.5.1    2020-10-15 [2] CRAN (R 4.1.0)
#  fs                     1.5.0    2020-07-31 [2] CRAN (R 4.1.0)
#  generics               0.1.0    2020-10-31 [2] CRAN (R 4.1.0)
#  GenomeInfoDb         * 1.28.1   2021-07-01 [2] Bioconductor
#  GenomeInfoDbData       1.2.6    2021-05-11 [2] Bioconductor
#  GenomicRanges        * 1.44.0   2021-05-19 [2] Bioconductor
#  GetoptLong             1.0.5    2020-12-15 [1] CRAN (R 4.1.0)
#  ggplot2                3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
#  ggrepel                0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
#  GlobalOptions          0.1.2    2020-06-10 [1] CRAN (R 4.1.0)
#  glue                   1.4.2    2020-08-27 [2] CRAN (R 4.1.0)
#  gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
#  here                 * 1.0.1    2020-12-13 [1] CRAN (R 4.1.0)
#  htmltools              0.5.1.1  2021-01-22 [2] CRAN (R 4.1.0)
#  htmlwidgets            1.5.3    2020-12-10 [2] CRAN (R 4.1.0)
#  httpuv                 1.6.1    2021-05-07 [2] CRAN (R 4.1.0)
#  igraph                 1.2.6    2020-10-06 [2] CRAN (R 4.1.0)
#  IRanges              * 2.26.0   2021-05-19 [2] Bioconductor
#  iSEE                 * 2.4.0    2021-05-19 [1] Bioconductor
#  iterators              1.0.13   2020-10-15 [2] CRAN (R 4.1.0)
#  jquerylib              0.1.4    2021-04-26 [2] CRAN (R 4.1.0)
#  jsonlite               1.7.2    2020-12-09 [2] CRAN (R 4.1.0)
#  later                  1.2.0    2021-04-23 [2] CRAN (R 4.1.0)
#  lattice                0.20-44  2021-05-02 [3] CRAN (R 4.1.0)
#  lifecycle              1.0.0    2021-02-15 [2] CRAN (R 4.1.0)
#  lobstr               * 1.1.1    2019-07-02 [2] CRAN (R 4.1.0)
#  magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
#  Matrix                 1.3-4    2021-06-01 [3] CRAN (R 4.1.0)
#  MatrixGenerics       * 1.4.0    2021-05-19 [2] Bioconductor
#  matrixStats          * 0.60.0   2021-07-26 [2] CRAN (R 4.1.0)
#  mgcv                   1.8-36   2021-06-01 [3] CRAN (R 4.1.0)
#  mime                   0.11     2021-06-23 [2] CRAN (R 4.1.0)
#  miniUI                 0.1.1.1  2018-05-18 [2] CRAN (R 4.1.0)
#  munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
#  nlme                   3.1-152  2021-02-04 [3] CRAN (R 4.1.0)
#  pillar                 1.6.1    2021-05-16 [2] CRAN (R 4.1.0)
#  pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
#  png                    0.1-7    2013-12-03 [2] CRAN (R 4.1.0)
#  promises               1.2.0.1  2021-02-11 [2] CRAN (R 4.1.0)
#  purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
#  R6                     2.5.0    2020-10-28 [2] CRAN (R 4.1.0)
#  rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.1.0)
#  RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
#  Rcpp                   1.0.7    2021-07-07 [2] CRAN (R 4.1.0)
#  RCurl                  1.98-1.3 2021-03-16 [2] CRAN (R 4.1.0)
#  rintrojs               0.3.0    2021-06-06 [1] CRAN (R 4.1.0)
#  rjson                  0.2.20   2018-06-08 [2] CRAN (R 4.1.0)
#  rlang                  0.4.11   2021-04-30 [2] CRAN (R 4.1.0)
#  rmote                  0.3.4    2021-05-25 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
#  rsconnect            * 0.8.18   2021-05-24 [2] CRAN (R 4.1.0)
#  S4Vectors            * 0.30.0   2021-05-19 [2] Bioconductor
#  sass                   0.4.0    2021-05-12 [2] CRAN (R 4.1.0)
#  scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
#  servr                  0.22     2021-04-14 [1] CRAN (R 4.1.0)
#  sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 4.1.0)
#  shape                  1.4.6    2021-05-19 [1] CRAN (R 4.1.0)
#  shiny                  1.6.0    2021-01-25 [2] CRAN (R 4.1.0)
#  shinyAce               0.4.1    2019-09-24 [1] CRAN (R 4.1.0)
#  shinydashboard         0.7.1    2018-10-17 [1] CRAN (R 4.1.0)
#  shinyjs                2.0.0    2020-09-09 [2] CRAN (R 4.1.0)
#  shinyWidgets           0.6.0    2021-03-15 [1] CRAN (R 4.1.0)
#  SingleCellExperiment * 1.14.1   2021-05-21 [2] Bioconductor
#  SummarizedExperiment * 1.22.0   2021-05-19 [2] Bioconductor
#  tibble                 3.1.3    2021-07-23 [2] CRAN (R 4.1.0)
#  tidyselect             1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
#  usethis              * 2.0.1    2021-02-10 [2] CRAN (R 4.1.0)
#  utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
#  vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
#  vipor                  0.4.5    2017-03-22 [1] CRAN (R 4.1.0)
#  viridisLite            0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
#  whisker              * 0.4      2019-08-28 [2] CRAN (R 4.1.0)
#  withr                * 2.4.2    2021-04-18 [2] CRAN (R 4.1.0)
#  xfun                   0.24     2021-06-15 [2] CRAN (R 4.1.0)
#  xtable                 1.8-4    2019-04-21 [2] CRAN (R 4.1.0)
#  XVector                0.32.0   2021-05-19 [2] Bioconductor
#  zlibbioc               1.38.0   2021-05-19 [2] Bioconductor
#
# [1] /users/lcollado/R/4.1
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/library
