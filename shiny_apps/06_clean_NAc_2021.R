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

load(here("rdas", "revision", "regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda"), verbose = TRUE)

source(here("shiny_apps", "00_clean_functions.R"))

explore_sce_original(sce.nac)
# [1] "Dimensions:"
# [1] 33538 20571
# [1] "Number of unique cell names:"
# [1] 20514
# [1] "Repeated cell names:"
#
# AAAGAACAGCGAGTCA-1 AAAGGTAGTTTCGTTT-1 AAAGTCCCAGCCTACG-1 AAATGGATCCCGTAAA-1
#                  2                  2                  2                  2
# AACACACGTGTGACCC-1 AAGTCGTAGTGGCAGT-1 ACCAACATCGGTCATA-1 ACGATGTGTAGTAAGT-1
#                  2                  2                  2                  2
# ACTGCAACAACATCGT-1 AGAAATGCACTTCCTG-1 AGCATCAGTCAAGTTC-1 AGGACTTTCTAGCAAC-1
#                  2                  2                  2                  2
# AGGCCACGTCCAGCAC-1 ATCAGGTTCATTCCTA-1 ATGACCAGTTAAGGAT-1 ATGAGGGCAGCTTTGA-1
#                  2                  2                  2                  2
# ATGAGTCTCGACATCA-1 ATTACCTAGACTCGAG-1 ATTGGGTTCCGTATAG-1 CACGGGTGTGAGCAGT-1
#                  2                  2                  2                  2
# CACTGTCGTCCTTTGC-1 CAGATACGTTGGCCGT-1 CAGGGCTTCAGACCCG-1 CAGTTCCAGTCGTTAC-1
#                  2                  2                  2                  2
# CATTCATGTCCTACGG-1 CCCTCAAAGCAATTCC-1 CCGTTCAGTCCATAGT-1 CCTCTCCGTTCAAGTC-1
#                  2                  2                  2                  2
# CCTGCATAGTGCTAGG-1 CGTTGGGTCTGTGCAA-1 CTATAGGTCACTGCTC-1 CTCACTGTCGCGTGCA-1
#                  2                  2                  2                  2
# CTCTCGACAGTTAAAG-1 GACCGTGGTTGTCTAG-1 GACCTTCGTCTGTAAC-1 GACTGATGTCCGATCG-1
#                  2                  2                  2                  2
# GATTTCTAGTGGCGAT-1 GCAGTTAAGTAAGACT-1 GCGTGCAGTATATGGA-1 GGACGTCCAGGTGACA-1
#                  2                  2                  2                  2
# GGCTTGGCATCGTGGC-1 GGGTTATGTGCGAGTA-1 GGTGATTAGCTGACAG-1 GTTCCGTGTTGGAGAC-1
#                  2                  2                  2                  2
# TAAGCACAGTCATGAA-1 TCACAAGAGGCTTAGG-1 TCCTAATTCGATAACC-1 TCGCAGGTCCGTTTCG-1
#                  2                  2                  2                  2
# TGAACGTGTGTCCTAA-1 TGATCAGGTAGGATAT-1 TGCTCGTTCCGCGAGT-1 TGGAACTGTGGTTCTA-1
#                  2                  2                  2                  2
# TGGGAAGAGAAGGGAT-1 TGGGATTGTGATGAAT-1 TGTGTGACACTGTGAT-1 TTGACCCGTTCCGCTT-1
#                  2                  2                  2                  2
# TTGTGTTCAATCGCGC-1
#                  2
# [1] "Number of unique genes names:"
# [1] 33538

sce_small <- create_small_sce_2021(sce.nac)
# * 2.120966 B
# * 1.093446 B
dim(sce_small)
# [1] 33538 19892
colData(sce_small)
rowData(sce_small)


## Test and get the "initial" code
# iSEE(sce_small)

save_sce_small(sce_small, "NAc", prefix = "tran2021_")
save_cell_colors(cell_colors.nac, "NAc")

create_app(sce_small, "NAc", cellmarkers = cellmarkers_fig_s7_2021, prefix = "tran2021_")

withr::with_dir(here("shiny_apps", "tran2021_NAc"), source("deploy.R"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.1.0 (2021-05-18)
#  os       macOS Big Sur 11.4
#  system   x86_64, darwin17.0
#  ui       RStudio
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2021-07-27
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version    date       lib source
#  assertthat             0.2.1      2019-03-21 [1] CRAN (R 4.1.0)
#  Biobase              * 2.52.0     2021-05-19 [1] Bioconductor
#  BiocGenerics         * 0.38.0     2021-05-19 [1] Bioconductor
#  bitops                 1.0-7      2021-04-24 [1] CRAN (R 4.1.0)
#  bslib                  0.2.5.1    2021-05-18 [1] CRAN (R 4.1.0)
#  cachem                 1.0.5      2021-05-15 [1] CRAN (R 4.1.0)
#  Cairo                  1.5-12.2   2020-07-07 [1] CRAN (R 4.1.0)
#  callr                  3.7.0      2021-04-20 [1] CRAN (R 4.1.0)
#  circlize               0.4.13     2021-06-09 [1] CRAN (R 4.1.0)
#  cli                    3.0.1      2021-07-17 [1] CRAN (R 4.1.0)
#  clue                   0.3-59     2021-04-16 [1] CRAN (R 4.1.0)
#  cluster                2.1.2      2021-04-17 [1] CRAN (R 4.1.0)
#  codetools              0.2-18     2020-11-04 [1] CRAN (R 4.1.0)
#  colorout               1.2-2      2020-11-03 [1] Github (jalvesaq/colorout@726d681)
#  colorspace             2.0-2      2021-06-24 [1] CRAN (R 4.1.0)
#  colourpicker           1.1.0      2020-09-14 [1] CRAN (R 4.1.0)
#  ComplexHeatmap         2.8.0      2021-05-19 [1] Bioconductor
#  crayon                 1.4.1      2021-02-08 [1] CRAN (R 4.1.0)
#  data.table             1.14.0     2021-02-21 [1] CRAN (R 4.1.0)
#  DBI                    1.1.1      2021-01-15 [1] CRAN (R 4.1.0)
#  DelayedArray           0.18.0     2021-05-19 [1] Bioconductor
#  desc                   1.3.0      2021-03-05 [1] CRAN (R 4.1.0)
#  devtools             * 2.4.2      2021-06-07 [1] CRAN (R 4.1.0)
#  digest                 0.6.27     2020-10-24 [1] CRAN (R 4.1.0)
#  doParallel             1.0.16     2020-10-16 [1] CRAN (R 4.1.0)
#  dplyr                  1.0.7      2021-06-18 [1] CRAN (R 4.1.0)
#  DT                     0.18       2021-04-14 [1] CRAN (R 4.1.0)
#  ellipsis               0.3.2      2021-04-29 [1] CRAN (R 4.1.0)
#  fansi                  0.5.0      2021-05-25 [1] CRAN (R 4.1.0)
#  fastmap                1.1.0      2021-01-25 [1] CRAN (R 4.1.0)
#  foreach                1.5.1      2020-10-15 [1] CRAN (R 4.1.0)
#  fs                     1.5.0      2020-07-31 [1] CRAN (R 4.1.0)
#  generics               0.1.0      2020-10-31 [1] CRAN (R 4.1.0)
#  GenomeInfoDb         * 1.28.1     2021-07-01 [1] Bioconductor
#  GenomeInfoDbData       1.2.6      2021-05-24 [1] Bioconductor
#  GenomicRanges        * 1.44.0     2021-05-19 [1] Bioconductor
#  GetoptLong             1.0.5      2020-12-15 [1] CRAN (R 4.1.0)
#  ggplot2                3.3.5      2021-06-25 [1] CRAN (R 4.1.0)
#  ggrepel                0.9.1      2021-01-15 [1] CRAN (R 4.1.0)
#  GlobalOptions          0.1.2      2020-06-10 [1] CRAN (R 4.1.0)
#  glue                   1.4.2      2020-08-27 [1] CRAN (R 4.1.0)
#  gtable                 0.3.0      2019-03-25 [1] CRAN (R 4.1.0)
#  here                 * 1.0.1      2020-12-13 [1] CRAN (R 4.1.0)
#  hms                    1.1.0      2021-05-17 [1] CRAN (R 4.1.0)
#  htmltools              0.5.1.1    2021-01-22 [1] CRAN (R 4.1.0)
#  htmlwidgets            1.5.3      2020-12-10 [1] CRAN (R 4.1.0)
#  httpuv                 1.6.1      2021-05-07 [1] CRAN (R 4.1.0)
#  igraph                 1.2.6      2020-10-06 [1] CRAN (R 4.1.0)
#  IRanges              * 2.26.0     2021-05-19 [1] Bioconductor
#  iSEE                 * 2.4.0      2021-05-19 [1] Bioconductor
#  iterators              1.0.13     2020-10-15 [1] CRAN (R 4.1.0)
#  jquerylib              0.1.4      2021-04-26 [1] CRAN (R 4.1.0)
#  jsonlite               1.7.2      2020-12-09 [1] CRAN (R 4.1.0)
#  later                  1.2.0      2021-04-23 [1] CRAN (R 4.1.0)
#  lattice                0.20-44    2021-05-02 [1] CRAN (R 4.1.0)
#  lifecycle              1.0.0      2021-02-15 [1] CRAN (R 4.1.0)
#  lobstr               * 1.1.1      2019-07-02 [1] CRAN (R 4.1.0)
#  lubridate              1.7.10     2021-02-26 [1] CRAN (R 4.1.0)
#  magrittr               2.0.1      2020-11-17 [1] CRAN (R 4.1.0)
#  Matrix                 1.3-4      2021-06-01 [1] CRAN (R 4.1.0)
#  MatrixGenerics       * 1.4.0      2021-05-19 [1] Bioconductor
#  matrixStats          * 0.60.0     2021-07-26 [1] CRAN (R 4.1.0)
#  memoise                2.0.0      2021-01-26 [1] CRAN (R 4.1.0)
#  mgcv                   1.8-36     2021-06-01 [1] CRAN (R 4.1.0)
#  mime                   0.11       2021-06-23 [1] CRAN (R 4.1.0)
#  miniUI                 0.1.1.1    2018-05-18 [1] CRAN (R 4.1.0)
#  munsell                0.5.0      2018-06-12 [1] CRAN (R 4.1.0)
#  nlme                   3.1-152    2021-02-04 [1] CRAN (R 4.1.0)
#  pillar                 1.6.1      2021-05-16 [1] CRAN (R 4.1.0)
#  pkgbuild               1.2.0      2020-12-15 [1] CRAN (R 4.1.0)
#  pkgconfig              2.0.3      2021-03-31 [1] Github (gaborcsardi/pkgconfig@b81ae03)
#  pkgload                1.2.1      2021-04-06 [1] CRAN (R 4.1.0)
#  png                    0.1-7      2013-12-03 [1] CRAN (R 4.1.0)
#  prettyunits            1.1.1      2020-01-24 [1] CRAN (R 4.1.0)
#  processx               3.5.2      2021-04-30 [1] CRAN (R 4.1.0)
#  promises               1.2.0.1    2021-02-11 [1] CRAN (R 4.1.0)
#  prompt                 1.0.1      2021-03-12 [1] CRAN (R 4.1.0)
#  pryr                   0.1.5      2021-07-26 [1] CRAN (R 4.1.0)
#  ps                     1.6.0      2021-02-28 [1] CRAN (R 4.1.0)
#  purrr                  0.3.4      2020-04-17 [1] CRAN (R 4.1.0)
#  R6                     2.5.0      2020-10-28 [1] CRAN (R 4.1.0)
#  rafalib              * 1.0.0      2015-08-09 [1] CRAN (R 4.1.0)
#  RColorBrewer           1.1-2      2014-12-07 [1] CRAN (R 4.1.0)
#  Rcpp                   1.0.7      2021-07-07 [1] CRAN (R 4.1.0)
#  RCurl                  1.98-1.3   2021-03-16 [1] CRAN (R 4.1.0)
#  remotes                2.4.0      2021-06-02 [1] CRAN (R 4.1.0)
#  rintrojs               0.3.0      2021-06-06 [1] CRAN (R 4.1.0)
#  rjson                  0.2.20     2018-06-08 [1] CRAN (R 4.1.0)
#  rlang                  0.4.11     2021-04-30 [1] CRAN (R 4.1.0)
#  rprojroot              2.0.2      2020-11-15 [1] CRAN (R 4.1.0)
#  rsconnect            * 0.8.18     2021-05-24 [1] CRAN (R 4.1.0)
#  rsthemes               0.2.1.9000 2021-02-12 [1] Github (gadenbuie/rsthemes@521572b)
#  rstudioapi             0.13       2020-11-12 [1] CRAN (R 4.1.0)
#  S4Vectors            * 0.30.0     2021-05-19 [1] Bioconductor
#  sass                   0.4.0      2021-05-12 [1] CRAN (R 4.1.0)
#  scales                 1.1.1      2020-05-11 [1] CRAN (R 4.1.0)
#  sessioninfo          * 1.1.1      2018-11-05 [1] CRAN (R 4.1.0)
#  shape                  1.4.6      2021-05-19 [1] CRAN (R 4.1.0)
#  shiny                  1.6.0      2021-01-25 [1] CRAN (R 4.1.0)
#  shinyAce               0.4.1      2019-09-24 [1] CRAN (R 4.1.0)
#  shinydashboard         0.7.1      2018-10-17 [1] CRAN (R 4.1.0)
#  shinyjs                2.0.0      2020-09-09 [1] CRAN (R 4.1.0)
#  shinyWidgets           0.6.0      2021-03-15 [1] CRAN (R 4.1.0)
#  SingleCellExperiment * 1.14.1     2021-05-21 [1] Bioconductor
#  stringi                1.7.3      2021-07-16 [1] CRAN (R 4.1.0)
#  stringr                1.4.0      2019-02-10 [1] CRAN (R 4.1.0)
#  SummarizedExperiment * 1.22.0     2021-05-19 [1] Bioconductor
#  suncalc                0.5.0      2019-04-03 [1] CRAN (R 4.1.0)
#  testthat             * 3.0.4      2021-07-01 [1] CRAN (R 4.1.0)
#  tibble                 3.1.3      2021-07-23 [1] CRAN (R 4.1.0)
#  tidyselect             1.1.1      2021-04-30 [1] CRAN (R 4.1.0)
#  usethis              * 2.0.1      2021-02-10 [1] CRAN (R 4.1.0)
#  utf8                   1.2.2      2021-07-24 [1] CRAN (R 4.1.0)
#  vctrs                  0.3.8      2021-04-29 [1] CRAN (R 4.1.0)
#  vipor                  0.4.5      2017-03-22 [1] CRAN (R 4.1.0)
#  viridisLite            0.4.0      2021-04-13 [1] CRAN (R 4.1.0)
#  whisker              * 0.4        2019-08-28 [1] CRAN (R 4.1.0)
#  withr                * 2.4.2      2021-04-18 [1] CRAN (R 4.1.0)
#  xtable                 1.8-4      2019-04-21 [1] CRAN (R 4.1.0)
#  XVector                0.32.0     2021-05-19 [1] Bioconductor
#  zlibbioc               1.38.0     2021-05-19 [1] Bioconductor
#
# [1] /Library/Frameworks/R.framework/Versions/4.1/Resources/library
