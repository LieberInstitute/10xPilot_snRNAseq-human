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

load(here("rdas", "revision", "regionSpecific_sACC-n5_cleaned-combined_SCE_MNT2021.rda"), verbose = TRUE)

source(here("shiny_apps", "00_clean_functions.R"))

explore_sce_original(sce.sacc)
# [1] "Dimensions:"
# [1] 33538 15669
# [1] "Number of unique cell names:"
# [1] 15649
# [1] "Repeated cell names:"
#
# ACGATCATCTGGTGGC-1 ATCCCTGCAAGCAATA-1 CAACGGCGTCCCTCAT-1 CACAGATAGAGCCCAA-1
#                  2                  2                  2                  2
# CCACACTGTGGACTGA-1 CCCTGATAGTAGTCAA-1 CCTACGTCACCACATA-1 CTACGGGAGGTAGGCT-1
#                  2                  2                  2                  2
# CTCATGCCACGATAGG-1 GAAGGACCAGCTATTG-1 GACCCTTAGTTTCGGT-1 GCAGGCTCACCCAATA-1
#                  2                  2                  2                  2
# GGGACCTCATGCGTGC-1 GTATTTCAGAGCAAGA-1 TAGACCAAGTGCTCAT-1 TCGCTCAAGCTCAGAG-1
#                  2                  2                  2                  2
# TGATTCTCATAACTCG-1 TGTAGACCAACCGTGC-1 TTCCTCTAGACAAGCC-1 TTTAGTCAGTGCAGGT-1
#                  2                  2                  2                  2
# [1] "Number of unique genes names:"
# [1] 33538

sce_small <- create_small_sce_2021(sce.sacc)
# * 1.59916 GB
# * 0.82885 GB
dim(sce_small)
# [1] 33538 15343
colData(sce_small)
rowData(sce_small)

save_sce_small(sce_small, "sACC", prefix = "tran2021_")
save_cell_colors(cell_colors.sacc, "sACC")

create_app(sce_small, "sACC", cellmarkers = cellmarkers_fig_s7_2021, prefix = "tran2021_")

withr::with_dir(here("shiny_apps", "tran2021_sACC"), source("deploy.R"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

## Same as 07_clean_AMY_2021.R
