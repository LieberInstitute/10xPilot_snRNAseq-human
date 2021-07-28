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

load(here("rdas", "revision", "regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda"), verbose = TRUE)

source(here("shiny_apps", "00_clean_functions.R"))

explore_sce_original(sce.dlpfc)
# [1] "Dimensions:"
# [1] 33538 11202
# [1] "Number of unique cell names:"
# [1] 11196
# [1] "Repeated cell names:"
#
# ACTTTGTCAGCTGAAG-1 CATCGTCCAATAGGGC-1 CGGCAGTCATTCACCC-1 CTAGACATCGCGGTAC-1
#                  2                  2                  2                  2
# GACCAATTCGTTAGAC-1 TGAGTCAAGACCATAA-1
#                  2                  2
# [1] "Number of unique genes names:"
# [1] 33538

sce_small <- create_small_sce_2021(sce.dlpfc)
# * 1.0002005 GB
# * 0.5226881 GB
dim(sce_small)
# [1] 33538 11202
colData(sce_small)
rowData(sce_small)

save_sce_small(sce_small, "DLPFC", prefix = "tran2021_")
save_cell_colors(cell_colors, "DLPFC")

create_app(sce_small, "DLPFC", cellmarkers = cellmarkers_fig_s7_2021, prefix = "tran2021_")

withr::with_dir(here("shiny_apps", "tran2021_DLPFC"), source("deploy.R"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

## Same as 07_clean_AMY_2021.R
