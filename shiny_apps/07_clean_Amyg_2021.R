library("SingleCellExperiment")
library("rafalib")
library("iSEE")
library("pryr")
library("here")
library("whisker")
library("usethis")
library("withr")
library("sessioninfo")

load(here("rdas", "revision", "regionSpecific_Amyg-n5_cleaned-combined_SCE_MNT2021.rda"), verbose = TRUE)

source(here("shiny_apps", "00_clean_functions.R"))

explore_sce_original(sce.amy.tsne.optb)
# [1] "Dimensions:"
# [1] 33538  6582
# [1] "Number of unique cell names:"
# [1] 6579
# [1] "Repeated cell names:"
#
# CAGTTCCTCTATTTCG-1 CGATCGGTCGGCATAT-1 GAAGCGATCGGTAGAG-1
# 2                  2                  2
# [1] "Number of unique genes names:"
# [1] 33538

sce_small <- create_small_sce(sce.amy.tsne.optb, cell_var = "cellType.split")
# 561 MB
# 298 MB
colData(sce_small)
rowData(sce_small)


## Test and get the "initial" code
# iSEE(sce_small)

save_sce_small(sce_small, "Amyg")

create_app(sce_small, "Amyg")

withr::with_dir(here("shiny_apps", "tran2020_Amyg"), source("deploy.R"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
