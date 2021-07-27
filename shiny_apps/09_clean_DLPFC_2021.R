library("SingleCellExperiment")
library("rafalib")
library("iSEE")
library("lobstr")
library("here")
library("whisker")
library("usethis")
library("withr")
library("sessioninfo")

load(here("rdas", "revision", "regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda"), verbose = TRUE)

source(here("shiny_apps", "00_clean_functions.R"))

explore_sce_original(sce.dlpfc.st)
# [1] "Dimensions:"
# [1] 33538  5399
# [1] "Number of unique cell names:"
# [1] 5398
# [1] "Repeated cell names:"
# CGGCAGTCATTCACCC-1
# 2
# [1] "Number of unique genes names:"
# [1] 33538

sce_small <- create_small_sce_2021(sce.dlpfc.st, cell_var = "cellType.split")
# 405 MB
# 220 MB
colData(sce_small)
rowData(sce_small)


## Test and get the "initial" code
# iSEE(sce_small)

save_sce_small(sce_small, "DLPFC", prefix = "tran2021_")

create_app(sce_small, "DLPFC", prefix = "tran2021_")

withr::with_dir(here("shiny_apps", "tran2021_DLPFC"), source("deploy.R"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
