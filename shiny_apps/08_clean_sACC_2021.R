library("SingleCellExperiment")
library("rafalib")
library("iSEE")
library("lobstr")
library("here")
library("whisker")
library("usethis")
library("withr")
library("sessioninfo")

load(here("rdas", "revision", "regionSpecific_sACC-n5_cleaned-combined_SCE_MNT2021.rda"), verbose = TRUE)

source(here("shiny_apps", "00_clean_functions.R"))

explore_sce_original(sce.sacc)
# [1] "Dimensions:"
# [1] 33538  7047
# [1] "Number of unique cell names:"
# [1] 7043
# [1] "Repeated cell names:"
#
# CACAGATAGAGCCCAA-1 CCTACGTCACCACATA-1 GTATTTCAGAGCAAGA-1 TGTAGACCAACCGTGC-1
# 2                  2                  2                  2
# [1] "Number of unique genes names:"
# [1] 33538

sce_small <- create_small_sce_2021(sce.sacc, cell_var = "cellType")
# 675 MB
# 353 MB
colData(sce_small)
rowData(sce_small)


## Test and get the "initial" code
# iSEE(sce_small)

save_sce_small(sce_small, "sACC", prefix = "tran2021_")

create_app(sce_small, "sACC", prefix = "tran2021_")

withr::with_dir(here("shiny_apps", "tran2021_sACC"), source("deploy.R"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
