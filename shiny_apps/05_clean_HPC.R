library("SingleCellExperiment")
library("rafalib")
library("iSEE")
library("pryr")
library("here")
library("whisker")
library("usethis")
library("withr")
library("sessioninfo")

load(here("rdas", "regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda"), verbose = TRUE)

source(here("shiny_apps", "00_clean_functions.R"))

explore_sce_original(sce.hpc)
# [1] "Dimensions:"
# [1] 33538 10444
# [1] "Number of unique cell names:"
# [1] 10431
# [1] "Repeated cell names:"
#
# AATGCCAGTGCACGCT-1 ACATCCCCAGGACTAG-1 ACATTTCTCTTCGTGC-1 AGCGCTGTCACTACTT-1 AGTGCCGCAGAACTAA-1 CCGAACGCACGAGAAC-1
# 2                  2                  2                  2                  2                  2
# CGGAATTTCAACTGAC-1 CTCCTTTGTAGAATGT-1 GCGATCGTCAGTAGGG-1 GTATTGGAGAGCAGAA-1 GTGTTAGTCGACCCAG-1 GTTGTCCTCTTAAGGC-1
# 2                  2                  2                  2                  2                  2
# TGAGTCATCTGGGCGT-1
# 2
# [1] "Number of unique genes names:"
# [1] 33538

sce_small <- create_small_sce(sce.hpc, cell_var = "cellType.split")
# 751 MB
# 395 MB
colData(sce_small)
rowData(sce_small)


## Test and get the "initial" code
# iSEE(sce_small)

save_sce_small(sce_small, "HPC")

create_app(sce_small, "HPC")

withr::with_dir(here("shiny_apps", "tran2020_HPC"), source("deploy.R"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
