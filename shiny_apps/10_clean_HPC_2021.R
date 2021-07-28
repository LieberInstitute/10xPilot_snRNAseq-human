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

load(here("rdas", "revision", "regionSpecific_HPC-n3_cleaned-combined_SCE_MNT2021.rda"), verbose = TRUE)

source(here("shiny_apps", "00_clean_functions.R"))

explore_sce_original(sce.hpc)
# [1] "Dimensions:"
# [1] 33538 10268
# [1] "Number of unique cell names:"
# [1] 10255
# [1] "Repeated cell names:"
#
# AATGCCAGTGCACGCT-1 ACATCCCCAGGACTAG-1 ACATTTCTCTTCGTGC-1 AGCGCTGTCACTACTT-1
#                  2                  2                  2                  2
# AGTGCCGCAGAACTAA-1 CCGAACGCACGAGAAC-1 CGGAATTTCAACTGAC-1 CTCCTTTGTAGAATGT-1
#                  2                  2                  2                  2
# GCGATCGTCAGTAGGG-1 GTATTGGAGAGCAGAA-1 GTGTTAGTCGACCCAG-1 GTTGTCCTCTTAAGGC-1
#                  2                  2                  2                  2
# TGAGTCATCTGGGCGT-1
#                  2
# [1] "Number of unique genes names:"
# [1] 33538

sce_small <- create_small_sce_2021(sce.hpc)
# * 0.7004298 GB
# * 0.3708309 GB
dim(sce_small)
# [1] 33538 10139
colData(sce_small)
rowData(sce_small)

save_sce_small(sce_small, "HPC", prefix = "tran2021_")
save_cell_colors(cell_colors.hpc, "HPC")

create_app(sce_small, "HPC", cellmarkers = cellmarkers_fig_s7_2021, prefix = "tran2021_")

withr::with_dir(here("shiny_apps", "tran2021_HPC"), source("deploy.R"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

## Same as 07_clean_AMY_2021.R
