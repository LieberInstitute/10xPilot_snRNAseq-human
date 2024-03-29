# 10xPilot_snRNAseq-human

## How to cite

* https://doi.org/10.1101/2020.10.07.329839

## Study design

TODO


## Explore the data interactively

We have provided 5 interactive websites that allow you to explore the data at single nucleus resolution for each of the brain regions. These interactive websites are powered by [`iSEE`](https://bioconductor.org/packages/iSEE) that allows you to add, hide, customize panels for visualizing the data. Please check the [`iSEE`](https://bioconductor.org/packages/iSEE) documentation for instructions on how to customize the panels. In particular, you might be interested in visualizing some of the marker genes from the lists provided below for the _region-specific analyses_.

* https://libd.shinyapps.io/tran2020_Amyg/
* https://libd.shinyapps.io/tran2020_DLPFC/
* https://libd.shinyapps.io/tran2020_HPC/
* https://libd.shinyapps.io/tran2020_NAc/
* https://libd.shinyapps.io/tran2020_sACC/

If you want to make these websites on your own computer, check the [`shiny_apps`](shiny_apps/) directory.


## Work with the data

The corresponding `SingleCellExperiment` objects (with `reducedDims`, annotations, etc.) for each of the five regions are publicly hosted at:

* https://libd-snrnaseq-pilot.s3.us-east-2.amazonaws.com/SCE_AMY_tran-etal.rda
* https://libd-snrnaseq-pilot.s3.us-east-2.amazonaws.com/SCE_DLPFC_tran-etal.rda
* https://libd-snrnaseq-pilot.s3.us-east-2.amazonaws.com/SCE_HPC_tran-etal.rda
* https://libd-snrnaseq-pilot.s3.us-east-2.amazonaws.com/SCE_NAc_tran-etal.rda
* https://libd-snrnaseq-pilot.s3.us-east-2.amazonaws.com/SCE_sACC_tran-etal.rda


## Marker lists and expression plots, top 40

### Region-specific analyses

Here, nuclei are clustered and annotated within each brain region, separately, with markers defined at that level.
-   **AMY**:

    https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/master/tables/top40genesLists-REFINED_Amyg-n2_cellType.split_Nov2020.csv

    https://github.com/LieberInstitute/10xPilot_snRNAseq-human/tree/master/pdfs/exploration/Amyg

-   **DLPFC**:

    https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/master/tables/top40genesLists_DLPFC-n2_cellType.split_SN-LEVEL-tests_May2020.csv

    https://github.com/LieberInstitute/10xPilot_snRNAseq-human/tree/master/pdfs/exploration/DLPFC

-   **HPC**:
    
    https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/master/tables/top40genesLists-REFINED_HPC-n3_cellType.split_Nov2020.csv
    
    https://github.com/LieberInstitute/10xPilot_snRNAseq-human/tree/master/pdfs/exploration/HPC

-   **NAc**:

    https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/master/tables/top40genesLists-REFINED_NAc-n5_cellType.final_Nov2020.csv

    https://github.com/LieberInstitute/10xPilot_snRNAseq-human/tree/master/pdfs/exploration/NAc-n5-markers

-   **sACC**:

    https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/master/tables/top40genesLists-REFINED_sACC-n2_cellType.split_Nov2020.csv

    https://github.com/LieberInstitute/10xPilot_snRNAseq-human/tree/master/pdfs/exploration/sACC


### Pan-brain-level analysis

Here, nuclei are clustered across all brain regions, together, and then annotated, with markers defined at _this_ level.

   https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/master/tables/top40genesLists_panBrain-n12_cellType_SN-LEVEL-tests_May2020.csv

   https://github.com/LieberInstitute/10xPilot_snRNAseq-human/tree/master/pdfs/exploration/panBrainMarkers


## LIBD internal

JHPCE location: `/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL`
