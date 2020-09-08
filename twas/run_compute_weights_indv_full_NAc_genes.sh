#!/bin/bash

## These mkdir steps + ln -s + "mkdir -p logs/NAc_genes" were typically done
## outside the loop at
## https://github.com/LieberInstitute/twas/blob/master/bsp2/compute_weights_indv.sh

## Required order for running this code:

# For the logs
mkdir -p logs/NAc_genes

## For output files
mkdir -p NAc_gene/tmp_files
mkdir -p NAc_gene/out_files

# For GEMMA
ln -s /dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/twas/NAc_gene/ NAc_gene/output

## For running the main script
qsub compute_weights_indv_full_NAc_genes.sh
