################################################################################
### LIBD pilot 10x snRNA-seq: DLPFC samples
###   **Region-specific analyses**
###   - R-batch job for detxn of optimal PC space with 'sce.dlpfc' object
###         -> see '10x-pilot_region-specific_DLPFC_step02_clust-annot_MNTJan2020.R'
###            for setup of the SCE
### MNT 07Feb2020
################################################################################

library(SingleCellExperiment)
library(scRNAseq)
library(batchelor)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(uwot)
library(DropletUtils)
library(jaffelab)
library(Rtsne)

# ===


load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_DLPFC-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=TRUE)
    # sce.dlpfc, chosen.hvgs.dlpfc

## PCA already done (interactively) - took top 100 PCs

## getClusteredPCs() to identify working PC space
pc.choice.dlpfc <- getClusteredPCs(reducedDim(sce.dlpfc))

# How many PCs should use in this space?
metadata(pc.choice.dlpfc)$chosen


    ## Plot n Clusters vs. d PCs
    pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/regionSpecific_DLPFC-n2_getClusteredPCs-results-w100pcs_MNTbFeb2020.pdf")
    plot(pc.choice.dlpfc$n.pcs, pc.choice.dlpfc$n.clusters,
         main=paste0("Combined DLPFC (n=2, Br5161-5212) samples (d PCs choice = ", metadata(pc.choice.dlpfc)$chosen, ")"))
    abline(v=metadata(pc.choice.dlpfc)$chosen, col="red", lty="dashed", lwd=0.8)  
    dev.off()


# Save
save(sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_DLPFC-n2_cleaned-combined_SCE_MNTFeb2020.rda")

rm(list=ls())
sessionInfo()
