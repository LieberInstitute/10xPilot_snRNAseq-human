################################################################################
### LIBD pilot 10x snRNA-seq: Amygdala samples
###   **Region-specific analyses**
###   - R-batch job for detxn of optimal PC space with 'sce.amy' object
###         -> see '10x-pilot_region-specific_Amyg_step02_clust-annot_MNTJan2020.R'
###            for setup of the SCE
### MNT 07Feb2020                ~29-30Jan2020~ (for iteration with MBNormalization)
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


load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=TRUE)
    # sce.amy, chosen.hvgs.amy

## PCA already done (interactively) - took top 100 PCs

## getClusteredPCs() to identify working PC space
pc.choice.amy <- getClusteredPCs(reducedDim(sce.amy))

# How many PCs should use in this space?
metadata(pc.choice.amy)$chosen


    ## Plot n Clusters vs. d PCs
    pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/regionSpecific_Amyg-n2_getClusteredPCs-results-w100pcs_MNTFeb2020.pdf")
    plot(pc.choice.amy$n.pcs, pc.choice.amy$n.clusters,
         main=paste0("Combined Amygdala (n=2, Br5161-5212) samples (d PCs choice = ", metadata(pc.choice.amy)$chosen, ")"))
    abline(v=metadata(pc.choice.amy)$chosen, col="red", lty="dashed", lwd=0.8)  
    dev.off()


# Save
save(sce.amy, chosen.hvgs.amy, pc.choice.amy,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda")

rm(list=ls())
sessionInfo()
