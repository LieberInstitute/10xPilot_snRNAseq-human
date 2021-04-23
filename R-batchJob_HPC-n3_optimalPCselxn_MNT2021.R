################################################################################
### LIBD pilot 10x snRNA-seq: HPC samples
###   **Region-specific analyses**
###   - R-batch job for detxn of optimal PC space with 'sce.hpc' object
###         -> see '10x_HPC-n3_step02_clust-annot_MNT.R'
###            for setup of the SCE
### MNT 23Apr2021
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

# ===


load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_HPC-n3_cleaned-combined_SCE_MNT2021.rda",
     verbose=TRUE)
    # sce.hpc, chosen.hvgs.hpc, ref.sampleInfo

## ** [corrected] PCA already done (interactively, using `fastMNN`) - took top 100 PCs

# getClusteredPCs() to identify working PC space
pc.choice.hpc <- getClusteredPCs(reducedDim(sce.hpc, "PCA_corrected"))

# How many PCs should use in this space?
metadata(pc.choice.hpc)$chosen


    ## Plot n Clusters vs. d PCs
    pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/revision/regionSpecific_HPC-n3_getClusteredPCs_MNT2021.pdf")
    plot(pc.choice.hpc$n.pcs, pc.choice.hpc$n.clusters,
         main=paste0("Combined HPC (n=3) samples (d PCs choice = ", metadata(pc.choice.hpc)$chosen, ")"))
    abline(v=metadata(pc.choice.hpc)$chosen, col="red", lty="dashed", lwd=0.8)  
    dev.off()


# Save
save(sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, ref.sampleInfo,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_HPC-n3_cleaned-combined_SCE_MNT2021.rda")

rm(list=ls())
sessionInfo()
