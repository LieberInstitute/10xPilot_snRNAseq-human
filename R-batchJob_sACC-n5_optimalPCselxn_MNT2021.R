################################################################################
### LIBD pilot 10x snRNA-seq: sACC samples
###   **Region-specific analyses**
###   - R-batch job for detxn of optimal PC space with 'sce.sacc' object
###         -> see '10x_sACC-n5_step02_clust-annot_MNT.R'
###            for setup of the SCE
### MNT 29Apr2021
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


load("/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_sACC-n5_cleaned-combined_SCE_MNT2021.rda",
     verbose=TRUE)
    # sce.sacc, chosen.hvgs.sacc, ref.sampleInfo, ref.sampleInfo.rev,

## ** [corrected] PCA already done (interactively, using `fastMNN`) - took top 100 PCs **

## getClusteredPCs() to identify working PC space
pc.choice.sacc <- getClusteredPCs(reducedDim(sce.sacc, "PCA_corrected"))

# How many PCs should use in this space?
metadata(pc.choice.sacc)$chosen


    ## Plot n Clusters vs. d PCs
    pdf("/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/revision/regionSpecific_sACC-n5_getClusteredPCs_MNT2021.pdf")
    plot(pc.choice.sacc$n.pcs, pc.choice.sacc$n.clusters,
         main=paste0("Combined sACC (n=5) samples (d PCs choice = ", metadata(pc.choice.sacc)$chosen, ")"))
    abline(v=metadata(pc.choice.sacc)$chosen, col="red", lty="dashed", lwd=0.8)  
    dev.off()


# Save
save(sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, ref.sampleInfo, ref.sampleInfo.rev,
     file="/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_sACC-n5_cleaned-combined_SCE_MNT2021.rda")

rm(list=ls())
sessionInfo()
