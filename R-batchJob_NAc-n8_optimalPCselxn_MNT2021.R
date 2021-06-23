################################################################################
### LIBD pilot 10x snRNA-seq: Nucleus accumbens samples - 8 samples (as of revision)
###   **Region-specific analyses**
###   - R-batch job for detxn of optimal PC space with 'sce.nac' object
###         -> see '10x_NAc-n8_step02_clust-annot_MNT.R'
###            for setup of the SCE
### MNT 23Jun2021
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


load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda",
     verbose=TRUE)
    # sce.nac, chosen.hvgs.nac, ref.sampleInfo

## ** [corrected] PCA already done (interactively, using `fastMNN`) - took top 100 PCs **

# getClusteredPCs() to identify working PC space
pc.choice.nac <- getClusteredPCs(reducedDim(sce.nac, "PCA_corrected"))

# How many PCs should use in this space?
metadata(pc.choice.nac)$chosen


## Plot n Clusters vs. d PCs
pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/revision/regionSpecific_NAc-n8_getClusteredPCs_MNT2021.pdf")
plot(pc.choice.nac$n.pcs, pc.choice.nac$n.clusters,
     main=paste0("Combined NAc (n=8) samples (d PCs choice = ", metadata(pc.choice.nac)$chosen, ")"))
abline(v=metadata(pc.choice.nac)$chosen, col="red", lty="dashed", lwd=0.8)  
dev.off()


# Save
save(sce.nac, chosen.hvgs.nac, ref.sampleInfo, pc.choice.nac,
     file="rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda")

rm(list=ls())
sessionInfo()
