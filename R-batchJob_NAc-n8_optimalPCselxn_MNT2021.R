################################################################################
### LIBD pilot 10x snRNA-seq: Nucleus accumbens samples - ALL FIVE
###   **Region-specific analyses**
###   - R-batch job for detxn of optimal PC space with 'sce.nac.all' object
###         -> see '10x-pilot_region-specific_NAc-n3_with_NeuN-sorted-n2_step02_clust-annot_MNTMar2020.R'
###            for setup of the SCE
### MNT 23Mar2020
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


load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/NeuN-sortedNAc_n2_cleaned-combined_MNTMar2020.rda",
     verbose=TRUE)
    # sce.nac.all, chosen.hvgs.nac.all

sce.nac.all
table(sce.nac.all$sample)   # 5 samples
dim(reducedDim(sce.nac.all, "PCA"))
    ## PCA already done (interactively) - took top 100 PCs

## getClusteredPCs() to identify working PC space
pc.choice.nac.all <- getClusteredPCs(reducedDim(sce.nac.all))

# How many PCs should use in this space?
metadata(pc.choice.nac.all)$chosen


    ## Plot n Clusters vs. d PCs
    pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/regionSpecific_NAc-ALL-n5_getClusteredPCs-results-w100pcs_MNTMar2020.pdf")
    plot(pc.choice.nac.all$n.pcs, pc.choice.nac.all$n.clusters,
         main=paste0("Combined NAc (n=3 homs, n=2 NeuN-sorts) samples (d PCs choice = ", metadata(pc.choice.nac.all)$chosen, ")"))
    abline(v=metadata(pc.choice.nac.all)$chosen, col="red", lty="dashed", lwd=0.8)  
    dev.off()


# Save
save(sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda")

rm(list=ls())
sessionInfo()
