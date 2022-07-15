################################################################################
### LIBD pilot 10x snRNA-seq: Amygdala samples
###   **Region-specific analyses**
###   - R-batch job for detxn of optimal PC space with 'sce.amy' object
###         -> see '10x_Amyg-n5_step02_clust-annot_MNT.R'
###            for setup of the SCE
### MNT 21Apr2021
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


load("/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_Amyg-n5_cleaned-combined_SCE_MNT2021.rda",
     verbose=TRUE)
    # sce.amy, chosen.hvgs.amy, ref.sampleInfo, ref.sampleInfo.rev

## ** [corrected] PCA already done (interactively, using `fastMNN`) - took top 100 PCs **

# getClusteredPCs() to identify working PC space
pc.choice.amy <- getClusteredPCs(reducedDim(sce.amy, "PCA_corrected"))

# How many PCs should use in this space?
metadata(pc.choice.amy)$chosen


    ## Plot n Clusters vs. d PCs
    pdf("/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/revision/regionSpecific_Amyg-n5_getClusteredPCs_MNT2021.pdf")
    plot(pc.choice.amy$n.pcs, pc.choice.amy$n.clusters,
         main=paste0("Combined Amygdala (n=5) samples (d PCs choice = ", metadata(pc.choice.amy)$chosen, ")"))
    abline(v=metadata(pc.choice.amy)$chosen, col="red", lty="dashed", lwd=0.8)  
    dev.off()


# Save
save(sce.amy, chosen.hvgs.amy, ref.sampleInfo, ref.sampleInfo.rev, pc.choice.amy,
     file="/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_Amyg-n5_cleaned-combined_SCE_MNT2021.rda")

rm(list=ls())
sessionInfo()
