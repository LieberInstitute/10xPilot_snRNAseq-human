################################################################################
### LIBD pilot 10x snRNA-seq: DLPFC samples
###   **Region-specific analyses**
###   - R-batch job for detxn of optimal PC space with 'sce.dlpfc' object
###         -> see '10x-pilot_region-specific_DLPFC_step02_clust-annot_MNTJan2020.R'
###            for setup of the SCE
### LAH 03May2021
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


load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda",
     verbose=TRUE)
    # sce.dlpfc, chosen.hvgs.dlpfc

## PCA already done (interactively) - took top 100 PCs

## getClusteredPCs() to identify working PC space
pc.choice.dlpfc <- getClusteredPCs(reducedDim(sce.dlpfc))

# How many PCs should use in this space?
metadata(pc.choice.dlpfc)$chosen


## Plot n Clusters vs. d PCs
pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/revision/regionSpecific_DLPFC-n3_getClusteredPCs-results-w100pcs_LAH2021.pdf")
plot(pc.choice.dlpfc$n.pcs, pc.choice.dlpfc$n.clusters,
     main=paste0("Combined DLPFC (n=3) samples (d PCs choice = ", metadata(pc.choice.dlpfc)$chosen, ")"))
abline(v=metadata(pc.choice.dlpfc)$chosen, col="red", lty="dashed", lwd=0.8)  
dev.off()


# Save
save(sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, ref.sampleInfo, ref.sampleInfo.rev,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda")

# sgejobs::job_single('R-batchJob_DLPFC-n3_optimalPCselxn_LAH2021', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript R-batchJob_DLPFC-n3_optimalPCselxn_LAH2021.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
