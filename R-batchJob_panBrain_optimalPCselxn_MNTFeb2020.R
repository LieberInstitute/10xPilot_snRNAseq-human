################################################################################
### LIBD pilot 10x-Frankenstein (n=12) snRNA-seq samples
###   ** PAN-BRAIN ANALYSIS **
###     - R-batch job for detxn of optimal PC space with ~'sce.pilot.n12'~ object
### MNT 07-08Feb2020 (now on 'sce.all.n12')              ~29-30Jan2020~
###     - change: using MBNorm to LSFs
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


load("rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda", verbose=TRUE)
    # sce.all.n12, chosen.hvgs.all.n12

## PCA already done (interactively) - took top 250 PCs

## getClusteredPCs() to identify working PC space
pc.choice.n12 <- getClusteredPCs(reducedDim(sce.all.n12))

# How many PCs should use at the pan-brain level?
metadata(pc.choice.n12)$chosen


    ## Plot $nClusters vs. nPCs
    pdf("pdfs/all-FACS-homogenates_n12_getClusteredPCs-results-w250pcs_MNTFeb2020.pdf")
    plot(pc.choice.n12$n.pcs, pc.choice.n12$n.clusters,
         main=paste0("Combined n=12 brain FACS-hom. samples (d PCs choice = ", metadata(pc.choice.n12)$chosen, ")"),
         cex=0.8)
    abline(v=metadata(pc.choice.n12)$chosen, col="red", lty="dashed", lwd=0.8)  
    dev.off()


# Save
save(sce.all.n12, chosen.hvgs.all.n12, pc.choice.n12, file="rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda")

rm(list=ls())
sessionInfo()