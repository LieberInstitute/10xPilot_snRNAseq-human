### MNT 10x snRNA-seq workflow: step 02
###   **Region-specific analyses**
###     - (2x) sACC samples from: Br5161 & Br5212
### Initiated MNT 29Jan2020
### MNT 29Apr2021: add expansion samples (n=3, incl'g 2 female)
#####################################################################

library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(batchelor)
library(DropletUtils)
library(jaffelab)
library(dendextend)
library(dynamicTreeCut)

### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

# ===


# Load 'pilot' samples
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/all-FACS-n14_preprint_SCEs_processing-QC_MNTMar2021.rda",
     verbose=T)
    # pilot.data, pilot.data.unfiltered, e.out, ref.sampleInfo
    rm(pilot.data.unfiltered, e.out)

# Load 2021 expansion set
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/all-FACS-n10_2021rev_SCEs_processing-QC_MNTMar2021.rda", verbose=T)
    # pilot.data.2, pilot.data.2.unfiltered, e.out.2, ref.sampleInfo.rev
    rm(pilot.data.2.unfiltered, e.out.2)


### MNT comment: At this point, each sample (which is a SCE object in the list, 'pilot.data') has been
  #              QC'd for cell/nucleus calling ('emptyDrops()' test) and mito rate thresholding
  #   Additionally, there have been a computed 'doubletScore', which will QC with, after
  #   clustering (e.g. that there are no doublet-driven clusters, etc.)

### Merging shared-region samples ============================================
    # Newest iterations for normalization: multiBatchNorm-alize
    
# Order as will do for `fastMNN()` (this shouldn't matter here)
sce.sacc <- cbind(pilot.data.2[["br5400.sacc"]], pilot.data[["br5212.sacc"]],
                  pilot.data[["br5161.sacc"]],
                  pilot.data.2[["br5701.sacc.neun"]], pilot.data.2[["br5276.sacc.neun"]]
)

sce.sacc
    #class: SingleCellExperiment 
    # dim: 33538 15669 
    # metadata(5): Samples Samples Samples Samples Samples
    # assays(1): counts
    # rownames(33538): MIR1302-2HG FAM138A ... AC213203.1 FAM231C
    # rowData names(6): gene_id gene_version ... gene_biotype Symbol.uniq
    # colnames(15669): AAACCCAAGAGTCTTC-1 AAACCCAAGGCCCGTT-1 ...
    #   TTTGACTGTATCGTGT-1 TTTGGTTAGCAGCACA-1
    # colData names(16): Sample Barcode ... protocol sequencer
    # reducedDimNames(0):
    # altExpNames(0):

# Use `multiBatchNorm()` to compute log-normalized counts, matching the scaling across samples
sce.sacc <- multiBatchNorm(sce.sacc, batch=sce.sacc$sampleID)

# Use the simple `modelGeneVar` - this makes more sense over `combineVar`, since the
#   cell composition is already known to be quite different (with NeuN selection)
geneVar.sacc <- modelGeneVar(sce.sacc)
chosen.hvgs.sacc <- geneVar.sacc$bio > 0
sum(chosen.hvgs.sacc)
    # [1] 11898


### Dimensionality reduction ================================================================

# Run `fastMNN` (internally uses `multiBatchPCA`), taking top 100 (instead of default 50 PCs)
set.seed(109)
mnn.hold <-  fastMNN(sce.sacc, batch=sce.sacc$sampleID,
                     merge.order=c("br5400.sacc","br5212.sacc","br5161.sacc",
                                   "br5701.sacc.neun","br5276.sacc.neun"),
                     subset.row=chosen.hvgs.sacc, d=100,
                     correct.all=TRUE, get.variance=TRUE,
                     BSPARAM=BiocSingular::IrlbaParam())
    # This temp file just used for getting batch-corrected components (drops a variety of entries)

table(colnames(mnn.hold) == colnames(sce.sacc))  # all TRUE
table(mnn.hold$batch == sce.sacc$sampleID) # all TRUE

# Add them to the SCE, as well as the metadata (though the latter might not be so usefl)
reducedDim(sce.sacc, "PCA_corrected") <- reducedDim(mnn.hold, "corrected") # 100 components
metadata(sce.sacc) <- metadata(mnn.hold)

# Save into a new region-specific SCE object/file
save(sce.sacc, chosen.hvgs.sacc, ref.sampleInfo, ref.sampleInfo.rev,
     file="rdas/revision/regionSpecific_sACC-n5_cleaned-combined_SCE_MNT2021.rda")


    ## 'getClusteredPCs()' evaluated in qsub mode (with 'R-batchJob_sACC-n2_optimalPCselxn_MNTFeb2020.R')
    #    --> saved into same .rda


### Picking up with optimally-defined PC space ===
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=TRUE)
    # sce.sacc, chosen.hvgs.sacc, pc.choice.sacc


# How many PCs is optimal?:
metadata(pc.choice.sacc)$chosen # [1] 96

## Assign this chosen ( PCs) to 'PCA_opt'
reducedDim(sce.sacc, "PCA_opt") <- reducedDim(sce.sacc, "PCA")[ ,1:(metadata(pc.choice.sacc)$chosen)]


## t-SNE
set.seed(109)
sce.sacc <- runTSNE(sce.sacc, dimred="PCA_opt")


## UMAP
set.seed(109)
sce.sacc <- runUMAP(sce.sacc, dimred="PCA_opt")

## Load in phenodata from pan-brain analysis -> colData for downstream use
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda",
     verbose=T)
# Want 'ref.sampleInfo'

sce.sacc$region <- ss(sce.sacc$sample,".5",1)
sce.sacc$donor <- paste0("Br",ss(sce.sacc$sample,"cc.",2))
sce.sacc$processDate <- ref.sampleInfo$realBatch[match(sce.sacc$sample, ref.sampleInfo$sampleID)]
sce.sacc$protocol <- ref.sampleInfo$protocol[match(sce.sacc$processDate, ref.sampleInfo$realBatch)]


# Save for now
save(sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, ref.sampleInfo, 
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda")


### Clustering: Two-step ======================================================
### Step 1: Perform graph-based clustering in this optimal PC space
#         - take k=20 NN to build graph
snn.gr <- buildSNNGraph(sce.sacc, k=20, use.dimred="PCA_opt")
clusters.k20 <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters.k20)
    ## 46 prelim clusters

# Assign as 'prelimCluster'
sce.sacc$prelimCluster <- factor(clusters.k20)

# Is sample driving this 'high-res' clustering at this level?
table(sce.sacc$prelimCluster, sce.sacc$sample)

table(sce.sacc$sample)

### Step 2: Hierarchical clustering of pseudo-bulked ("PB'd") counts with most robust normalization
  #         (as determined in: 'side-Rscript_testingStep2_HC-normalizn-approaches_wAmygData_MNTJan2020.R')
#           ** That is, to pseudo-bulk (aka 'cluster-bulk') on raw counts, on all [non-zero] genes,
#              normalize with `librarySizeFactors()`, log2-transform, then perform HC'ing


# Preliminary cluster index for pseudo-bulking
clusIndexes = splitit(sce.sacc$prelimCluster)
prelimCluster.PBcounts <- sapply(clusIndexes, function(ii){
  rowSums(assays(sce.sacc)$counts[ ,ii])
  }
)

# And btw...
table(rowSums(prelimCluster.PBcounts)==0)
#FALSE  TRUE
#27623  5915

# Compute LSFs at this level
sizeFactors.PB.all  <- librarySizeFactors(prelimCluster.PBcounts)

# Normalize with these LSFs
geneExprs.temp <- t(apply(prelimCluster.PBcounts, 1, function(x) {log2(x/sizeFactors.PB.all + 1)}))


## Perform hierarchical clustering
dist.clusCollapsed <- dist(t(geneExprs.temp))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")
#tree.clusCollapsed$labels

dend <- as.dendrogram(tree.clusCollapsed, hang=0.2)

# Just for observation
par(cex=.6)
myplclust(tree.clusCollapsed, main="2x sACC prelim-kNN-cluster relationships", cex.main=2, cex.lab=1.5, cex=1.8)


clust.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                               minClusterSize=2, deepSplit=1, cutHeight=550)


table(clust.treeCut)
unname(clust.treeCut[order.dendrogram(dend)])
    ## Cutting at 550 looks good - go ahead and proceed with this

# Add new labels to those prelimClusters cut off
clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)]==0)] <- max(clust.treeCut)+c(1)

# A posteriori: re-assign prelimCluster 3 as to not collapse with 23, 29
clust.treeCut[order.dendrogram(dend)][which(order.dendrogram(dend)==3)] <- max(clust.treeCut)+c(1)

labels_colors(dend) <- tableau20[clust.treeCut[order.dendrogram(dend)]]

# Print for future reference
pdf("pdfs/regionSpecific_sACC-n2_HC-prelimCluster-relationships_Feb2020.pdf")
par(cex=1, font=2)
plot(dend, main="2x sACC prelim-kNN-cluster relationships")
dev.off()


# Make reference for new cluster assignment
clusterRefTab.sacc <- data.frame(origClust=order.dendrogram(dend),
                                merged=clust.treeCut[order.dendrogram(dend)])

# Assign as 'collapsedCluster'
sce.sacc$collapsedCluster <- factor(clusterRefTab.sacc$merged[match(sce.sacc$prelimCluster, clusterRefTab.sacc$origClust)])


# Print some visualizations:
pdf("pdfs/regionSpecific_sACC-n2_reducedDims-with-collapsedClusters_Feb2020.pdf")
plotReducedDim(sce.sacc, dimred="PCA", ncomponents=5, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.sacc, colour_by="sample", point_alpha=0.5)
plotTSNE(sce.sacc, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.sacc, colour_by="sum", point_alpha=0.5)
plotUMAP(sce.sacc, colour_by="collapsedCluster", point_alpha=0.5)
dev.off()


## Print marker genes for annotation
markers.mathys.custom = list(
  'neurons' = c('SYT1', 'SNAP25', 'GRIN1'),
  'excitatory_neuron' = c('CAMK2A', 'NRGN','SLC17A7', 'SLC17A6', 'SLC17A8'),
  'inhibitory_neuron' = c('GAD1', 'GAD2', 'SLC32A1'),
  'oligodendrocyte' = c('MBP', 'MOBP', 'PLP1'),
  'oligodendrocyte_precursor' = c('PDGFRA', 'VCAN', 'CSPG4'),
  'microglia' = c('CD74', 'CSF1R', 'C3'),
  'astrocytes' = c('GFAP', 'TNC', 'AQP4', 'SLC1A2'),
  'endothelial' = c('CLDN5', 'FLT1', 'VTN')
)

pdf("pdfs/regionSpecific_sACC-n2_marker-logExprs_collapsedClusters_Mar2020.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.sacc, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="collapsedCluster", colour_by="collapsedCluster", point_alpha=0.5, point_size=.7,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:11], length(markers.mathys.custom[[i]])))
  )
}
dev.off()

    ## Observation: microglia not really seen at this level of collapsing...

    # Looking at collapsedCluster 7
    sce.sacc.cc7 <- sce.sacc[ ,sce.sacc$collapsedCluster==7]
    unique(sce.sacc.cc7$prelimCluster)
    
    plotExpression(sce.sacc.cc7, exprs_values = "logcounts", features=c(markers.mathys.custom[["microglia"]]),
                   x="prelimCluster", colour_by="prelimCluster", point_alpha=0.5, point_size=.7,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:3], length(markers.mathys.custom[["microglia"]])))
    
    plotExpression(sce.sacc.cc7, exprs_values = "logcounts", features=c(markers.mathys.custom[["oligodendrocyte"]]),
                   x="prelimCluster", colour_by="prelimCluster", point_alpha=0.5, point_size=.7,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:3], length(markers.mathys.custom[["microglia"]])))

        ## So it looks like prelimCluster 3 is the one.  However at any cutHeights from 400:550, this
         # would never be separated from the oligos, so manually assign, going back to line 176
         #      (and re-name these plots with prefix 'zold_')
    
    
    

## Add annotations, looking at marker gene expression
annotationTab.sacc <- data.frame(cluster=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
                                 cellType=c("Excit.1", "Inhib.1", "Excit.2", "Inhib.2", "Excit.3",
                                            "Astro", "Oligo", "OPC", "Excit.4", "Ambig.lowNtrxts",
                                            "Micro")
                                 )

sce.sacc$cellType <- annotationTab.sacc$cellType[match(sce.sacc$collapsedCluster,
                                                       annotationTab.sacc$cluster)]




## Save for now MNT 21Feb2020
save(sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc, ref.sampleInfo,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda")



# Observation: cluster 6 not an obvious 'cell type'...
newClusIndex <- splitit(sce.sacc$collapsedCluster)
sapply(newClusIndex, function(x) {quantile(sce.sacc[,x]$sum)})
    #             1      2         3       4      5       6        7       8     9
    #0%     3985.00   2448   2175.00  3595.0   4625   464.0   514.00  1276.0  5475
    #25%   36862.00  25268  35064.75 24582.5  34369  4807.5  4368.50  8600.5 36320
    #50%   56701.00  34905  46186.50 30895.0  53997  7014.0  6131.50 11122.0 47976
    #75%   74135.75  44749  62363.00 38684.0  68638  9473.5  8041.25 13572.5 58484
    #100% 196431.00 121477 127499.00 78788.0 118374 27467.0 25828.00 34023.0 97538
    #         10      11
    #0%    103.0  543.00
    #25%   248.5 3072.25
    #50%   455.0 4669.50
    #75%   817.5 5761.50
    #100% 4828.0 9973.00



### MNT 25Mar2020 === === ===
# Re-print marker expression with cell type labels and dropping 'ambig.lowNtrxts' cluster
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)

table(sce.sacc$cellType)

# First drop "Ambig.lowNtrxts" (43 nuclei)
sce.sacc <- sce.sacc[ ,sce.sacc$cellType != "Ambig.lowNtrxts"]
sce.sacc$cellType <- droplevels(sce.sacc$cellType)

pdf("pdfs/regionSpecific_sACC-n2_marker-logExprs_collapsedClusters_Mar2020.pdf", height=6, width=12)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.sacc, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,
                   add_legend=F) +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
                   width = 0.3, colour=rep(tableau10medium[1:10], length(markers.mathys.custom[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(markers.mathys.custom)[i], " markers"))
  )
}
dev.off()


      ## -> proceed to 'step03_markerDetxn-analyses[...].R'



### Session info for 29Apr2021 ==========================================================
sessionInfo()
# R version 4.0.4 RC (2021-02-08 r79975)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS:   /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/lib/libRblas.so
# LAPACK: /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/lib/libRlapack.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices datasets  utils    
# [8] methods   base     
# 
# other attached packages:
#   [1] dynamicTreeCut_1.63-1       dendextend_1.14.0          
# [3] jaffelab_0.99.30            rafalib_1.0.0              
# [5] DropletUtils_1.10.3         batchelor_1.6.2            
# [7] scran_1.18.5                scater_1.18.6              
# [9] ggplot2_3.3.3               EnsDb.Hsapiens.v86_2.99.0  
# [11] ensembldb_2.14.1            AnnotationFilter_1.14.0    
# [13] GenomicFeatures_1.42.3      AnnotationDbi_1.52.0       
# [15] SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0
# [17] Biobase_2.50.0              GenomicRanges_1.42.0       
# [19] GenomeInfoDb_1.26.7         IRanges_2.24.1             
# [21] S4Vectors_0.28.1            BiocGenerics_0.36.1        
# [23] MatrixGenerics_1.2.1        matrixStats_0.58.0         
# 
# loaded via a namespace (and not attached):
#   [1] googledrive_1.0.1         ggbeeswarm_0.6.0         
# [3] colorspace_2.0-0          ellipsis_0.3.1           
# [5] scuttle_1.0.4             bluster_1.0.0            
# [7] XVector_0.30.0            BiocNeighbors_1.8.2      
# [9] rstudioapi_0.13           bit64_4.0.5              
# [11] fansi_0.4.2               xml2_1.3.2               
# [13] splines_4.0.4             R.methodsS3_1.8.1        
# [15] sparseMatrixStats_1.2.1   cachem_1.0.4             
# [17] Rsamtools_2.6.0           ResidualMatrix_1.0.0     
# [19] dbplyr_2.1.1              R.oo_1.24.0              
# [21] HDF5Array_1.18.1          compiler_4.0.4           
# [23] httr_1.4.2                dqrng_0.2.1              
# [25] assertthat_0.2.1          Matrix_1.3-2             
# [27] fastmap_1.1.0             lazyeval_0.2.2           
# [29] limma_3.46.0              BiocSingular_1.6.0       
# [31] prettyunits_1.1.1         tools_4.0.4              
# [33] rsvd_1.0.3                igraph_1.2.6             
# [35] gtable_0.3.0              glue_1.4.2               
# [37] GenomeInfoDbData_1.2.4    dplyr_1.0.5              
# [39] rappdirs_0.3.3            Rcpp_1.0.6               
# [41] vctrs_0.3.6               Biostrings_2.58.0        
# [43] rhdf5filters_1.2.0        rtracklayer_1.50.0       
# [45] DelayedMatrixStats_1.12.3 stringr_1.4.0            
# [47] beachmat_2.6.4            lifecycle_1.0.0          
# [49] irlba_2.3.3               statmod_1.4.35           
# [51] XML_3.99-0.6              edgeR_3.32.1             
# [53] zlibbioc_1.36.0           scales_1.1.1             
# [55] hms_1.0.0                 ProtGenerics_1.22.0      
# [57] rhdf5_2.34.0              RColorBrewer_1.1-2       
# [59] curl_4.3                  memoise_2.0.0            
# [61] gridExtra_2.3             segmented_1.3-3          
# [63] biomaRt_2.46.3            stringi_1.5.3            
# [65] RSQLite_2.2.7             BiocParallel_1.24.1      
# [67] rlang_0.4.10              pkgconfig_2.0.3          
# [69] bitops_1.0-7              lattice_0.20-41          
# [71] purrr_0.3.4               Rhdf5lib_1.12.1          
# [73] GenomicAlignments_1.26.0  bit_4.0.4                
# [75] tidyselect_1.1.0          magrittr_2.0.1           
# [77] R6_2.5.0                  generics_0.1.0           
# [79] DelayedArray_0.16.3       DBI_1.1.1                
# [81] pillar_1.6.0              withr_2.4.2              
# [83] RCurl_1.98-1.3            tibble_3.1.1             
# [85] crayon_1.4.1              utf8_1.2.1               
# [87] BiocFileCache_1.14.0      viridis_0.6.0            
# [89] progress_1.2.2            locfit_1.5-9.4           
# [91] grid_4.0.4                blob_1.2.1               
# [93] R.utils_2.10.1            openssl_1.4.3            
# [95] munsell_0.5.0             beeswarm_0.3.1           
# [97] viridisLite_0.4.0         vipor_0.4.5              
# [99] askpass_1.1  



