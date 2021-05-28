### MNT 10x snRNA-seq workflow: step 02
###   **Region-specific analyses**
###     - (3x) HPC samples from: Br5161 & Br5212 & Br5287
### Initiated MNT 07Feb2020
### MNT 23Apr2021: Updated QC'd SCE (no add'l donors)
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

source("plotExpressionCustom.R")

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

    ### MNT comment: At this point, each sample (which is a SCE object in the list, 'pilot.data') has been
    #               QC'd for cell/nucleus calling ('emptyDrops()' test) and mito rate thresholding
    #   Additionally, there have been a computed 'doubletScore', which will QC with, after
    #   clustering (e.g. that there are no doublet-driven clusters, etc.)

    
### Merging shared-region samples ============================================
  # Newest iterations for normalization: multiBatchNorm-alize

sce.hpc <- cbind(pilot.data[["br5161.hpc"]], pilot.data[["br5212.hpc"]], pilot.data[["br5287.hpc"]])

sce.hpc
    # class: SingleCellExperiment 
    # dim: 33538 10268 
    # metadata(3): Samples Samples Samples
    # assays(1): counts
    # rownames(33538): MIR1302-2HG FAM138A ... AC213203.1 FAM231C
    # rowData names(6): gene_id gene_version ... gene_biotype Symbol.uniq
    # colnames(10268): AAACCCATCTGTCAGA-1 AAACCCATCTGTCGCT-1 ...
    #   TTTGGTTGTGGTCCGT-1 TTTGTTGCAGAAACCG-1
    # colData names(16): Sample Barcode ... protocol sequencer
    # reducedDimNames(0):
    # altExpNames(0):    

# Use `multiBatchNorm()` to compute log-normalized counts, matching the scaling across samples
sce.hpc <- multiBatchNorm(sce.hpc, batch=sce.hpc$sampleID)

# Use the simple `modelGeneVar` - this makes more sense over `combineVar`, since the
#   cell composition is already known to be quite different (with NeuN selection)
geneVar.hpc <- modelGeneVar(sce.hpc)
chosen.hvgs.hpc <- geneVar.hpc$bio > 0
sum(chosen.hvgs.hpc)
    # [1] 8696


### Dimensionality reduction ================================================================

# Run `fastMNN` (internally uses `multiBatchPCA`), taking top 100 (instead of default 50 PCs)
set.seed(109)
mnn.hold <-  fastMNN(sce.hpc, batch=sce.hpc$sampleID,
                     merge.order=c("br5161.hpc","br5212.hpc","br5287.hpc"),
                     subset.row=chosen.hvgs.hpc, d=100,
                     correct.all=TRUE, get.variance=TRUE,
                     BSPARAM=BiocSingular::IrlbaParam())
    # This temp file just used for getting batch-corrected components (drops a variety of entries)

table(colnames(mnn.hold) == colnames(sce.hpc))  # all TRUE
table(mnn.hold$batch == sce.hpc$sampleID) # all TRUE

# Add them to the SCE, as well as the metadata (though the latter might not be so usefl)
reducedDim(sce.hpc, "PCA_corrected") <- reducedDim(mnn.hold, "corrected") # 100 components
metadata(sce.hpc) <- metadata(mnn.hold)

# Save into a new region-specific SCE object/flie
save(sce.hpc, chosen.hvgs.hpc, ref.sampleInfo,
     file="rdas/revision/regionSpecific_HPC-n3_cleaned-combined_SCE_MNT2021.rda")


    ## 'getClusteredPCs()' evaluated in qsub mode (with 'R-batchJob_HPC-n3_optimalPCselxn_MNT2021.R')
    #    --> saved into same .rda


### Picking up with optimally-defined PC space ===
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_HPC-n3_cleaned-combined_SCE_MNT2021.rda",
     verbose=TRUE)
    # sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, ref.sampleInfo


# How many PCs is optimal?:
metadata(pc.choice.hpc)$chosen
    ## 59

## Assign this chosen ( PCs) to 'PCA_opt'
reducedDim(sce.hpc, "PCA_opt") <- reducedDim(sce.hpc, "PCA_corrected")[ ,1:(metadata(pc.choice.hpc)$chosen)]


## t-SNE
set.seed(109)
sce.hpc <- runTSNE(sce.hpc, dimred="PCA_opt")


## UMAP
set.seed(109)
sce.hpc <- runUMAP(sce.hpc, dimred="PCA_opt")

# How do these look?
plotReducedDim(sce.hpc, dimred="TSNE", colour_by="sampleID")
plotReducedDim(sce.hpc, dimred="TSNE", colour_by="prelimCluster")
plotReducedDim(sce.hpc, dimred="UMAP", colour_by="sampleID")


### Clustering: Two-step ======================================================
### Step 1: Perform graph-based clustering in this optimal PC space
#         - take k=20 NN to build graph
snn.gr <- buildSNNGraph(sce.hpc, k=20, use.dimred="PCA_opt")
clusters.k20 <- igraph::cluster_walktrap(snn.gr)$membership

table(clusters.k20)
    ##    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
    #   236  302  117  767  111   69  112   35   35 1830   26  123  117   26 
    #    15   16   17   18   19   20   21   22   23   24   25   26   27   28 
    #   168   56 1512   84   51   27  892   43   68   34    6  269    5   36 
    #    29   30   31   32   33   34   35   36   37   38   39   40   41   42 
    #    32   46  231   33   15 1703   91   23   30  117   50    6   20    8 
    #    43   44   45   46   47   48   49   50 
    #     5    2   26  567   54    6   17   29

# Assign as 'prelimCluster'
sce.hpc$prelimCluster <- factor(clusters.k20)

# Is sample driving this 'high-res' clustering at this level?
table(sce.hpc$prelimCluster, sce.hpc$sampleID)  

table(sce.hpc$sampleID)
    #br5161.hpc br5212.hpc br5287.hpc 
    #      4421       3977       1870

### Step 2: Hierarchical clustering of pseudo-bulked ("PB'd") counts
#           ** That is, to pseudo-bulk (aka 'cluster-bulk') on raw counts, on all [non-zero] genes,
#              normalize with `librarySizeFactors()`, log2-transform, then perform HC'ing


# Preliminary cluster index for pseudo-bulking
clusIndexes = splitit(sce.hpc$prelimCluster)
prelimCluster.PBcounts <- sapply(clusIndexes, function(ii){
  rowSums(assays(sce.hpc)$counts[ ,ii])
  }
)

# Compute LSFs at this level
sizeFactors.PB.all  <- librarySizeFactors(prelimCluster.PBcounts)

# Normalize with these LSFs
geneExprs.temp <- t(apply(prelimCluster.PBcounts, 1, function(x) {log2(x/sizeFactors.PB.all + 1)}))

## Perform hierarchical clustering
dist.clusCollapsed <- dist(t(geneExprs.temp))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

dend <- as.dendrogram(tree.clusCollapsed, hang=0.2)

# Just for observation
par(cex=.6)
myplclust(tree.clusCollapsed, main="3x HPC prelim-kNN-cluster relationships", cex.main=2, cex.lab=1.5, cex=1.8)


clust.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                               minClusterSize=2, deepSplit=1, cutHeight=325)


table(clust.treeCut)
unname(clust.treeCut[order.dendrogram(dend)])
    ## Cutting at 325 looks the best - go ahead and proceed with this

    # (post-hoc:) All/most the inhibs are collapsed to one cluster... let's cut shallower for these
    neuron.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                                  minClusterSize=1, deepSplit=1, cutHeight=200)
    unname(neuron.treeCut[order.dendrogram(dend)])
    
    # Take those and re-assign to the first assignments; leave 20,24,45 together
    clust.treeCut[order.dendrogram(dend)][c(6:21)] <- ifelse(neuron.treeCut[order.dendrogram(dend)][c(6:21)] == 0,
                                                              0, neuron.treeCut[order.dendrogram(dend)][c(6:21)] + 8)
    
    unname(clust.treeCut[order.dendrogram(dend)])


# Add new labels to those (18x) prelimClusters cut off
clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)]==0)] <-
  max(clust.treeCut)+c(1:8, 9,9, 10, 11,11, 12:14, 15,15)

# 'Re-write', since there are missing numbers
clust.treeCut[order.dendrogram(dend)] <- as.numeric(as.factor(clust.treeCut[order.dendrogram(dend)]))

## Define color pallet
cluster_colors <- unique(c(tableau10medium, tableau20)[clust.treeCut[order.dendrogram(dend)]])
names(cluster_colors) <- unique(clust.treeCut[order.dendrogram(dend)])
labels_colors(dend) <- cluster_colors[as.character(clust.treeCut[order.dendrogram(dend)])]

# Print for future reference
pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/revision/regionSpecific_HPC-n3_HC-prelimCluster-relationships_MNT2021.pdf",width=9)
par(cex=1.1, font=2)
plot(dend, main="3x HPC prelim-kNN-cluster relationships")
dev.off()


# Make reference for new cluster assignment
clusterRefTab.hpc <- data.frame(origClust=order.dendrogram(dend),
                                  merged=clust.treeCut[order.dendrogram(dend)])


# Assign as 'collapsedCluster'
sce.hpc$collapsedCluster <- factor(clusterRefTab.hpc$merged[match(sce.hpc$prelimCluster, clusterRefTab.hpc$origClust)])

# Print some visualizations:
pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/revision/regionSpecific_HPC-n3_reducedDims-with-collapsedClusters_MNT2021.pdf")
plotReducedDim(sce.hpc, dimred="PCA_corrected", ncomponents=5, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.hpc, colour_by="sampleID", point_alpha=0.5)
plotTSNE(sce.hpc, colour_by="protocol", point_alpha=0.5)
plotTSNE(sce.hpc, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.hpc, colour_by="sum", point_alpha=0.5)
plotTSNE(sce.hpc, colour_by="doubletScore", point_alpha=0.5)
# And some more informative UMAPs
plotUMAP(sce.hpc, colour_by="sampleID", point_alpha=0.5)
plotUMAP(sce.hpc, colour_by="collapsedCluster", point_alpha=0.5)
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
  'endothelial' = c('CLDN5', 'FLT1', 'VTN', 'PECAM1'),
  # Updated T cell markers - these seem to be more stable (across regions)
  'Tcell' = c('SKAP1', 'ITK', 'CD247')
)

pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/revision/regionSpecific_HPC-n3_marker-logExprs_collapsedClusters_MNT2021.pdf",
    height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.hpc, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="collapsedCluster", colour_by="collapsedCluster", point_alpha=0.5, point_size=.7,
                   add_legend=F) + stat_summary(fun = median, fun.min = median, fun.max = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:18], length(markers.mathys.custom[[i]])))
  )
}
dev.off()


## QC - How do the total UMI distribution look?
newClusIndex <- splitit(sce.hpc$collapsedCluster)
sapply(newClusIndex, function(x) {quantile(sce.hpc[,x]$sum)})
    #            1        2       3         4     5       6      7     8        9
    # 0%     2841.0   475.00  1940.0    981.00   454   535.0   7674  1127  1693.00
    # 25%   20273.5  4540.25  5631.5  25939.00  3580  6381.5  13944  3609  7598.50
    # 50%   30830.0  6347.50  7996.5  33836.50  4661  9330.0  35400  5767 12681.00
    # 75%   41428.0  8634.50 11802.5  44510.25  5966 12623.5  53016  8420 14749.25
    # 100% 111453.0 29556.00 30085.0 114615.00 12463 27854.0 150461 20088 52661.00
    #           10       11       12      13   14      15    16      17      18
    # 0%    2665.0  1651.00  3641.00 1105.00  451   621.0 11237  1380.0   434.0
    # 25%   4124.5  2478.50  7335.00 2472.25  915  3720.5 13892  4535.5   545.5
    # 50%   4814.0  6049.00  9693.00 3215.00 1459  5352.0 14277  9119.0   913.0
    # 75%   7011.0 19584.75 15610.25 3756.50 1806  7971.0 18065 12166.5  1097.0
    # 100% 10390.0 58784.00 48371.00 6038.00 4235 17296.0 30332 20135.0 54810.0

table(sce.hpc$collapsedCluster)
    #   1    2    3    4    5    6    7    8    9   10   11   12   13   14 
    # 331 5912  936  486 1161  823   87  234    6   35    6   38   26  105 
    #  15   16   17   18 
    #  43    5   15   19


## doublet score?
sapply(newClusIndex, function(x) {round(quantile(sce.hpc$doubletScore[x]),2)})
    #         1     2     3    4    5     6    7    8    9   10   11   12
    # 0%   0.03  0.00  0.00 0.02 0.00  0.00 0.04 0.01 0.72 0.03 0.80 0.27
    # 25%  0.77  0.20  0.08 0.09 0.03  0.04 0.07 0.06 1.39 0.04 0.86 0.74
    # 50%  1.14  0.51  0.15 0.16 0.07  0.12 0.20 0.08 1.51 0.04 1.01 0.92
    # 75%  1.38  1.12  0.25 0.54 0.14  0.22 0.68 0.12 1.54 0.05 1.12 0.99
    # 100% 6.93 15.58 16.04 7.74 4.90 12.03 1.44 2.48 1.89 0.67 6.98 6.50
    #        13   14   15    16   17   18
    # 0%   0.05 0.00 0.05  6.77 0.15 0.00
    # 25%  0.06 0.02 0.27  6.80 0.56 0.01
    # 50%  0.09 0.04 0.34  6.93 0.69 0.02
    # 75%  0.11 0.05 0.41  7.84 1.35 0.04
    # 100% 0.12 0.66 1.98 10.81 2.46 1.17

# At 'prelimCluster' level?
clusIndex <- splitit(sce.hpc$prelimCluster)
sapply(clusIndex, function(x) {round(quantile(sce.hpc$doubletScore[x]),2)})
    # 29 & 47 for sure - 41 & 43 relatively high--43 being 'collapsedCluster' 16


## Add annotations, looking at marker gene expression
 #    (canonical, above, and from markers looked at in 'step3')
annotationTab.hpc <- data.frame(collapsedCluster=c(1:23))
annotationTab.hpc$cellType <- NA
annotationTab.hpc$cellType[c(1)] <- paste0("Oligo")
annotationTab.hpc$cellType[c(3,4)] <- c("Micro","OPC")
annotationTab.hpc$cellType[c(2,6)] <- paste0("Astro_", c("A","B"))
annotationTab.hpc$cellType[c(5,8:12,14,15)] <- paste0("Excit_", c("A","B","C","D","E","F","G","H"))
annotationTab.hpc$cellType[c(7,13,16,17)] <- paste0("Inhib_", c("A","B","C","D"))
annotationTab.hpc$cellType[18] <- "Tcell"
annotationTab.hpc$cellType[c(19,23)] <- paste0("drop.lowNTx_", c("A","B"))
annotationTab.hpc$cellType[20] <- "Mural"
annotationTab.hpc$cellType[21] <- "drop.doublet"
annotationTab.hpc$cellType[c(22)] <- "OPC_COP"


sce.hpc$cellType <- annotationTab.hpc$cellType[match(sce.hpc$collapsedCluster,
                                                         annotationTab.hpc$collapsedCluster)]
sce.hpc$cellType <- factor(sce.hpc$cellType)




## Re-print marker expression with cell type labels ===
# load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_sACC-n5_cleaned-combined_SCE_MNT2021.rda",
#      verbose=T)
cell_colors.hpc <- cluster_colors[order(as.integer(names(cluster_colors)))]
names(cell_colors.hpc) <- annotationTab.hpc$cellType
cell_colors.hpc
    #         Oligo       Astro_A         Micro           OPC       Excit_A       Astro_B 
    #     "#729ECE"     "#FF9E4A"     "#67BF5C"     "#ED665D"     "#AD8BC9"     "#A8786E" 
    #       Inhib_A       Excit_B       Excit_C       Excit_D       Excit_E       Excit_F 
    #     "#ED97CA"     "#A2A2A2"     "#CDCC5D"     "#6DCCDA"     "#1F77B4"     "#AEC7E8" 
    #       Inhib_B       Excit_G       Excit_H       Inhib_C       Inhib_D         Tcell 
    #     "#FF7F0E"     "#FFBB78"     "#2CA02C"     "#98DF8A"     "#D62728"     "#FF9896" 
    # drop.lowNTx_A         Mural  drop.doublet       OPC_COP drop.lowNTx_B 
    #     "#9467BD"     "#C5B0D5"     "#8C564B"     "#C49C94"     "#E377C2"


## Save
save(sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, ref.sampleInfo, clusterRefTab.hpc, annotationTab.hpc, cell_colors.hpc,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_HPC-n3_cleaned-combined_SCE_MNT2021.rda")


pdf("pdfs/revision/regionSpecific_HPC-n3_marker-logExprs_collapsedClusters_MNT2021.pdf", height=6, width=10)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpressionCustom(sce = sce.hpc,
                         features = markers.mathys.custom[[i]], 
                         features_name = names(markers.mathys.custom)[[i]], 
                         anno_name = "cellType") +
      scale_color_manual(values = cell_colors.hpc)
  )
}
dev.off()

    # Optionally, drop doublet cluster & 'drop.lowNTx' and re-print
    sce.hpc <- sce.hpc[ ,-grep("drop.",sce.hpc$cellType)]
    sce.hpc$cellType <- droplevels(sce.hpc$cellType)

# Final comment - the Micro cluster--its prelim clusters were checked for endothelial
#                 markers, but doesn't seem they were merged together.  Perhaps this dataset
#                 just didn't capture any, as with sACC (unlike for the AMY)

    
## Re-print reducedDims with these annotations ===
pdf("pdfs/revision/regionSpecific_HPC-n3_reducedDims-with-collapsedClusters_MNT2021.pdf")
plotReducedDim(sce.hpc, dimred="PCA_corrected", ncomponents=5, colour_by="cellType", point_alpha=0.5)
plotTSNE(sce.hpc, colour_by="sampleID", point_alpha=0.5)
plotTSNE(sce.hpc, colour_by="protocol", point_alpha=0.5)
plotTSNE(sce.hpc, colour_by="prelimCluster", text_by="prelimCluster",
         text_size=3, point_alpha=0.5)
plotTSNE(sce.hpc, colour_by="cellType", text_by="cellType",
         text_size=3, point_alpha=0.5)
plotTSNE(sce.hpc, colour_by="sum", point_alpha=0.5)
plotTSNE(sce.hpc, colour_by="doubletScore", point_alpha=0.5)
# And some more informative UMAPs
plotUMAP(sce.hpc, colour_by="prelimCluster", text_by="prelimCluster",
         text_size=3, point_alpha=0.5)
plotUMAP(sce.hpc, colour_by="cellType", text_by="cellType",
         text_size=3, point_alpha=0.5)
dev.off()
    
    
      ## -> proceed to 'step03_markerDetxn-analyses[...].R'



## And finally, for reference:
table(sce.hpc$cellType.split, sce.hpc$sample)
    #                br5161.hpc br5212.hpc br5287.hpc
    # Astro_A              424        375        137
    # Astro_B               83        125         26
    # drop.doublet           4          1          0
    # drop.lowNTx_A         42         54          9
    # drop.lowNTx_B          9          5          5
    # Excit_A              155        314         17
    # Excit_B                4          9         74
    # Excit_C                1          0          5
    # Excit_D                2          1         32
    # Excit_E                4          2          0
    # Excit_F               35          2          1
    # Inhib                170         87         74
    # Micro                487        481        193
    # Mural                 20         19          4
    # Oligo               2586       2235       1091
    # OPC                  374        255        194
    # OPC_COP                7          3          5
    # Tcell                 14          9          3


### Session info for 07May2021 ==========================================================
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
#   [1] parallel  stats4    stats     graphics  grDevices datasets 
# [7] utils     methods   base     
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
# [3] colorspace_2.0-0          ellipsis_0.3.2           
# [5] scuttle_1.0.4             bluster_1.0.0            
# [7] XVector_0.30.0            BiocNeighbors_1.8.2      
# [9] rstudioapi_0.13           farver_2.1.0             
# [11] bit64_4.0.5               fansi_0.4.2              
# [13] xml2_1.3.2                splines_4.0.4            
# [15] R.methodsS3_1.8.1         sparseMatrixStats_1.2.1  
# [17] cachem_1.0.4              Rsamtools_2.6.0          
# [19] ResidualMatrix_1.0.0      dbplyr_2.1.1             
# [21] R.oo_1.24.0               HDF5Array_1.18.1         
# [23] compiler_4.0.4            httr_1.4.2               
# [25] dqrng_0.2.1               assertthat_0.2.1         
# [27] Matrix_1.3-2              fastmap_1.1.0            
# [29] lazyeval_0.2.2            limma_3.46.0             
# [31] BiocSingular_1.6.0        prettyunits_1.1.1        
# [33] tools_4.0.4               rsvd_1.0.3               
# [35] igraph_1.2.6              gtable_0.3.0             
# [37] glue_1.4.2                GenomeInfoDbData_1.2.4   
# [39] dplyr_1.0.5               rappdirs_0.3.3           
# [41] Rcpp_1.0.6                vctrs_0.3.6              
# [43] Biostrings_2.58.0         rhdf5filters_1.2.0       
# [45] rtracklayer_1.50.0        DelayedMatrixStats_1.12.3
# [47] stringr_1.4.0             beachmat_2.6.4           
# [49] lifecycle_1.0.0           irlba_2.3.3              
# [51] statmod_1.4.35            XML_3.99-0.6             
# [53] edgeR_3.32.1              zlibbioc_1.36.0          
# [55] scales_1.1.1              hms_1.0.0                
# [57] ProtGenerics_1.22.0       rhdf5_2.34.0             
# [59] RColorBrewer_1.1-2        curl_4.3                 
# [61] memoise_2.0.0             gridExtra_2.3            
# [63] segmented_1.3-3           biomaRt_2.46.3           
# [65] stringi_1.5.3             RSQLite_2.2.7            
# [67] BiocParallel_1.24.1       rlang_0.4.10             
# [69] pkgconfig_2.0.3           bitops_1.0-7             
# [71] lattice_0.20-41           purrr_0.3.4              
# [73] Rhdf5lib_1.12.1           labeling_0.4.2           
# [75] GenomicAlignments_1.26.0  cowplot_1.1.1            
# [77] bit_4.0.4                 tidyselect_1.1.1         
# [79] magrittr_2.0.1            R6_2.5.0                 
# [81] generics_0.1.0            DelayedArray_0.16.3      
# [83] DBI_1.1.1                 pillar_1.6.0             
# [85] withr_2.4.2               RCurl_1.98-1.3           
# [87] tibble_3.1.1              crayon_1.4.1             
# [89] utf8_1.2.1                BiocFileCache_1.14.0     
# [91] viridis_0.6.0             progress_1.2.2           
# [93] locfit_1.5-9.4            grid_4.0.4               
# [95] blob_1.2.1                digest_0.6.27            
# [97] R.utils_2.10.1            openssl_1.4.3            
# [99] munsell_0.5.0             beeswarm_0.3.1           
# [101] viridisLite_0.4.0         vipor_0.4.5              
# [103] askpass_1.1


