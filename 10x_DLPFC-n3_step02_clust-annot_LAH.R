
### MNT 10x snRNA-seq workflow: step 02
###   **Region-specific analyses**
###     - (2x) DLPFC samples from: Br5161 & Br5212
### Initiated MNT 12Feb2020
### LAH 27Apr2021: add expansion samples (n=3)
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
library(purrr)

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

# Load 2021 expansion set
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/all-FACS-n10_2021rev_SCEs_processing-QC_MNTMar2021.rda", verbose=T)
    # pilot.data.2, pilot.data.2.unfiltered, e.out.2, ref.sampleInfo.rev
    rm(pilot.data.2.unfiltered, e.out.2)

## filter for dlpfc samples
pilot.data <- pilot.data[grep("dlpfc", names(pilot.data))]
pilot.data.2 <- pilot.data.2[grep("dlpfc", names(pilot.data.2))]

map_int(c(pilot.data, pilot.data.2), ncol)
# br5161.dlpfc br5212.dlpfc br5207.dlpfc 
# 4215         1693         5294 

### MNT comment: At this point, each sample (which is a SCE object in the list, 'pilot.data') has been
  #              QC'd for cell/nucleus calling ('emptyDrops()' test) and mito rate thresholding
  #   Additionally, there have been a computed 'doubletScore', which will QC with, after
  #   clustering (e.g. that there are no doublet-driven clusters, etc.)

### Merging shared-region samples ============================================
  # Newest iterations for normalization: multiBatchNorm-alize

# Order as will do for `fastMNN()` (this shouldn't matter here)
sce.dlpfc <- cbind(pilot.data.2[["br5207.dlpfc"]], 
                   pilot.data[["br5161.dlpfc"]],
                   pilot.data[["br5212.dlpfc"]]
)

sce.dlpfc
# class: SingleCellExperiment 
# dim: 33538 11202 
# metadata(3): Samples Samples Samples
# assays(1): counts
# rownames(33538): MIR1302-2HG FAM138A ... AC213203.1 FAM231C
# rowData names(6): gene_id gene_version ... gene_biotype Symbol.uniq
# colnames(11202): AAACCCACAGTCGTTA-1 AAACCCAGTGCATACT-1 ... TTTGTTGGTGACCGAA-1 TTTGTTGTCTCGAACA-1
# colData names(16): Sample Barcode ... protocol sequencer
# reducedDimNames(0):
#   altExpNames(0):

table(sce.dlpfc$sampleID)
# br5161.dlpfc br5207.dlpfc br5212.dlpfc 
# 4215         5294         1693 

# Use `multiBatchNorm()` to compute log-normalized counts, matching the scaling across samples
sce.dlpfc <- multiBatchNorm(sce.dlpfc, batch=sce.dlpfc$sampleID)

# Use the simple `modelGeneVar` - this makes more sense over `combineVar`, since the
#   cell composition is already known to be quite different (with NeuN selection)
geneVar <- modelGeneVar(sce.dlpfc)
chosen.hvgs.dlpfc <- geneVar$bio > 0
sum(chosen.hvgs.dlpfc)
# [1] 11317


### Dimensionality reduction ================================================================

# Run `fastMNN` (internally uses `multiBatchPCA`), taking top 100 (instead of default 50 PCs)
set.seed(109)
mnn.hold <-  fastMNN(sce.dlpfc, batch=sce.dlpfc$sampleID,
                     merge.order=c("br5207.dlpfc","br5161.dlpfc","br5212.dlpfc"),
                     subset.row=chosen.hvgs.dlpfc, d=100,
                     correct.all=TRUE, get.variance=TRUE,
                     BSPARAM=BiocSingular::IrlbaParam())
    # This temp file just used for getting batch-corrected components (drops a variety of entries)

table(colnames(mnn.hold) == colnames(sce.dlpfc))  # all TRUE
table(mnn.hold$batch == sce.dlpfc$sampleID) # all TRUE

# Add them to the SCE, as well as the metadata (though the latter might not be so usefl)
reducedDim(sce.dlpfc, "PCA_corrected") <- reducedDim(mnn.hold, "corrected") # 100 components
metadata(sce.dlpfc) <- metadata(mnn.hold)

# Save into a new region-specific SCE object/flie
save(sce.dlpfc, chosen.hvgs.dlpfc, ref.sampleInfo, ref.sampleInfo.rev,
     file="rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda")


    ## 'getClusteredPCs()' evaluated in qsub mode (with 'R-batchJob_Amyg-n2_optimalPCselxn_MNTFeb2020.R')
    #    --> saved into same .rda


### Picking up with optimally-defined PC space ===
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda",
     verbose=TRUE)
    # sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, ref.sampleInfo, ref.sampleInfo.rev

# How many PCs is optimal?:
metadata(pc.choice.dlpfc)$chosen
    # [1] 99

## Assign this chosen (99 PCs) to 'PCA_opt'
reducedDim(sce.dlpfc, "PCA_opt") <- reducedDim(sce.dlpfc, "PCA_corrected")[ ,1:(metadata(pc.choice.dlpfc)$chosen)]


## t-SNE
set.seed(109)
sce.dlpfc <- runTSNE(sce.dlpfc, dimred="PCA_opt")


## UMAP
set.seed(109)
sce.dlpfc <- runUMAP(sce.dlpfc, dimred="PCA_opt")


# How do these look?
plotReducedDim(sce.dlpfc, dimred="TSNE", colour_by="sampleID")
plotReducedDim(sce.dlpfc, dimred="UMAP", colour_by="sampleID")


### Clustering: Two-step ======================================================
### Step 1: Perform graph-based clustering in this optimal PC space
#         - take k=20 NN to build graph
snn.gr <- buildSNNGraph(sce.dlpfc, k=20, use.dimred="PCA_opt")
clusters.k20 <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters.k20)
    ## 107 prelim clusters

# Assign as 'prelimCluster'
sce.dlpfc$prelimCluster <- factor(clusters.k20)
plotReducedDim(sce.dlpfc, dimred="TSNE", colour_by="prelimCluster")

# Is sample driving this 'high-res' clustering at this level?
(sample_prelimClusters <- table(sce.dlpfc$prelimCluster, sce.dlpfc$sampleID))  # (a little bit, but is typical)
sample_prelimClusters[which(rowSums(sample_prelimClusters == 0) == 2),]
# 39 - only 4 samples all from Br5207

# rbind the ref.sampleInfo[.rev]
ref.sampleInfo <- rbind(ref.sampleInfo, ref.sampleInfo.rev)

## check doublet score for each prelim clust
clusIndexes = splitit(sce.dlpfc$prelimCluster)
prelimCluster.medianDoublet <- sapply(clusIndexes, function(ii){
  median(sce.dlpfc$doubletScore[ii])
}
)

summary(prelimCluster.medianDoublet)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.01059  0.07083  0.14823  0.53264  0.30064 14.79144 

hist(prelimCluster.medianDoublet)

## watch in clustering
prelimCluster.medianDoublet[prelimCluster.medianDoublet > 5]
# 19       32       73 
# 14.79144  7.98099 10.20462 

table(sce.dlpfc$prelimCluster)[c(19, 32, 73)]
# 19 32 73 
# 27 32  8 

# Save for now
save(sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, ref.sampleInfo,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda")


### Step 2: Hierarchical clustering of pseudo-bulked ("PB'd") counts with most robust normalization
  #         (as determined in: 'side-Rscript_testingStep2_HC-normalizn-approaches_wAmygData_MNTJan2020.R')
  #           ** That is, to pseudo-bulk (aka 'cluster-bulk') on raw counts, on all [non-zero] genes,
  #              normalize with `librarySizeFactors()`, log2-transform, then perform HC'ing

# Preliminary cluster index for pseudo-bulking
clusIndexes = splitit(sce.dlpfc$prelimCluster)
prelimCluster.PBcounts <- sapply(clusIndexes, function(ii){
  rowSums(assays(sce.dlpfc)$counts[ ,ii])
  }
)
    
    # And btw...
    table(rowSums(prelimCluster.PBcounts)==0)
    # FALSE  TRUE 
    # 29310  4228

# Compute LSFs at this level
sizeFactors.PB.all  <- librarySizeFactors(prelimCluster.PBcounts)

# Normalize with these LSFs
geneExprs.temp <- t(apply(prelimCluster.PBcounts, 1, function(x) {log2(x/sizeFactors.PB.all + 1)}))

## Perform hierarchical clustering
dist.clusCollapsed <- dist(t(geneExprs.temp))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

dend <- as.dendrogram(tree.clusCollapsed, hang=0.2)

# labels(dend)[grep(c("19|32|73"),labels(dend))] <- paste0(labels(dend)[grep(c("19|32|73"),labels(dend))], "*")

# Just for observation
par(cex=.6)
myplclust(tree.clusCollapsed, cex.main=2, cex.lab=1.5, cex=1.8)

dend %>% 
  set("labels_cex", 0.8) %>%
  plot(horiz = TRUE)
abline(v = 325, lty = 2)

clust.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                               minClusterSize=2, deepSplit=1, cutHeight=325)

table(clust.treeCut)
unname(clust.treeCut[order.dendrogram(dend)])
    ## Cutting at 250 looks good for the main neuronal branch, but a lot of glial
     #    prelim clusters are dropped off (0's)

    # Cut at 400 for broad glia branch (will manually merge remaining dropped off)    
    glia.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                                  minClusterSize=2, deepSplit=1, cutHeight=400)
    unname(glia.treeCut[order.dendrogram(dend)])
    
    # Take those and re-assign to the first assignments
    
# clust <- clust.treeCut[order.dendrogram(dend)]
# clust2 <- name_zeros(clust, list(c(1,2), c(106,107)))
# unname(clust2)

# Add new labels to those prelimClusters cut off
clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)]==0)] <- max(clust.treeCut)+c(1, 1, 2, 3, 4, 5, 5)

# 'Re-write', since there are missing numbers
# clust.treeCut[order.dendrogram(dend)] <- as.numeric(as.factor(clust2))
clust.treeCut[order.dendrogram(dend)] <- as.numeric(as.factor(clust.treeCut[order.dendrogram(dend)]))

## Define color pallet
cluster_colors <- unique(tableau20[clust.treeCut[order.dendrogram(dend)]])
names(cluster_colors) <- unique(clust.treeCut[order.dendrogram(dend)])
labels_colors(dend) <- cluster_colors[clust.treeCut[order.dendrogram(dend)]]

# Print for future reference
pdf("pdfs/revision/regionSpecific_DLPFC-n3_HC-prelimCluster-relationships_LAH2021.pdf", height = 9)
par(cex=0.6, font=2)
plot(dend, main="3x DLPFC prelim-kNN-cluster relationships with collapsed assignments", horiz = TRUE)
abline(v = 325, lty = 2)
dev.off()


# Make reference for new cluster assignment
clusterRefTab.dlpfc <- data.frame(origClust=order.dendrogram(dend),
                                merged=clust.treeCut[order.dendrogram(dend)])


# Assign as 'collapsedCluster'
sce.dlpfc$collapsedCluster <- factor(clusterRefTab.dlpfc$merged[match(sce.dlpfc$prelimCluster, clusterRefTab.dlpfc$origClust)])
n_clusters <- length(levels(sce.dlpfc$collapsedCluster))
# Print some visualizations:
pdf("pdfs/revision/regionSpecific_DLPFC-n3_reducedDims-with-collapsedClusters_LAH2021.pdf")
plotReducedDim(sce.dlpfc, dimred="PCA_corrected", ncomponents=5, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="sampleID", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="protocol", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="sum", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="doubletScore", point_alpha=0.5)
# And some more informative UMAPs
plotUMAP(sce.dlpfc, colour_by="sampleID", point_alpha=0.5)
plotUMAP(sce.dlpfc, colour_by="collapsedCluster", point_alpha=0.5)
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
  'endothelial' = c('CLDN5', 'FLT1', 'VTN', 'PECAM1')
)


pdf("pdfs/revision/regionSpecific_DLPFC-n3_marker-logExprs_collapsedClusters_LAH2021.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpressionCustom(sce = sce.dlpfc,
                         features = markers.mathys.custom[[i]], 
                         features_name = names(markers.mathys.custom)[[i]], 
                         anno_name = "collapsedCluster")
  )
}
dev.off()


## Add annotations, looking at marker gene expression
annotationTab.dlpfc <- data.frame(collapsedCluster=c(1:n_clusters))
annotationTab.dlpfc$cellType <- NA
annotationTab.dlpfc$cellType[c(1:3, 6, 17,18)] <- paste0("Inhib_", c("A","B","C","D","E","F"))
annotationTab.dlpfc$cellType[c(4,5,8,10:13)] <- paste0("Excit_", c("A","B","C","D","E","F","G"))
annotationTab.dlpfc$cellType[c(7, 9, 13, 16)] <- c("Astro", "Oligo", "OPC", "Micro")
annotationTab.dlpfc$cellType[c(14,15)] <- paste0("ambig.glial_", c("A","B"))


sce.dlpfc$cellType <- annotationTab.dlpfc$cellType[match(sce.dlpfc$collapsedCluster,
                                                         annotationTab.dlpfc$collapsedCluster)]
sce.dlpfc$cellType <- factor(sce.dlpfc$cellType)

## QC - How do the total UMI distribution look?
# newClusIndex <- splitit(sce.dlpfc$collapsedCluster)
newClusIndex <- splitit(sce.dlpfc$cellType)
sapply(newClusIndex, function(x) {quantile(sce.dlpfc$sum[x])})
#          1        2     3     4     5     6        7         8       9       10       11    12       13   14     15
# 0%    3045  1708.00  1803  1201  1838  1749   884.00   2295.00   850.0  9100.00   3029.0  1099  2084.00 1774 1979.0
# 25%  15330 11474.25 18791 28468 19205 23015  4008.75  39673.25  4940.5 26131.75  37044.5 27809  7371.75 2332 3210.5
# 50%  18968 16423.50 24138 34997 25241 29623  5737.00  49345.00  6385.0 32065.50  49414.0 37430  9112.00 2762 4692.5
# 75%  23036 22355.75 28996 43590 33292 35364  7953.00  59114.25  7986.0 39545.00  61162.0 47099 11052.75 4379 5218.0
# 100% 55574 66556.00 81134 87392 69309 67503 26618.00 115449.00 25379.0 69202.00 101625.0 83826 23492.00 6919 8043.0
#            16      17       18
# 0%     879.00  7966.0  8833.00
# 25%   3019.25  9056.5 19428.75
# 50%   3883.50 14241.0 20281.50
# 75%   4911.75 17472.0 22677.00
# 100% 11137.00 42588.0 27504.00

sapply(newClusIndex, function(x) {quantile(sce.dlpfc[,x]$doubletScore)})


table(sce.dlpfc$collapsedCluster)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18 
# 333  454  365  529  773  413  782  524 5455  132  187  243  572   19   18  388    7    8 


# Save
save(sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc, ref.sampleInfo, annotationTab.dlpfc,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda")


## Re-print marker expression plots with annotated cluster names ===
cell_colors <- cluster_colors[order(as.integer(names(cluster_colors)))]
names(cell_colors) <- annotationTab.dlpfc$cellType
cell_colors
# Inhib_A       Inhib_B       Inhib_C       Excit_A       Excit_B       Inhib_D         Astro       Excit_C 
# "#1F77B4"     "#AEC7E8"     "#FF7F0E"     "#FFBB78"     "#2CA02C"     "#98DF8A"     "#D62728"     "#FF9896" 
# Oligo       Excit_D       Excit_E       Excit_F           OPC ambig.glial_A ambig.glial_B         Micro 
# "#9467BD"     "#C5B0D5"     "#8C564B"     "#C49C94"     "#E377C2"     "#F7B6D2"     "#7F7F7F"     "#C7C7C7" 
# Inhib_E       Inhib_F 
# "#BCBD22"     "#DBDB8D" 

pdf("pdfs/revision/regionSpecific_DLPFC-n3_marker-logExprs_collapsedClusters_LAH2021.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpressionCustom(sce = sce.dlpfc,
                         features = markers.mathys.custom[[i]], 
                         features_name = names(markers.mathys.custom)[[i]], 
                         anno_name = "cellType") +
      scale_color_manual(values = cell_colors)
  )
}
dev.off()

    ## Optionally, for a cleaner version, drop those 'drop.lowNTx_'s and re-print
    sce.dlpfc <- sce.dlpfc[ ,-grep("drop.", sce.dlpfc$cellType.prelim)]
    sce.dlpfc$cellType.prelim <- droplevels(sce.dlpfc$cellType.prelim)

      ## -> proceed to 'step03_markerDetxn-analyses[...].R'



### For reference === == === == ===
table(sce.dlpfc$cellType, sce.dlpfc$sampleID)
#               br5161.dlpfc br5207.dlpfc br5212.dlpfc
# ambig.glial_A            3           11            5
# ambig.glial_B            3           13            2
# Astro                  371          274          137
# Excit_A                111          298          120
# Excit_B                 75          544          154
# Excit_C                 44          325          155
# Excit_D                 22           83           27
# Excit_E                 77           85           25
# Excit_F                102          105           36
# Inhib_A                 39          205           89
# Inhib_B                 98          250          106
# Inhib_C                 47          262           56
# Inhib_D                119          216           78
# Inhib_E                  2            3            2
# Inhib_F                  0            7            1
# Micro                  152          144           92
# Oligo                 2754         2184          517
# OPC                    196          285           91

table(sce.dlpfc$prelimCluster, sce.dlpfc$sampleID)

table(sce.dlpfc$cellType, sce.dlpfc$collapsedCluster)
#                  1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18
# ambig.glial_A    0    0    0    0    0    0    0    0    0    0    0    0    0   19    0    0    0    0
# ambig.glial_B    0    0    0    0    0    0    0    0    0    0    0    0    0    0   18    0    0    0
# Astro            0    0    0    0    0    0  782    0    0    0    0    0    0    0    0    0    0    0
# Excit_A          0    0    0  529    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# Excit_B          0    0    0    0  773    0    0    0    0    0    0    0    0    0    0    0    0    0
# Excit_C          0    0    0    0    0    0    0  524    0    0    0    0    0    0    0    0    0    0
# Excit_D          0    0    0    0    0    0    0    0    0  132    0    0    0    0    0    0    0    0
# Excit_E          0    0    0    0    0    0    0    0    0    0  187    0    0    0    0    0    0    0
# Excit_F          0    0    0    0    0    0    0    0    0    0    0  243    0    0    0    0    0    0
# Inhib_A        333    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# Inhib_B          0  454    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# Inhib_C          0    0  365    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# Inhib_D          0    0    0    0    0  413    0    0    0    0    0    0    0    0    0    0    0    0
# Inhib_E          0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    7    0
# Inhib_F          0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    8
# Micro            0    0    0    0    0    0    0    0    0    0    0    0    0    0    0  388    0    0
# Oligo            0    0    0    0    0    0    0    0 5455    0    0    0    0    0    0    0    0    0
# OPC              0    0    0    0    0    0    0    0    0    0    0    0  572    0    0    0    0    0



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

