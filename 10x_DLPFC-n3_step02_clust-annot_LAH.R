
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

## Define DLPFC pallet
# labels_colors(dend) <- tableau20[clust2]
dlpfc_pallet <- unique(tableau20[clust.treeCut[order.dendrogram(dend)]])
names(dlpfc_pallet) <- unique(clust.treeCut[order.dendrogram(dend)])

# labels_colors(dend) <- tableau20[clust.treeCut[order.dendrogram(dend)]]
labels_colors(dend) <- dlpfc_pallet[clust.treeCut[order.dendrogram(dend)]]

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
    plotExpression(sce.dlpfc, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]), ncol=2,
                   x="collapsedCluster", colour_by="collapsedCluster", point_alpha=0.4, point_size=.7,
                   add_legend=F) + stat_summary(fun = median, fun.min = median, fun.max = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:n_clusters], length(markers.mathys.custom[[i]]))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
    ggtitle(label=paste0(names(markers.mathys.custom)[i], " markers"))
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
dlpfc_pallet <- dlpfc_pallet[order(as.integer(names(dlpfc_pallet)))]
names(dlpfc_pallet) <- annotationTab.dlpfc$cellType

pdf("pdfs/revision/regionSpecific_DLPFC-n3_marker-logExprs_collapsedClusters_LAH2021.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.dlpfc, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="cellType.prelim", colour_by="cellType.prelim", point_alpha=0.2, point_size=.7,
                   add_legend=F) +
      stat_summary(fun = median, fun.min = median, fun.max = median,
                   geom = "crossbar", width = 0.3,
                   #colour=rep(tableau20[1:17], length(markers.mathys.custom[[i]]))) +
                   # colour=rep(tableau20[1:n_clusters], length(markers.mathys.custom[[i]]))
      )+
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(markers.mathys.custom)[i], " markers"))+
      scale_color_manual(values = dlpfc_pallet)
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



### MNT 06May2020: subcluster-level annotations ========================================
# Some motivations:
#   i) to see if the high-VCAN-neuronal cluster exists in prelimCluster
#      as seen in the pan-brain level
#   ii) To now formally define subcluster-level populations

load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda",
     verbose=T)
    # sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc, ref.sampleInfo

# Look at some marker expression at the prelimCluster level
pdf("pdfs/revision/ztemp_amyg-prelimCluster-neuronalMarkerExpression.pdf", width=8, height=7)
plotExpression(sce.dlpfc, exprs_values="logcounts", features=c("SNAP25","GAD1","GAD2","SLC17A6","SLC17A7","VCAN"),
               x="prelimCluster", colour_by="prelimCluster", ncol=2)
dev.off()


## More-manual annotations for neuronal subpops:
clusterRefTab.dlpfc$cellType <- sce.dlpfc$cellType[match(clusterRefTab.dlpfc$merged, sce.dlpfc$collapsedCluster)]
clusterRefTab.dlpfc$cellType <- as.character(clusterRefTab.dlpfc$cellType)

# Make new column for subclusers
clusterRefTab.dlpfc$manual <- clusterRefTab.dlpfc$cellType

    ## Excit subclusters: ====
    clusterRefTab.dlpfc$manual <- ifelse(clusterRefTab.dlpfc$origClust %in% c(3,9,16,20),
                                       paste0(clusterRefTab.dlpfc$cellType, ".1"),
                                       as.character(clusterRefTab.dlpfc$manual))
    
    # 10 and 12 would cut off at different heights; 12 expresses SLC17A7 and 10 doesn't really any:
    clusterRefTab.dlpfc$manual <- ifelse(clusterRefTab.dlpfc$origClust %in% c(10),
                                       paste0(clusterRefTab.dlpfc$cellType, ".2"),
                                       as.character(clusterRefTab.dlpfc$manual))
    
    clusterRefTab.dlpfc$manual <- ifelse(clusterRefTab.dlpfc$origClust %in% c(12),
                                       paste0(clusterRefTab.dlpfc$cellType, ".3"),
                                       as.character(clusterRefTab.dlpfc$manual))
    
    ## Inhib subclusters: merge 2/6 and 7/14 pairs, then split the rest
    clusterRefTab.dlpfc$manual <- ifelse(clusterRefTab.dlpfc$origClust %in% c(2,6),
                                       paste0(clusterRefTab.dlpfc$cellType, ".1"),
                                       as.character(clusterRefTab.dlpfc$manual))
    
    clusterRefTab.dlpfc$manual <- ifelse(clusterRefTab.dlpfc$origClust %in% c(7,14),
                                       paste0(clusterRefTab.dlpfc$cellType, ".2"),
                                       as.character(clusterRefTab.dlpfc$manual))
    
    clusterRefTab.dlpfc$manual <- ifelse(clusterRefTab.dlpfc$origClust %in% c(23),
                                       paste0(clusterRefTab.dlpfc$cellType, ".3"),
                                       as.character(clusterRefTab.dlpfc$manual))
    
    clusterRefTab.dlpfc$manual <- ifelse(clusterRefTab.dlpfc$origClust %in% c(18),
                                       paste0(clusterRefTab.dlpfc$cellType, ".4"),
                                       as.character(clusterRefTab.dlpfc$manual))
    
    clusterRefTab.dlpfc$manual <- ifelse(clusterRefTab.dlpfc$origClust %in% c(21),
                                       paste0(clusterRefTab.dlpfc$cellType, ".5"),
                                       as.character(clusterRefTab.dlpfc$manual))
    
    ## All other glial types will be kept the same
    
        ## Post-hoc: Looks like Excit.2 truly inhibitory, and Inhib.5 truly [or more so] excitatory
        #           (the latter is the ~39-40 high-VCAN nuclei ID'd in this sample at pan-brain) 
        #  -> swap them
        
        clusterRefTab.dlpfc$manual <- ifelse(clusterRefTab.dlpfc$origClust == 10,  # ("Excit.2")
                                           "Inhib.5",
                                           as.character(clusterRefTab.dlpfc$manual))
        clusterRefTab.dlpfc$manual <- ifelse(clusterRefTab.dlpfc$origClust == 21,  # ("Inhib.5")
                                           "Excit.2",
                                           as.character(clusterRefTab.dlpfc$manual))
    
        # --> THEN re-run the below
        
    ## end cluster splitting chunk ====

clusterRefTab.dlpfc
    
## Add new annotations
sce.dlpfc$cellType.split <- clusterRefTab.dlpfc$manual[match(sce.dlpfc$prelimCluster,
                                                         clusterRefTab.dlpfc$origClust)]
sce.dlpfc$cellType.split <- factor(sce.dlpfc$cellType.split)

table(sce.dlpfc$cellType.split, sce.dlpfc$cellType)
    # good

table(sce.dlpfc$cellType.split) # (printing post-hoc-corrected annotations)
    #Ambig.lowNtrxts           Astro         Excit.1         Excit.2         Excit.3
    #             50             852             334              40              55
    #        Inhib.1         Inhib.2         Inhib.3         Inhib.4         Inhib.5
    #            171             109              35              24              98
    #          Micro           Oligo             OPC
    #            764            3473             627

## Save these
save(sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc, ref.sampleInfo,
     file="rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda")


## Also print expression at this level of partitioning ===

# First remove "Ambig.lowNtrxts" (50 nuclei):
sce.dlpfc <- sce.dlpfc[ ,sce.dlpfc$cellType.split != "Ambig.lowNtrxts"]
sce.dlpfc$cellType.split <- droplevels(sce.dlpfc$cellType.split)

pdf("pdfs/revision/regionSpecific_Amyg-n2_marker-logExprs_cellTypesSplit_May2020.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.dlpfc, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7,
                   add_legend=F) + stat_summary(fun = median, fun.min = median, fun.max = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:12], length(markers.mathys.custom[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(markers.mathys.custom)[i], " markers"))
  )
}
dev.off()


## Let's also re-plot reducedDims with new [broad & split] cell type annotations
#        (and rename old file with prefix 'zold_')
pdf("pdfs/revision/regionSpecific_Amyg-n2_reducedDims-with-collapsedClusters_May2020.pdf")
plotReducedDim(sce.dlpfc, dimred="PCA", ncomponents=5, colour_by="cellType", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="sample", point_size=3.5, point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="prelimCluster", point_size=3.5, point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="cellType", point_size=3.5, point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="cellType.split", point_size=3.5, point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="sum", point_size=3.5, point_alpha=0.5)
plotUMAP(sce.dlpfc, colour_by="cellType", point_size=3.5, point_alpha=0.5)
plotUMAP(sce.dlpfc, colour_by="cellType.split", point_size=3.5, point_alpha=0.5)
dev.off()


## And finally, for reference:
table(sce.dlpfc$cellType.split, sce.dlpfc$sample)
    #          amy.5161 amy.5212
    # Astro        489      363
    # Excit.1      141      193
    # Excit.2        0       40
    # Excit.3        0       55
    # Inhib.1       16      155
    # Inhib.2       33       76
    # Inhib.3       11       24
    # Inhib.4       24        0
    # Inhib.5       85       13
    # Micro        425      339
    # Oligo       1697     1776
    # OPC          335      292


## Some diggings-into for paper =======
# BLA or central amygdala-specific (Cartpt in all subregions... in mouse at least)
plotExpression(sce.dlpfc, exprs_values = "logcounts", features=toupper(c("Cyp26b1", "Bmp3", "Cartpt")),
               x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7,
               add_legend=F) + stat_summary(fun = median, fun.min = median, fun.max = median,
                                            geom = "crossbar", width = 0.3,
                                            colour=rep(tableau20[1:12], 3)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

# Split by sample
sce.dlpfc$sample.cellType.split <- paste0(ss(sce.dlpfc$sample,"\\.",2), "_", sce.dlpfc$cellType.split)
plotExpression(sce.dlpfc, exprs_values = "logcounts", features=toupper(c("Cyp26b1", "Bmp3")),
               x="sample.cellType.split", colour_by="sample.cellType.split", point_alpha=0.5, point_size=.7,
               add_legend=F, show_median=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
    ## Looks like some BLA

# Actually faceting is nicer and maintains colors but have to print genes, separately
# From bulk-RNA-seq (AnJa shared): Look at PART1 (BLA >) & GABRQ/SYTL5 (MeA >)
pdf("pdfs/revision/exploration/zGenesByDonor_PART1-GABRQ-SYTL5_AMY-subclusters_MNT.pdf", height=3, width=7)
# BLA-enriched
plotExpression(sce.dlpfc, exprs_values = "logcounts", features="PART1",
               x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7,
               add_legend=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7)) + facet_grid(~ sce.dlpfc$donor) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3,
               colour=rep(tableau20[1:12][c(1:2,5:12, 1:7,9:12)])) 
# MeA-enriched
plotExpression(sce.dlpfc, exprs_values = "logcounts", features="GABRQ",
               x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7,
               add_legend=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7)) + facet_grid(~ sce.dlpfc$donor) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3,
               colour=rep(tableau20[1:12][c(1:2,5:12, 1:7,9:12)])) 
plotExpression(sce.dlpfc, exprs_values = "logcounts", features="SYTL5",
               x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7,
               add_legend=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7)) + facet_grid(~ sce.dlpfc$donor) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3,
               colour=rep(tableau20[1:12][c(1:2,5:12, 1:7,9:12)])) 
dev.off()






## Added MNT 24May2020: tSNE in lower dims =======================================================
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc, ref.sampleInfo


# How many PCs?
head(attr(reducedDim(sce.dlpfc, "PCA"), "percentVar"), n=50)
    # [1] 20.56604550 10.39668224  6.66923256  3.20654266  1.88015247  0.90499720
    # [7]  0.63195463  0.30153787  0.28714069  0.27722945  0.25484711  0.23154164
    # [13]  0.21200112  0.19865766  0.19058467  0.17696622  0.15895906  0.14111373
    # [19]  0.12953405  0.12661346  0.11187693  0.10874758  0.09528034  0.08876223
    # [25]  0.08503460  0.08069491  0.07432969  0.07332327  0.07270196  0.07085366
    # [31]  0.06895861  0.06840260  0.06628307  0.06566241  0.06499767  0.06456945
    # [37]  0.06311090  0.06122234  0.05977570  0.05910767  0.05887602  0.05824736
    # [43]  0.05790745  0.05736661  0.05646818  0.05604596  0.05573424  0.05512331
    # [49]  0.05473937  0.05442114

# Btw: metadata(pc.choice.dlpfc)$chosen == 45 [PCs]

# 0.1% var or greater
reducedDim(sce.dlpfc, "PCA_22") <- reducedDim(sce.dlpfc, "PCA")[ ,c(1:22)]
# 0.2% var or greater
reducedDim(sce.dlpfc, "PCA_13") <- reducedDim(sce.dlpfc, "PCA")[ ,c(1:13)]
# Top 10 (as Mathys, et al)
reducedDim(sce.dlpfc, "PCA_10") <- reducedDim(sce.dlpfc, "PCA")[ ,c(1:10)]

# First remove this reducedDim bc this has caused trouble previously
reducedDim(sce.dlpfc, "TSNE") <- NULL

# This is interesting:
sapply(c(1:45), function(x){round(cor(reducedDim(sce.dlpfc,"PCA")[ ,x], as.numeric(as.factor(sce.dlpfc$donor))),3)})
    #[1] -0.03913056  0.18884743 -0.08465435  0.08171014  0.70429054  0.03941123
    #[7] -0.07619896  0.02419195  0.08840512 -0.05000060
boxplot(reducedDim(sce.dlpfc,"PCA")[ ,5] ~ sce.dlpfc$donor)
    # Two distributions for sure, though 'close' to one another
# Alternatively
plotReducedDim(sce.dlpfc, dimred="PCA", ncomponents=5, colour_by="donor", point_alpha=0.5)
    # Can't really appreciate as well actually

    ## -> let's get rid of PC5 and try on "PCA_optb"
reducedDim(sce.dlpfc, "PCA_optb") <- reducedDim(sce.dlpfc, "PCA")[ ,c(1:4, 6:(metadata(pc.choice.dlpfc)$chosen))]




# ## 22 PCs tSNE === 
# set.seed(109)
# sce.dlpfc.tsne.22pcs <- runTSNE(sce.dlpfc, dimred="PCA_22")
# 
# ## 13 PCs tSNE ===
# set.seed(109)
# sce.dlpfc.tsne.13pcs <- runTSNE(sce.dlpfc, dimred="PCA_13")
# 
# ## 10 PCs tSNE ===
# set.seed(109)
# sce.dlpfc.tsne.10pcs <- runTSNE(sce.dlpfc, dimred="PCA_10")

## "optimal-b" PCs tSNE ===
set.seed(109)
sce.dlpfc.tsne.optb <- runTSNE(sce.dlpfc, dimred="PCA_optb")
    ## Overall this is the best.  Still very strong donor/batch effects, but this _could_ also be
     #   because they seem to be quite different subdivisions of the amygdala...


# Drop "Ambig.lowNtrxts" cluster as always
# sce.dlpfc.tsne.22pcs <- sce.dlpfc.tsne.22pcs[ ,sce.dlpfc.tsne.22pcs$cellType.split != "Ambig.lowNtrxts"] # 50
# sce.dlpfc.tsne.22pcs$cellType.split <- droplevels(sce.dlpfc.tsne.22pcs$cellType.split)
# 
# sce.dlpfc.tsne.13pcs <- sce.dlpfc.tsne.13pcs[ ,sce.dlpfc.tsne.13pcs$cellType.split != "Ambig.lowNtrxts"] # 50
# sce.dlpfc.tsne.13pcs$cellType.split <- droplevels(sce.dlpfc.tsne.13pcs$cellType.split)
# 
# sce.dlpfc.tsne.10pcs <- sce.dlpfc.tsne.10pcs[ ,sce.dlpfc.tsne.10pcs$cellType.split != "Ambig.lowNtrxts"] # 50
# sce.dlpfc.tsne.10pcs$cellType.split <- droplevels(sce.dlpfc.tsne.10pcs$cellType.split)

sce.dlpfc.tsne.optb <- sce.dlpfc.tsne.optb[ ,sce.dlpfc.tsne.optb$cellType.split != "Ambig.lowNtrxts"] # 50
sce.dlpfc.tsne.optb$cellType.split <- droplevels(sce.dlpfc.tsne.optb$cellType.split)


pdf("pdfs/revision/exploration/zExplore_Amyg-n2_tSNE_22-13-10-optb-PCs_MNTMay2020.pdf", width=8)
# 22 PCs
plotTSNE(sce.dlpfc.tsne.22pcs, colour_by="cellType.split", point_alpha=0.5, point_size=4.0,
         text_size=8, theme_size=18) +
  ggtitle("t-SNE on top 22 PCs (>= 0.1% var)") + theme(plot.title = element_text(size=19))
# 13 PCs
plotTSNE(sce.dlpfc.tsne.13pcs, colour_by="cellType.split", point_alpha=0.5, point_size=4.0,
         text_size=8, theme_size=18) +
  ggtitle("t-SNE on top 13 PCs (>= 0.2% var)") + theme(plot.title = element_text(size=19))
# 10 PCs
plotTSNE(sce.dlpfc.tsne.10pcs, colour_by="cellType.split", point_alpha=0.5, point_size=4.0,
         text_size=8, theme_size=18) +
  ggtitle("t-SNE on top 10 PCs") + theme(plot.title = element_text(size=19))
# optimal PCs, version b (PC 5 removed)
plotTSNE(sce.dlpfc.tsne.optb, colour_by="cellType.split", point_alpha=0.5, point_size=4.0, text_by="cellType.split",
         text_size=8, theme_size=18) +
  ggtitle("t-SNE on optimal PCs (45), donor-correlated PC[5] removed") + theme(plot.title = element_text(size=15))
# and color by sample
plotTSNE(sce.dlpfc.tsne.optb, colour_by="sample", point_size=4.5, point_alpha=0.5,
         text_size=8, theme_size=18) +
  ggtitle("t-SNE on optimal PCs (45), donor-correlated PC[5] removed") + theme(plot.title = element_text(size=15))
dev.off()



# Save the candidates
Readme <- "This AMY SCE already has 50 'ambig.lowNtrxts' nuclei removed, and tSNE is on optimal 45 [-PC5] PCs"
save(sce.dlpfc.tsne.optb, Readme, file="rdas/ztemp_Amyg-n2_SCE-with-tSNEonOptPCs-minus-PC5_MNT.rda")



# Adapted from scater::plotReducedDim():
# Hidden function needed for this chunk ====
    .coerce_to_factor <- function(x, level.limit, msg) {
      if (!is.null(x)) {
        x <- as.factor(x)
        if (nlevels(x) > level.limit) {
          stop(sprintf("more than %i levels for '%s'", level.limit, msg))
        }
      }
      x
    }
    # ====

text_by <- "cellType.split"
text_out <- retrieveCellInfo(sce.dlpfc.tsne.optb, text_by, search="colData")
text_out$val <- .coerce_to_factor(text_out$val, level.limit=Inf)
## actually not necessary if the colData chosen (usually cellType[.etc] is factorized)
df_to_plot <- data.frame(reducedDim(sce.dlpfc.tsne.optb, "TSNE"))
by_text_x <- vapply(split(df_to_plot$X1, text_out$val), median, FUN.VALUE=0)
by_text_y <- vapply(split(df_to_plot$X2, text_out$val), median, FUN.VALUE=0)
# plot_out <- plot_out + annotate("text", x=by_text_x, y=by_text_y, 
#                             label=names(by_text_x), size=text_size, colour=text_colour)



plotTSNE(sce.dlpfc.tsne.optb, colour_by="cellType.split", point_size=4.5, point_alpha=0.5,
         text_size=8, theme_size=18) +
  annotate("text", x=by_text_x, y=by_text_y, 
           label=names(by_text_x), size=6) +
  ggtitle("t-SNE on optimal PCs (45), donor-correlated PC[5] removed")

# OR
sce.dlpfc.tsne.optb$labels <- ifelse(!duplicated(sce.dlpfc.tsne.optb$cellType.split), as.character(sce.dlpfc.tsne.optb$cellType.split), NA)
Labs.df <- data.frame(by_text_x, by_text_y, labs=names(by_text_x))

colDF <- data.frame(colData(sce.dlpfc.tsne.optb))
DFforLabs <- cbind(reducedDim(sce.dlpfc.tsne.optb,"TSNE"), data.frame(colDF$labels))
colnames(DFforLabs) <- c("X","Y","labels")

# plotTSNE(sce.dlpfc.tsne.optb, colour_by="cellType.split", point_size=4.5, point_alpha=0.5,
#          text_size=8, theme_size=18) +
#   geom_text_repel(data=DFforLabs,
#                   aes(label=labels)) +
#   ggtitle("t-SNE on top 15 PCs (>= 0.25% var)")
## OK THIS FINALLY WORKS
# -> can replace those X,Y with the median positions for those labels?

DFforLabs.edit <- DFforLabs
DFforLabs.edit$X[!is.na(DFforLabs$labels)] <- by_text_x[match(as.character(DFforLabs$labels[!is.na(DFforLabs$labels)]),
                                                              names(by_text_x))]
DFforLabs.edit$Y[!is.na(DFforLabs$labels)] <- by_text_y[match(as.character(DFforLabs$labels[!is.na(DFforLabs$labels)]),
                                                              names(by_text_y))]

## Finally print
library(ggrepel)

pdf("pdfs/revision/pubFigures/FINAL_pilotPaper_Amyg-n2_tSNE_optPCs-PC5_MNTMay2020.pdf", width=8)
set.seed(109)
plotTSNE(sce.dlpfc.tsne.optb, colour_by="cellType.split", point_size=6, point_alpha=0.5,
         theme_size=18) +
  geom_text_repel(data=DFforLabs.edit, size=6.0,
                  aes(label=labels)) +
  ggtitle("t-SNE on optimal PCs (45), donor-correlated PC[5] removed")
dev.off()









apply(table(sce.dlpfc$cellType.split, sce.dlpfc$donor),2,function(x){round(prop.table(x),3)})
    #                 Br5161 Br5212
    # Ambig.lowNtrxts  0.010  0.005
    # Astro            0.149  0.109
    # Excit.1          0.043  0.058
    # Excit.2          0.000  0.012
    # Excit.3          0.000  0.016
    # Inhib.1          0.005  0.046
    # Inhib.2          0.010  0.023
    # Inhib.3          0.003  0.007
    # Inhib.4          0.007  0.000
    # Inhib.5          0.026  0.004
    # Micro            0.129  0.101
    # Oligo            0.516  0.531
    # OPC              0.102  0.087

round(apply(apply(table(sce.dlpfc$cellType.split, sce.dlpfc$donor),2,prop.table),1,mean),3)
    # Ambig.lowNtrxts           Astro         Excit.1         Excit.2         Excit.3
    #           0.008           0.129           0.050           0.006           0.008
    #         Inhib.1         Inhib.2         Inhib.3         Inhib.4         Inhib.5
    #           0.026           0.016           0.005           0.004           0.015
    #           Micro           Oligo             OPC
    #           0.115           0.524           0.095


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

