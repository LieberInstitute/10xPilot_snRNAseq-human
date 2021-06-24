### MNT 10x snRNA-seq workflow: step 02
###   **Region-specific analyses: nucleus accumbens (NAc)**
###     - Preprint: (3x) un-selected samples + (2x) NeuN-sorted samples
###     - Revision: (3x) samples (2 female, 1 NeuN-sorted)
### Initiated MNT 04Mar2020
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

source('plotExpressionCustom.R')

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
    # pilot.data, pilot.data.unfiltered, e.out, ref.sampleInfo, pilot.data.alt, e.out.alt
    rm(pilot.data.unfiltered, e.out, e.out.alt)

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
sce.nac <- cbind(pilot.data[["br5161.nac"]], pilot.data[["br5212.nac"]], pilot.data[["br5287.nac"]],
                 # Revision un-selected samples:
                 pilot.data.2[["br5400.nac"]], pilot.data.2[["br5276.nac"]],
                 # NeuN-enriched:
                 pilot.data[["br5207.nac.neun"]],
                 # Re-sequenced Br5182-NAc sample:
                 pilot.data.alt[["br5182.nac.neun"]],
                 pilot.data.2[["br5701.nac.neun"]]  # (revision, NeuN-enriched)
                 )

sce.nac
    # class: SingleCellExperiment 
    # dim: 33538 20571 
    # metadata(8): Samples Samples ... Samples Samples
    # assays(1): counts
    # rownames(33538): MIR1302-2HG FAM138A ... AC213203.1 FAM231C
    # rowData names(6): gene_id gene_version ... gene_biotype Symbol.uniq
    # colnames(20571): AAACCCACATCGAACT-1 AAACCCATCCAACCAA-1 ...
    #   TTTGGTTAGCAGCACA-1 TTTGGTTAGGTCACAG-1
    # colData names(16): Sample Barcode ... protocol sequencer
    # reducedDimNames(0):
    # altExpNames(0):


# Use `multiBatchNorm()` to compute log-normalized counts, matching the scaling across samples
sce.nac <- multiBatchNorm(sce.nac, batch=sce.nac$sampleID)

# Use the simple `modelGeneVar` - this makes more sense over `combineVar`, since the
#   cell composition is already known to be quite different (with NeuN selection)
geneVar.nac <- modelGeneVar(sce.nac)
chosen.hvgs.nac <- geneVar.nac$bio > 0
sum(chosen.hvgs.nac)
    # [1] 11372


### Dimensionality reduction ================================================================

# Run `fastMNN` (internally uses `multiBatchPCA`), taking top 100 (instead of default 50 PCs)
set.seed(109)
mnn.hold <-  fastMNN(sce.nac, batch=sce.nac$sampleID,
                     merge.order=c("br5161.nac","br5212.nac","br5287.nac",
                                   "br5400.nac","br5276.nac",
                                   "br5207.nac.neun","br5182.nac.neun",
                                   "br5701.nac.neun"),
                     subset.row=chosen.hvgs.nac, d=100,
                     correct.all=TRUE, get.variance=TRUE,
                     BSPARAM=BiocSingular::IrlbaParam())
    # This temp file just used for getting batch-corrected components (drops a variety of entries)

table(colnames(mnn.hold) == colnames(sce.nac))  # all TRUE
table(mnn.hold$batch == sce.nac$sampleID) # all TRUE

# Add them to the SCE, as well as the metadata (though the latter might not be so usefl)
reducedDim(sce.nac, "PCA_corrected") <- reducedDim(mnn.hold, "corrected") # 100 components
metadata(sce.nac) <- metadata(mnn.hold)

# rbind the ref.sampleInfo[.rev]
ref.sampleInfo <- rbind(ref.sampleInfo, ref.sampleInfo.rev)

# Save into a new data file, which will dedicate for pan-brain-analyses
save(sce.nac, chosen.hvgs.nac, ref.sampleInfo,
     file="rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda")


    ## 'getClusteredPCs()' evaluated in qsub mode (with 'R-batchJob_NAc-n8_optimalPCselxn_MNT2021.R')
    #    --> saved into same .rda


### Resume work in optimal PC space
load("rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda", verbose=T)
    # sce.nac, chosen.hvgs.nac, ref.sampleInfo, pc.choice.nac

## How many PCs is optimal?:
metadata(pc.choice.nac)$chosen
    # 94

# Assign this chosen (94 PCs) to 'PCA_opt'
reducedDim(sce.nac, "PCA_opt") <- reducedDim(sce.nac, "PCA_corrected")[ ,1:(metadata(pc.choice.nac)$chosen)]


## t-SNE
set.seed(109)
sce.nac <- runTSNE(sce.nac, dimred="PCA_opt")


## UMAP
set.seed(109)
sce.nac <- runUMAP(sce.nac, dimred="PCA_opt")


# How do these look?
plotReducedDim(sce.nac, dimred="TSNE", colour_by="sampleID")
plotReducedDim(sce.nac, dimred="UMAP", colour_by="sampleID")




### Clustering: Two-step ======================================================
### Step 1: Perform graph-based clustering in this optimal PC space
#         - take k=20 NN to build graph
snn.gr <- buildSNNGraph(sce.nac, k=20, use.dimred="PCA_opt")
clusters.k20 <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters.k20)
    #    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
    # 3927 4262  239  251  988  283  285  718   40 5146  314  638  529   63   98   99 
    #   17   18   19   20   21   22   23   24   25   26   27   28   29 
    # 1000  651   22   58   86  429   36   52  240   41   37   21   18

# Assign as 'prelimCluster'
sce.nac$prelimCluster <- factor(clusters.k20)
table(sce.nac$prelimCluster, sce.nac$donor)

# Save for now
save(sce.nac, chosen.hvgs.nac, ref.sampleInfo, pc.choice.nac,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda")


### Step 2: Hierarchical clustering of pseudo-bulked ("PB'd") counts
#           ** That is, to pseudo-bulk (aka 'cluster-bulk') on raw counts, on all [non-zero] genes,
#              normalize with `librarySizeFactors()`, log2-transform, then perform HC'ing


## Preliminary cluster index for pseudo-bulking
clusIndexes = splitit(sce.nac$prelimCluster)
prelimCluster.PBcounts <- sapply(clusIndexes, function(ii){
  rowSums(assays(sce.nac)$counts[ ,ii])
 }
)

    # Btw: prelimCluster 23, 24, 26 may be doublet-driven:
    round(sapply(clusIndexes, function(x){quantile(sce.nac$doubletScore[x])}), 3)
        #           1      2     3     4     5     6      7     8     9     10    11     12    13
        # 0%    0.000  0.000 0.017 0.000 0.000 0.000  0.027 0.000 0.366  0.000 0.000  0.000 0.000
        # 25%   0.230  0.248 0.094 0.027 0.194 0.060  0.213 0.205 0.530  0.125 0.018  0.060 0.095
        # 50%   0.435  0.443 0.119 0.043 0.436 0.121  0.305 0.332 0.890  0.358 0.043  0.102 0.302
        # 75%   0.843  0.823 0.182 1.330 0.781 0.204  0.613 0.587 1.413  0.830 0.075  0.204 0.594
        # 100% 20.479 11.620 2.213 3.736 4.856 2.315 11.178 8.403 3.119 11.355 5.496 14.674 9.045
        #         14    15    16    17    18    19    20     21    22     23     24    25     26
        # 0%   0.000 0.017 0.008 0.000 0.000 0.008 0.034  0.009 0.000  2.269  2.028 0.000  3.550
        # 25%  0.021 0.053 0.058 0.079 0.047 0.025 0.131  0.036 0.025  4.380  8.185 0.026 13.219
        # 50%  0.028 0.080 0.095 0.145 0.089 0.059 0.234  0.060 0.053  7.026 12.780 0.044 16.139
        # 75%  0.049 1.136 0.127 0.294 0.131 0.102 0.406  0.478 0.105  8.306 15.285 0.112 18.726
        # 100% 0.089 2.919 1.871 6.723 4.215 1.328 1.856 17.153 1.517 12.332 25.141 4.006 24.444
        #         27     28    29
        # 0%   0.362  0.825 1.832
        # 25%  0.936  2.701 3.478
        # 50%  1.166  6.252 4.889   - tested findMarkers interactively on these and there are indeed
        # 75%  1.956  6.918 5.006     no real PW markers for prelimClusters 23, 24, 26, 28; can drop
        # 100% 6.708 10.558 5.225   - 29 actually has real markers, including SOX4, BCAN, GPR17, TNS3 (OPC_COP)
    
# Compute LSFs at this level
sizeFactors.PB.all  <- librarySizeFactors(prelimCluster.PBcounts)

# Normalize with these LSFs
geneExprs.temp <- t(apply(prelimCluster.PBcounts, 1, function(x) {log2(x/sizeFactors.PB.all + 1)}))


## Perform hierarchical clustering
dist.clusCollapsed <- dist(t(geneExprs.temp))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

dend <- as.dendrogram(tree.clusCollapsed, hang=0.2)

# Just for observation
myplclust(tree.clusCollapsed, main="all NAc sample clusters (n=8) prelim-kNN-cluster relationships",
          cex.main=1, cex.lab=0.8, cex=0.6)

sapply(clusIndexes, function(x) {quantile(sce.nac$sum[x])})
    # prelimCluster 13 is definitely a 'drop.lowNTx_'

clust.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                               minClusterSize=2, deepSplit=1, cutHeight=220)

table(clust.treeCut)
unname(clust.treeCut[order.dendrogram(dend)])

    ## Cut at 220, but manually merge others based on other 'cutHeight's
     #    (225, 300, 425--for glial)

## Add new labels to those prelimClusters cut off
clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)]==0)] <-
  max(clust.treeCut)+c(1:3, 4,4,4, 5:6, 7:9, 10,10,10, 11,11,12,12,13,rep(14,4))

## Define color pallet
cluster_colors <- unique(c(tableau20, tableau10medium)[clust.treeCut[order.dendrogram(dend)]])
names(cluster_colors) <- unique(clust.treeCut[order.dendrogram(dend)])
labels_colors(dend) <- cluster_colors[as.character(clust.treeCut[order.dendrogram(dend)])]

# Print for future reference
pdf("pdfs/revision/regionSpecific_NAc-n8_HC-prelimCluster-relationships_MNT2021.pdf")
par(cex=1.2, font=2)
plot(dend, main="All NAc (n=8) prelim-kNN-cluster relationships")
dev.off()


# Make reference for new cluster assignment
clusterRefTab.nac <- data.frame(origClust=order.dendrogram(dend),
                                merged=clust.treeCut[order.dendrogram(dend)])


# Assign as 'collapsedCluster'
sce.nac$collapsedCluster <- factor(clusterRefTab.nac$merged[match(sce.nac$prelimCluster, clusterRefTab.nac$origClust)])

# Print some visualizations:
pdf("pdfs/revision/regionSpecific_NAc-n8_reducedDims-with-collapsedClusters_MNT2021.pdf")
plotReducedDim(sce.nac, dimred="PCA_corrected", ncomponents=5, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.nac, colour_by="sampleID", point_alpha=0.5)
plotTSNE(sce.nac, colour_by="protocol", point_alpha=0.5)
plotTSNE(sce.nac, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.nac, colour_by="sum", point_alpha=0.5)
plotTSNE(sce.nac, colour_by="doubletScore", point_alpha=0.5)
# And some more informative UMAPs
plotUMAP(sce.nac, colour_by="sampleID", point_alpha=0.5)
plotUMAP(sce.nac, colour_by="collapsedCluster", point_alpha=0.5)
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
  'endothelial' = c('CLDN5', 'FLT1', 'VTN'),
  'Tcell' = c('SKAP1', 'ITK', 'CD247'),
  # Post-hoc: Some intersecting markers from mural cells in HPC, AMY
  'Mural' = c('COL1A2', 'TBX18', 'RBPMS'),
  # Kristen's MSN markers - not printed for the broader collapsed clusters
  'MSNs.pan' = c("PPP1R1B","BCL11B"),# "CTIP2")
  'MSNs.D1' = c("DRD1", "PDYN", "TAC1"),
  'MSNs.D2' = c("DRD2", "PENK")
)

## Print broad cell type markers
pdf("pdfs/revision/regionSpecific_NAc-n8_marker-logExprs_collapsedClusters_MNT2021.pdf", height=5, width=9)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpressionCustom(sce = sce.nac,
                         features = markers.mathys.custom[[i]], 
                         features_name = names(markers.mathys.custom)[[i]], 
                         anno_name = "collapsedCluster")# +
      #scale_color_manual(values = cell_colors.nac)
  )
}
dev.off()

# Do this for the prelimClusters too, since there's just 29
    # First remove the will-be 'drop.lowNTx_' prelimCluster 13
    # Post-hoc: also drop 24 & 26
    sce.temp <- sce.nac
    sce.temp <- sce.temp[ ,!(sce.temp$prelimCluster == 13)]
    sce.temp <- sce.temp[ ,!(sce.temp$prelimCluster %in% c(24,26))]
    sce.temp$prelimCluster <- droplevels(sce.temp$prelimCluster)
    
    pdf("pdfs/revision/regionSpecific_NAc-n8_marker-logExprs_prelimClusters_MNT2021.pdf", height=5, width=12)
    for(i in 1:length(markers.mathys.custom)){
      print(
        plotExpressionCustom(sce = sce.temp,
                             features = markers.mathys.custom[[i]], 
                             features_name = names(markers.mathys.custom)[[i]], 
                             anno_name = "prelimCluster") +
        scale_color_manual(values = c(tableau10medium, tableau20))
      )
    }
    dev.off()

# Check nuclear library sizes for annotation of 'lessCollapsed' cluster 6:
clusIndexes <- splitit(sce.nac$lessCollapsed)
sapply(clusIndexes, function(x) {quantile(sce.nac[ ,x]$sum)})
    #             1        2       3     4       5    6     7     8        9      10
    # 0%    1422.00  2471.00   365.0   635  1379.0  101  1247  2371  1197.00  2739.0
    # 25%  12988.75 19138.00  3618.5  6305  6937.5  133  3205 12217 19541.50 16635.5
    # 50%  16667.50 27300.00  5569.0  8380  9145.0  284  4525 15150 23778.00 20483.0
    # 75%  21967.75 37641.75  8010.5 10767 11354.0  608  6323 18557 30245.25 24723.5
    # 100% 44444.00 68589.00 31312.0 32867 31898.0 1648 12537 39311 69559.00 67355.0
    #         11       12      13    14       15       16
    # 0%    6603  6143.00  3473.0  1388  3895.00  1801.00
    # 25%  15968 11429.75 14588.5 10748 20340.00 18638.50
    # 50%  20199 15010.50 18914.0 13472 24847.50 22881.50
    # 75%  33506 24642.25 24444.5 17118 30074.75 28108.75
    # 100% 60522 68219.00 67288.0 43590 85460.00 84761.00


## Add annotations for 'lessCollapsed' clusters, looking at marker gene expression
annotationTab.nac <- data.frame(cluster=c(1, 2, 3, 4, 5,
                                              6, 7, 8, 9, 10,
                                              11, 12, 13, 14, 15,
                                              16),
                                    cellType.split=c("MSN.D1.1", "MSN.broad", "Oligo", "Astro", "OPC",
                                               "ambig.lowNtrxts", "Micro", "MSN.D1.2", "MSN.D2.1","MSN.D1.3",
                                               "Inhib.1","Inhib.2","Inhib.3","Inhib.4","MSN.D1.4",
                                               "MSN.D2.2"),
                                    cellType.broad=c("MSN.D1", "MSN.broad", "Oligo", "Astro", "OPC",
                                                     "ambig.lowNtrxts", "Micro", "MSN.D1", "MSN.D2", "MSN.D1",
                                                     rep("Inhib", 4), "MSN.D1",
                                                     "MSN.D2")
)

# Add info to clusterRefTab.nac
clusterRefTab.nac$cellType.broad <- annotationTab.nac$cellType.broad[match(clusterRefTab.nac$merged.man.b,
                                                                                   annotationTab.nac$cluster)]
clusterRefTab.nac$cellType.split <- annotationTab.nac$cellType.split[match(clusterRefTab.nac$merged.man.b,
                                                                                   annotationTab.nac$cluster)]


## Then add to SCE - like for DLPFC, give 'cellType' & 'cellType.split'
sce.nac$cellType <- annotationTab.nac$cellType.broad[match(sce.nac$lessCollapsed,
                                                                annotationTab.nac$cluster)]

sce.nac$cellType.split <- annotationTab.nac$cellType.split[match(sce.nac$lessCollapsed,
                                                                 annotationTab.nac$cluster)]


## Save
save(sce.nac, chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac, ref.sampleInfo, 
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda")


# Re-print marker expression with cell type labels and dropping 'ambig.lowNtrxts' cluster
table(sce.dlpfc$cellType.split)

# First drop "Ambig.lowNtrxts" (93 nuclei)
sce.nac <- sce.nac[ ,sce.nac$cellType.split != "ambig.lowNtrxts"]
sce.nac$cellType.split <- droplevels(sce.nac$cellType.split)

pdf("pdfs/regionSpecific_NAc-ALL-n5_marker-logExprs_collapsedClusters_Mar2020.pdf", height=6, width=12)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.nac, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7,
                   add_legend=F) +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
                   width = 0.3, colour=rep(tableau20[1:15], length(markers.mathys.custom[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(markers.mathys.custom)[i], " markers"))
  )
}
dev.off()


## Re-print some reducedDims with this information
#pdf("pdfs/zold_regionSpecific_NAc-ALL-n5_reducedDims-with-collapsedClusters_Apr2020.pdf")
pdf("pdfs/regionSpecific_NAc-ALL-n5_reducedDims-with-cellType.final_Apr2020.pdf")
plotReducedDim(sce.nac, dimred="PCA", ncomponents=5, colour_by="cellType", point_alpha=0.5)
plotTSNE(sce.nac, colour_by="processDate", point_size=3.5, point_alpha=0.5) + ggtitle("t-SNE on opt PCs")
plotTSNE(sce.nac, colour_by="sample", point_size=3.5, point_alpha=0.5) + ggtitle("t-SNE on opt PCs")
plotTSNE(sce.nac, colour_by="cellType", point_size=3.5, point_alpha=0.5) + ggtitle("t-SNE on opt PCs")
#plotTSNE(sce.nac, colour_by="cellType.split", point_size=3.5, point_alpha=0.5) + ggtitle("t-SNE on opt PCs")
plotTSNE(sce.nac, colour_by="cellType.final", point_size=3.5, point_alpha=0.5) + ggtitle("t-SNE on opt PCs")
plotTSNE(sce.nac, colour_by="sum", point_size=3.5, point_alpha=0.5) + ggtitle("t-SNE on opt PCs")
# UMAP
#plotUMAP(sce.nac, colour_by="cellType.split", point_size=3.5, point_alpha=0.5) + ggtitle("UMAP on opt PCs")
plotUMAP(sce.nac, colour_by="cellType.final", point_size=3.5, point_alpha=0.5) + ggtitle("UMAP on opt PCs")
dev.off()

      ## -> proceed to 'step03_markerDetxn-analyses[...].R'





### MNT 06Apr2020: Further dive into 'MSN.broad' cluster, which is entirely Br5212 =============
  # -> revisit SNN clsutering - see if `cut_at()` can tease these differences out
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda",
     verbose=T)
    # sce.nac, chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac, ref.sampleInfo


table(sce.nac$cellType.split, sce.nac$sample)
    #                 nac.5161 nac.5212 nac.5287 nac.neun.5182 nac.neun.5207
    # ambig.lowNtrxts       19       42       22             7             3
    # Astro                149      384       12             0             0
    # Inhib.1                1        3        0            16             5
    # Inhib.2                1        1        1            42            11
    # Inhib.3                7        7        9            86           167
    # Inhib.4                9        8        4           104            58
    # Micro                 72       72       37             0             0
    # MSN.broad              0      266        0             0             0
    # MSN.D1.1               2        0        0           117            13
    # MSN.D1.2              10        3        0           285             3
    # MSN.D1.3              17        8        6           369           319
    # MSN.D1.4             178        2       72          1505          1829
    # MSN.D2.1               9        6        3           134           148
    # MSN.D2.2              41       14        5          1602          1870
    # Oligo               1454      854      499             0             0
    # OPC                   98      104       37             0             0




## Trying `igraph::cut_at` ===============
library(igraph)
snn.gr <- buildSNNGraph(sce.nac, k=20, use.dimred="PCA_opt")
#clusters.k20 <- igraph::cluster_walktrap(snn.gr)$membership
clusters.k20.comms <- cluster_walktrap(snn.gr)
# Store the initial assignment
og.membership.nac <- clusters.k20.comms$membership
table(og.membership.nac)
    #    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16
    # 1847  276  266  300 1925 1739  719  161  301 1619 1913  181  132   56  882  183
    #   17   18   19   20   21   22
    #  384   31  139   62  100   25

is_hierarchical(clusters.k20.comms) # TRUE
    
    # Btw: it's 'prelimCluster' 3 (== "MSN.broad"; 266 nuclei) that we want to break up
    table(sce.nac$prelimCluster, sce.nac$cellType)

# "`cut_at()` returns a numeric vector, the membership vector of the vertices."
test.membership <- cut_at(clusters.k20.comms, no=23)

table(test.membership, og.membership.nac)
table(test.membership[which(og.membership.nac==3)])

for(i in c(23:30)){
  print(paste0("no. of output clusters: ", i))
  test.membership <- cut_at(clusters.k20.comms, no=i)
  print(table(test.membership[which(og.membership.nac==3)]))
}
    ## [...] [1] "no. of output clusters: 25"
            #  16  23
            # 167  99
    ## so no = 25 is where this 'MSN.broad' cluster is split
    rm(i, test.membership)
    new.membership.5212msn <- cut_at(clusters.k20.comms, no=25)
    
    # Check
    table(new.membership.5212msn, og.membership.nac)  # good
    
    # Check order
    table(sce.nac$prelimCluster == og.membership.nac) # all TRUE

# Add this new level, first pushing to 'high' factor so isn't accidentally added to another
sce.nac$prelimCluster.split <- ifelse(sce.nac$prelimCluster==3,
                                          new.membership.5212msn[og.membership.nac==3]+22,
                                          sce.nac$prelimCluster)
#    # 38, then 45
#
#    # Actually when tabulate it's unexpected...
#    table(sce.nac$prelimCluster.split)  # 165 + 101, instead of 167 + 99 ......
#                                            # also when took this approach (proceeded with
#                                            # 'Reassign to 23, 24' chunk), it doesn't look
#                                            # properly separated by D1-vs-D2...
    
        # Just add new info as is and check 'new clusters' 16 + 23
        sce.nac$prelimCluster.temp <- factor(new.membership.5212msn)
            # Then plotted.  This showed this `cut_at()` approach works actually.
        sce.nac$prelimCluster.temp <- NULL
        
        
# Try this way:
new.membership.5212msn[og.membership.nac==3] <- new.membership.5212msn[og.membership.nac==3] + 22
# Now apply as'prelimCluster.split'
sce.nac$prelimCluster.split <- ifelse(sce.nac$prelimCluster==3,
                                          new.membership.5212msn,
                                          sce.nac$prelimCluster)
    
    
# Reassign to 23, 24
sce.nac$prelimCluster.split[sce.nac$prelimCluster.split==38] <- 23
sce.nac$prelimCluster.split[sce.nac$prelimCluster.split==45] <- 24
table(sce.nac$prelimCluster.split)

sce.nac$cellType.moreSplit <- ifelse(sce.nac$prelimCluster.split %in% c(23,24),
                                       paste0(sce.nac$cellType.split,".",sce.nac$prelimCluster.split),
                                       as.character(sce.nac$cellType.split))

sce.nac$cellType.moreSplit <- factor(sce.nac$cellType.moreSplit)

table(sce.nac$cellType.moreSplit, sce.nac$cellType.split)
    # good.


## Check
plotExpression(sce.nac, exprs_values = "logcounts", features=c(markers.mathys.custom[["MSNs.D1"]]),
               x="cellType.moreSplit", colour_by="cellType.moreSplit", point_alpha=0.5, point_size=.7,
               add_legend=F) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
               width = 0.3, colour=rep(tableau20[1:17], length(markers.mathys.custom[["MSNs.D1"]]))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ## good.


      # Driven by low n transcripts?
      cc_idx <- splitit(sce.nac$cellType.moreSplit)
      sapply(cc_idx, function(x){quantile(sce.nac[ ,x]$sum)})
          ## doesn't seem to be
      
      
      sce.nac.5212msns <- sce.nac[ ,sce.nac$cellType.split=="MSN.broad"]
      
      pheatmap::pheatmap(assay(sce.nac.5212msns,"logcounts")[c(markers.mathys.custom[["MSNs.D1"]],
                                                               markers.mathys.custom[["MSNs.D2"]],
                                                               markers.mathys.custom[["MSNs.pan"]]), ],
                         show_colnames=FALSE)
          ## Alright dichotomy is definitely there...
           # - revisiting OSCA Ch. 10, there is a subclustering approach


## OSCA Ch. 10.7 Subclustering approach =============
# Repeating modelling and PCA on the subset.
  # use 'sce.nac.5212msns' from above
dec.5212 <- modelGeneVar(sce.nac.5212msns)
sce.nac.5212msns <- denoisePCA(sce.nac.5212msns, technical=dec.5212,
                         subset.row=getTopHVGs(dec.5212, prop=0.1))
    # reduced "PCA" entry to 30 PCs
g.5212 <- buildSNNGraph(sce.nac.5212msns, use.dimred="PCA")
clust.5212 <- igraph::cluster_walktrap(g.5212)$membership
table(clust.5212)
#  1   2
# 98 168

plotExpression(sce.nac.5212msns, features=c("DRD1", "DRD2", "BCL11B"),
               x=I(factor(clust.5212)))
    # Dope this works too.

  ## end test for this approach - go with `cut_at()` results ====


# First save
save(sce.nac, chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac, ref.sampleInfo,
     file="rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda")

# First drop "Ambig.lowNtrxts" (93 nuclei)
sce.nac <- sce.nac[ ,sce.nac$cellType.moreSplit != "ambig.lowNtrxts"]
sce.nac$cellType.moreSplit <- droplevels(sce.nac$cellType.moreSplit)

## Print broad cell type markers
#pdf("pdfs/ztemp_regionSpecific_NAc-ALL-n5_marker-logExprs_cellType.moreSplit_MNTApr2020.pdf", height=6, width=12)
pdf("pdfs/regionSpecific_NAc-ALL-n5_marker-logExprs_cellType.final_MNTApr2020.pdf", height=6, width=12)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.nac, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
#                   x="cellType.moreSplit", colour_by="cellType.moreSplit", point_alpha=0.5, point_size=.7,
                   x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7,
                   add_legend=F) +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
                   width = 0.3, colour=rep(tableau20[1:14], length(markers.mathys.custom[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(markers.mathys.custom)[i], " markers"))
  )
}
dev.off()


## Notes from this insight:
 #    - MSN.broad.23 D1-expressing for sure; most like D1.3 or D1.4
 #    - MSN.broad.24 like MSN.D2.2 in that expresses DRD2 & PENK

# Pseudo-bulk across JUST this cellType.moreSplit and run some correlation to re-assign

sce.nac.PB <- aggregateAcrossCells(sce.nac, ids=sce.nac$cellType.moreSplit,
                                   use_exprs_values="counts")
# Drop genes with all 0's
sce.nac.PB <- sce.nac.PB[!rowSums(assay(sce.nac.PB, "counts"))==0, ]
    ## keeps 29236 genes

# Remove stored `sizeFactors()` because this will mess you up
#     * Also, to be safe, can always provide manually-computed SFs:
sizeFactors(sce.nac.PB) <- NULL
LSFvec <- librarySizeFactors(sce.nac.PB)
# Re-compute logcounts
sce.nac.PB <- logNormCounts(sce.nac.PB, size_factors=LSFvec)

pdf("pdfs/ztemp_regionSpecific_NAc-ALL-n5_corExprs_PB-cellType.moreSplit_MNTApr2020.pdf")
pheatmap::pheatmap(cor(assay(sce.nac.PB, "logcounts")),display_numbers=T, fontsize_number=7)
dev.off()
    # so MSN.broad.23 most ~ MSN.D1.4
    # and MSN.broad.24 most ~ MSN.D2.2  (yep, see below)
cor(assay(sce.nac.PB, "logcounts"))[grep("broad", colnames(sce.nac.PB)),grep("MSN",colnames(sce.nac.PB))]
    #              MSN.broad.23 MSN.broad.24  MSN.D1.1  MSN.D1.2  MSN.D1.3  MSN.D1.4
    # MSN.broad.23    1.0000000    0.9552372 0.9056483 0.9276173 0.9334208 0.9559685*
    # MSN.broad.24    0.9552372    1.0000000 0.8983732 0.9162969 0.9253921 0.9479631
    #               MSN.D2.1  MSN.D2.2
    # MSN.broad.23 0.9212994 0.9539614
    # MSN.broad.24 0.9171996 0.9504078*

# Re-load actually, bc had dropped "ambig.lowNtrxts" actually
load("rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda", verbose=T)

## We'll reassign as such then, calling a new 'cellType.final':
sce.nac$cellType.final <- sce.nac$cellType.moreSplit

sce.nac$cellType.final[sce.nac$cellType.moreSplit=="MSN.broad.23"] <- "MSN.D1.4"
sce.nac$cellType.final[sce.nac$cellType.moreSplit=="MSN.broad.24"] <- "MSN.D2.2"
sce.nac$cellType.final <- droplevels(sce.nac$cellType.final)

table(sce.nac$cellType.moreSplit, sce.nac$cellType.final) # dope.

# Save
save(sce.nac, chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac, ref.sampleInfo,
     file="rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda")



### Some checks
plotTSNE(sce.nac, colour_by="cellType.split", point_alpha=0.5, point_size=3.5, text_by="cellType.split")

## MSN.D1.3-specific top [pt-coding] genes: 
plotExpression(sce.nac, exprs_values = "logcounts", features=c("DRD1", "CASZ1", "CRHR2", "RXFP1"),
               x="cellType.moreSplit", colour_by="cellType.moreSplit", point_alpha=0.5, point_size=.7,
               add_legend=F) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
               width = 0.3, colour=rep(tableau20[1:17], 4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


    ##  -> then went ahead and re-printed broad markers with 'cellType.final'




## Added MNT 12May2020: tSNE in lower dims ===
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda",
     verbose=T)
    # sce.nac, chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac, ref.sampleInfo


# How many PCs?
head(attr(reducedDim(sce.nac, "PCA"), "percentVar"), n=50)
    # [1] 27.87823928  4.31556377  2.28210852  1.53394544  1.34230365  0.84148529
    # [7]  0.80660587  0.68569439  0.60321278  0.46073112  0.43747455  0.32713745
    # [13]  0.30032199  0.29464166  0.28570101  0.21158092  0.18735573  0.16784829
    # [19]  0.15606074  0.13760112  0.12960231  0.11897204  0.11822540  0.10666966
    # [25]  0.09385858  0.09166851  0.08853377  0.08355615  0.08057045  0.07368351
    # [31]  0.06726170  0.06580073  0.06347227  0.06081381  0.05676892  0.05429146
    # [37]  0.05142856  0.05020298  0.04900668  0.04745873  0.04702622  0.04585275


# 0.05% var or greater
reducedDim(sce.nac, "PCA_38") <- reducedDim(sce.nac, "PCA")[ ,c(1:38)]
# 0.1% var or greater
reducedDim(sce.nac, "PCA_24") <- reducedDim(sce.nac, "PCA")[ ,c(1:24)]
# 0.25% var or greater
reducedDim(sce.nac, "PCA_15") <- reducedDim(sce.nac, "PCA")[ ,c(1:15)]
# Top 10 (as Mathys, et al)
reducedDim(sce.nac, "PCA_10") <- reducedDim(sce.nac, "PCA")[ ,c(1:10)]

# First remove this reducedDim bc this has caused trouble previously
reducedDim(sce.nac, "TSNE") <- NULL

## 38 PCs tsNE === 
set.seed(109)
sce.all.tsne.38pcs <- runTSNE(sce.nac, dimred="PCA_38")

## 24 PCs tSNE ===    (not much different)
# set.seed(109)
# sce.all.tsne.24pcs <- runTSNE(sce.nac, dimred="PCA_24")

## 15 PCs tsNE ===
set.seed(109)
sce.all.tsne.15pcs <- runTSNE(sce.nac, dimred="PCA_15")

## 10 PCs tsNE ===
set.seed(109)
sce.all.tsne.10pcs <- runTSNE(sce.nac, dimred="PCA_10")



# Drop "Ambig.lowNtrxts" cluster as always
sce.all.tsne.38pcs <- sce.all.tsne.38pcs[ ,sce.all.tsne.38pcs$cellType.final != "ambig.lowNtrxts"] # 93
sce.all.tsne.38pcs$cellType.final <- droplevels(sce.all.tsne.38pcs$cellType.final)

sce.all.tsne.15pcs <- sce.all.tsne.15pcs[ ,sce.all.tsne.15pcs$cellType.final != "ambig.lowNtrxts"] # 93
sce.all.tsne.15pcs$cellType.final <- droplevels(sce.all.tsne.15pcs$cellType.final)

sce.all.tsne.10pcs <- sce.all.tsne.10pcs[ ,sce.all.tsne.10pcs$cellType.final != "ambig.lowNtrxts"] # 93
sce.all.tsne.10pcs$cellType.final <- droplevels(sce.all.tsne.10pcs$cellType.final)


pdf("pdfs/pubFigures/NAc-ALL-n5_tSNE_38-15-10-PCs_MNTMay2020.pdf", width=9)
# 38 PCs
plotTSNE(sce.all.tsne.38pcs, colour_by="cellType.final", point_alpha=0.5, point_size=4.0,
         text_size=8, theme_size=18) +
  ggtitle("t-SNE on top 38 PCs (>= 0.05% var)") + theme(plot.title = element_text(size=19))
# # 24 PCs
# plotTSNE(sce.all.tsne.24pcs, colour_by="cellType.final", point_alpha=0.5, point_size=4.0,
#          text_size=8, theme_size=18) +
#   ggtitle("t-SNE on top 24 PCs (>= 0.10% var)") + theme(plot.title = element_text(size=19))
# 15 PCs
plotTSNE(sce.all.tsne.15pcs, colour_by="cellType.final", point_alpha=0.5, point_size=4.0,
         text_size=8, theme_size=18) +
  ggtitle("t-SNE on top 15 PCs (>= 0.25% var)") + theme(plot.title = element_text(size=19))
# 10 PCs
plotTSNE(sce.all.tsne.10pcs, colour_by="cellType.final", point_alpha=0.5, point_size=4.0,
         text_size=8, theme_size=18) +
  ggtitle("t-SNE on top 10 PCs") + theme(plot.title = element_text(size=19))
dev.off()



## Then just some reducedDims once again with these - just sample/processDate & cellType.final
pdf("pdfs/zforRef_regionSpecific_NAc-ALL-n5_reducedDims-with-cellType.final_15PCs_May2020.pdf", width=8)
plotTSNE(sce.all.tsne.15pcs, colour_by="processDate", point_size=4.5, point_alpha=0.5,
         text_size=8, theme_size=18) +
  ggtitle("t-SNE on top 15 PCs (>= 0.25% var)")
plotTSNE(sce.all.tsne.15pcs, colour_by="sample", point_size=4.5, point_alpha=0.5,
         text_size=8, theme_size=18) +
  ggtitle("t-SNE on top 15 PCs (>= 0.25% var)")
plotTSNE(sce.all.tsne.15pcs, colour_by="cellType.final", point_size=4.5, point_alpha=0.5,
         text_size=8, theme_size=18, text_by="cellType.final") +
  ggtitle("t-SNE on top 15 PCs (>= 0.25% var)")
dev.off()


# Save the candidates
save(sce.all.tsne.38pcs, file="rdas/ztemp_NAc-ALL-n5_SCE-with-tSNEon38PCs_MNT.rda")

save(sce.all.tsne.10pcs, sce.all.tsne.15pcs,
     file="rdas/ztemp_NAc-ALL-n5_SCE-with-tSNEon15-10PCs_MNT.rda")


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
      
        text_by <- "cellType.final"
        text_out <- retrieveCellInfo(sce.all.tsne.15pcs, text_by, search="colData")
        text_out$val <- .coerce_to_factor(text_out$val, level.limit=Inf)
            ## actually not necessary if the colData chosen (usually cellType[.etc] is factorized)
            df_to_plot <- data.frame(reducedDim(sce.all.tsne.15pcs, "TSNE"))
        by_text_x <- vapply(split(df_to_plot$X1, text_out$val), median, FUN.VALUE=0)
        by_text_y <- vapply(split(df_to_plot$X2, text_out$val), median, FUN.VALUE=0)
        # plot_out <- plot_out + annotate("text", x=by_text_x, y=by_text_y, 
        #                             label=names(by_text_x), size=text_size, colour=text_colour)
    
    

## TSNE with [repelled] labels (ggrepel::geom_text_repel()) ===
 # Added 21May2020 - to finalize-finalize (& for grants)
load("rdas/ztemp_NAc-ALL-n5_SCE-with-tSNEon15-10PCs_MNT.rda", verbose=T)
    # sce.all.tsne.10pcs, sce.all.tsne.15pcs
    rm(sce.all.tsne.10pcs)




plotTSNE(sce.all.tsne.15pcs, colour_by="cellType.final", point_size=4.5, point_alpha=0.5,
         text_size=8, theme_size=18) +
  annotate("text", x=by_text_x, y=by_text_y, 
           label=names(by_text_x), size=6) +
  ggtitle("t-SNE on top 15 PCs (>= 0.25% var)")

# OR
sce.all.tsne.15pcs$labels <- ifelse(!duplicated(sce.all.tsne.15pcs$cellType.final), as.character(sce.all.tsne.15pcs$cellType.final), NA)
Labs.df <- data.frame(by_text_x, by_text_y, labs=names(by_text_x))

colDF <- data.frame(colData(sce.all.tsne.15pcs))
DFforLabs <- cbind(reducedDim(sce.all.tsne.15pcs,"TSNE"), data.frame(colDF$labels))
colnames(DFforLabs) <- c("X","Y","labels")

# plotTSNE(sce.all.tsne.15pcs, colour_by="cellType.final", point_size=4.5, point_alpha=0.5,
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
pdf("pdfs/pubFigures/FINAL_pilotPaper_NAc-ALL-n5_tSNE_15PCs_MNTMay2020.pdf", width=8)
set.seed(109)
plotTSNE(sce.all.tsne.15pcs, colour_by="cellType.final", point_size=5, point_alpha=0.5,
         theme_size=18) +
  geom_text_repel(data=DFforLabs.edit, size=5.5,
                  aes(label=labels)) +
  ggtitle("t-SNE on top 15 PCs (>= 0.25% var)")
dev.off()

# Grant version - different seeds still overlay D1.1 & D1.3 w/ bigger font, so just adjust those positions
#               - crank up the axes label sizes too - will just need to trim off the title & legend
#                 bc these are linked
DFforLabs.edit[!is.na(DFforLabs$labels), ]
DFforLabs.edit$Y[which(DFforLabs$labels=="MSN.D1.3")] <- DFforLabs.edit$Y[which(DFforLabs$labels=="MSN.D1.3")]-2

pdf("pdfs/pubFigures/FINAL_forGrant_NAc-ALL-n5_tSNE_15PCs_MNTMay2020.pdf", width=9)
set.seed(109)
plotTSNE(sce.all.tsne.15pcs, colour_by="cellType.final", point_size=6.5, point_alpha=0.5,
         theme_size=25) +
  geom_text_repel(data=DFforLabs.edit, size=7,
                  aes(label=labels)) +
  ggtitle("t-SNE on top 15 PCs (>= 0.25% var)")
dev.off()



## Heatmap of broad marker genes (leaving out defined markers/subcluster for now) ===
load("rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda", verbose=T)
    # sce.nac, chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac, ref.sampleInfo

# First drop "Ambig.lowNtrxts" (93 nuclei)
sce.nac <- sce.nac[ ,sce.nac$cellType.final != "ambig.lowNtrxts"]
sce.nac$cellType.final <- droplevels(sce.nac$cellType.final)

cell.idx <- splitit(sce.nac$cellType.final)
dat <- as.matrix(assay(sce.nac, "logcounts"))

pdf('pdfs/pubFigures/heatmap-geneExprs_all-NAc-n5_cellType.final_mean-broadMarkers_MNT21May2020.pdf', useDingbats=TRUE, height=6, width=6)
#pdf('pdfs/pubFigures/heatmap-geneExprs_all-NAc-n5_cellType.final_median-broadMarkers_MNT21May2020.pdf', height=6, width=6)
genes <- c('DRD1','TAC1','DRD2','PENK','PPP1R1B','GAD1','SNAP25','CAMK2A','MBP','PDGFRA','AQP4','CD74')
current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMeans(dat[genes, ii])))
# # or medians:
# current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMedians(dat[genes, ii])))
    # # for some reason `rowMedians()` doesn't keep row names...
    # rownames(current_dat) <- genes
current_dat <- current_dat[ ,c(7:12, 2:5, 13:14,1,6)]
pheatmap(current_dat, cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(0.02, 4, length.out = 101),
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "OrRd"))(100), fontsize = 19,
         legend=T)
grid::grid.text(label="Expression \n (log2)",x=0.93,y=0.60)
dev.off()




## For Day Lab ( & Keri's requests) - MNT 03Apr2020 ========================
# First drop "Ambig.lowNtrxts" (93 nuclei)
sce.nac <- sce.nac[ ,sce.nac$cellType.final != "ambig.lowNtrxts"]
sce.nac$cellType.final <- droplevels(sce.nac$cellType.final)

## AFTER re-assigning MSN.broad (and added all DRD genes for kicks)
pdf("pdfs/exploration/zref_logExprs_gene-requests_Martinowich.pdf", width=7, height=9)
plotExpression(sce.nac, exprs_values = "logcounts", features=c("PVALB", "KIT", "CHAT", "RELN", "TH",
                                                                   paste0("DRD", 1:5)),
               x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7,
               add_legend=F, ncol=3) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
               width = 0.3, colour=rep(tableau20[1:14], 10)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle(label="Miscellaneous genes of interest")
pdf.options(height=4.5)
dev.off()

# Added 27Apr2020
pdf("pdfs/exploration/zref_logExprs_gene-requests_Martinowich.2.pdf", width=7, height=5)
plotExpression(sce.nac, exprs_values = "logcounts", features=c("NGFR","NTRK1","NTRK2","NTRK3","BDNF"),
               x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7,
               add_legend=F, ncol=3) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
               width = 0.3, colour=rep(tableau20[1:14], 5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


c("NGFR", "NTRK1", "NTRK2", "NTRK3", "BDNF") %in% rownames(sce.nac)

# Other checks by request - not to print
plotExpression(sce.nac, exprs_values = "logcounts", features="HTR7",
               x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7,
               add_legend=F, ncol=3) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
               width = 0.3, colour=rep(tableau20[1:14], 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



## Session info for 23Jun2021 =============================================
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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices datasets  utils     methods  
# [9] base     
# 
# other attached packages:
#   [1] scDblFinder_1.4.0           gridExtra_2.3               dynamicTreeCut_1.63-1      
# [4] dendextend_1.14.0           jaffelab_0.99.30            rafalib_1.0.0              
# [7] DropletUtils_1.10.3         batchelor_1.6.2             scran_1.18.5               
# [10] scater_1.18.6               ggplot2_3.3.3               EnsDb.Hsapiens.v86_2.99.0  
# [13] ensembldb_2.14.1            AnnotationFilter_1.14.0     GenomicFeatures_1.42.3     
# [16] AnnotationDbi_1.52.0        SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0
# [19] Biobase_2.50.0              GenomicRanges_1.42.0        GenomeInfoDb_1.26.7        
# [22] IRanges_2.24.1              S4Vectors_0.28.1            BiocGenerics_0.36.1        
# [25] MatrixGenerics_1.2.1        matrixStats_0.58.0         
# 
# loaded via a namespace (and not attached):
#   [1] googledrive_1.0.1         ggbeeswarm_0.6.0          colorspace_2.0-0         
# [4] ellipsis_0.3.2            scuttle_1.0.4             bluster_1.0.0            
# [7] XVector_0.30.0            BiocNeighbors_1.8.2       rstudioapi_0.13          
# [10] farver_2.1.0              bit64_4.0.5               fansi_0.4.2              
# [13] xml2_1.3.2                splines_4.0.4             R.methodsS3_1.8.1        
# [16] sparseMatrixStats_1.2.1   cachem_1.0.4              Rsamtools_2.6.0          
# [19] ResidualMatrix_1.0.0      dbplyr_2.1.1              R.oo_1.24.0              
# [22] HDF5Array_1.18.1          compiler_4.0.4            httr_1.4.2               
# [25] dqrng_0.2.1               assertthat_0.2.1          Matrix_1.3-2             
# [28] fastmap_1.1.0             lazyeval_0.2.2            limma_3.46.0             
# [31] BiocSingular_1.6.0        prettyunits_1.1.1         tools_4.0.4              
# [34] rsvd_1.0.3                igraph_1.2.6              gtable_0.3.0             
# [37] glue_1.4.2                GenomeInfoDbData_1.2.4    dplyr_1.0.5              
# [40] rappdirs_0.3.3            Rcpp_1.0.6                vctrs_0.3.6              
# [43] Biostrings_2.58.0         rhdf5filters_1.2.0        rtracklayer_1.50.0       
# [46] DelayedMatrixStats_1.12.3 stringr_1.4.0             beachmat_2.6.4           
# [49] lifecycle_1.0.0           irlba_2.3.3               statmod_1.4.35           
# [52] XML_3.99-0.6              edgeR_3.32.1              zlibbioc_1.36.0          
# [55] scales_1.1.1              hms_1.0.0                 ProtGenerics_1.22.0      
# [58] rhdf5_2.34.0              RColorBrewer_1.1-2        curl_4.3                 
# [61] memoise_2.0.0             segmented_1.3-3           biomaRt_2.46.3           
# [64] stringi_1.5.3             RSQLite_2.2.7             BiocParallel_1.24.1      
# [67] rlang_0.4.10              pkgconfig_2.0.3           bitops_1.0-7             
# [70] lattice_0.20-41           purrr_0.3.4               Rhdf5lib_1.12.1          
# [73] labeling_0.4.2            GenomicAlignments_1.26.0  cowplot_1.1.1            
# [76] bit_4.0.4                 tidyselect_1.1.1          magrittr_2.0.1           
# [79] R6_2.5.0                  generics_0.1.0            DelayedArray_0.16.3      
# [82] DBI_1.1.1                 pillar_1.6.0              withr_2.4.2              
# [85] RCurl_1.98-1.3            tibble_3.1.1              crayon_1.4.1             
# [88] xgboost_1.3.2.1           utf8_1.2.1                BiocFileCache_1.14.0     
# [91] viridis_0.6.0             progress_1.2.2            locfit_1.5-9.4           
# [94] grid_4.0.4                data.table_1.14.0         blob_1.2.1               
# [97] digest_0.6.27             R.utils_2.10.1            openssl_1.4.3            
# [100] munsell_0.5.0             beeswarm_0.3.1            viridisLite_0.4.0        
# [103] vipor_0.4.5               askpass_1.1 

