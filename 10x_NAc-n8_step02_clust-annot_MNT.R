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

    # Btw: prelimCluster 24 & 26, especially, may be doublet-driven, but also 23, 28 & 29:
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
    
    # And total N transcripts / cluster:
    sapply(clusIndexes, function(x){quantile(sce.nac$sum[x])})
        #             1         2       3       4        5     6      7        8        9       10
        # 0%     1774.0   2680.00  2362.0  1951.0  1828.00  2531   3807   4290.0  5246.00   565.00
        # 25%   21349.0  20138.00 16942.5 12086.0  7261.25 13269  21212  23099.0 12017.75  4231.75
        # 50%   26944.0  25311.00 20716.0 16269.0  9707.50 16479  27878  29959.5 16953.50  5894.50
        # 75%   35404.5  33227.75 25313.0 21893.5 12737.75 20451  34620  40725.5 24388.50  7531.75
        # 100% 114562.0 171929.00 72727.0 70890.0 29548.00 77527 109993 124105.0 65162.00 31312.00
        #            11        12    13     14        15      16       17      18     19       20
        # 0%    1197.00   2500.00   101  105.0   7257.00   798.0   635.00  1379.0 1164.0  4534.00
        # 25%  20204.75  17329.75   289  469.5  21798.25  4049.5  7457.75  8630.0 2566.5 11219.00
        # 50%  25540.50  22217.00   461  866.0  28172.50  6524.0 10811.50 11284.0 4229.5 13969.50
        # 75%  33380.75  27756.00   881 1327.5  38033.25  9175.5 15545.00 14231.5 6374.0 17356.75
        # 100% 92007.00 113864.00 71212 5239.0 118689.00 17601.0 37265.00 29825.0 8855.0 21624.00
        #            21    22       23       24       25    26    27    28       29
        # 0%     751.00  1504  7539.00 18719.00  2579.00  3430  4273  3388  9916.00
        # 25%  12951.25  4268 13242.25 40880.25 14446.50 28996 16622  9667 14691.75
        # 50%  15630.50  5545 17571.00 49361.00 18832.00 37389 24143 12934 18738.50
        # 75%  22642.50  6789 25108.25 63095.75 25553.25 45372 35972 15563 24131.25
        # 100% 73538.00 12537 36220.00 95998.00 80360.00 71127 69839 26000 42806.00
    
### MNT comment 25Jun: skip 'collapsing' of prelimClusters, since several dropped from above QC:
# # Compute LSFs at this level
# sizeFactors.PB.all  <- librarySizeFactors(prelimCluster.PBcounts)
# 
# # Normalize with these LSFs
# geneExprs.temp <- t(apply(prelimCluster.PBcounts, 1, function(x) {log2(x/sizeFactors.PB.all + 1)}))
# 
# 
# ## Perform hierarchical clustering
# dist.clusCollapsed <- dist(t(geneExprs.temp))
# tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")
# 
# dend <- as.dendrogram(tree.clusCollapsed, hang=0.2)
# 
# # Just for observation
# myplclust(tree.clusCollapsed, main="all NAc sample clusters (n=8) prelim-kNN-cluster relationships",
#           cex.main=1, cex.lab=0.8, cex=0.6)
# 
# sapply(clusIndexes, function(x) {quantile(sce.nac$sum[x])})
#     # prelimCluster 13 is definitely a 'drop.lowNTx_'
# 
# clust.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
#                                minClusterSize=2, deepSplit=1, cutHeight=220)
# 
# table(clust.treeCut)
# unname(clust.treeCut[order.dendrogram(dend)])
# 
#     ## Cut at 220, but manually merge others based on other 'cutHeight's
#      #    (225, 300, 425--for glial)
# 
# ## Add new labels to those prelimClusters cut off
# clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)]==0)] <-
#   max(clust.treeCut)+c(1:3, 4,4,4, 5:6, 7:9, 10,10,10, 11,11,12,12,13,rep(14,4))
# 
# ## Define color pallet
# cluster_colors <- unique(c(tableau20, tableau10medium)[clust.treeCut[order.dendrogram(dend)]])
# names(cluster_colors) <- unique(clust.treeCut[order.dendrogram(dend)])
# labels_colors(dend) <- cluster_colors[as.character(clust.treeCut[order.dendrogram(dend)])]
# 
# # Print for future reference
# pdf("pdfs/revision/regionSpecific_NAc-n8_HC-prelimCluster-relationships_MNT2021.pdf")
# par(cex=1.2, font=2)
# plot(dend, main="All NAc (n=8) prelim-kNN-cluster relationships")
# dev.off()
# 
# 
# # Make reference for new cluster assignment
# clusterRefTab.nac <- data.frame(origClust=order.dendrogram(dend),
#                                 merged=clust.treeCut[order.dendrogram(dend)])
# 
# 
# # Assign as 'collapsedCluster'
# sce.nac$collapsedCluster <- factor(clusterRefTab.nac$merged[match(sce.nac$prelimCluster, clusterRefTab.nac$origClust)])

# # Print some visualizations:    (just print after annotating)
# pdf("pdfs/revision/regionSpecific_NAc-n8_reducedDims-with-collapsedClusters_MNT2021.pdf")
# plotReducedDim(sce.nac, dimred="PCA_corrected", ncomponents=5, colour_by="collapsedCluster", point_alpha=0.5)
# plotTSNE(sce.nac, colour_by="sampleID", point_alpha=0.5)
# plotTSNE(sce.nac, colour_by="protocol", point_alpha=0.5)
# plotTSNE(sce.nac, colour_by="prelimCluster", point_alpha=0.5)
# plotTSNE(sce.nac, colour_by="sum", point_alpha=0.5)
# plotTSNE(sce.nac, colour_by="doubletScore", point_alpha=0.5)
# # And some more informative UMAPs
# plotUMAP(sce.nac, colour_by="sampleID", point_alpha=0.5)
# plotUMAP(sce.nac, colour_by="prelimCluster", point_alpha=0.5)
# dev.off()

## Print marker genes for annotation
markers.mathys.custom = list(
  'neurons' = c('SYT1', 'SNAP25', 'GRIN1'),
  'excitatory_neuron' = c('CAMK2A', 'NRGN','SLC17A7', 'SLC17A6', 'SLC17A8'),
  'inhibitory_neuron' = c('GAD1', 'GAD2', 'SLC32A1'),
  'oligodendrocyte' = c('MBP', 'MOBP', 'PLP1'),
  'oligodendrocyte_precursor' = c('PDGFRA', 'VCAN', 'CSPG4'),
  # Post-hoc from interactive marker analysis
  'differn_committed_OPC' = c("SOX4", "BCAN", "GPR17", "TNS3"),
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
                         anno_name = "collapsedCluster")
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


## Save
save(sce.nac, chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac, ref.sampleInfo, 
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda")


### MNT Note 25Jun2021 ===
  #     Since there were only 29 'prelimCluster' from kNN graph-based clustering, and 5 of these
  #     are being dropped due to doublet ID/low Tx-driven clustering, we are not not utilizing any
  #     hierarchical clustering like with the other regions for a 'collapsedCluster' column ===


## Annotate 'prelimCluster's, based on broad cell type markers ===
#    (canonical, above, and from markers looked at in 'step3')
annotationTab.nac <- data.frame(cluster=c(1:29))
annotationTab.nac$cellType <- NA
annotationTab.nac$cellType[c(4,9,15,25,27)] <- paste0("Inhib_", c("A","B","C","D","E"))
annotationTab.nac$cellType[c(1,3,6,8,12,21)] <- paste0("MSN.D1_", c("A","B","C","D","E","F"))
annotationTab.nac$cellType[c(2,7,11,20)] <- paste0("MSN.D2_", c("A","B","C","D"))
annotationTab.nac$cellType[c(5,10)] <- c("Oligo_A","Oligo_B")
annotationTab.nac$cellType[c(16,17,22)] <- c("Astro_A","Astro_B","Micro")
annotationTab.nac$cellType[c(18,29)] <- c("OPC","OPC_COP")
annotationTab.nac$cellType[c(14,19)] <- c("Micro_resting","Macrophage")
annotationTab.nac$cellType[c(13, 23,24,26,28)] <- c("drop.lowNTx",
                                                    paste0("drop.doublet_", c("A","B","C","D")))


sce.nac$cellType <- annotationTab.nac$cellType[match(sce.nac$prelimCluster,
                                                     annotationTab.nac$cluster)]
sce.nac$cellType <- factor(sce.nac$cellType)

cell_colors.nac <- c(tableau20, tableau10medium)[1:29]
names(cell_colors.nac) <- annotationTab.nac$cellType
cell_colors.nac
    #      MSN.D1_A       MSN.D2_A       MSN.D1_B        Inhib_A        Oligo_A       MSN.D1_C 
    #     "#1F77B4"      "#AEC7E8"      "#FF7F0E"      "#FFBB78"      "#2CA02C"      "#98DF8A" 
    #      MSN.D2_B       MSN.D1_D        Inhib_B        Oligo_B       MSN.D2_C       MSN.D1_E 
    #     "#D62728"      "#FF9896"      "#9467BD"      "#C5B0D5"      "#8C564B"      "#C49C94" 
    #   drop.lowNTx  Micro_resting        Inhib_C        Astro_A        Astro_B            OPC 
    #     "#E377C2"      "#F7B6D2"      "#7F7F7F"      "#C7C7C7"      "#BCBD22"      "#DBDB8D" 
    #    Macrophage       MSN.D2_D       MSN.D1_F          Micro drop.doublet_A drop.doublet_B 
    #     "#17BECF"      "#9EDAE5"      "#729ECE"      "#FF9E4A"      "#67BF5C"      "#ED665D" 
    #       Inhib_D drop.doublet_C        Inhib_E drop.doublet_D        OPC_COP 
    #     "#AD8BC9"      "#A8786E"      "#ED97CA"      "#A2A2A2"      "#CDCC5D"


# Save
save(sce.nac, chosen.hvgs.nac, pc.choice.nac, ref.sampleInfo, annotationTab.nac, cell_colors.nac,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda")



## Re-print marker expression plots with annotated cluster names ===
pdf("pdfs/revision/regionSpecific_NAc-n8_marker-logExprs_annotatedClusters_MNT2021.pdf", height=5, width=9)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpressionCustom(sce = sce.nac,
                         features = markers.mathys.custom[[i]], 
                         features_name = names(markers.mathys.custom)[[i]], 
                         anno_name = "cellType") +
      scale_color_manual(values = cell_colors.nac)
  )
}
dev.off()

    ## Optionally, for a cleaner version, drop those 'drop.[lowNTx/doublet]_'s and re-print
    sce.nac <- sce.nac[ ,-grep("drop.", sce.nac$cellType)]
    sce.nac$cellType <- droplevels(sce.nac$cellType)



## Re-print reducedDims with these annotations (keep the 'drop.' clusters here) ===
pdf("pdfs/revision/regionSpecific_NAc-n8_reducedDims-with-annotatedClusters_MNT2021.pdf",width=8)
plotReducedDim(sce.nac, dimred="PCA_corrected", ncomponents=5, colour_by="cellType", point_alpha=0.5) +
  scale_color_manual(values = cell_colors.nac) + labs(colour="Cell type")
plotTSNE(sce.nac, colour_by="sampleID", point_alpha=0.5, point_size=2)
plotTSNE(sce.nac, colour_by="protocol", point_alpha=0.5, point_size=2)
plotTSNE(sce.nac, colour_by="prelimCluster", text_by="prelimCluster",
         text_size=3, point_alpha=0.5, point_size=2)
plotTSNE(sce.nac, colour_by="cellType", text_by="cellType",
         text_size=3, point_alpha=0.5, point_size=2) +
  scale_color_manual(values = cell_colors.nac,
                     labels=paste0(levels(sce.nac$cellType)," (",table(sce.nac$cellType),")")) +
  labs(colour="Cell type")
plotTSNE(sce.nac, colour_by="sum", point_alpha=0.5, point_size=2)
plotTSNE(sce.nac, colour_by="doubletScore", point_alpha=0.5, point_size=2)
# And some more informative UMAPs
plotUMAP(sce.nac, colour_by="prelimCluster", text_by="prelimCluster",
         text_size=3, point_alpha=0.5, point_size=2)
plotUMAP(sce.nac, colour_by="cellType", text_by="cellType",
         text_size=3, point_alpha=0.5, point_size=2) +
  scale_color_manual(values = cell_colors.nac,
                     labels=paste0(levels(sce.nac$cellType)," (",table(sce.nac$cellType),")")) +
  labs(colour="Cell type")
dev.off()

      ## -> proceed to 'step03_markerDetxn-analyses[...].R'




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

    
    
### For main Fig (revision): TSNE with [repelled] labels (ggrepel::geom_text_repel()) ===
load("rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda", verbose=T)
    # sce.nac, chosen.hvgs.nac, pc.choice.nac, ref.sampleInfo, annotationTab.nac, cell_colors.nac
        
# Drop those 'drop.[lowNTx/doublet]_'s and re-print
sce.nac <- sce.nac[ ,-grep("drop.", sce.nac$cellType)]

# A priori: Also remove rare / not-homeostatic cell pops - just for main figure/section
sce.nac <- sce.nac[ ,-c(grep("Macro_infilt", sce.nac$cellType),
                        grep("Micro_resting",sce.nac$cellType),
                        grep("OPC_COP", sce.nac$cellType))]
sce.nac$cellType <- droplevels(sce.nac$cellType)

    text_by <- "cellType"
    text_out <- retrieveCellInfo(sce.nac, text_by, search="colData")
    text_out$val <- .coerce_to_factor(text_out$val, level.limit=Inf)
        ## actually not necessary if the colData chosen (usually cellType[.etc] is factorized)
        df_to_plot <- data.frame(reducedDim(sce.nac, "TSNE"))
    by_text_x <- vapply(split(df_to_plot$X1, text_out$val), median, FUN.VALUE=0)
    by_text_y <- vapply(split(df_to_plot$X2, text_out$val), median, FUN.VALUE=0)
            # plot_out <- plot_out + annotate("text", x=by_text_x, y=by_text_y, 
            #                             label=names(by_text_x), size=text_size, colour=text_colour)

plotTSNE(sce.nac, colour_by="cellType", point_size=4.5, point_alpha=0.5,
         text_size=8, theme_size=18) +
  annotate("text", x=by_text_x, y=by_text_y, 
           label=names(by_text_x), size=6) +
  scale_color_manual(values = cell_colors.nac) +
  ggtitle("t-SNE on optimal fastMNN-corrected PCs (94)")

# OR
sce.nac$labels <- ifelse(!duplicated(sce.nac$cellType), as.character(sce.nac$cellType), NA)
Labs.df <- data.frame(by_text_x, by_text_y, labs=names(by_text_x))

colDF <- data.frame(colData(sce.nac))
DFforLabs <- cbind(reducedDim(sce.nac,"TSNE"), data.frame(colDF$labels))
colnames(DFforLabs) <- c("X","Y","labels")

# -> can replace those X,Y with the median positions for those labels
DFforLabs.edit <- DFforLabs
DFforLabs.edit$X[!is.na(DFforLabs$labels)] <- by_text_x[match(as.character(DFforLabs$labels[!is.na(DFforLabs$labels)]),
                                                         names(by_text_x))]
DFforLabs.edit$Y[!is.na(DFforLabs$labels)] <- by_text_y[match(as.character(DFforLabs$labels[!is.na(DFforLabs$labels)]),
                                                              names(by_text_y))]

## Finally print
library(ggrepel)
pdf("pdfs/revision/pubFigures/regionSpecific_NAc-n8_main-fig-TSNE_MNT2021.pdf", height=6.5, width=10.5)
set.seed(109)
plotTSNE(sce.nac, colour_by="cellType", point_size=4.5, point_alpha=0.5,
         theme_size=16) + labs(colour="Cell type") +
  geom_text_repel(data=DFforLabs.edit, size=4.5,
                  aes(label=labels)) +
  scale_color_manual(values = cell_colors.nac,
                     labels=paste0(levels(sce.nac$cellType)," (",table(sce.nac$cellType),")")) +
  ggtitle("t-SNE on optimal fastMNN-corrected PCs (94)")
dev.off()

# # Grant version - different seeds still overlay D1.1 & D1.3 w/ bigger font, so just adjust those positions
# #               - crank up the axes label sizes too - will just need to trim off the title & legend
# #                 bc these are linked
# DFforLabs.edit[!is.na(DFforLabs$labels), ]
# DFforLabs.edit$Y[which(DFforLabs$labels=="MSN.D1.3")] <- DFforLabs.edit$Y[which(DFforLabs$labels=="MSN.D1.3")]-2
# 
# pdf("pdfs/pubFigures/FINAL_forGrant_NAc-ALL-n5_tSNE_15PCs_MNTMay2020.pdf", width=9)
# set.seed(109)
# plotTSNE(sce.nac, colour_by="cellType", point_size=6.5, point_alpha=0.5,
#          theme_size=25) +
#   geom_text_repel(data=DFforLabs.edit, size=7,
#                   aes(label=labels)) +
#   ggtitle("t-SNE on top 15 PCs (>= 0.25% var)")
# dev.off()



## Heatmap of broad marker genes (leaving out defined markers/subcluster for now) ===
cell.idx <- splitit(sce.nac$cellType)
dat <- as.matrix(assay(sce.nac, "logcounts"))

pdf('pdfs/revision/pubFigures/heatmap-geneExprs_NAc-n8_mean-broadMarkers_MNT2021.pdf', useDingbats=TRUE, height=6, width=7)
#pdf('pdfs/revision/pubFigures/heatmap-geneExprs_NAc-n8_median-broadMarkers_MNT2021.pdf', useDingbats=TRUE, height=6, width=7)
genes <- c('DRD1','TAC1','DRD2','PENK','PPP1R1B','GAD1','SNAP25','CAMK2A','MBP','PDGFRA','AQP4','CD74')
current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMeans(dat[genes, ii])))
# # or medians:
    # current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMedians(dat[genes, ii])))
    # # for some reason `rowMedians()` doesn't keep row names...
    # rownames(current_dat) <- genes
current_dat <- current_dat[ ,c(9:18, 3:7, 19:21,1:2,8)]
pheatmap(current_dat, cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(0.02, 4, length.out = 101),
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "OrRd"))(100), fontsize = 17,
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



### Session info for 25Jun2021 ============
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
#   [1] dynamicTreeCut_1.63-1       dendextend_1.14.0           jaffelab_0.99.30           
# [4] rafalib_1.0.0               DropletUtils_1.10.3         batchelor_1.6.2            
# [7] scran_1.18.5                scater_1.18.6               ggplot2_3.3.3              
# [10] EnsDb.Hsapiens.v86_2.99.0   ensembldb_2.14.1            AnnotationFilter_1.14.0    
# [13] GenomicFeatures_1.42.3      AnnotationDbi_1.52.0        SingleCellExperiment_1.12.0
# [16] SummarizedExperiment_1.20.0 Biobase_2.50.0              GenomicRanges_1.42.0       
# [19] GenomeInfoDb_1.26.7         IRanges_2.24.1              S4Vectors_0.28.1           
# [22] BiocGenerics_0.36.1         MatrixGenerics_1.2.1        matrixStats_0.58.0         
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
# [61] memoise_2.0.0             gridExtra_2.3             segmented_1.3-3          
# [64] biomaRt_2.46.3            stringi_1.5.3             RSQLite_2.2.7            
# [67] BiocParallel_1.24.1       rlang_0.4.10              pkgconfig_2.0.3          
# [70] bitops_1.0-7              lattice_0.20-41           purrr_0.3.4              
# [73] Rhdf5lib_1.12.1           labeling_0.4.2            GenomicAlignments_1.26.0 
# [76] cowplot_1.1.1             bit_4.0.4                 tidyselect_1.1.1         
# [79] magrittr_2.0.1            R6_2.5.0                  generics_0.1.0           
# [82] DelayedArray_0.16.3       DBI_1.1.1                 pillar_1.6.0             
# [85] withr_2.4.2               RCurl_1.98-1.3            tibble_3.1.1             
# [88] crayon_1.4.1              utf8_1.2.1                BiocFileCache_1.14.0     
# [91] viridis_0.6.0             progress_1.2.2            locfit_1.5-9.4           
# [94] grid_4.0.4                blob_1.2.1                digest_0.6.27            
# [97] R.utils_2.10.1            openssl_1.4.3             munsell_0.5.0            
# [100] beeswarm_0.3.1            viridisLite_0.4.0         vipor_0.4.5              
# [103] askpass_1.1 

