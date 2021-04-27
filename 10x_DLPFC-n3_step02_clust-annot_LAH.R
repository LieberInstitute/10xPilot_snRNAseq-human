
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
     file="rdas/revision/regionSpecific_DLPFC_n3_cleaned-combined_SCE_LAH2021.rda")


    ## 'getClusteredPCs()' evaluated in qsub mode (with 'R-batchJob_Amyg-n2_optimalPCselxn_MNTFeb2020.R')
    #    --> saved into same .rda


### Picking up with optimally-defined PC space ===
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_DLPFC_n3_cleaned-combined_SCE_LAH2021.rda",
     verbose=TRUE)
    # sce.dlpfc, chosen.hvgs.amy, pc.choice.amy, ref.sampleInfo, ref.sampleInfo.rev

# How many PCs is optimal?:
metadata(pc.choice.amy)$chosen
    # [1] 87

## Assign this chosen (87 PCs) to 'PCA_opt'
reducedDim(sce.dlpfc, "PCA_opt") <- reducedDim(sce.dlpfc, "PCA_corrected")[ ,1:(metadata(pc.choice.amy)$chosen)]


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
    ## 59 prelim clusters

# Assign as 'prelimCluster'
sce.dlpfc$prelimCluster <- factor(clusters.k20)
plotReducedDim(sce.dlpfc, dimred="TSNE", colour_by="prelimCluster")

# Is sample driving this 'high-res' clustering at this level?
table(sce.dlpfc$prelimCluster, sce.dlpfc$sampleID)  # (a little bit, but is typical)

# rbind the ref.sampleInfo[.rev]
ref.sampleInfo <- rbind(ref.sampleInfo, ref.sampleInfo.rev)

# Save for now
save(sce.dlpfc, chosen.hvgs.amy, pc.choice.amy, ref.sampleInfo,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_DLPFC_n3_cleaned-combined_SCE_LAH2021.rda")


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
        #FALSE  TRUE
        #29381  4157

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
myplclust(tree.clusCollapsed, main="5x Amyg prelim-kNN-cluster relationships", cex.main=2, cex.lab=1.5, cex=1.8)


clust.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                               minClusterSize=2, deepSplit=1, cutHeight=250)


table(clust.treeCut)
unname(clust.treeCut[order.dendrogram(dend)])
    ## Cutting at 250 looks good for the main neuronal branch, but a lot of glial
     #    prelim clusters are dropped off (0's)

    # Cut at 400 for broad glia branch (will manually merge remaining dropped off)    
    glia.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                                  minClusterSize=2, deepSplit=1, cutHeight=400)
    unname(glia.treeCut[order.dendrogram(dend)])
    
    # Take those and re-assign to the first assignments
    clust.treeCut[order.dendrogram(dend)][c(38:59)] <- ifelse(glia.treeCut[order.dendrogram(dend)][c(38:59)] == 0,
                                                                      0, glia.treeCut[order.dendrogram(dend)][c(38:59)] + 10)
                                                                      
    unname(clust.treeCut[order.dendrogram(dend)])

# Add new labels to those prelimClusters cut off
clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)]==0)] <- max(clust.treeCut)+c(1:2, 3,3, 4:6)

# 'Re-write', since there are missing numbers
clust.treeCut[order.dendrogram(dend)] <- as.numeric(as.factor(clust.treeCut[order.dendrogram(dend)]))

labels_colors(dend) <- tableau20[clust.treeCut[order.dendrogram(dend)]]

# Print for future reference
pdf("pdfs/revision/regionSpecific_DLPFC_n3_HC-prelimCluster-relationships_MNT2021.pdf")
par(cex=0.8, font=2)
plot(dend, main="5x Amyg prelim-kNN-cluster relationships with collapsed assignments")
dev.off()


# Make reference for new cluster assignment
clusterRefTab.amy <- data.frame(origClust=order.dendrogram(dend),
                                merged=clust.treeCut[order.dendrogram(dend)])


# Assign as 'collapsedCluster'
sce.dlpfc$collapsedCluster <- factor(clusterRefTab.amy$merged[match(sce.dlpfc$prelimCluster, clusterRefTab.amy$origClust)])

# Print some visualizations:
pdf("pdfs/revision/regionSpecific_DLPFC_n3_reducedDims-with-collapsedClusters_MNT2021.pdf")
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

pdf("pdfs/revision/regionSpecific_DLPFC_n3_marker-logExprs_collapsedClusters_MNT2021.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.dlpfc, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]), ncol=2,
                   x="collapsedCluster", colour_by="collapsedCluster", point_alpha=0.4, point_size=.7,
                   add_legend=F) + stat_summary(fun = median, fun.min = median, fun.max = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:17], length(markers.mathys.custom[[i]]))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
    ggtitle(label=paste0(names(markers.mathys.custom)[i], " markers"))
  )
}
dev.off()




## QC - How do the total UMI distribution look?
newClusIndex <- splitit(sce.dlpfc$collapsedCluster)
sapply(newClusIndex, function(x) {quantile(sce.dlpfc[,x]$sum)})
    #             1      2        3      4        5       6         7       8       9
    # 0%     2876.0    831   7358.0   1685   566.00  2219.0   1238.00  1182.0   641.0
    # 25%   25828.5  22550  31449.0   9983  6921.25 12195.5  20472.75  7471.5  4634.5
    # 50%   48907.5  34143  50313.0  18481 11103.50 36275.0  45897.50 11275.0  6510.5
    # 75%   66821.5  43717  70976.5  39510 17152.25 49672.5  67085.00 16469.5  8964.0
    # 100% 156945.0 115493 165583.0 121273 95032.00 94599.0 171608.00 35728.0 30028.0

    #           10    11       12    *13       14    *15     16     17
    # 0%    1209.0   435  4036.00  121.0  1711.00  101.0  438.0  188.0
    # 25%   9787.5  3952  8556.00  245.0  4021.75  373.5 2382.0  850.5
    # 50%  13243.0  5381 10896.00  313.0  5210.00  743.0 3657.0 1156.0
    # 75%  16967.5  7055 13478.25  468.5  7050.25 1018.0 5029.5 1842.5
    # 100% 39830.0 20500 27167.00 2770.0 15214.00 1556.0 6237.0 7739.0


    ## Obser.: looks like collapsedCluster's 13 & 15 should be dropped - their marker
    ##         expression is noisy too;;  16 is unclear, but it may just be dropped
    ##         from consideration in more downstream analyses
    ##       - DOES look like 14 are endothelials...! (weirdly they came from the NeuN-enriched...)
    ## === === === === === === ==

table(sce.dlpfc$collapsedCluster)
    #   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17 
    # 780  541  399  525  500  555  216 1555 6080 1459 1201   44 1067   70   71   31   83
        #  * 13 is quite overrepresented in Br5400, but this sample also had a mode
        #    of <1000 total UMIs that were still kept after QC - thus probably a low trxts-driven cluster
        #  * 15 is pretty evenly-ish distributed (16 too, but will keep 16)

## Add annotations, looking at marker gene expression
annotationTab.amy <- data.frame(collapsedCluster=c(1:17))
annotationTab.amy$cellType <- NA
annotationTab.amy$cellType[c(1:2,4:7)] <- paste0("Inhib_", c("A","B","C","D","E","F"))
annotationTab.amy$cellType[c(3,12)] <- paste0("Excit_", c("A","B"))
annotationTab.amy$cellType[c(8:11)] <- c("Astro_A", "Oligo", "OPC", "Micro")
annotationTab.amy$cellType[c(13,15)] <- paste0("drop.lowNTx_", c("A","B"))
annotationTab.amy$cellType[c(14,16:17)] <- c("Endo", "ambig.glial", "Astro_B")


sce.dlpfc$cellType.prelim <- annotationTab.amy$cellType[match(sce.dlpfc$collapsedCluster,
                                                         annotationTab.amy$collapsedCluster)]
sce.dlpfc$cellType.prelim <- factor(sce.dlpfc$cellType.prelim)

# Save
save(sce.dlpfc, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo, annotationTab.amy,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_DLPFC_n3_cleaned-combined_SCE_LAH2021.rda")



## Re-print marker expression plots with annotated cluster names ===
pdf("pdfs/revision/regionSpecific_DLPFC_n3_marker-logExprs_collapsedClusters_MNT2021.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.dlpfc, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="cellType.prelim", colour_by="cellType.prelim", point_alpha=0.2, point_size=.7,
                   add_legend=F) +
      stat_summary(fun = median, fun.min = median, fun.max = median,
                   geom = "crossbar", width = 0.3,
                   #colour=rep(tableau20[1:17], length(markers.mathys.custom[[i]]))) +
                   colour=rep(tableau20[1:15], length(markers.mathys.custom[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(markers.mathys.custom)[i], " markers"))
  )
}
dev.off()

    ## Optionally, for a cleaner version, drop those 'drop.lowNTx_'s and re-print
    sce.dlpfc <- sce.dlpfc[ ,-grep("drop.", sce.dlpfc$cellType.prelim)]
    sce.dlpfc$cellType.prelim <- droplevels(sce.dlpfc$cellType.prelim)

      ## -> proceed to 'step03_markerDetxn-analyses[...].R'



### For reference === == === == ===
table(sce.dlpfc$cellType, sce.dlpfc$sample)
    #                 amy.5161 amy.5212
    # Ambig.lowNtrxts       34       16
    # Astro                489      363
    # Excit                226      261
    # Inhib                 84      295
    # Micro                425      339
    # Oligo               1697     1776
    # OPC                  335      292

table(sce.dlpfc$prelimCluster, sce.dlpfc$sample)
    #    amy.5161 amy.5212
    # 1       425      339
    # 2         8      113
    # 3         0      123
    # 4         2     1712
    # 5      1695       64
    # 6         8       42
    # 7        33       16
    # 8       473       20
    # 9        73        1
    # 10       85       13
    # 11        0      272
    # 12        0       55
    # 13        2       62
    # 14        0       60
    # 15       34       16
    # 16       68        0
    # 17        0      255
    # 18       24        0
    # 19      335       37
    # 20        0       69
    # 21        0       40
    # 22       14        9
    # 23       11       24

table(sce.dlpfc$cellType, sce.dlpfc$collapsedCluster)
#                    1    2    3    4    5    6    7
# Ambig.lowNtrxts    0    0    0    0    0   50    0
# Astro              0    0  852    0    0    0    0
# Excit              0  487    0    0    0    0    0
# Inhib            379    0    0    0    0    0    0
# Micro              0    0    0    0    0    0  764
# Oligo              0    0    0 3473    0    0    0
# OPC                0    0    0    0  627    0    0





### MNT 06May2020: subcluster-level annotations ========================================
# Some motivations:
#   i) to see if the high-VCAN-neuronal cluster exists in prelimCluster
#      as seen in the pan-brain level
#   ii) To now formally define subcluster-level populations

load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.dlpfc, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo

# Look at some marker expression at the prelimCluster level
pdf("pdfs/revision/ztemp_amyg-prelimCluster-neuronalMarkerExpression.pdf", width=8, height=7)
plotExpression(sce.dlpfc, exprs_values="logcounts", features=c("SNAP25","GAD1","GAD2","SLC17A6","SLC17A7","VCAN"),
               x="prelimCluster", colour_by="prelimCluster", ncol=2)
dev.off()


## More-manual annotations for neuronal subpops:
clusterRefTab.amy$cellType <- sce.dlpfc$cellType[match(clusterRefTab.amy$merged, sce.dlpfc$collapsedCluster)]
clusterRefTab.amy$cellType <- as.character(clusterRefTab.amy$cellType)

# Make new column for subclusers
clusterRefTab.amy$manual <- clusterRefTab.amy$cellType

    ## Excit subclusters: ====
    clusterRefTab.amy$manual <- ifelse(clusterRefTab.amy$origClust %in% c(3,9,16,20),
                                       paste0(clusterRefTab.amy$cellType, ".1"),
                                       as.character(clusterRefTab.amy$manual))
    
    # 10 and 12 would cut off at different heights; 12 expresses SLC17A7 and 10 doesn't really any:
    clusterRefTab.amy$manual <- ifelse(clusterRefTab.amy$origClust %in% c(10),
                                       paste0(clusterRefTab.amy$cellType, ".2"),
                                       as.character(clusterRefTab.amy$manual))
    
    clusterRefTab.amy$manual <- ifelse(clusterRefTab.amy$origClust %in% c(12),
                                       paste0(clusterRefTab.amy$cellType, ".3"),
                                       as.character(clusterRefTab.amy$manual))
    
    ## Inhib subclusters: merge 2/6 and 7/14 pairs, then split the rest
    clusterRefTab.amy$manual <- ifelse(clusterRefTab.amy$origClust %in% c(2,6),
                                       paste0(clusterRefTab.amy$cellType, ".1"),
                                       as.character(clusterRefTab.amy$manual))
    
    clusterRefTab.amy$manual <- ifelse(clusterRefTab.amy$origClust %in% c(7,14),
                                       paste0(clusterRefTab.amy$cellType, ".2"),
                                       as.character(clusterRefTab.amy$manual))
    
    clusterRefTab.amy$manual <- ifelse(clusterRefTab.amy$origClust %in% c(23),
                                       paste0(clusterRefTab.amy$cellType, ".3"),
                                       as.character(clusterRefTab.amy$manual))
    
    clusterRefTab.amy$manual <- ifelse(clusterRefTab.amy$origClust %in% c(18),
                                       paste0(clusterRefTab.amy$cellType, ".4"),
                                       as.character(clusterRefTab.amy$manual))
    
    clusterRefTab.amy$manual <- ifelse(clusterRefTab.amy$origClust %in% c(21),
                                       paste0(clusterRefTab.amy$cellType, ".5"),
                                       as.character(clusterRefTab.amy$manual))
    
    ## All other glial types will be kept the same
    
        ## Post-hoc: Looks like Excit.2 truly inhibitory, and Inhib.5 truly [or more so] excitatory
        #           (the latter is the ~39-40 high-VCAN nuclei ID'd in this sample at pan-brain) 
        #  -> swap them
        
        clusterRefTab.amy$manual <- ifelse(clusterRefTab.amy$origClust == 10,  # ("Excit.2")
                                           "Inhib.5",
                                           as.character(clusterRefTab.amy$manual))
        clusterRefTab.amy$manual <- ifelse(clusterRefTab.amy$origClust == 21,  # ("Inhib.5")
                                           "Excit.2",
                                           as.character(clusterRefTab.amy$manual))
    
        # --> THEN re-run the below
        
    ## end cluster splitting chunk ====

clusterRefTab.amy
    
## Add new annotations
sce.dlpfc$cellType.split <- clusterRefTab.amy$manual[match(sce.dlpfc$prelimCluster,
                                                         clusterRefTab.amy$origClust)]
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
save(sce.dlpfc, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo,
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
    # sce.dlpfc, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo


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

# Btw: metadata(pc.choice.amy)$chosen == 45 [PCs]

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
reducedDim(sce.dlpfc, "PCA_optb") <- reducedDim(sce.dlpfc, "PCA")[ ,c(1:4, 6:(metadata(pc.choice.amy)$chosen))]




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



### Session info for 21Apr2021 ==============================
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
# [7] scran_1.18.5                EnsDb.Hsapiens.v86_2.99.0  
# [9] ensembldb_2.14.1            AnnotationFilter_1.14.0    
# [11] GenomicFeatures_1.42.3      AnnotationDbi_1.52.0       
# [13] scater_1.18.6               ggplot2_3.3.3              
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
# [65] stringi_1.5.3             RSQLite_2.2.6            
# [67] BiocParallel_1.24.1       rlang_0.4.10             
# [69] pkgconfig_2.0.3           bitops_1.0-6             
# [71] lattice_0.20-41           purrr_0.3.4              
# [73] Rhdf5lib_1.12.1           GenomicAlignments_1.26.0 
# [75] labeling_0.4.2            cowplot_1.1.1            
# [77] bit_4.0.4                 tidyselect_1.1.0         
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

    
