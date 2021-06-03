### MNT 10x snRNA-seq workflow: step 02
###   **Region-specific analyses**
###     - (2x) amygdala samples from: Br5161 & Br5212
### Initiated MNT 29Jan2020
### MNT 21Apr2021: add expansion samples (n=3, incl'g 2 female)
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
sce.amy <- cbind(pilot.data[["br5161.amy"]], pilot.data[["br5212.amy"]],
                 pilot.data.2[["br5701.amy"]],
                 pilot.data.2[["br5276.amy.neun"]], pilot.data.2[["br5400.amy.neun"]]
                 )

sce.amy
    #class: SingleCellExperiment 
    # dim: 33538 15177 
    # metadata(5): Samples Samples Samples Samples Samples
    # assays(1): counts
    # rownames(33538): MIR1302-2HG FAM138A ... AC213203.1 FAM231C
    # rowData names(6): gene_id gene_version ... gene_biotype Symbol.uniq
    # colnames(15177): AAACCCAAGCACGATG-1 AAACCCACAGCGGTCT-1 ...
    #   TTTGGTTGTTGTTTGG-1 TTTGGTTTCAGACCGC-1
    # colData names(16): Sample Barcode ... protocol sequencer
    # reducedDimNames(0):
    # altExpNames(0):

# Use `multiBatchNorm()` to compute log-normalized counts, matching the scaling across samples
sce.amy <- multiBatchNorm(sce.amy, batch=sce.amy$sampleID)

# Use the simple `modelGeneVar` - this makes more sense over `combineVar`, since the
#   cell composition is already known to be quite different (with NeuN selection)
geneVar.amy <- modelGeneVar(sce.amy)
chosen.hvgs.amy <- geneVar.amy$bio > 0
sum(chosen.hvgs.amy)
    # [1] 10109


### Dimensionality reduction ================================================================

# Run `fastMNN` (internally uses `multiBatchPCA`), taking top 100 (instead of default 50 PCs)
set.seed(109)
mnn.hold <-  fastMNN(sce.amy, batch=sce.amy$sampleID,
                     merge.order=c("br5161.amy","br5212.amy","br5701.amy",
                                   "br5276.amy.neun","br5400.amy.neun"),
                     subset.row=chosen.hvgs.amy, d=100,
                     correct.all=TRUE, get.variance=TRUE,
                     BSPARAM=BiocSingular::IrlbaParam())
    # This temp file just used for getting batch-corrected components (drops a variety of entries)

table(colnames(mnn.hold) == colnames(sce.amy))  # all TRUE
table(mnn.hold$batch == sce.amy$sampleID) # all TRUE

# Add them to the SCE, as well as the metadata (though the latter might not be so usefl)
reducedDim(sce.amy, "PCA_corrected") <- reducedDim(mnn.hold, "corrected") # 100 components
metadata(sce.amy) <- metadata(mnn.hold)

# Save into a new region-specific SCE object/flie
save(sce.amy, chosen.hvgs.amy, ref.sampleInfo, ref.sampleInfo.rev,
     file="rdas/revision/regionSpecific_Amyg-n5_cleaned-combined_SCE_MNT2021.rda")


    ## 'getClusteredPCs()' evaluated in qsub mode (with 'R-batchJob_Amyg-n2_optimalPCselxn_MNTFeb2020.R')
    #    --> saved into same .rda


### Picking up with optimally-defined PC space ===
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_Amyg-n5_cleaned-combined_SCE_MNT2021.rda",
     verbose=TRUE)
    # sce.amy, chosen.hvgs.amy, pc.choice.amy, ref.sampleInfo, ref.sampleInfo.rev

# How many PCs is optimal?:
metadata(pc.choice.amy)$chosen
    # [1] 87

## Assign this chosen (87 PCs) to 'PCA_opt'
reducedDim(sce.amy, "PCA_opt") <- reducedDim(sce.amy, "PCA_corrected")[ ,1:(metadata(pc.choice.amy)$chosen)]


## t-SNE
set.seed(109)
sce.amy <- runTSNE(sce.amy, dimred="PCA_opt")


## UMAP
set.seed(109)
sce.amy <- runUMAP(sce.amy, dimred="PCA_opt")


# How do these look?
plotReducedDim(sce.amy, dimred="TSNE", colour_by="sampleID")
plotReducedDim(sce.amy, dimred="UMAP", colour_by="sampleID")


### Clustering: Two-step ======================================================
### Step 1: Perform graph-based clustering in this optimal PC space
#         - take k=20 NN to build graph
snn.gr <- buildSNNGraph(sce.amy, k=20, use.dimred="PCA_opt")
clusters.k20 <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters.k20)
    ## 59 prelim clusters

# Assign as 'prelimCluster'
sce.amy$prelimCluster <- factor(clusters.k20)
plotReducedDim(sce.amy, dimred="TSNE", colour_by="prelimCluster")

# Is sample driving this 'high-res' clustering at this level?
table(sce.amy$prelimCluster, sce.amy$sampleID)  # (a little bit, but is typical)

# rbind the ref.sampleInfo[.rev]
ref.sampleInfo <- rbind(ref.sampleInfo, ref.sampleInfo.rev)

# Save for now
save(sce.amy, chosen.hvgs.amy, pc.choice.amy, ref.sampleInfo,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_Amyg-n5_cleaned-combined_SCE_MNT2021.rda")


### Step 2: Hierarchical clustering of pseudo-bulked ("PB'd") counts
  #           ** That is, to pseudo-bulk (aka 'cluster-bulk') on raw counts, on all [non-zero] genes,
  #              normalize with `librarySizeFactors()`, log2-transform, then perform HC'ing


# Preliminary cluster index for pseudo-bulking
clusIndexes = splitit(sce.amy$prelimCluster)
prelimCluster.PBcounts <- sapply(clusIndexes, function(ii){
  rowSums(assays(sce.amy)$counts[ ,ii])
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
myplclust(tree.clusCollapsed, main="5x Amyg prelim-kNN-cluster relationships", cex.main=2, cex.lab=1.5, cex=1.8)


clust.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                               minClusterSize=2, deepSplit=1, cutHeight=225)


table(clust.treeCut)
unname(clust.treeCut[order.dendrogram(dend)])
    ## Cutting at 225 looks good for the main neuronal branch, but a lot of glial
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
clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)]==0)] <- max(clust.treeCut)+c(1:10)

# 'Re-write', since there are missing numbers
clust.treeCut[order.dendrogram(dend)] <- as.numeric(as.factor(clust.treeCut[order.dendrogram(dend)]))

## Define color pallet
cluster_colors <- unique(c(tableau20, tableau10medium)[clust.treeCut[order.dendrogram(dend)]])
names(cluster_colors) <- unique(clust.treeCut[order.dendrogram(dend)])
labels_colors(dend) <- cluster_colors[as.character(clust.treeCut[order.dendrogram(dend)])]

# Print for future reference
pdf("pdfs/revision/regionSpecific_Amyg-n5_HC-prelimCluster-relationships_MNT2021.pdf")
par(cex=0.8, font=2)
plot(dend, main="5x Amyg prelim-kNN-cluster relationships with collapsed assignments")
dev.off()


# Make reference for new cluster assignment
clusterRefTab.amy <- data.frame(origClust=order.dendrogram(dend),
                                merged=clust.treeCut[order.dendrogram(dend)])


# Assign as 'collapsedCluster'
sce.amy$collapsedCluster <- factor(clusterRefTab.amy$merged[match(sce.amy$prelimCluster, clusterRefTab.amy$origClust)])

# Print some visualizations:
pdf("pdfs/revision/regionSpecific_Amyg-n5_reducedDims-with-collapsedClusters_MNT2021.pdf")
plotReducedDim(sce.amy, dimred="PCA_corrected", ncomponents=5, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.amy, colour_by="sampleID", point_alpha=0.5)
plotTSNE(sce.amy, colour_by="protocol", point_alpha=0.5)
plotTSNE(sce.amy, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.amy, colour_by="sum", point_alpha=0.5)
plotTSNE(sce.amy, colour_by="doubletScore", point_alpha=0.5)
# And some more informative UMAPs
plotUMAP(sce.amy, colour_by="sampleID", point_alpha=0.5)
plotUMAP(sce.amy, colour_by="collapsedCluster", point_alpha=0.5)
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
  # Added post-hoc, looking at markers (step3) - MNT 25May2021
  'Tcell' = c('SKAP1', 'ITK', 'CD247')
)

pdf("pdfs/revision/regionSpecific_Amyg-n5_marker-logExprs_collapsedClusters_MNT2021.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.amy, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]), ncol=2,
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
newClusIndex <- splitit(sce.amy$collapsedCluster)
sapply(newClusIndex, function(x) {quantile(sce.amy[,x]$sum)})
    #            1      2      3       4         5        6         7       8
    # 0%     2876.0    831   1685  2219.0   7358.00   566.00   1238.00  1182.0
    # 25%   29547.5  22550   9983 12195.5  37322.50  7653.50  20472.75  7471.5
    # 50%   51398.0  34143  18481 36275.0  55357.00 12241.00  45897.50 11275.0
    # 75%   67826.0  43717  39510 49672.5  73867.75 18813.75  67085.00 16469.5
    # 100% 156945.0 115493 121273 94599.0 165583.00 82358.00 171608.00 35728.0
    #            9      10    11       12     13       14      15       16      17
    # 0%     641.0  1209.0   435  4036.00  121.0   756.00  7739.0  3346.00  1711.0
    # 25%   4634.5  9787.5  3952  8556.00  245.0  4527.75 16961.0 11318.50  3453.0
    # 50%   6510.5 13243.0  5381 10896.00  313.0  7092.00 20019.0 14411.50  4584.0
    # 75%   8964.0 16967.5  7055 13478.25  468.5  9291.50 30862.5 26971.75  5476.5
    # 100% 30028.0 39830.0 20500 27167.00 2770.0 95032.00 58639.0 92477.00 15214.0
    #           18     19     20     21
    # 0%    1927.0  101.0  438.0  188.0
    # 25%   4826.5  373.5 2382.0  850.5
    # 50%   5850.0  743.0 3657.0 1156.0
    # 75%   7433.5 1018.0 5029.5 1842.5
    # 100% 11453.0 1556.0 6237.0 7739.0


    ## Obser.: looks like collapsedCluster's 13 & 19 should be dropped - their marker
    ##         expression is noisy too


table(sce.amy$collapsedCluster)
    #   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
    # 728  541  525  555  344  414  216 1555 6080 1459 1201   44 1067   86   55   52 
    #  17   18   19   20   21 
    #  31   39   71   31   83
        #  * 13 is quite overrepresented in Br5400, but this sample also had a median
        #    of <1000 total UMIs that were still kept after QC - thus probably a low trxts-driven cluster
        #  * 19 is pretty evenly-ish distributed

## doublet score?
sapply(newClusIndex, function(x) {round(quantile(sce.amy$doubletScore[x]),2)})
    #        1    2     3    4    5    6    7     8    9    10    11   12   13   14
    # 0%   0.22 0.15  0.16 0.11 0.01 0.13 0.20  0.00 0.00  0.00  0.00 0.66 0.00 0.16
    # 25%  0.55 0.37  0.36 0.31 0.10 0.28 0.97  0.07 0.16  0.07  0.03 1.02 0.02 0.28
    # 50%  0.84 0.84  0.69 0.76 0.16 0.40 1.40  0.16 0.39  0.19  0.08 1.19 0.06 0.36
    # 75%  1.02 1.09  1.12 1.00 0.33 0.64 1.54  0.29 0.86  0.35  0.19 1.31 0.33 0.68
    # 100% 4.94 6.94 17.29 3.33 5.61 3.40 7.24 12.15 7.72 13.11 10.28 9.84 2.45 1.19
    # 1       5   16   17   18   19   20   21
    # 0%   0.01 0.28 0.17 0.20 0.00 0.01 0.01
    # 25%  0.04 0.40 0.26 0.33 0.01 0.05 0.04
    # 50%  0.05 0.52 0.33 0.52 0.02 0.07 0.16
    # 75%  0.08 0.63 0.49 0.64 0.11 0.12 0.21
    # 100% 1.18 2.11 1.01 0.83 1.00 0.85 0.77

# At 'prelimCluster' level?
clusIndex <- splitit(sce.amy$prelimCluster)
sapply(clusIndex, function(x) {round(quantile(sce.amy$doubletScore[x]),2)})
    # 47 for sure (maybe 25 too) - 55 relatively high
    #     - none of these are single-donor-specific though (where doublets would arise)


## Add annotations, looking at marker gene expression
#    (canonical, above, and from markers looked at in 'step3')
annotationTab.amy <- data.frame(collapsedCluster=c(1:21))
annotationTab.amy$cellType <- NA
annotationTab.amy$cellType[c(1:4,6:7,14,16)] <- paste0("Inhib_", c("A","B","C","D","E","F","G","H"))
annotationTab.amy$cellType[c(5,12,15)] <- paste0("Excit_", c("A","B","C"))
annotationTab.amy$cellType[c(8:11)] <- c("Astro_A", "Oligo", "OPC", "Micro")
annotationTab.amy$cellType[c(13,19)] <- paste0("drop.lowNTx_", c("A","B"))
annotationTab.amy$cellType[c(17:18,20:21)] <- c("Endo","Mural", "Tcell", "Astro_B")


sce.amy$cellType <- annotationTab.amy$cellType[match(sce.amy$collapsedCluster,
                                                         annotationTab.amy$collapsedCluster)]
sce.amy$cellType <- factor(sce.amy$cellType)

cell_colors.amy <- cluster_colors[order(as.integer(names(cluster_colors)))]
names(cell_colors.amy) <- annotationTab.amy$cellType
cell_colors.amy
    #       Inhib_A       Inhib_B       Inhib_C       Inhib_D       Excit_A       Inhib_E 
    #     "#1F77B4"     "#AEC7E8"     "#FF7F0E"     "#FFBB78"     "#2CA02C"     "#98DF8A" 
    #       Inhib_F       Astro_A         Oligo           OPC         Micro       Excit_B 
    #     "#D62728"     "#FF9896"     "#9467BD"     "#C5B0D5"     "#8C564B"     "#C49C94" 
    # drop.lowNTx_A       Inhib_G       Excit_C       Inhib_H          Endo         Mural 
    #     "#E377C2"     "#F7B6D2"     "#7F7F7F"     "#C7C7C7"     "#BCBD22"     "#DBDB8D" 
    # drop.lowNTx_B         Tcell       Astro_B 
    #     "#17BECF"     "#9EDAE5"     "#729ECE" 

# Save
save(sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo, annotationTab.amy, cell_colors.amy,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_Amyg-n5_cleaned-combined_SCE_MNT2021.rda")



## Re-print marker expression plots with annotated cluster names ===
pdf("pdfs/revision/regionSpecific_Amyg-n5_marker-logExprs_collapsedClusters_MNT2021.pdf", height=5, width=9)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpressionCustom(sce = sce.amy,
                         features = markers.mathys.custom[[i]], 
                         features_name = names(markers.mathys.custom)[[i]], 
                         anno_name = "cellType") +
      scale_color_manual(values = cell_colors.amy)
  )
}
dev.off()

    ## Optionally, for a cleaner version, drop those 'drop.lowNTx_'s and re-print
    sce.amy <- sce.amy[ ,-grep("drop.", sce.amy$cellType)]
    sce.amy$cellType <- droplevels(sce.amy$cellType)


    
## Re-print reducedDims with these annotations (keep the 'drop.' clusters here) ===
pdf("pdfs/revision/regionSpecific_Amyg-n5_reducedDims-with-collapsedClusters_MNT2021.pdf",width=8)
plotReducedDim(sce.amy, dimred="PCA_corrected", ncomponents=5, colour_by="cellType", point_alpha=0.5) +
  scale_color_manual(values = cell_colors.amy) + labs(colour="Cell type")
plotTSNE(sce.amy, colour_by="sampleID", point_alpha=0.5, point_size=2)
plotTSNE(sce.amy, colour_by="protocol", point_alpha=0.5, point_size=2)
plotTSNE(sce.amy, colour_by="prelimCluster", text_by="prelimCluster",
         text_size=3, point_alpha=0.5, point_size=2)
plotTSNE(sce.amy, colour_by="cellType", text_by="cellType",
         text_size=3, point_alpha=0.5, point_size=2) +
  scale_color_manual(values = cell_colors.amy,
                     labels=paste0(levels(sce.amy$cellType)," (",table(sce.amy$cellType),")")) +
  labs(colour="Cell type")
plotTSNE(sce.amy, colour_by="sum", point_alpha=0.5, point_size=2)
plotTSNE(sce.amy, colour_by="doubletScore", point_alpha=0.5, point_size=2)
# And some more informative UMAPs
plotUMAP(sce.amy, colour_by="prelimCluster", text_by="prelimCluster",
         text_size=3, point_alpha=0.5, point_size=2)
plotUMAP(sce.amy, colour_by="cellType", text_by="cellType",
         text_size=3, point_alpha=0.5, point_size=2) +
  scale_color_manual(values = cell_colors.amy,
                     labels=paste0(levels(sce.amy$cellType)," (",table(sce.amy$cellType),")")) +
  labs(colour="Cell type")
dev.off()



      ## -> proceed to 'step03_markerDetxn-analyses[...].R'



## Some diggings-into for paper =======
# BLA or central amygdala-specific (Cartpt in all subregions... in mouse at least)
plotExpression(sce.amy, exprs_values = "logcounts", features=toupper(c("Cyp26b1", "Bmp3", "Cartpt")),
               x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7,
               add_legend=F) + stat_summary(fun = median, fun.min = median, fun.max = median,
                                            geom = "crossbar", width = 0.3,
                                            colour=rep(tableau20[1:12], 3)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

# Split by sample
sce.amy$sample.cellType.split <- paste0(ss(sce.amy$sample,"\\.",2), "_", sce.amy$cellType.split)
plotExpression(sce.amy, exprs_values = "logcounts", features=toupper(c("Cyp26b1", "Bmp3")),
               x="sample.cellType.split", colour_by="sample.cellType.split", point_alpha=0.5, point_size=.7,
               add_legend=F, show_median=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
    ## Looks like some BLA

# Actually faceting is nicer and maintains colors but have to print genes, separately
# From bulk-RNA-seq (AnJa shared): Look at PART1 (BLA >) & GABRQ/SYTL5 (MeA >)
pdf("pdfs/revision/exploration/zGenesByDonor_PART1-GABRQ-SYTL5_AMY-subclusters_MNT.pdf", height=3, width=7)
# BLA-enriched
plotExpression(sce.amy, exprs_values = "logcounts", features="PART1",
               x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7,
               add_legend=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7)) + facet_grid(~ sce.amy$donor) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3,
               colour=rep(tableau20[1:12][c(1:2,5:12, 1:7,9:12)])) 
# MeA-enriched
plotExpression(sce.amy, exprs_values = "logcounts", features="GABRQ",
               x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7,
               add_legend=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7)) + facet_grid(~ sce.amy$donor) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3,
               colour=rep(tableau20[1:12][c(1:2,5:12, 1:7,9:12)])) 
plotExpression(sce.amy, exprs_values = "logcounts", features="SYTL5",
               x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7,
               add_legend=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7)) + facet_grid(~ sce.amy$donor) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3,
               colour=rep(tableau20[1:12][c(1:2,5:12, 1:7,9:12)])) 
dev.off()





## Reference chunk: If wanting to manually edit 'text_by=' coordinates =======
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
text_out <- retrieveCellInfo(sce.amy.tsne.optb, text_by, search="colData")
text_out$val <- .coerce_to_factor(text_out$val, level.limit=Inf)
## actually not necessary if the colData chosen (usually cellType[.etc] is factorized)
df_to_plot <- data.frame(reducedDim(sce.amy.tsne.optb, "TSNE"))
by_text_x <- vapply(split(df_to_plot$X1, text_out$val), median, FUN.VALUE=0)
by_text_y <- vapply(split(df_to_plot$X2, text_out$val), median, FUN.VALUE=0)
# plot_out <- plot_out + annotate("text", x=by_text_x, y=by_text_y, 
#                             label=names(by_text_x), size=text_size, colour=text_colour)



plotTSNE(sce.amy.tsne.optb, colour_by="cellType.split", point_size=4.5, point_alpha=0.5,
         text_size=8, theme_size=18) +
  annotate("text", x=by_text_x, y=by_text_y, 
           label=names(by_text_x), size=6) +
  ggtitle("t-SNE on optimal PCs (45), donor-correlated PC[5] removed")

# OR
sce.amy.tsne.optb$labels <- ifelse(!duplicated(sce.amy.tsne.optb$cellType.split), as.character(sce.amy.tsne.optb$cellType.split), NA)
Labs.df <- data.frame(by_text_x, by_text_y, labs=names(by_text_x))

colDF <- data.frame(colData(sce.amy.tsne.optb))
DFforLabs <- cbind(reducedDim(sce.amy.tsne.optb,"TSNE"), data.frame(colDF$labels))
colnames(DFforLabs) <- c("X","Y","labels")

# plotTSNE(sce.amy.tsne.optb, colour_by="cellType.split", point_size=4.5, point_alpha=0.5,
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

#pdf("pdfs/revision/pubFigures/FINAL_pilotPaper_Amyg-n2_tSNE_optPCs-PC5_MNTMay2020.pdf", width=8)
set.seed(109)
plotTSNE(sce.amy.tsne.optb, colour_by="cellType.split", point_size=6, point_alpha=0.5,
         theme_size=18) +
  geom_text_repel(data=DFforLabs.edit, size=6.0,
                  aes(label=labels)) +
  ggtitle("t-SNE on optimal PCs (45), donor-correlated PC[5] removed")
#dev.off()

### end reference chunk ========








apply(table(sce.amy$cellType, sce.amy$donor),2,function(x){round(prop.table(x),3)})
    #          br5161 br5212 br5276 br5400 br5701
    # Astro_A  0.148  0.107  0.102  0.064  0.108
    # Astro_B  0.002  0.003  0.022  0.007  0.001
    # Endo     0.000  0.000  0.009  0.004  0.001
    # Excit_A  0.032  0.062  0.002  0.008  0.005
    # Excit_B  0.000  0.012  0.000  0.000  0.001
    # Excit_C  0.002  0.013  0.000  0.004  0.000
    # Inhib_A  0.000  0.000  0.162  0.208  0.000
    # Inhib_B  0.011  0.035  0.031  0.141  0.021
    # Inhib_C  0.039  0.005  0.126  0.049  0.003
    # Inhib_D  0.011  0.023  0.055  0.156  0.014
    # Inhib_E  0.000  0.000  0.179  0.004  0.001
    # Inhib_F  0.007  0.021  0.016  0.046  0.002
    # Inhib_G  0.000  0.000  0.034  0.005  0.000
    # Inhib_H  0.000  0.000  0.022  0.001  0.000
    # Micro    0.126  0.093  0.006  0.067  0.101
    # Mural    0.001  0.000  0.011  0.004  0.002
    # Oligo    0.516  0.533  0.134  0.177  0.582
    # OPC      0.104  0.089  0.088  0.053  0.153
    # Tcell    0.001  0.002  0.001  0.002  0.004

round(apply(apply(table(sce.amy$cellType, sce.amy$donor),2,prop.table),1,mean),3)
    #Astro_A Astro_B    Endo Excit_A Excit_B Excit_C Inhib_A Inhib_B Inhib_C Inhib_D 
    #  0.106   0.007   0.003   0.022   0.003   0.004   0.074   0.048   0.044   0.052 
    #Inhib_E Inhib_F Inhib_G Inhib_H   Micro   Mural   Oligo     OPC   Tcell 
    #  0.037   0.019   0.008   0.005   0.079   0.003   0.389   0.097   0.002



### Session info for 03Jun2021 ==============================
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

    