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
load("/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/all-FACS-n14_preprint_SCEs_processing-QC_MNTMar2021.rda",
     verbose=T)
    # pilot.data, pilot.data.unfiltered, e.out, ref.sampleInfo
    rm(pilot.data.unfiltered, e.out)

# Load 2021 expansion set
load("/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/all-FACS-n10_2021rev_SCEs_processing-QC_MNTMar2021.rda", verbose=T)
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


    ## 'getClusteredPCs()' evaluated in qsub mode (with 'R-batchJob_sACC-n5_optimalPCselxn_MNT2021.R')
    #    --> saved into same .rda


### Picking up with optimally-defined PC space ===
load("/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_sACC-n5_cleaned-combined_SCE_MNT2021.rda",
     verbose=TRUE)
    # sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, ref.sampleInfo, ref.sampleInfo.rev


# How many PCs is optimal?:
metadata(pc.choice.sacc)$chosen # [1] 97

## Assign this chosen (97 PCs) to 'PCA_opt'
reducedDim(sce.sacc, "PCA_opt") <- reducedDim(sce.sacc, "PCA_corrected")[ ,1:(metadata(pc.choice.sacc)$chosen)]


## t-SNE
set.seed(109)
sce.sacc <- runTSNE(sce.sacc, dimred="PCA_opt")


## UMAP
set.seed(109)
sce.sacc <- runUMAP(sce.sacc, dimred="PCA_opt")


# How do these look?
plotReducedDim(sce.sacc, dimred="TSNE", colour_by="sampleID")
plotReducedDim(sce.sacc, dimred="UMAP", colour_by="sampleID")


### Clustering: Two-step ======================================================
### Step 1: Perform graph-based clustering in this optimal PC space
#         - take k=20 NN to build graph
snn.gr <- buildSNNGraph(sce.sacc, k=20, use.dimred="PCA_opt")
clusters.k20 <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters.k20)
    ## 67 prelim clusters

# Assign as 'prelimCluster'
sce.sacc$prelimCluster <- factor(clusters.k20)
plotReducedDim(sce.sacc, dimred="TSNE", colour_by="prelimCluster")

# Is sample driving this 'high-res' clustering at this level?
table(sce.sacc$prelimCluster, sce.sacc$sampleID)  # (not much!)

# rbind the ref.sampleInfo[.rev]
ref.sampleInfo <- rbind(ref.sampleInfo, ref.sampleInfo.rev)

## Save for now
save(sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, ref.sampleInfo, 
     file="/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_sACC-n5_cleaned-combined_SCE_MNT2021.rda")



### Step 2: Hierarchical clustering of pseudo-bulked ("PB'd") counts
#           ** That is, to pseudo-bulk (aka 'cluster-bulk') on raw counts, on all [non-zero] genes,
#              normalize with `librarySizeFactors()`, log2-transform, then perform HC'ing


# Preliminary cluster index for pseudo-bulking
clusIndexes = splitit(sce.sacc$prelimCluster)
prelimCluster.PBcounts <- sapply(clusIndexes, function(ii){
  rowSums(assays(sce.sacc)$counts[ ,ii])
  }
)

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
myplclust(tree.clusCollapsed, main="5x sACC prelim-kNN-cluster relationships", cex.main=2, cex.lab=1.5, cex=1.4)


clust.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                               minClusterSize=2, deepSplit=1, cutHeight=225)


table(clust.treeCut)
unname(clust.treeCut[order.dendrogram(dend)])
    #  [1]  0  0 16 16  0  0  0  0  0 12 12  0 15 15  4  4  4  4  9  9  9  8  8  8 13
    # [26] 13  0  5  5  5  5  0  2  2  2  2  2  0 10 10 10  3  3  3  3  3 14 14  6  6
    # [51]  6  6  0  1  1  1  1  1  1  0 11 11 11  7  7  7  7

        ## Cutting at 225 looks good for the main neuronal branch, but a lot of glial
         # prelim clusters are dropped off (0's)
    

## Cut at 325 shows that prelimClust's 32, 41 should be merged - the rest will leave, due to more
 # extreme height differences:    
    glia.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                                   minClusterSize=2, deepSplit=1, cutHeight=325)
    unname(glia.treeCut[order.dendrogram(dend)])
        #  [1]  9  9  9  9  0  0 12 12 10 10 10  0  4  4  4  4  4  4  2  2  2  2  2  2  2
        # [26]  2  7  7  7  7  7  0  6  6  6  6  6  1  1  1  1  1  1  1  1  1  5  5  5  5
        # [51]  5  5  3  3  3  3  3  3  3  8 11 11 11  8  8  8  8
    
    clust.treeCut[order.dendrogram(dend)][c(7:8)] <- max(clust.treeCut)+c(1)  # '17'


# Add new labels to rest of the prelimClusters cut off
clust.treeCut[order.dendrogram(dend)][which(
  clust.treeCut[order.dendrogram(dend)]==0)] <- max(clust.treeCut) + c(1,2, 3,3, 4:10)

# ( 'Re-write', in case there are missing numbers )
clust.treeCut[order.dendrogram(dend)] <- as.numeric(as.factor(clust.treeCut[order.dendrogram(dend)]))

## Define color pallet
cluster_colors <- unique(c(tableau10medium, tableau20)[clust.treeCut[order.dendrogram(dend)]])
names(cluster_colors) <- unique(clust.treeCut[order.dendrogram(dend)])
labels_colors(dend) <- cluster_colors[as.character(clust.treeCut[order.dendrogram(dend)])]

# Print for future reference
pdf("pdfs/revision/regionSpecific_sACC-n5_HC-prelimCluster-relationships_MNT2021.pdf", width=9)
par(cex=0.8, font=2)
plot(dend, main="5x sACC prelim-kNN-cluster relationships with collapsed assignments")
dev.off()


# Make reference for new cluster assignment
clusterRefTab.sacc <- data.frame(origClust=order.dendrogram(dend),
                                merged=clust.treeCut[order.dendrogram(dend)])

# Assign as 'collapsedCluster'
sce.sacc$collapsedCluster <- factor(clusterRefTab.sacc$merged[match(sce.sacc$prelimCluster, clusterRefTab.sacc$origClust)])


# Print some visualizations:
pdf("pdfs/revision/regionSpecific_sACC-n5_reducedDims-with-collapsedClusters_MNT2021.pdf")
plotReducedDim(sce.sacc, dimred="PCA_corrected", ncomponents=5, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.sacc, colour_by="sampleID", point_alpha=0.5)
plotTSNE(sce.sacc, colour_by="protocol", point_alpha=0.5)
plotTSNE(sce.sacc, colour_by="collapsedCluster", text_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.sacc, colour_by="sum", point_alpha=0.5)
plotTSNE(sce.sacc, colour_by="doubletScore", point_alpha=0.5)
# And some more informative UMAPs
plotUMAP(sce.sacc, colour_by="sampleID", point_alpha=0.5)
plotUMAP(sce.sacc, colour_by="collapsedCluster", text_by="collapsedCluster", point_alpha=0.5)
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
  # T cell markers - these seem to be somewhat common (across regions)
  'Tcell' = c('SKAP1', 'ITK', 'CD247')
)

pdf("pdfs/revision/regionSpecific_sACC-n5_marker-logExprs_collapsedClusters_MNT2021.pdf", height=6, width=10)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.sacc, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]), ncol=2,
                   x="collapsedCluster", colour_by="collapsedCluster", point_alpha=0.3, point_size=.7,
                   add_legend=F, show_median=T) +
      # (since there are more than 20 collapsed clusters)
      # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
      #              geom = "crossbar", width = 0.3,
      #              colour=rep(tableau20[1:11], length(markers.mathys.custom[[i]])))
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(markers.mathys.custom)[i], " markers"))
  )
}
dev.off()



## QC - How do the total UMI distribution look?
newClusIndex <- splitit(sce.sacc$collapsedCluster)
sapply(newClusIndex, function(x) {quantile(sce.sacc[,x]$sum)})
    #             1        2     3         4        5       6        7        8        9    10       11
    # 0%    2909.00   2448.0  3360   2101.00   2390.0  2465.0  2313.00   2754.0   2984.0  2728  1472.00
    # 25%  16433.75  22017.5 13669  26971.75  25059.5 20242.0 14973.75  30164.5  29526.5 13140 20121.75
    # 50%  21474.50  30746.0 20617  40460.00  39576.0 26094.5 20633.00  43945.0  48494.0 19108 29225.00
    # 75%  27697.50  41323.5 27494  56412.75  56040.5 32440.5 27588.75  60824.0  67895.0 25248 38971.00
    # 100% 63006.00 121477.0 93363 127499.00 123619.0 69549.0 67423.00 172041.0 196431.0 73742 84097.00
    #         12        13       14       15      16      17       18       19       20     21     *22
    # 0%    1218   2156.00  2412.00  3238.00   174.0  1276.0  4614.00   102.00   113.00  111.0   103.0
    # 25%   4046  30561.25 22631.00 23977.00  5781.0  7510.5  9990.25  1977.25  3198.50 1540.5   215.0
    # 50%   5643  43868.50 32647.00 32559.00  7575.0 10068.0 12587.50  3872.50  4378.00 1706.0   321.0
    # 75%   7539  58960.75 40748.75 39867.75 10077.5 12682.5 17470.25  6176.25  5458.25 1926.5   681.5
    # 100% 25828 113463.00 73499.00 70217.00 23974.0 34023.0 27467.00 14206.00 11896.00 2649.0 14701.0
    #             23       24      25       26    27
    # 0%     3436.00  2333.00  2740.0  2874.00 12801
    # 25%   19197.00  5914.25 21647.0 13851.50 23846
    # 50%   47066.50  7686.00 33647.0 17916.50 27375
    # 75%   58822.25  9597.50 45194.5 26535.75 30262
    # 100% 102363.00 13985.00 69890.0 40784.00 48721

## doublet score?
sapply(newClusIndex, function(x) {round(quantile(sce.sacc$doubletScore[x]),2)})
    #         1     2    3     4    5    6    7     8    9    10    11   12    13   14   15   16   17
    # 0%   0.02  0.00 0.01  0.01 0.00 0.00 0.05  0.00 0.04  0.02  0.04 0.00  0.02 0.02 0.03 0.00 0.00
    # 25%  0.11  0.08 0.05  0.08 0.17 0.08 0.12  0.11 0.28  0.08  0.28 0.03  0.18 0.06 0.17 0.01 0.02
    # 50%  0.18  0.16 0.08  0.17 0.99 0.13 0.23  0.23 0.51  0.13  0.46 0.13  0.31 0.28 0.29 0.05 0.06
    # 75%  0.31  0.30 0.18  0.39 1.71 0.48 0.54  0.41 0.83  0.21  0.70 0.45  0.54 0.98 0.47 0.12 0.13
    # 100% 8.45 14.34 9.20 10.08 3.66 3.81 5.88 23.05 5.30 11.64 17.43 7.19 12.44 4.31 5.08 2.00 5.66
    
    #        *18   19   20   21   22    23   24   25   26   27
    # 0%    2.37 0.00 0.00 0.00 0.01  0.05 0.39 0.02 0.04 0.10
    # 25%   3.73 0.01 0.00 0.00 0.05  1.24 0.44 0.18 0.10 0.21
    # 50%   5.34 0.01 0.02 0.01 0.12  2.09 0.71 0.39 0.14 0.28
    # 75%   6.43 0.04 0.05 0.02 0.36  2.41 1.08 0.74 0.21 1.05
    # 100% 20.33 3.21 5.84 0.30 1.68 13.00 1.20 1.59 4.74 2.00


# At 'prelimCluster' level?
clusIndex <- splitit(sce.sacc$prelimCluster)
sapply(clusIndex, function(x) {round(quantile(sce.sacc$doubletScore[x]),2)})
    # 62 [moderately, but with median > 5], which is == collapsedCluster 18 (see above); 32, 53 (meds 4.1, 3.9)

table(sce.sacc$collapsedCluster)
    #   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
    # 842  912  465  856  575  384  330 1735  311  521  206 4389  428  208  228  747 
    #  17   18   19   20   21   22   23   24   25   26   27 
    # 911   28  160  784  195  298   30   20   39   42   25 
    

## Add annotations, looking at marker gene expression
annotationTab.sacc <- data.frame(collapsedCluster=c(1:27))
annotationTab.sacc$cellType <- NA
annotationTab.sacc$cellType[c(4:5,8:9,13,15,23)] <- paste0("Excit_", c("A","B","C","D","E","F","G"))
annotationTab.sacc$cellType[c(1:3,6:7,10:11,14,25:27)] <- paste0("Inhib_",
                                                                 c("A","B","C","D","E","F","G","H","I","J","K"))

annotationTab.sacc$cellType[c(12,21)] <- paste0("Oligo_", c("A","B"))
annotationTab.sacc$cellType[17] <- c("OPC")
annotationTab.sacc$cellType[18] <- "drop.doublet"  # ahhh, these are doublets (see above)
    # they express both astro & oligo markers - btw it's prelimCluster 62
annotationTab.sacc$cellType[c(16,19)] <- paste0("Astro_", c("A","B"))
annotationTab.sacc$cellType[20] <- "Micro"
annotationTab.sacc$cellType[22] <- "drop.lowNTx"
annotationTab.sacc$cellType[24] <- "Neu_FAT2.CDH15"
    # perhaps this one is axonal debris??  It was prelimCluster 67 (stuck out a ton)


# Then add these to the SCE
sce.sacc$cellType <- annotationTab.sacc$cellType[match(sce.sacc$collapsedCluster,
                                                       annotationTab.sacc$collapsedCluster)]
sce.sacc$cellType <- factor(sce.sacc$cellType)

# Set colors to be stable for this region/these cell classes
cell_colors.sacc <- cluster_colors[order(as.integer(names(cluster_colors)))]
names(cell_colors.sacc) <- annotationTab.sacc$cellType
cell_colors.sacc
    #   Inhib_A        Inhib_B        Inhib_C        Excit_A        Excit_B        Inhib_D 
    # "#729ECE"      "#FF9E4A"      "#67BF5C"      "#ED665D"      "#AD8BC9"      "#A8786E" 
    #   Inhib_E        Excit_C        Excit_D        Inhib_F        Inhib_G        Oligo_A 
    # "#ED97CA"      "#A2A2A2"      "#CDCC5D"      "#6DCCDA"      "#1F77B4"      "#AEC7E8" 
    #   Excit_E        Inhib_H        Excit_F        Astro_A            OPC   drop.doublet 
    # "#FF7F0E"      "#FFBB78"      "#2CA02C"      "#98DF8A"      "#D62728"      "#FF9896" 
    #   Astro_B          Micro        Oligo_B    drop.lowNTx        Excit_G Neu_FAT2.CDH15 
    # "#9467BD"      "#C5B0D5"      "#8C564B"      "#C49C94"      "#E377C2"      "#F7B6D2" 
    #   Inhib_I        Inhib_J        Inhib_K 
    # "#7F7F7F"      "#C7C7C7"      "#BCBD22"


## Save
save(sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc, ref.sampleInfo, annotationTab.sacc, cell_colors.sacc,
     file="/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_sACC-n5_cleaned-combined_SCE_MNT2021.rda")



## Re-print marker expression with cell type labels ===
# load("/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_sACC-n5_cleaned-combined_SCE_MNT2021.rda",
#      verbose=T)

pdf("pdfs/revision/regionSpecific_sACC-n5_marker-logExprs_collapsedClusters_MNT2021.pdf", height=4.5, width=10)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpressionCustom(sce = sce.sacc,
                         features = markers.mathys.custom[[i]], 
                         features_name = names(markers.mathys.custom)[[i]], 
                         anno_name = "cellType") +
      scale_color_manual(values = cell_colors.sacc)
  )
}
dev.off()

    # Optionally, drop doublet cluster & 'drop.lowNTx' and re-print
    sce.sacc <- sce.sacc[ ,-grep("drop.",sce.sacc$cellType)]
    sce.sacc$cellType <- droplevels(sce.sacc$cellType)

    # Final comment - the Micro cluster--its prelim clusters were checked for endothelial
    #                 markers, but doesn't seem they were merged together.  Perhaps this dataset
    #                 just didn't capture any (unlike for the AMY)
    

    
## Re-print reducedDims with these annotations ===
pdf("pdfs/revision/regionSpecific_sACC-n5_reducedDims-with-collapsedClusters_MNT2021.pdf", width=8)
plotReducedDim(sce.sacc, dimred="PCA_corrected", ncomponents=5, colour_by="cellType", point_alpha=0.5) +
  scale_color_manual(values = cell_colors.sacc) + labs(colour="Cell type")
plotTSNE(sce.sacc, colour_by="sampleID", point_alpha=0.5, point_size=2)
plotTSNE(sce.sacc, colour_by="protocol", point_alpha=0.5, point_size=2)
plotTSNE(sce.sacc, colour_by="prelimCluster", text_by="prelimCluster",
         text_size=3, point_alpha=0.5, point_size=2)
plotTSNE(sce.sacc, colour_by="cellType", text_by="cellType",
         text_size=3, point_alpha=0.5, point_size=2) +
  scale_color_manual(values = cell_colors.sacc,
                     labels=paste0(levels(sce.sacc$cellType)," (",table(sce.sacc$cellType),")")) +
  labs(colour="Cell type")
plotTSNE(sce.sacc, colour_by="sum", point_alpha=0.5, point_size=2)
plotTSNE(sce.sacc, colour_by="doubletScore", point_alpha=0.5, point_size=2)
# And some more informative UMAPs
plotUMAP(sce.sacc, colour_by="prelimCluster", text_by="prelimCluster",
         text_size=3, point_alpha=0.5, point_size=2)
plotUMAP(sce.sacc, colour_by="cellType", text_by="cellType",
         text_size=3, point_alpha=0.5, point_size=2) +
  scale_color_manual(values = cell_colors.sacc,
                     labels=paste0(levels(sce.sacc$cellType)," (",table(sce.sacc$cellType),")")) +
  labs(colour="Cell type")
dev.off()



  
      ## -> proceed to 'step03_markerDetxn-analyses[...].R'



### Session info for 06Jun2021 ==========================================================
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
# [10] bit64_4.0.5               fansi_0.4.2               xml2_1.3.2               
# [13] splines_4.0.4             R.methodsS3_1.8.1         sparseMatrixStats_1.2.1  
# [16] cachem_1.0.4              Rsamtools_2.6.0           ResidualMatrix_1.0.0     
# [19] dbplyr_2.1.1              R.oo_1.24.0               HDF5Array_1.18.1         
# [22] compiler_4.0.4            httr_1.4.2                dqrng_0.2.1              
# [25] assertthat_0.2.1          Matrix_1.3-2              fastmap_1.1.0            
# [28] lazyeval_0.2.2            limma_3.46.0              BiocSingular_1.6.0       
# [31] prettyunits_1.1.1         tools_4.0.4               rsvd_1.0.3               
# [34] igraph_1.2.6              gtable_0.3.0              glue_1.4.2               
# [37] GenomeInfoDbData_1.2.4    dplyr_1.0.5               rappdirs_0.3.3           
# [40] Rcpp_1.0.6                vctrs_0.3.6               Biostrings_2.58.0        
# [43] rhdf5filters_1.2.0        rtracklayer_1.50.0        DelayedMatrixStats_1.12.3
# [46] stringr_1.4.0             beachmat_2.6.4            lifecycle_1.0.0          
# [49] irlba_2.3.3               statmod_1.4.35            XML_3.99-0.6             
# [52] edgeR_3.32.1              zlibbioc_1.36.0           scales_1.1.1             
# [55] hms_1.0.0                 ProtGenerics_1.22.0       rhdf5_2.34.0             
# [58] RColorBrewer_1.1-2        curl_4.3                  memoise_2.0.0            
# [61] gridExtra_2.3             segmented_1.3-3           biomaRt_2.46.3           
# [64] stringi_1.5.3             RSQLite_2.2.7             BiocParallel_1.24.1      
# [67] rlang_0.4.10              pkgconfig_2.0.3           bitops_1.0-7             
# [70] lattice_0.20-41           purrr_0.3.4               Rhdf5lib_1.12.1          
# [73] GenomicAlignments_1.26.0  bit_4.0.4                 tidyselect_1.1.1         
# [76] magrittr_2.0.1            R6_2.5.0                  generics_0.1.0           
# [79] DelayedArray_0.16.3       DBI_1.1.1                 pillar_1.6.0             
# [82] withr_2.4.2               RCurl_1.98-1.3            tibble_3.1.1             
# [85] crayon_1.4.1              utf8_1.2.1                BiocFileCache_1.14.0     
# [88] viridis_0.6.0             progress_1.2.2            locfit_1.5-9.4           
# [91] grid_4.0.4                blob_1.2.1                R.utils_2.10.1           
# [94] openssl_1.4.3             munsell_0.5.0             beeswarm_0.3.1           
# [97] viridisLite_0.4.0         vipor_0.4.5               askpass_1.1 


