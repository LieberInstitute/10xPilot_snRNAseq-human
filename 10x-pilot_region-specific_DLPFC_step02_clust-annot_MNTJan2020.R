### MNT 10x snRNA-seq workflow: step 02
###   **Region-specific analyses**
###     - (2x) DLPFC samples from: Br5161 & Br5212
### Initiated MNT 29Jan2020
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


load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_processing-QC_MNTJan2020.rda",
     verbose=T)
    # pilot.data, pilot.data.unfiltered, e.out

### MNT comment: At this point, each sample (which is a SCE object in the list, 'pilot.data') has been
#              QC'd for cell/nucleus calling ('emptyDrops()' test) and mito rate thresholding


### Merging shared-region samples ============================================
# Take those DLPFC samples
pilot.dlpfc <- list(pilot.data[["dlpfc.5161"]],
                  pilot.data[["dlpfc.5212"]])
names(pilot.dlpfc) <- c("dlpfc.5161","dlpfc.5212")

### Newest iterations for normalization: cbind, THEN take scaled LSFs computed on all nuclei
# Add $sample identity
for(i in 1:length(pilot.dlpfc)){
  pilot.dlpfc[[i]]$sample <- names(pilot.dlpfc)[i]
}

sce.dlpfc <- cbind(pilot.dlpfc[[1]], pilot.dlpfc[[2]])

# Remove $logcounts
assay(sce.dlpfc, "logcounts") <- NULL
# Re-generate log-normalized counts
sce.dlpfc <- logNormCounts(sce.dlpfc)

geneVar.dlpfc <- modelGeneVar(sce.dlpfc)
chosen.hvgs.dlpfc <- geneVar.dlpfc$bio > 0
sum(chosen.hvgs.dlpfc)
    # [1] 9313


### Dimensionality reduction ================================================================

# Run PCA, taking top 100 (instead of default 50 PCs)
set.seed(109)
sce.dlpfc <- runPCA(sce.dlpfc, subset_row=chosen.hvgs.dlpfc, ncomponents=100,
                  BSPARAM=BiocSingular::RandomParam())

# Save into a new data file, which will dedicate for pan-brain-analyses
save(sce.dlpfc, chosen.hvgs.dlpfc, file="rdas/regionSpecific_DLPFC-n2_cleaned-combined_SCE_MNTFeb2020.rda")


## 'getClusteredPCs()' evaluated in qsub mode (with 'R-batchJob_DLPFC-n2_optimalPCselxn_MNTFeb2020.R')
#    --> saved into same .rda


### Picking up with optimally-defined PC space ===
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_DLPFC-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=TRUE)
    # sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc


# How many PCs is optimal?:
metadata(pc.choice.dlpfc)$chosen
    # [1] 92

## Assign this chosen ( PCs) to 'PCA_opt'
#reducedDim(sce.dlpfc, "PCA_opt") <- reducedDim(sce.dlpfc, "PCA")[ ,1:(metadata(pc.choice.dlpfc)$chosen)]
reducedDim(sce.dlpfc, "PCA_opt") <- reducedDim(sce.dlpfc, "PCA")[ ,c(1:92)]
    ## 21Feb2020 comment: this is what was taken the first time; data accidentally written over, so
     #                    though a different 'metadata(pc.choice.dlpfc)$chosen' the second time,
     #                    keeping this initial choice to recapitulate work


## t-SNE
set.seed(109)
sce.dlpfc <- runTSNE(sce.dlpfc, dimred="PCA_opt")


## UMAP
set.seed(109)
sce.dlpfc <- runUMAP(sce.dlpfc, dimred="PCA_opt")


## Load in phenodata from pan-brain analysis -> colData for downstream use
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda",
     verbose=T)
    # Want 'ref.sampleInfo'

sce.dlpfc$region <- ss(sce.dlpfc$sample,".5",1)
sce.dlpfc$donor <- paste0("Br",ss(sce.dlpfc$sample,"c.",2))
sce.dlpfc$processDate <- ref.sampleInfo$realBatch[match(sce.dlpfc$sample, ref.sampleInfo$sampleID)]
sce.dlpfc$protocol <- ref.sampleInfo$protocol[match(sce.dlpfc$processDate, ref.sampleInfo$realBatch)]



### Clustering: Two-step ======================================================
### Step 1: Perform graph-based clustering in this optimal PC space
#         - take k=20 NN to build graph
snn.gr <- buildSNNGraph(sce.dlpfc, k=20, use.dimred="PCA_opt")
clusters.k20 <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters.k20)
    ##

# Assign as 'prelimCluster'
sce.dlpfc$prelimCluster <- factor(clusters.k20)

# Is sample driving this 'high-res' clustering at this level?
table(sce.dlpfc$prelimCluster, sce.dlpfc$sample)  # (a little bit, but is typical)

# Save for now
save(sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, ref.sampleInfo,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_DLPFC-n2_cleaned-combined_SCE_MNTFeb2020.rda")


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
        #28128  5410

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
myplclust(tree.clusCollapsed, main="2x DLPFC prelim-kNN-cluster relationships", cex.main=2, cex.lab=1.5, cex=1.8)


clust.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                               minClusterSize=2, deepSplit=1, cutHeight=365)


table(clust.treeCut)
unname(clust.treeCut[order.dendrogram(dend)])
    ## Cutting at 365 looks the best - go ahead and proceed with this

# The first cut-off prelimCluster only has 5 nuclei - re-merge with its originalmembers
clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)]==0)[1]] <- 2

# Add new labels to those remaining (2x) prelimClusters cut off
clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)]==0)] <- max(clust.treeCut)+c(1,2)


labels_colors(dend) <- tableau10medium[clust.treeCut[order.dendrogram(dend)]]

# Print for future reference
pdf("pdfs/regionSpecific_DLPFC-n2_HC-prelimCluster-relationships_Feb2020.pdf")
par(cex=1.1, font=2)
plot(dend, main="2x DLPFC prelim-kNN-cluster relationships")
dev.off()


# Make reference for new cluster assignment
clusterRefTab.dlpfc <- data.frame(origClust=order.dendrogram(dend),
                                   merged=clust.treeCut[order.dendrogram(dend)])


# Assign as 'collapsedCluster'
sce.dlpfc$collapsedCluster <- factor(clusterRefTab.dlpfc$merged[match(sce.dlpfc$prelimCluster, clusterRefTab.dlpfc$origClust)])

# Print some visualizations:
pdf("pdfs/regionSpecific_DLPFC-n2_reducedDims-with-collapsedClusters_Feb2020.pdf")
plotReducedDim(sce.dlpfc, dimred="PCA", ncomponents=5, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="sample", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="sum", point_alpha=0.5)
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
  'endothelial' = c('CLDN5', 'FLT1', 'VTN')
)

pdf("pdfs/zold_regionSpecific_DLPFC-n2_marker-logExprs_collapsedClusters_Feb2020.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.dlpfc, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="collapsedCluster", colour_by="collapsedCluster", point_alpha=0.5, point_size=.7,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:7], length(markers.mathys.custom[[i]])))
  )
}
dev.off()


# Observation: cluster 6 not an obvious 'cell type'...
newClusIndex <- splitit(sce.dlpfc$collapsedCluster)
sapply(newClusIndex, function(x) {quantile(sce.dlpfc[,x]$sum)})
    #           1        2     3     4       5       6        7
    #0%    1099.0  1509.00   523   432   981.0  101.00   499.00
    #25%  25255.5 14841.25  4372  3215  7471.0  131.75  2763.50
    #50%  36050.0 22512.00  6179  4717  9564.0  185.50  3962.50
    #75%  46973.0 31195.00  8002  6463 11972.5  296.50  5267.25
    #100% 89414.0 65971.00 25379 26618 23306.0 5789.00 11137.00

    # Looks like that collapsedCluster 6 is just driven by low # transcripts...


## Add annotations, looking at marker gene expression
annotationTab.dlpfc <- data.frame(cluster=c(1, 2, 3, 4, 5, 6, 7),
                                 cellType=c("Excit", "Inhib", "Oligo",
                                             "Astro", "OPC", "Ambig.lowNtrxts", "Micro")
)

sce.dlpfc$cellType <- annotationTab.dlpfc$cellType[match(sce.dlpfc$collapsedCluster,
                                                         annotationTab.dlpfc$cluster)]

## Save for now MNT 14Feb2020
save(sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc, ref.sampleInfo,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_DLPFC-n2_cleaned-combined_SCE_MNTFeb2020.rda")



### MNT 25Mar2020 === === ===
  # Re-print marker expression plots with annotated cluster names, after dropping 'Ambig.lowNtrxts'
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_DLPFC-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
table(sce.dlpfc$cellType)

# First drop "Ambig.lowNtrxts" (168 nuclei)
sce.dlpfc <- sce.dlpfc[ ,sce.dlpfc$cellType != "Ambig.lowNtrxts"]
sce.dlpfc$cellType <- droplevels(sce.dlpfc$cellType)


pdf("pdfs/regionSpecific_DLPFC-n2_marker-logExprs_collapsedClusters_Mar2020.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.dlpfc, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:6], length(markers.mathys.custom[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(markers.mathys.custom)[i], " markers"))
  )
}
dev.off()


    ## -> proceed to 'step03_markerDetxn-analyses[...].R'





