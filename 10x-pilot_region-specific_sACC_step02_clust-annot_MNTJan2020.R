### MNT 10x snRNA-seq workflow: step 02
###   **Region-specific analyses**
###     - (2x) sACC samples from: Br5161 & Br5212
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
# Take those sACC samples
pilot.sacc <- list(pilot.data[["sacc.5161"]],
                  pilot.data[["sacc.5212"]])
names(pilot.sacc) <- c("sacc.5161","sacc.5212")

### Newest iterations for normalization: cbind, THEN take scaled LSFs computed on all nuclei
# Add $sample identity
for(i in 1:length(pilot.sacc)){
  pilot.sacc[[i]]$sample <- names(pilot.sacc)[i]
}

sce.sacc <- cbind(pilot.sacc[[1]], pilot.sacc[[2]])

# Remove $logcounts
assay(sce.sacc, "logcounts") <- NULL
# Re-generate log-normalized counts
sce.sacc <- logNormCounts(sce.sacc)

geneVar.sacc <- modelGeneVar(sce.sacc)
chosen.hvgs.sacc <- geneVar.sacc$bio > 0
sum(chosen.hvgs.sacc)
    # [1] 10791


### Dimensionality reduction ================================================================

# Run PCA, taking top 100 (instead of default 50 PCs)
set.seed(109)
sce.sacc <- runPCA(sce.sacc, subset_row=chosen.hvgs.sacc, ncomponents=100,
                  BSPARAM=BiocSingular::RandomParam())

# Save into a new data file, which will dedicate for pan-brain-analyses
save(sce.sacc, chosen.hvgs.sacc, file="rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda")


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



      ## -> proceed to 'step03_markerDetxn-analyses[...].R'



