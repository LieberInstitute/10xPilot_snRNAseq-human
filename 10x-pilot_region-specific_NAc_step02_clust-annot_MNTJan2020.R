### MNT 10x snRNA-seq workflow: step 02
###   **Region-specific analyses**
###     - (2x) NAc samples from: Br5161 & Br5212 & Br5287
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
# Take those NAc samples
pilot.nac <- list(pilot.data[["nac.5161"]],
                  pilot.data[["nac.5212"]],
                  pilot.data[["nac.5287"]])
names(pilot.nac) <- c("nac.5161","nac.5212","nac.5287")

### Newest iterations for normalization: cbind, THEN take scaled LSFs computed on all nuclei
# Add $sample identity
for(i in 1:length(pilot.nac)){
  pilot.nac[[i]]$sample <- names(pilot.nac)[i]
}

sce.nac <- cbind(pilot.nac[[1]], pilot.nac[[2]], pilot.nac[[3]])

# Remove $logcounts
assay(sce.nac, "logcounts") <- NULL
# Re-generate log-normalized counts
sce.nac <- logNormCounts(sce.nac)

geneVar.nac <- modelGeneVar(sce.nac)
chosen.hvgs.nac <- geneVar.nac$bio > 0
sum(chosen.hvgs.nac)
    # [1] 9326


### Dimensionality reduction ================================================================

# Run PCA, taking top 100 (instead of default 50 PCs)
set.seed(109)
sce.nac <- runPCA(sce.nac, subset_row=chosen.hvgs.nac, ncomponents=100,
                  BSPARAM=BiocSingular::RandomParam())

# Save into a new data file, which will dedicate for pan-brain-analyses
save(sce.nac, chosen.hvgs.nac, file="rdas/regionSpecific_NAc-n3_cleaned-combined_SCE_MNTFeb2020.rda")


## 'getClusteredPCs()' evaluated in qsub mode (with 'R-batchJob_NAc-n3_optimalPCselxn_MNTFeb2020.R')
#    --> saved into same .rda


### Picking up with optimally-defined PC space ===
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_NAc-n3_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=TRUE)
    # sce.nac, chosen.hvgs.nac, pc.choice.nac


# How many PCs is optimal?:
metadata(pc.choice.nac)$chosen  # [1] 42

## Assign this chosen ( PCs) to 'PCA_opt'
reducedDim(sce.nac, "PCA_opt") <- reducedDim(sce.nac, "PCA")[ ,1:(metadata(pc.choice.nac)$chosen)]


## t-SNE
set.seed(109)
sce.nac <- runTSNE(sce.nac, dimred="PCA_opt")


## UMAP
set.seed(109)
sce.nac <- runUMAP(sce.nac, dimred="PCA_opt")

## Load in phenodata from pan-brain analysis -> colData for downstream use
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda",
     verbose=T)
# Want 'ref.sampleInfo'

sce.nac$region <- ss(sce.nac$sample,".5",1)
sce.nac$donor <- paste0("Br",ss(sce.nac$sample,"c.",2))
sce.nac$processDate <- ref.sampleInfo$realBatch[match(sce.nac$sample, ref.sampleInfo$sampleID)]
sce.nac$protocol <- ref.sampleInfo$protocol[match(sce.nac$processDate, ref.sampleInfo$realBatch)]


# Save for now
save(sce.nac, chosen.hvgs.nac, pc.choice.nac, ref.sampleInfo,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_NAc-n3_cleaned-combined_SCE_MNTFeb2020.rda")


### Clustering: Two-step ======================================================
### Step 1: Perform graph-based clustering in this optimal PC space
#         - take k=20 NN to build graph
snn.gr <- buildSNNGraph(sce.nac, k=20, use.dimred="PCA_opt")
clusters.k20 <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters.k20)
    ## 18 prelim clusters

# Assign as 'prelimCluster'
sce.nac$prelimCluster <- factor(clusters.k20)

# Is sample driving this 'high-res' clustering at this level?
table(sce.nac$prelimCluster, sce.nac$sample)

table(sce.nac$sample)

### Step 2: Hierarchical clustering of pseudo-bulked ("PB'd") counts with most robust normalization
  #         (as determined in: 'side-Rscript_testingStep2_HC-normalizn-approaches_wAmygData_MNTJan2020.R')
#           ** That is, to pseudo-bulk (aka 'cluster-bulk') on raw counts, on all [non-zero] genes,
#              normalize with `librarySizeFactors()`, log2-transform, then perform HC'ing


# Preliminary cluster index for pseudo-bulking
clusIndexes = splitit(sce.nac$prelimCluster)
prelimCluster.PBcounts <- sapply(clusIndexes, function(ii){
  rowSums(assays(sce.nac)$counts[ ,ii])
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
myplclust(tree.clusCollapsed, main="2x Amyg prelim-kNN-cluster relationships", cex.main=2, cex.lab=1.5, cex=1.8)


clust.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                               minClusterSize=2, deepSplit=1, cutHeight=400)


table(clust.treeCut)
unname(clust.treeCut[order.dendrogram(dend)])
  ## Cutting at 400 looks good - go ahead and proceed with this

# Add new labels to those prelimClusters cut off
clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)]==0)] <- max(clust.treeCut)+c(1)

labels_colors(dend) <- tableau10medium[clust.treeCut[order.dendrogram(dend)]]

# Print for future reference
pdf("pdfs/regionSpecific_NAc-n3_HC-prelimCluster-relationships_Feb2020.pdf")
par(cex=1.1, font=2)
plot(dend)
dev.off()


# Make reference for new cluster assignment
clusterRefTab.nac <- data.frame(origClust=order.dendrogram(dend),
                                merged=clust.treeCut[order.dendrogram(dend)])

# Assign as 'collapsedCluster'
sce.nac$collapsedCluster <- factor(clusterRefTab.nac$merged[match(sce.nac$prelimCluster, clusterRefTab.nac$origClust)])


# Print some visualizations:
pdf("pdfs/regionSpecific_NAc-n3_reducedDims-with-collapsedClusters_Feb2020.pdf")
plotReducedDim(sce.nac, dimred="PCA", ncomponents=5, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.nac, colour_by="sample", point_alpha=0.5)
plotTSNE(sce.nac, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.nac, colour_by="sum", point_alpha=0.5)
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
  'endothelial' = c('CLDN5', 'FLT1', 'VTN')
)

pdf("pdfs/regionSpecific_NAc-n3_marker-logExprs_collapsedClusters_Feb2020.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.nac, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="collapsedCluster", colour_by="collapsedCluster", point_alpha=0.5, point_size=.7,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:7], length(markers.mathys.custom[[i]])))
  )
}
dev.off()




# Save for now
save(sce.nac, chosen.hvgs.nac, pc.choice.nac, ref.sampleInfo, clusterRefTab.nac,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_NAc-n3_cleaned-combined_SCE_MNTFeb2020.rda")


# Observation: cluster 6 not an obvious 'cell type'...
newClusIndex <- splitit(sce.nac$collapsedCluster)
sapply(newClusIndex, function(x) {quantile(sce.nac[,x]$sum)})
    #           1     2        3       4     5        6      7
    #0%     365.0  1197   635.00  1379.0   938   218.00  101.0
    #25%   3619.0 23523  6301.50  6937.5 13104  2335.25  118.0
    #50%   5573.0 31839  8379.00  9145.0 20825  3994.00  144.5
    #75%   8011.5 41978 10755.75 11354.0 34874  6160.50  267.5
    #100% 31312.0 78013 32867.00 31898.0 67288 12537.00 1970.0

    # Looks like that collapsedCluster 7 is just driven by low # transcripts...
    library(pheatmap)
    cc7 <- assay(sce.nac, "logcounts")[ ,sce.nac$collapsedCluster==7]
    cor.cc7 <- cor(as.matrix(cc7))
    pheatmap(cor.cc7) # preeeety poor
    cc5 <- assay(sce.nac, "logcounts")[ ,sce.nac$collapsedCluster==5]
    cor.cc5 <- cor(as.matrix(cc5))
    pheatmap(cor.cc5) # oh yeah.  Then these are definitely just very poorly-sequenced nuclei
                      #           or had a bunch of random debris partitioned into a GEM

## Add annotations, looking at marker gene expression
annotationTab.nac <- data.frame(cluster=c(1, 2, 3, 4, 5, 6, 7),
                                  cellType=c("Oligo", "Inhib.NRGNpos", "Astro",
                                             "OPC", "Inhib.NRGNneg", "Micro", "Ambig.lowNtrxts")
)

sce.nac$cellType <- annotationTab.nac$cellType[match(sce.nac$collapsedCluster,
                                                         annotationTab.nac$cluster)]

## Save for now MNT 14Feb2020
save(sce.nac, chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac, ref.sampleInfo,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_NAc-n3_cleaned-combined_SCE_MNTFeb2020.rda")



      ## -> proceed to 'step03_markerDetxn-analyses[...].R'





