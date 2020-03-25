### MNT 10x snRNA-seq workflow: step 02
###   **Region-specific analyses**
###     - (3x) HPC samples from: Br5161 & Br5212 & Br5287
### Initiated MNT 07Feb2020
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
# Take those HPC samples
pilot.hpc <- list(pilot.data[["hpc.5161"]],
                  pilot.data[["hpc.5212"]],
                  pilot.data[["hpc.5287"]])
names(pilot.hpc) <- c("hpc.5161","hpc.5212","hpc.5287")

### Newest iterations for normalization: cbind, THEN take scaled LSFs computed on all nuclei
# Add $sample identity
for(i in 1:length(pilot.hpc)){
  pilot.hpc[[i]]$sample <- names(pilot.hpc)[i]
}

sce.hpc <- cbind(pilot.hpc[[1]], pilot.hpc[[2]], pilot.hpc[[3]])

# Remove $logcounts
assay(sce.hpc, "logcounts") <- NULL
# Re-generate log-normalized counts
sce.hpc <- logNormCounts(sce.hpc)

geneVar.hpc <- modelGeneVar(sce.hpc)
chosen.hvgs.hpc <- geneVar.hpc$bio > 0
sum(chosen.hvgs.hpc)
    # [1] 8997


### Dimensionality reduction ================================================================

# Run PCA, taking top 100 (instead of default 50 PCs)
set.seed(109)
sce.hpc <- runPCA(sce.hpc, subset_row=chosen.hvgs.hpc, ncomponents=100,
                  BSPARAM=BiocSingular::RandomParam())

# Save into a new data file, which will dedicate for pan-brain-analyses
save(sce.hpc, chosen.hvgs.hpc, file="rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda")


## 'getClusteredPCs()' evaluated in qsub mode (with 'R-batchJob_HPC-n3_optimalPCselxn_MNTFeb2020.R')
#    --> saved into same .rda


### Picking up with optimally-defined PC space ===
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=TRUE)
    # sce.hpc, chosen.hvgs.hpc, pc.choice.hpc


# How many PCs is optimal?:
metadata(pc.choice.hpc)$chosen
    ## 49

## Assign this chosen ( PCs) to 'PCA_opt'
reducedDim(sce.hpc, "PCA_opt") <- reducedDim(sce.hpc, "PCA")[ ,1:(metadata(pc.choice.hpc)$chosen)]


## t-SNE
set.seed(109)
sce.hpc <- runTSNE(sce.hpc, dimred="PCA_opt")


## UMAP
set.seed(109)
sce.hpc <- runUMAP(sce.hpc, dimred="PCA_opt")

## Load in phenodata from pan-brain analysis -> colData for downstream use
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda",
     verbose=T)
    # Want 'ref.sampleInfo'

sce.hpc$region <- ss(sce.hpc$sample,".5",1)
sce.hpc$donor <- paste0("Br",ss(sce.hpc$sample,"c.",2))
sce.hpc$processDate <- ref.sampleInfo$realBatch[match(sce.hpc$sample, ref.sampleInfo$sampleID)]
sce.hpc$protocol <- ref.sampleInfo$protocol[match(sce.hpc$processDate, ref.sampleInfo$realBatch)]

# Save for now
save(sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, ref.sampleInfo,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda")


### Clustering: Two-step ======================================================
### Step 1: Perform graph-based clustering in this optimal PC space
#         - take k=20 NN to build graph
snn.gr <- buildSNNGraph(sce.hpc, k=20, use.dimred="PCA_opt")
clusters.k20 <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters.k20)
    ## 

# Assign as 'prelimCluster'
sce.hpc$prelimCluster <- factor(clusters.k20)

# Is sample driving this 'high-res' clustering at this level?
table(sce.hpc$prelimCluster, sce.hpc$sample)  

table(sce.hpc$sample)

### Step 2: Hierarchical clustering of pseudo-bulked ("PB'd") counts with most robust normalization
#         (as determined in: 'side-Rscript_testingStep2_HC-normalizn-approaches_wAmygData_MNTJan2020.R')
#           ** That is, to pseudo-bulk (aka 'cluster-bulk') on raw counts, on all [non-zero] genes,
#              normalize with `librarySizeFactors()`, log2-transform, then perform HC'ing


# Preliminary cluster index for pseudo-bulking
clusIndexes = splitit(sce.hpc$prelimCluster)
prelimCluster.PBcounts <- sapply(clusIndexes, function(ii){
  rowSums(assays(sce.hpc)$counts[ ,ii])
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
myplclust(tree.clusCollapsed, main="3x HPC prelim-kNN-cluster relationships", cex.main=2, cex.lab=1.5, cex=1.8)


clust.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                               minClusterSize=2, deepSplit=1, cutHeight=475)


table(clust.treeCut)
unname(clust.treeCut[order.dendrogram(dend)])
## Cutting at 475 looks the best - go ahead and proceed with this

    ## The first cut-off prelimCluster only has 5 nuclei - re-merge with its originalmembers
    #clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)]==0)[1]] <- 2

# Add new labels to those (2x) prelimClusters cut off
clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)]==0)] <- max(clust.treeCut)+c(1,2)


labels_colors(dend) <- tableau10medium[clust.treeCut[order.dendrogram(dend)]]

# Print for future reference
pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/regionSpecific_HPC-n3_HC-prelimCluster-relationships_Feb2020.pdf")
par(cex=1.1, font=2)
plot(dend, main="3x HPC prelim-kNN-cluster relationships")
dev.off()


# Make reference for new cluster assignment
clusterRefTab.hpc <- data.frame(origClust=order.dendrogram(dend),
                                  merged=clust.treeCut[order.dendrogram(dend)])


# Assign as 'collapsedCluster'
sce.hpc$collapsedCluster <- factor(clusterRefTab.hpc$merged[match(sce.hpc$prelimCluster, clusterRefTab.hpc$origClust)])

# Print some visualizations:
pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/regionSpecific_HPC-n3_reducedDims-with-collapsedClusters_Feb2020.pdf")
plotReducedDim(sce.hpc, dimred="PCA", ncomponents=5, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.hpc, colour_by="sample", point_alpha=0.5)
plotTSNE(sce.hpc, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.hpc, colour_by="sum", point_alpha=0.5)
plotUMAP(sce.hpc, colour_by="collapsedCluster", point_alpha=0.5)
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
  # Added MNT 20Mar2020
  'Tcell' = c('TRAC','SKAP1','CCL5')
)

pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/zold_regionSpecific_HPC-n3_marker-logExprs_collapsedClusters_Feb2020.pdf",
    height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.hpc, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="collapsedCluster", colour_by="collapsedCluster", point_alpha=0.5, point_size=.7,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:8], length(markers.mathys.custom[[i]])))
  )
}
dev.off()


# Observation: cluster 7 or 8 not an obvious 'cell type'...
newClusIndex <- splitit(sce.hpc$collapsedCluster)
sapply(newClusIndex, function(x) {quantile(sce.hpc[,x]$sum)})
    #            1       2      3       4     5     6    7       8
    #0%      981.0   231.0   1152   343.0   328   168  102  783.00
    #25%   22657.0  4215.5  18243  6145.0  4488  3481  127 2314.75
    #50%   32994.0  7025.0  28322  9248.0  6266  4661  156 3215.00
    #75%   46054.5 10610.0  40240 12553.5  8558  6103  236 3756.50
    #100% 150461.0 30085.0 111453 30332.0 29556 22692 1739 6038.00

    # Looks like that collapsedCluster 6 is just driven by low # transcripts...
    library(pheatmap)
    cc3 <- assay(sce.hpc, "logcounts")[ ,sce.hpc$collapsedCluster==3]
    cor.cc3 <- cor(as.matrix(cc3))
    pheatmap(cor.cc3,show_rownames=F,show_colnames=F)
    
    cc7 <- assay(sce.hpc, "logcounts")[ ,sce.hpc$collapsedCluster==7]
    cor.cc7 <- cor(as.matrix(cc7))  # pretty poor, as expected
    pheatmap(cor.cc7,show_rownames=F,show_colnames=F)
    quantile(cor.cc7)
        #          0%          25%          50%          75%         100%
        #-0.005752025  0.021315889  0.033826675  0.049805456  1.000000000
    
    cc8 <- assay(sce.hpc, "logcounts")[ ,sce.hpc$collapsedCluster==8]
    cor.cc8 <- cor(as.matrix(cc8))  # pretty poor - less bad than 7
    pheatmap(cor.cc8,show_rownames=F,show_colnames=F)
    quantile(cor.cc8)
        #       0%       25%       50%       75%      100%
        #0.2003143 0.3150381 0.3570639 0.3946995 1.0000000

## Add annotations, looking at marker gene expression
annotationTab.hpc <- data.frame(cluster=c(1, 2, 3, 4, 5, 6, 7, 8),
                                  cellType=c("Excit", "Astro", "Inhib", "OPC",
                                             "Oligo", "Micro", "Ambig.lowNtrxts", "Ambig.glial")
)

sce.hpc$cellType <- annotationTab.hpc$cellType[match(sce.hpc$collapsedCluster,
                                                         annotationTab.hpc$cluster)]

## Save for now MNT 20Feb2020
save(sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda")



### MNT 20Mar2020 === === ===
# Re-print marker expression plots with annotated cluster names, after dropping 'Ambig.lowNtrxts'
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
table(sce.hpc$cellType)

# First drop "Ambig.lowNtrxts" (168 nuclei)
sce.hpc <- sce.hpc[ ,sce.hpc$cellType != "Ambig.lowNtrxts"]
# Then rename "Ambig.glial" to "Tcell" (26 nuclei)
#     (A posteriori - from downstream marker exploration)
sce.hpc.temp <- sce.hpc
sce.hpc.temp$cellType <- droplevels(sce.hpc.temp$cellType)
sce.hpc.temp$cellType <- factor(gsub(pattern="Ambig.glial", "Tcell", sce.hpc.temp$cellType))


pdf("pdfs/regionSpecific_HPC-n3_marker-logExprs_collapsedClusters_Mar2020.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.hpc.temp, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:7], length(markers.mathys.custom[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(markers.mathys.custom)[i], " markers"))
  )
}
dev.off()

      ## -> proceed to 'step03_markerDetxn-analyses[...].R'










