### MNT 10x snRNA-seq workflow: step 02
###   **Region-specific analyses**
###     - (2x) amygdala samples from: Br5161 & Br5212
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
# Take those amygdala samples
pilot.amy <- list(pilot.data[["amy.5161"]],
                  pilot.data[["amy.5212"]])
names(pilot.amy) <- c("amy.5161","amy.5212")

### Newest iterations for normalization: cbind, THEN take scaled LSFs computed on all nuclei
# Add $sample identity
for(i in 1:length(pilot.amy)){
  pilot.amy[[i]]$sample <- names(pilot.amy)[i]
}

sce.amy <- cbind(pilot.amy[[1]], pilot.amy[[2]])

# Remove $logcounts
assay(sce.amy, "logcounts") <- NULL
# Re-generate log-normalized counts
sce.amy <- logNormCounts(sce.amy)

geneVar.amy <- modelGeneVar(sce.amy)
chosen.hvgs.amy <- geneVar.amy$bio > 0
sum(chosen.hvgs.amy)
    # [1] 7819


### Dimensionality reduction ================================================================

# Run PCA, taking top 100 (instead of default 50 PCs)
set.seed(109)
sce.amy <- runPCA(sce.amy, subset_row=chosen.hvgs.amy, ncomponents=100,
                  BSPARAM=BiocSingular::RandomParam())

# Save into a new data file, which will dedicate for pan-brain-analyses
#save(sce.amy, chosen.hvgs.amy, file="rdas/zold_regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTJan2020.rda")
save(sce.amy, chosen.hvgs.amy, file="rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda")


## 'getClusteredPCs()' evaluated in qsub mode (with 'R-batchJob_Amyg-n2_optimalPCselxn_MNTFeb2020.R')
#    --> saved into same .rda


### Picking up with optimally-defined PC space ===
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=TRUE)
    # sce.amy, chosen.hvgs.amy, pc.choice.amy

# How many PCs is optimal?:
metadata(pc.choice.amy)$chosen
    # [1] 45

## Assign this chosen (40 PCs) to 'PCA_opt'
reducedDim(sce.amy, "PCA_opt") <- reducedDim(sce.amy, "PCA")[ ,1:(metadata(pc.choice.amy)$chosen)]


## t-SNE
set.seed(109)
sce.amy <- runTSNE(sce.amy, dimred="PCA_opt")
plotTSNE(sce.amy, colour_by="sample")


## UMAP
set.seed(109)
sce.amy <- runUMAP(sce.amy, dimred="PCA_opt")


## Load in phenodata from pan-brain analysis -> colData for downstream use
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda",
     verbose=T)
    # Want 'ref.sampleInfo'

sce.amy$region <- ss(sce.amy$sample,".5",1)
sce.amy$donor <- paste0("Br",ss(sce.amy$sample,"y.",2))
sce.amy$processDate <- ref.sampleInfo$realBatch[match(sce.amy$sample, ref.sampleInfo$sampleID)]
sce.amy$protocol <- ref.sampleInfo$protocol[match(sce.amy$processDate, ref.sampleInfo$realBatch)]



### Clustering: Two-step ======================================================
### Step 1: Perform graph-based clustering in this optimal PC space
#         - take k=20 NN to build graph
snn.gr <- buildSNNGraph(sce.amy, k=20, use.dimred="PCA_opt")
clusters.k20 <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters.k20)
    ## 23 prelim clusters

# Assign as 'prelimCluster'
sce.amy$prelimCluster <- factor(clusters.k20)

# Is sample driving this 'high-res' clustering at this level?
table(sce.amy$prelimCluster, sce.amy$sample)  # (a little bit, but is typical)

# Save for now
save(sce.amy, chosen.hvgs.amy, pc.choice.amy, ref.sampleInfo,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda")


### Step 2: Hierarchical clustering of pseudo-bulked ("PB'd") counts with most robust normalization
  #         (as determined in: 'side-Rscript_testingStep2_HC-normalizn-approaches_wAmygData_MNTJan2020.R')
  #           ** That is, to pseudo-bulk (aka 'cluster-bulk') on raw counts, on all [non-zero] genes,
  #              normalize with `librarySizeFactors()`, log2-transform, then perform HC'ing


# Preliminary cluster index for pseudo-bulking
clusIndexes = splitit(sce.amy$prelimCluster)
prelimCluster.PBcounts <- sapply(clusIndexes, function(ii){
  rowSums(assays(sce.amy)$counts[ ,ii])
  }
)
    
    # And btw...
    table(rowSums(prelimCluster.PBcounts)==0)
        #FALSE  TRUE
        #28470  5068

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
                               minClusterSize=2, deepSplit=1, cutHeight=375)


table(clust.treeCut)
unname(clust.treeCut[order.dendrogram(dend)])
    ## Cutting at 375 looks good - go ahead and proceed with this

# Add new labels to those prelimClusters cut off
clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)]==0)] <- max(clust.treeCut)+c(1,2)

labels_colors(dend) <- tableau10medium[clust.treeCut[order.dendrogram(dend)]]

# Print for future reference
pdf("pdfs/regionSpecific_Amyg-n2_HC-prelimCluster-relationships_Feb2020.pdf")
par(cex=1.1, font=2)
plot(dend)
dev.off()


# Make reference for new cluster assignment
clusterRefTab.amy <- data.frame(origClust=order.dendrogram(dend),
                                merged=clust.treeCut[order.dendrogram(dend)])


# Assign as 'collapsedCluster'
sce.amy$collapsedCluster <- factor(clusterRefTab.amy$merged[match(sce.amy$prelimCluster, clusterRefTab.amy$origClust)])

# Print some visualizations:
pdf("pdfs/regionSpecific_Amyg-n2_reducedDims-with-collapsedClusters_Feb2020.pdf")
plotReducedDim(sce.amy, dimred="PCA", ncomponents=5, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.amy, colour_by="sample", point_alpha=0.5)
plotTSNE(sce.amy, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.amy, colour_by="sum", point_alpha=0.5)
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
  'endothelial' = c('CLDN5', 'FLT1', 'VTN')
)

pdf("pdfs/zold_regionSpecific_Amyg-n2_marker-logExprs_collapsedClusters_Feb2020.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.amy, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="collapsedCluster", colour_by="collapsedCluster", point_alpha=0.5, point_size=.7,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:7], length(markers.mathys.custom[[i]])))
  )
}
dev.off()




# Observation: cluster 6 not an obvious 'cell type'
newClusIndex <- splitit(sce.amy$collapsedCluster)
sapply(newClusIndex, function(x) {quantile(sce.amy[,x]$sum)})
    #          1        2       3     4       5       6        7
    #0%     1974    787.0   597.0   756  2005.0  104.00   241.00
    #25%   23879  30981.5  6979.5  4952  9529.0  126.50  3486.75
    #50%   38317  48113.0 10514.0  7185 12947.0  177.50  5277.00
    #75%   49701  67138.5 15754.5  9996 16867.5  256.75  6874.75
    #100% 110271 165583.0 34702.0 32011 32885.0 7465.00 14531.00

## Add annotations, looking at marker gene expression
annotationTab.amy <- data.frame(cluster=c(1, 2, 3, 4, 5, 6, 7),
                                  cellType=c("Inhib", "Excit", "Astro",
                                             "Oligo", "OPC", "Ambig.lowNtrxts", "Micro")
)

sce.amy$cellType <- annotationTab.amy$cellType[match(sce.amy$collapsedCluster,
                                                         annotationTab.amy$cluster)]

# Save for now MNT 14Feb2020
save(sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo, 
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda")



### MNT 20Mar2020 === === ===
  # Re-print marker expression plots with annotated cluster names, after dropping 'Ambig.lowNtrxts'
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
table(sce.amy$cellType)

# First drop "Ambig.lowNtrxts" (50 nuclei)
sce.amy <- sce.amy[ ,sce.amy$cellType != "Ambig.lowNtrxts"]
sce.amy$cellType <- droplevels(sce.amy$cellType)


pdf("pdfs/regionSpecific_Amyg-n2_marker-logExprs_collapsedClusters_Mar2020.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.amy, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
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



### For reference === == === == ===
table(sce.amy$cellType, sce.amy$sample)
    #                 amy.5161 amy.5212
    # Ambig.lowNtrxts       34       16
    # Astro                489      363
    # Excit                226      261
    # Inhib                 84      295
    # Micro                425      339
    # Oligo               1697     1776
    # OPC                  335      292

table(sce.amy$prelimCluster, sce.amy$sample)
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

table(sce.amy$cellType, sce.amy$collapsedCluster)
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
    # sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo

# Look at some marker expression at the prelimCluster level
pdf("pdfs/ztemp_amyg-prelimCluster-neuronalMarkerExpression.pdf", width=8, height=7)
plotExpression(sce.amy, exprs_values="logcounts", features=c("SNAP25","GAD1","GAD2","SLC17A6","SLC17A7","VCAN"),
               x="prelimCluster", colour_by="prelimCluster", ncol=2)
dev.off()


## More-manual annotations for neuronal subpops:
clusterRefTab.amy$cellType <- sce.amy$cellType[match(clusterRefTab.amy$merged, sce.amy$collapsedCluster)]
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
sce.amy$cellType.split <- clusterRefTab.amy$manual[match(sce.amy$prelimCluster,
                                                         clusterRefTab.amy$origClust)]
sce.amy$cellType.split <- factor(sce.amy$cellType.split)

table(sce.amy$cellType.split, sce.amy$cellType)
    # good

table(sce.amy$cellType.split) # (printing post-hoc-corrected annotations)
    #Ambig.lowNtrxts           Astro         Excit.1         Excit.2         Excit.3
    #             50             852             334              40              55
    #        Inhib.1         Inhib.2         Inhib.3         Inhib.4         Inhib.5
    #            171             109              35              24              98
    #          Micro           Oligo             OPC
    #            764            3473             627

## Save these
save(sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo,
     file="rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda")


## Also print expression at this level of partitioning ===

# First remove "Ambig.lowNtrxts" (50 nuclei):
sce.amy <- sce.amy[ ,sce.amy$cellType.split != "Ambig.lowNtrxts"]
sce.amy$cellType.split <- droplevels(sce.amy$cellType.split)

pdf("pdfs/regionSpecific_Amyg-n2_marker-logExprs_cellTypesSplit_May2020.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.amy, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:12], length(markers.mathys.custom[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(markers.mathys.custom)[i], " markers"))
  )
}
dev.off()


## Let's also re-plot reducedDims with new [broad & split] cell type annotations
#        (and rename old file with prefix 'zold_')
pdf("pdfs/regionSpecific_Amyg-n2_reducedDims-with-collapsedClusters_May2020.pdf")
plotReducedDim(sce.amy, dimred="PCA", ncomponents=5, colour_by="cellType", point_alpha=0.5)
plotTSNE(sce.amy, colour_by="sample", point_size=3.5, point_alpha=0.5)
plotTSNE(sce.amy, colour_by="prelimCluster", point_size=3.5, point_alpha=0.5)
plotTSNE(sce.amy, colour_by="cellType", point_size=3.5, point_alpha=0.5)
plotTSNE(sce.amy, colour_by="cellType.split", point_size=3.5, point_alpha=0.5)
plotTSNE(sce.amy, colour_by="sum", point_size=3.5, point_alpha=0.5)
plotUMAP(sce.amy, colour_by="cellType", point_size=3.5, point_alpha=0.5)
plotUMAP(sce.amy, colour_by="cellType.split", point_size=3.5, point_alpha=0.5)
dev.off()


## And finally, for reference:
table(sce.amy$cellType.split, sce.amy$sample)
    #         amy.5161 amy.5212
    # Astro        489      363
    # Excit.1      141      193
    # Excit.2       85       13
    # Excit.3        0       55
    # Inhib.1       16      155
    # Inhib.2       33       76
    # Inhib.3       11       24
    # Inhib.4       24        0
    # Inhib.5        0       40
    # Micro        425      339
    # Oligo       1697     1776
    # OPC          335      292

    