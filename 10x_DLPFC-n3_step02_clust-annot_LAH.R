
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
library(here)

source("plotExpressionCustom.R")

load(here("rdas","revision","tableau_colors.rda"), verbose = TRUE)

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
     file="rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda")


    ## 'getClusteredPCs()' evaluated in qsub mode (with 'R-batchJob_Amyg-n2_optimalPCselxn_MNTFeb2020.R')
    #    --> saved into same .rda


### Picking up with optimally-defined PC space ===
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda",
     verbose=TRUE)
    # sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, ref.sampleInfo, ref.sampleInfo.rev

# How many PCs is optimal?:
metadata(pc.choice.dlpfc)$chosen
    # [1] 99

## Assign this chosen (99 PCs) to 'PCA_opt'
reducedDim(sce.dlpfc, "PCA_opt") <- reducedDim(sce.dlpfc, "PCA_corrected")[ ,1:(metadata(pc.choice.dlpfc)$chosen)]


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
    ## 107 prelim clusters

# Assign as 'prelimCluster'
sce.dlpfc$prelimCluster <- factor(clusters.k20)
plotReducedDim(sce.dlpfc, dimred="TSNE", colour_by="prelimCluster")

# Is sample driving this 'high-res' clustering at this level?
(sample_prelimClusters <- table(sce.dlpfc$prelimCluster, sce.dlpfc$sampleID))  # (a little bit, but is typical)
sample_prelimClusters[which(rowSums(sample_prelimClusters == 0) == 2),]
# 39 - only 4 samples all from Br5207

# rbind the ref.sampleInfo[.rev]
ref.sampleInfo <- rbind(ref.sampleInfo, ref.sampleInfo.rev)

## check doublet score for each prelim clust
clusIndexes = splitit(sce.dlpfc$prelimCluster)
prelimCluster.medianDoublet <- sapply(clusIndexes, function(ii){
  median(sce.dlpfc$doubletScore[ii])
}
)

summary(prelimCluster.medianDoublet)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.01059  0.07083  0.14823  0.53264  0.30064 14.79144

hist(prelimCluster.medianDoublet)

## watch in clustering
prelimCluster.medianDoublet[prelimCluster.medianDoublet > 5]
# 19       32       73
# 14.79144  7.98099 10.20462

table(sce.dlpfc$prelimCluster)[c(19, 32, 73)]
# 19 32 73
# 27 32  8

# Save for now
save(sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, ref.sampleInfo,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda")


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
    # FALSE  TRUE
    # 29310  4228

# Compute LSFs at this level
sizeFactors.PB.all  <- librarySizeFactors(prelimCluster.PBcounts)

# Normalize with these LSFs
geneExprs.temp <- t(apply(prelimCluster.PBcounts, 1, function(x) {log2(x/sizeFactors.PB.all + 1)}))

## Perform hierarchical clustering
dist.clusCollapsed <- dist(t(geneExprs.temp))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

dend <- as.dendrogram(tree.clusCollapsed, hang=0.2)

# labels(dend)[grep(c("19|32|73"),labels(dend))] <- paste0(labels(dend)[grep(c("19|32|73"),labels(dend))], "*")

# Just for observation
par(cex=.6)
myplclust(tree.clusCollapsed, cex.main=2, cex.lab=1.5, cex=1.8)

dend %>%
  set("labels_cex", 0.8) %>%
  plot(horiz = TRUE)
abline(v = 325, lty = 2)

clust.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                               minClusterSize=2, deepSplit=1, cutHeight=325)

table(clust.treeCut)
unname(clust.treeCut[order.dendrogram(dend)])
    ## Cutting at 250 looks good for the main neuronal branch, but a lot of glial
     #    prelim clusters are dropped off (0's)

    # # Cut at 400 for broad glia branch (will manually merge remaining dropped off)
    # glia.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
    #                               minClusterSize=2, deepSplit=1, cutHeight=400)
    # unname(glia.treeCut[order.dendrogram(dend)])

    # Take those and re-assign to the first assignments

# clust <- clust.treeCut[order.dendrogram(dend)]
# clust2 <- name_zeros(clust, list(c(1,2), c(106,107)))
# unname(clust2)

# Add new labels to those prelimClusters cut off
clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)]==0)] <- max(clust.treeCut)+c(1, 6, 2, 3, 4, 5, 5)

# 'Re-write', since there are missing numbers
# clust.treeCut[order.dendrogram(dend)] <- as.numeric(as.factor(clust2))
clust.treeCut[order.dendrogram(dend)] <- as.numeric(as.factor(clust.treeCut[order.dendrogram(dend)]))

## Define color pallet
cluster_colors <- unique(tableau20[clust.treeCut[order.dendrogram(dend)]])
names(cluster_colors) <- unique(clust.treeCut[order.dendrogram(dend)])
labels_colors(dend) <- cluster_colors[clust.treeCut[order.dendrogram(dend)]]

# Print for future reference
pdf("pdfs/revision/regionSpecific_DLPFC-n3_HC-prelimCluster-relationships_LAH2021.pdf", height = 9)
par(cex=0.6, font=2)
plot(dend, main="3x DLPFC prelim-kNN-cluster relationships with collapsed assignments", horiz = TRUE)
abline(v = 325, lty = 2)
dev.off()


# Make reference for new cluster assignment
clusterRefTab.dlpfc <- data.frame(origClust=order.dendrogram(dend),
                                merged=clust.treeCut[order.dendrogram(dend)])


# Assign as 'collapsedCluster'
sce.dlpfc$collapsedCluster <- factor(clusterRefTab.dlpfc$merged[match(sce.dlpfc$prelimCluster, clusterRefTab.dlpfc$origClust)])
n_clusters <- length(levels(sce.dlpfc$collapsedCluster))
# Print some visualizations:
pdf("pdfs/revision/regionSpecific_DLPFC-n3_reducedDims-with-collapsedClusters_LAH2021.pdf")
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

tail(table(sce.dlpfc$prelimCluster, sce.dlpfc$collapsedCluster),20)

## Print marker genes for annotation
load(here("rdas","revision","markers.rda"), verbose = TRUE)

pdf("pdfs/revision/regionSpecific_DLPFC-n3_marker-logExprs_collapsedClusters_LAH2021.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpressionCustom(sce = sce.dlpfc,
                         features = markers.mathys.custom[[i]],
                         features_name = names(markers.mathys.custom)[[i]],
                         anno_name = "collapsedCluster")
  )
}
dev.off()


## Add annotations, looking at marker gene expression
annotationTab.dlpfc <- data.frame(collapsedCluster=c(1:n_clusters))
annotationTab.dlpfc$cellType <- NA
annotationTab.dlpfc$cellType[c(1:3, 6, 17,18)] <- paste0("Inhib_", c("A","B","C","D","E","F"))
annotationTab.dlpfc$cellType[c(4,5,8,10:13)] <- paste0("Excit_", c("A","B","C","D","E","F","G"))
annotationTab.dlpfc$cellType[c(7, 9, 13, 16)] <- c("Astro", "Oligo", "OPC", "Micro")
annotationTab.dlpfc$cellType[c(14,15,19)] <- c("Macrophage","Mural","Tcell")


sce.dlpfc$cellType <- annotationTab.dlpfc$cellType[match(sce.dlpfc$collapsedCluster,
                                                         annotationTab.dlpfc$collapsedCluster)]
sce.dlpfc$cellType <- factor(sce.dlpfc$cellType)

table(sce.dlpfc$cellType)
# Astro    Excit_A    Excit_B    Excit_C    Excit_D    Excit_E    Excit_F    Inhib_A    Inhib_B    Inhib_C 
# 782        529        773        524        132        187        243        333        454        365 
# Inhib_D    Inhib_E    Inhib_F Macrophage      Micro      Mural      Oligo        OPC      Tcell 
# 413          7          8         10        388         18       5455        572          9

## QC - How do the total UMI distribution look?
# newClusIndex <- splitit(sce.dlpfc$collapsedCluster)
newClusIndex <- splitit(sce.dlpfc$cellType)
sapply(newClusIndex, function(x) {quantile(sce.dlpfc$sum[x])})
# Astro Excit_A Excit_B   Excit_C  Excit_D  Excit_E Excit_F Inhib_A  Inhib_B Inhib_C Inhib_D Inhib_E  Inhib_F
# 0%     884.00    1201    1838   2295.00  9100.00   3029.0    1099    3045  1708.00    1803    1749  7966.0  8833.00
# 25%   4008.75   28468   19205  39673.25 26131.75  37044.5   27809   15330 11474.25   18791   23015  9056.5 19428.75
# 50%   5737.00   34997   25241  49345.00 32065.50  49414.0   37430   18968 16423.50   24138   29623 14241.0 20281.50
# 75%   7953.00   43590   33292  59114.25 39545.00  61162.0   47099   23036 22355.75   28996   35364 17472.0 22677.00
# 100% 26618.00   87392   69309 115449.00 69202.00 101625.0   83826   55574 66556.00   81134   67503 42588.0 27504.00
# Macrophage    Micro  Mural   Oligo      OPC Tcell
# 0%      1800.00   879.00 1979.0   850.0  2084.00  1774
# 25%     2113.75  3019.25 3210.5  4940.5  7371.75  2493
# 50%     3655.00  3883.50 4692.5  6385.0  9112.00  2668
# 75%     4413.50  4911.75 5218.0  7986.0 11052.75  3031
# 100%    5076.00 11137.00 8043.0 25379.0 23492.00  6919

sapply(newClusIndex, function(x) {quantile(sce.dlpfc$doubletScore[x])}) 
#          Astro  Excit_A   Excit_B  Excit_C    Excit_D  Excit_E  Excit_F  Inhib_A   Inhib_B  Inhib_C   Inhib_D  Inhib_E
# 0%    0.000000 0.000000  0.010588 0.021176  0.0252900  0.09273 0.074116 0.021176  0.000000 0.000000  0.000000 0.021176
# 25%   0.025290 0.042150  0.077878 0.179996  0.1058800  0.21918 0.635280 0.063528  0.052940 0.095292  0.042352 0.031764
# 50%   0.063528 0.109590  0.137644 0.292889  0.2117600  0.37058 1.230780 0.128668  0.118020 0.243524  0.075870 0.042150
# 75%   0.116468 0.370920  0.232936 0.441309  0.4515585  0.56481 1.461144 0.490970  0.243524 0.885150  0.201172 0.107278
# 100% 11.565960 7.972764 26.671172 9.793900 11.0962240 10.52447 7.232940 6.119864 15.780960 8.788040 10.270360 0.145598
#       Inhib_F Macrophage    Micro    Mural     Oligo      OPC    Tcell
# 0%   0.063528   0.006772 0.000000 0.021176  0.000000 0.000000 0.006772
# 25%  0.074116   0.007726 0.003386 0.031764  0.042352 0.016860 0.010588
# 50%  0.153526   0.010588 0.016860 0.052940  0.189616 0.044018 0.010588
# 75%  0.275288   0.010588 0.044018 0.119115  0.584085 0.092730 0.031764
# 100% 1.147854   0.025290 4.324590 0.556380 11.001150 8.502164 0.067440

sapply(newClusIndex, function(x) {median(sce.dlpfc$doubletScore[x])})

table(sce.dlpfc$collapsedCluster)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19 
# 333  454  365  529  773  413  782  524 5455  132  187  243  572   10   18  388    7    8    9 

cell_colors <- cluster_colors[order(as.integer(names(cluster_colors)))]
names(cell_colors) <- annotationTab.dlpfc$cellType
cell_colors
# Inhib_A    Inhib_B    Inhib_C    Excit_A    Excit_B    Inhib_D      Astro    Excit_C      Oligo    Excit_D 
# "#1F77B4"  "#AEC7E8"  "#FF7F0E"  "#FFBB78"  "#2CA02C"  "#98DF8A"  "#D62728"  "#FF9896"  "#9467BD"  "#C5B0D5" 
# Excit_E    Excit_F        OPC Macrophage      Mural      Micro    Inhib_E    Inhib_F      Tcell 
# "#8C564B"  "#C49C94"  "#E377C2"  "#F7B6D2"  "#7F7F7F"  "#C7C7C7"  "#BCBD22"  "#DBDB8D"  "#17BECF" 

# Save
save(sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc, ref.sampleInfo, annotationTab.dlpfc, cell_colors,
     file=here("rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda"))

## Re-print marker expression plots with annotated cluster names ===
pdf("pdfs/revision/regionSpecific_DLPFC-n3_marker-logExprs_collapsedClusters_LAH2021.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpressionCustom(sce = sce.dlpfc,
                         features = markers.mathys.custom[[i]],
                         features_name = names(markers.mathys.custom)[[i]],
                         anno_name = "cellType") +
      scale_color_manual(values = cell_colors)
  )
}
dev.off()

# Reprint some visualizations:
pdf("pdfs/revision/regionSpecific_DLPFC-n3_reducedDims-with-collapsedClusters_LAH2021.pdf")
plotReducedDim(sce.dlpfc, dimred="PCA_corrected", ncomponents=5, colour_by="cellType", point_alpha=0.5)+
  scale_color_manual(values = cell_colors) + labs(colour="Cell type")
plotTSNE(sce.dlpfc, colour_by="sampleID", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="protocol", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="cellType", point_alpha=0.5)+
  scale_color_manual(values = cell_colors) + labs(colour="Cell type")
plotTSNE(sce.dlpfc, colour_by="sum", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="doubletScore", point_alpha=0.5)
# And some more informative UMAPs
plotUMAP(sce.dlpfc, colour_by="sampleID", point_alpha=0.5)
plotUMAP(sce.dlpfc, colour_by="cellType", point_alpha=0.5)+
  scale_color_manual(values = cell_colors) + labs(colour="Cell type")
dev.off()



 ## -> proceed to 'step03_markerDetxn-analyses[...].R'

### For reference === == === == ===
table(sce.dlpfc$cellType, sce.dlpfc$sampleID)
#            br5161.dlpfc br5207.dlpfc br5212.dlpfc
# Astro               371          274          137
# Excit_A             111          298          120
# Excit_B              75          544          154
# Excit_C              44          325          155
# Excit_D              22           83           27
# Excit_E              77           85           25
# Excit_F             102          105           36
# Inhib_A              39          205           89
# Inhib_B              98          250          106
# Inhib_C              47          262           56
# Inhib_D             119          216           78
# Inhib_E               2            3            2
# Inhib_F               0            7            1
# Macrophage            1            6            3
# Micro               152          144           92
# Mural                 3           13            2
# Oligo              2754         2184          517
# OPC                 196          285           91
# Tcell                 2            5            2

table(sce.dlpfc$prelimCluster, sce.dlpfc$sampleID)

annotationTab.dlpfc
# collapsedCluster   cellType
# 1                 1    Inhib_A
# 2                 2    Inhib_B
# 3                 3    Inhib_C
# 4                 4    Excit_A
# 5                 5    Excit_B
# 6                 6    Inhib_D
# 7                 7      Astro
# 8                 8    Excit_C
# 9                 9      Oligo
# 10               10    Excit_D
# 11               11    Excit_E
# 12               12    Excit_F
# 13               13        OPC
# 14               14 Macrophage
# 15               15      Mural
# 16               16      Micro
# 17               17    Inhib_E
# 18               18    Inhib_F
# 19               19      Tcell

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

