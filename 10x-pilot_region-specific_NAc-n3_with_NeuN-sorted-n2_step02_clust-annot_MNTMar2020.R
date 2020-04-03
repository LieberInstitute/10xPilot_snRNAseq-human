### MNT 10x snRNA-seq workflow: step 02
###   **Region-specific analyses**
###     - (3x) NAc samples from: Br5161 & Br5212 & Br5287
###     - (2x) NeuN-sorted samples from: Br5207 & Br5182
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

### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

# ===


## Bring in (2x) NeuN-sorted NAc samples
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/NeuN-sortedNAc_n2_processing-QC_MNTMar2020.rda",
     verbose=T)
    # pilot.neun, pilot.neun.unfiltered, e.out
    rm(e.out, pilot.neun.unfiltered)

## Also load in homogenate NAc samples
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_processing-QC_MNTJan2020.rda",
     verbose=T)
    # pilot.data, pilot.data.unfiltered, e.out
    rm(e.out, pilot.data.unfiltered)



### MNT comment: At this point, each sample (which is a SCE object in the list, 'pilot.data') has been
#              QC'd for cell/nucleus calling ('emptyDrops()' test) and mito rate thresholding


### Merging shared-region samples ============================================
  # Newest iterations for normalization: cbind, THEN take scaled LSFs computed on all nuclei
  # (i.e. no more MBN, because batch is so confounded with sample)


# Add $sample identity
names(pilot.data)
for(i in 1:length(pilot.data)){
  pilot.data[[i]]$sample <- names(pilot.data)[i]
}
# Neun couple
for(i in 1:length(pilot.neun)){
  pilot.neun[[i]]$sample <- names(pilot.neun)[i]
}

# Remove $logcounts in pilot.data
for(i in 1:length(pilot.data)){
  assay(pilot.data[[i]], "logcounts") <- NULL
}

# Since these were processed separately, make sure rownames are same
table(rownames(pilot.data[["nac.5161"]]) == rownames(pilot.neun[["nac.neun.5207"]]))  # all TRUE

# Also remove internal colData, bc 'pilot.data' had size factors previously generated
for(i in 1:length(pilot.data)){
  int_colData(pilot.data[[i]])$size_factor <- NULL
}


sce.nac.all <- cbind(pilot.data[["nac.5161"]], pilot.data[["nac.5212"]], pilot.data[["nac.5287"]], 
                     pilot.neun[["nac.neun.5182"]], pilot.neun[["nac.neun.5207"]])

    # For reference
    table(sce.nac.all$sample)
        #nac.5161      nac.5212      nac.5287 nac.neun.5182 nac.neun.5207
        #    2067          1774           707          4267          4426


# Remove any potential sizeFactors first, kept from previous processing (of 'pilot.data' set)
# (but removing that column from the internal colData should have done this job)
sizeFactors(sce.nac.all) <- NULL

# Generate log-normalized counts
sce.nac.all <- logNormCounts(sce.nac.all)

geneVar.nac.all <- modelGeneVar(sce.nac.all)
chosen.hvgs.nac.all <- geneVar.nac.all$bio > 0
sum(chosen.hvgs.nac.all)
    # [1] 12407


### Dimensionality reduction ================================================================

# Run PCA, taking top 100 (instead of default 50 PCs)
set.seed(109)
sce.nac.all <- runPCA(sce.nac.all, subset_row=chosen.hvgs.nac.all, ncomponents=100,
                  BSPARAM=BiocSingular::RandomParam())

# Save into a new data file, which will dedicate for pan-brain-analyses
save(sce.nac.all, chosen.hvgs.nac.all, file="rdas/NeuN-sortedNAc_n2_cleaned-combined_MNTMar2020.rda")



    # === === === === === === === ===
    ## 'getClusteredPCs()' evaluated in qsub mode (with 'R-batchJob_panBrain_optimalPCselxn_MNTMar2020.R')
    #    --> saved into same .rda

    ## MNT 06Mar2020: skipping this step for now - just work in top 100 PCs
    sum(attr(reducedDim(sce.nac.all), "percentVar")[1:50])  # 45.24% what would be by default
    sum(attr(reducedDim(sce.nac.all), "percentVar")[1:100]) # 46.92% - not much more

    ## *** FIND THIS PCA_opt-skip STEP IN "side-Rscript_exploringALL-NAc-samples_100PCs_MNT06Mar.R"
    ## === === === === === === === ==


    ### skip ==============
    ## Explore that $chosen PC space: "rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda" MNT 23Mar2020
    load("rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda", verbose=T)
        # sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all
    
    metadata(pc.choice.nac.all) # 87
    rm(list=setdiff(ls(), "pc.choice.nac.all"))
    
    # Bring in explored object(s) with clustering, etc. at 100 PCs
    load("rdas/zref_NeuN-sortedNAc-with-homs_exploreTop100pcs_MNT_06Mar2020.rda", verbose=T)
        # sce.nac.all, sce.nac.10pcs, chosen.hvgs.nac.all, ref.sampleInfo, clusterRefTab.nac.all
    
    
    # Add into new reducedDim entry 
    reducedDim(sce.nac.all, "PCA_opt") <- reducedDim(sce.nac.all, "PCA")[ ,1:(metadata(pc.choice.nac.all)$chosen)]
    
    # Run clustering to compare prelimClusters-on-100 PCs
    snn.gr.87 <- buildSNNGraph(sce.nac.all, k=20, use.dimred="PCA_opt")
    clusters.k20.87 <- igraph::cluster_walktrap(snn.gr.87)$membership
    table(clusters.k20.87)
        ## Previously (w/ 100 PCs):
        #    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16
        #  279 1732  275 1948 1948  718  300  293  161  188 1588 1845  133   27   56  857
        #   17   18   19   20   21   22
        #  184  384  136  103   25   61
    
        ## in this "PCA_opt"
        #clusters.k20.87
        #   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16
        #1847  276  266  300 1925 1739  719  161  301 1619 1913  181  132   56  882  183
        #  17   18   19   20   21   22
        # 384   31  139   62  100   25
    
    # To test if same/similar to the above, may have to pseudo-bulk-HC, then compare assignments
    
    sce.nac.all$prelimCluster.87 <- factor(clusters.k20.87)
    table(sce.nac.all$prelimCluster, sce.nac.all$prelimCluster.87)
        ## Kind of an identity line
    table(sce.nac.all$prelimCluster.87, sce.nac.all$lessCollapsed)
        ## For the most part they assign to a unique 'lessCollapsed' cluster
    table(sce.nac.all$prelimCluster.87, sce.nac.all$collapsedCluster)
        ## same here
    
        ## Seems like 'prelimCluster.87' 18 & 20 pertain to the 'ambig.lowNtrxts' cluster in the
         # 100-PC analysis (collapsedCluster 6, or prelimCluster 14 & 22)
    
         # --> Nevertheless, proceed taking the $chosen dimensions for due diligence
         # === === === ===
    
    rm(list=ls())
    # end skip/explore ========================


### Resume work in optimal PC space
load("rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda", verbose=T)
    # sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all

## How many PCs is optimal?:
metadata(pc.choice.nac.all)$chosen
    # 87

# Assign this chosen ( PCs) to 'PCA_opt'
reducedDim(sce.nac.all, "PCA_opt") <- reducedDim(sce.nac.all, "PCA")[ ,1:(metadata(pc.choice.nac.all)$chosen)]


## t-SNE
set.seed(109)
sce.nac.all <- runTSNE(sce.nac.all, dimred="PCA_opt")


## UMAP
set.seed(109)
sce.nac.all <- runUMAP(sce.nac.all, dimred="PCA_opt")


### Add some colData - use "rdas/REFERENCE_sampleInfo_n14.rda" previously generated in explore phase
load("rdas/REFERENCE_sampleInfo_n14.rda", verbose=T)
    # ref.sampleInfo 

# Add to sce.nac.all colData
sce.nac.all$region <- ss(sce.nac.all$sample,"\\.", 1)
sce.nac.all$donor <- paste0("Br",substr(sce.nac.all$sample, start=nchar(sce.nac.all$sample)-3, stop=nchar(sce.nac.all$sample)))
sce.nac.all$processDate <- ref.sampleInfo$realBatch[match(sce.nac.all$sample, ref.sampleInfo$sampleID)]
sce.nac.all$protocol <- ref.sampleInfo$protocol[match(sce.nac.all$sample, ref.sampleInfo$sampleID)]

table(sce.nac.all$protocol, sce.nac.all$donor)
    #           Br5161 Br5182 Br5207 Br5212 Br5287
    #Frank           0      0      0      0    707
    #Frank.NeuN      0   4267   4426      0      0
    #pseudoSort   2067      0      0   1774      0

# Save for now
save(sce.nac.all, chosen.hvgs.nac.all, ref.sampleInfo, pc.choice.nac.all,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda")



### Clustering: Two-step ======================================================
### Step 1: Perform graph-based clustering in this optimal PC space
#         - take k=20 NN to build graph
snn.gr <- buildSNNGraph(sce.nac.all, k=20, use.dimred="PCA_opt")
clusters.k20 <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters.k20)
    ##   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16
    # 1847  276  266  300 1925 1739  719  161  301 1619 1913  181  132   56  882  183
    #   17   18   19   20   21   22
    #  384   31  139   62  100   25   (as above)

# Assign as 'prelimCluster'
sce.nac.all$prelimCluster <- factor(clusters.k20)


### Step 2: Hierarchical clustering of pseudo-bulked ("PB'd") counts with most robust normalization
#         (as determined in: 'side-Rscript_testingStep2_HC-normalizn-approaches_wAmygData_MNTJan2020.R')
#           ** That is, to pseudo-bulk (aka 'cluster-bulk') on raw counts, on all [non-zero] genes,
#              normalize with `librarySizeFactors()`, log2-transform, then perform HC'ing


## Preliminary cluster index for pseudo-bulking
# Will need to split out by single-vs-multiple-nuclei-containing prelimClusters bc
#`rowSums()` doesn't like the former
clusIndexes = splitit(sce.nac.all$prelimCluster)

prelimCluster.PBcounts <- sapply(clusIndexes, function(ii){
  rowSums(assays(sce.nac.all)$counts[ ,ii])
 }
)
    ## of dim 33538    22


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
myplclust(tree.clusCollapsed, main="all NAc sample clusters (n=5) prelim-kNN-cluster relationships",
          cex.main=1, cex.lab=0.8, cex=0.6)

# 18 and 20 clear outliers (as predicted earlier)
sapply(clusIndexes, function(x) {quantile(sce.nac.all[ ,x]$sum)})
  # yep they're driven by low transcript capture


clust.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                               minClusterSize=2, deepSplit=1, cutHeight=400)


table(clust.treeCut)
unname(clust.treeCut[order.dendrogram(dend)])
    ## cutHeight at [400] looks best; proceed with this

## Add new labels to those prelimClusters cut off
#    - going to just put 18 & 20 together (see above)
clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)]==0)] <- max(clust.treeCut)+c(1,1,2)

labels_colors(dend) <- tableau10medium[clust.treeCut[order.dendrogram(dend)]]

# Print for future reference
pdf("pdfs/regionSpecific_NAc-ALL-n5_HC-prelimCluster-relationships_Mar2020.pdf")
par(cex=1.2, font=2)
myplclust(tree.clusCollapsed, lab.col=tableau10medium[clust.treeCut],
          main="All NAc (n=5) prelim-kNN-cluster relationships")
dev.off()


## After deciding to make a lessCollapsed partitioning of these (below):
clust.treeCut[order.dendrogram(dend)] <- clusterRefTab.nac.all$merged.man.b

pdf("pdfs/regionSpecific_NAc-ALL-n5_HC-prelimCluster-relationships_lessCollapsedCols_Mar2020.pdf")
par(cex=1.2, font=2)
myplclust(tree.clusCollapsed, lab.col=tableau20[clust.treeCut],
          main="All NAc (n=5) prelim-kNN-cluster relationships")
dev.off()

# Make reference for new cluster assignment
clusterRefTab.nac.all <- data.frame(origClust=order.dendrogram(dend),
                                    merged=clust.treeCut[order.dendrogram(dend)])

# Add more manual merging/partitioning of these prelim clusters - mainly in collapsedCluster 1/2...
clusterRefTab.nac.all$merged.man <- clusterRefTab.nac.all$merged
clusterRefTab.nac.all$merged.man[clusterRefTab.nac.all$merged %in% c(1:2)] <-  c(1, rep(2,5), 8:14)
    ## Justification for these, generally, is looking at how do samples separate by
     # prelimCluster - if sample-specific, merging is more liberal.

# A less-less-collapsed version to separate D1/D2 prelim clusters
clusterRefTab.nac.all$merged.man.b <- clusterRefTab.nac.all$merged.man
clusterRefTab.nac.all$merged.man.b[5:8] <- c(15,16,15,16)


## Assign as 'collapsedCluster' or 'lessCollapsed'
sce.nac.all$collapsedCluster <- factor(clusterRefTab.nac.all$merged[match(sce.nac.all$prelimCluster,
                                                                          clusterRefTab.nac.all$origClust)])

#sce.nac.all$lessCollapsed <- factor(clusterRefTab.nac.all$merged.man[match(sce.nac.all$prelimCluster,
#                                                                           clusterRefTab.nac.all$origClust)])

# Go with this:
sce.nac.all$lessCollapsed <- factor(clusterRefTab.nac.all$merged.man.b[match(sce.nac.all$prelimCluster,
                                                                               clusterRefTab.nac.all$origClust)])


table(sce.nac.all$collapsedCluster, sce.nac.all$sample)
    #   nac.5161 nac.5212 nac.5287 nac.neun.5182 nac.neun.5207
    # 1      257      299       86          4012          4182
    # 2       18       19       14           248           241
    # 3     1454      854      499             0             0
    # 4      149      384       12             0             0
    # 5       98      104       37             0             0
    # 6       19       42       22             7             3
    # 7       72       72       37             0             0

table(sce.nac.all$lessCollapsed, sce.nac.all$sample)
    # nac.5161 nac.5212 nac.5287 nac.neun.5182 nac.neun.5207
    # 1         2        0        0           117            13
    # 2         0      266        0             0             0
    # 3      1454      854      499             0             0
    # 4       149      384       12             0             0
    # 5        98      104       37             0             0
    # 6        19       42       22             7             3
    # 7        72       72       37             0             0
    # 8        10        3        0           285             3
    # 9         9        6        3           134           148
    # 10       17        8        6           369           319
    # 11        1        3        0            16             5
    # 12        1        1        1            42            11
    # 13        7        7        9            86           167
    # 14        9        8        4           104            58
    # 15      178        2       72          1505          1829
    # 16       41       14        5          1602          1870

## original prelim cluster
table(sce.nac.all$prelimCluster, sce.nac.all$sample)
    #    nac.5161 nac.5212 nac.5287 nac.neun.5182 nac.neun.5207
    # 1       153        0       66          1395           233
    # 2         7        7        9            86           167
    # 3         0      266        0             0             0
    # 4         9        6        3           134           148
    # 5      1417       15      493             0             0
    # 6        25        2        6           110          1596
    # 7        17        8        6           369           319
    # 8       149        0       12             0             0
    # 9        10        3        0           285             3
    # 10       32       13        5          1538            31
    # 11        9        1        0            64          1839
    # 12       72       72       37             0             0
    # 13        2        0        0           117            13
    # 14        1        1        1            42            11
    # 15       37      839        6             0             0
    # 16        9        8        4           104            58
    # 17        0      384        0             0             0
    # 18        1       30        0             0             0
    # 19       98        4       37             0             0
    # 20       18       12       22             7             3
    # 21        0      100        0             0             0
    # 22        1        3        0            16             5



# Print some visualizations:
pdf("pdfs/regionSpecific_NAc-ALL-n5_reducedDims-with-collapsedClusters_Mar2020.pdf")
plotReducedDim(sce.nac.all, dimred="PCA", ncomponents=5, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.nac.all, colour_by="processDate", point_alpha=0.5) + ggtitle("t-SNE on opt PCs", )
plotTSNE(sce.nac.all, colour_by="sample", point_alpha=0.5) + ggtitle("t-SNE on opt PCs", )
plotTSNE(sce.nac.all, colour_by="collapsedCluster", point_alpha=0.5) + ggtitle("t-SNE on opt PCs", )
plotTSNE(sce.nac.all, colour_by="lessCollapsed", point_alpha=0.5) + ggtitle("t-SNE on opt PCs", )
plotTSNE(sce.nac.all, colour_by="sum", point_alpha=0.5) + ggtitle("t-SNE on opt PCs", )
# UMAP
plotUMAP(sce.nac.all, colour_by="collapsedCluster", point_alpha=0.5) + ggtitle("UMAP on opt PCs", )
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
  # Kristen's MSN markers - not printed for the broader collapsed clusters
  'MSNs.D1' = c("DRD1", "PDYN", "TAC1"),
  'MSNs.D2' = c("DRD2", "PENK"),
  'MSNs.pan' = c("PPP1R1B","BCL11B")# "CTIP2")
)

## Print broad cell type markers
#pdf("pdfs/regionSpecific_NAc-ALL-n5_marker-logExprs_collapsedClusters_Mar2020.pdf", height=6, width=12)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.nac.all, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="lessCollapsed", colour_by="lessCollapsed", point_alpha=0.5, point_size=.7,
                   add_legend=F) +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
                   width = 0.3, colour=rep(tableau20[1:16], length(markers.mathys.custom[[i]])))
  )
}
#dev.off()


# Check nuclear library sizes for annotation of 'lessCollapsed' cluster 6:
clusIndexes <- splitit(sce.nac.all$lessCollapsed)
sapply(clusIndexes, function(x) {quantile(sce.nac.all[ ,x]$sum)})
    #             1        2       3     4       5    6     7     8        9      10
    # 0%    1422.00  2471.00   365.0   635  1379.0  101  1247  2371  1197.00  2739.0
    # 25%  12988.75 19138.00  3618.5  6305  6937.5  133  3205 12217 19541.50 16635.5
    # 50%  16667.50 27300.00  5569.0  8380  9145.0  284  4525 15150 23778.00 20483.0
    # 75%  21967.75 37641.75  8010.5 10767 11354.0  608  6323 18557 30245.25 24723.5
    # 100% 44444.00 68589.00 31312.0 32867 31898.0 1648 12537 39311 69559.00 67355.0
    #         11       12      13    14       15       16
    # 0%    6603  6143.00  3473.0  1388  3895.00  1801.00
    # 25%  15968 11429.75 14588.5 10748 20340.00 18638.50
    # 50%  20199 15010.50 18914.0 13472 24847.50 22881.50
    # 75%  33506 24642.25 24444.5 17118 30074.75 28108.75
    # 100% 60522 68219.00 67288.0 43590 85460.00 84761.00


## Add annotations for 'lessCollapsed' clusters, looking at marker gene expression
annotationTab.nac.all <- data.frame(cluster=c(1, 2, 3, 4, 5,
                                              6, 7, 8, 9, 10,
                                              11, 12, 13, 14, 15,
                                              16),
                                    cellType.split=c("MSN.D1.1", "MSN.broad", "Oligo", "Astro", "OPC",
                                               "ambig.lowNtrxts", "Micro", "MSN.D1.2", "MSN.D2.1","MSN.D1.3",
                                               "Inhib.1","Inhib.2","Inhib.3","Inhib.4","MSN.D1.4",
                                               "MSN.D2.2"),
                                    cellType.broad=c("MSN.D1", "MSN.broad", "Oligo", "Astro", "OPC",
                                                     "ambig.lowNtrxts", "Micro", "MSN.D1", "MSN.D2", "MSN.D1",
                                                     rep("Inhib", 4), "MSN.D1",
                                                     "MSN.D2")
)

# Add info to clusterRefTab.nac.all
clusterRefTab.nac.all$cellType.broad <- annotationTab.nac.all$cellType.broad[match(clusterRefTab.nac.all$merged.man.b,
                                                                                   annotationTab.nac.all$cluster)]
clusterRefTab.nac.all$cellType.split <- annotationTab.nac.all$cellType.split[match(clusterRefTab.nac.all$merged.man.b,
                                                                                   annotationTab.nac.all$cluster)]


## Then add to SCE - like for DLPFC, give 'cellType' & 'cellType.split'
sce.nac.all$cellType <- annotationTab.nac.all$cellType.broad[match(sce.nac.all$lessCollapsed,
                                                                annotationTab.nac.all$cluster)]

sce.nac.all$cellType.split <- annotationTab.nac.all$cellType.split[match(sce.nac.all$lessCollapsed,
                                                                 annotationTab.nac.all$cluster)]


## Save
save(sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all, ref.sampleInfo, 
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda")


# Re-print marker expression with cell type labels and dropping 'ambig.lowNtrxts' cluster
table(sce.dlpfc$cellType.split)

# First drop "Ambig.lowNtrxts" (93 nuclei)
sce.nac.all <- sce.nac.all[ ,sce.nac.all$cellType.split != "ambig.lowNtrxts"]
sce.nac.all$cellType.split <- droplevels(sce.nac.all$cellType.split)

pdf("pdfs/regionSpecific_NAc-ALL-n5_marker-logExprs_collapsedClusters_Mar2020.pdf", height=6, width=12)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.nac.all, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7,
                   add_legend=F) +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
                   width = 0.3, colour=rep(tableau20[1:15], length(markers.mathys.custom[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(markers.mathys.custom)[i], " markers"))
  )
}
dev.off()


## Re-print some reducedDims with this information
pdf("pdfs/regionSpecific_NAc-ALL-n5_reducedDims-with-collapsedClusters_Apr2020.pdf")
plotReducedDim(sce.nac.all, dimred="PCA", ncomponents=5, colour_by="cellType", point_alpha=0.5)
plotTSNE(sce.nac.all, colour_by="processDate", point_size=3.5, point_alpha=0.5) + ggtitle("t-SNE on opt PCs")
plotTSNE(sce.nac.all, colour_by="sample", point_size=3.5, point_alpha=0.5) + ggtitle("t-SNE on opt PCs")
plotTSNE(sce.nac.all, colour_by="cellType", point_size=3.5, point_alpha=0.5) + ggtitle("t-SNE on opt PCs")
plotTSNE(sce.nac.all, colour_by="cellType.split", point_size=3.5, point_alpha=0.5) + ggtitle("t-SNE on opt PCs")
plotTSNE(sce.nac.all, colour_by="sum", point_size=3.5, point_alpha=0.5) + ggtitle("t-SNE on opt PCs")
# UMAP
plotUMAP(sce.nac.all, colour_by="cellType.split", point_size=3.5, point_alpha=0.5) + ggtitle("UMAP on opt PCs")
dev.off()

      ## -> proceed to 'step03_markerDetxn-analyses[...].R'




## For Day Lab - MNT 03Apr2020
plotExpression(sce.nac.all, exprs_values = "logcounts", features="RELN",
               x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7,
               add_legend=F) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
               width = 0.3, colour=rep(tableau20[1:15], 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


## CHAT? - for cholinergic interneurons
grep("CHAT", rownames(sce.nac.all))
plotExpression(sce.nac.all, exprs_values = "logcounts", features="CHAT",
               x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7,
               add_legend=F) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
               width = 0.3, colour=rep(tableau20[1:15], 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


plotExpression(sce.nac.all, exprs_values = "logcounts", features="PVALB",
               x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7,
               add_legend=F) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", 
               width = 0.3, colour=rep(tableau20[1:15], 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



### Further dive into 'MSN.broad' cluster, which is entirely Br5212
  # -> see if 



