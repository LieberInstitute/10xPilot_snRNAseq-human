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


# Load 'pilot' samples
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/all-FACS-n14_preprint_SCEs_processing-QC_MNTMar2021.rda",
     verbose=T)
    # pilot.data, pilot.data.unfiltered, e.out, ref.sampleInfo
    rm(pilot.data.unfiltered, e.out)

    ### MNT comment: At this point, each sample (which is a SCE object in the list, 'pilot.data') has been
    #               QC'd for cell/nucleus calling ('emptyDrops()' test) and mito rate thresholding
    #   Additionally, there have been a computed 'doubletScore', which will QC with, after
    #   clustering (e.g. that there are no doublet-driven clusters, etc.)

    
### Merging shared-region samples ============================================
  # Newest iterations for normalization: multiBatchNorm-alize

sce.hpc <- cbind(pilot.data[["br5161.hpc"]], pilot.data[["br5212.hpc"]], pilot.data[["br5287.hpc"]])

sce.hpc
    # class: SingleCellExperiment 
    # dim: 33538 10268 
    # metadata(3): Samples Samples Samples
    # assays(1): counts
    # rownames(33538): MIR1302-2HG FAM138A ... AC213203.1 FAM231C
    # rowData names(6): gene_id gene_version ... gene_biotype Symbol.uniq
    # colnames(10268): AAACCCATCTGTCAGA-1 AAACCCATCTGTCGCT-1 ...
    #   TTTGGTTGTGGTCCGT-1 TTTGTTGCAGAAACCG-1
    # colData names(16): Sample Barcode ... protocol sequencer
    # reducedDimNames(0):
    # altExpNames(0):    

# Use `multiBatchNorm()` to compute log-normalized counts, matching the scaling across samples
sce.hpc <- multiBatchNorm(sce.hpc, batch=sce.hpc$sampleID)

# Use the simple `modelGeneVar` - this makes more sense over `combineVar`, since the
#   cell composition is already known to be quite different (with NeuN selection)
geneVar.hpc <- modelGeneVar(sce.hpc)
chosen.hvgs.hpc <- geneVar.hpc$bio > 0
sum(chosen.hvgs.hpc)
    # [1] 8696


### Dimensionality reduction ================================================================

# Run `fastMNN` (internally uses `multiBatchPCA`), taking top 100 (instead of default 50 PCs)
set.seed(109)
mnn.hold <-  fastMNN(sce.hpc, batch=sce.hpc$sampleID,
                     merge.order=c("br5161.hpc","br5212.hpc","br5287.hpc"),
                     subset.row=chosen.hvgs.hpc, d=100,
                     correct.all=TRUE, get.variance=TRUE,
                     BSPARAM=BiocSingular::IrlbaParam())
    # This temp file just used for getting batch-corrected components (drops a variety of entries)

table(colnames(mnn.hold) == colnames(sce.hpc))  # all TRUE
table(mnn.hold$batch == sce.hpc$sampleID) # all TRUE

# Add them to the SCE, as well as the metadata (though the latter might not be so usefl)
reducedDim(sce.hpc, "PCA_corrected") <- reducedDim(mnn.hold, "corrected") # 100 components
metadata(sce.hpc) <- metadata(mnn.hold)

# Save into a new region-specific SCE object/flie
save(sce.hpc, chosen.hvgs.hpc, ref.sampleInfo,
     file="rdas/revision/regionSpecific_HPC-n3_cleaned-combined_SCE_MNT2021.rda")


    ## 'getClusteredPCs()' evaluated in qsub mode (with 'R-batchJob_HPC-n3_optimalPCselxn_MNTFeb2020.R')
    #    --> saved into same .rda


### Picking up with optimally-defined PC space ===
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_HPC-n3_cleaned-combined_SCE_MNT2021.rda",
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
                               minClusterSize=2, deepSplit=1, cutHeight=214)


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

# First drop "Ambig.lowNtrxts" (101 nuclei)
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



### MNT 30Apr2020: subcluster-level annotations ========================================
  # Some motivations:
  #   i) to see if the high-VCAN-neuronal cluster exists in prelimCluster
  #      as seen in the pan-brain level
  #   ii) Tcell cluster ('Ambig.glial') is probably its own prelimCluster

load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo

table(sce.hpc$cellType)
table(sce.hpc$prelimCluster, sce.hpc$cellType)
    # yes 'Tcell' (26 nuclei) cluster is entirely prelimCluster 31
table(sce.hpc$cellType, sce.hpc$donor)
    #     and comes from all three donors

# Look at:   (& the HC dendrogram, printed before)
table(sce.hpc$prelimCluster, sce.hpc$sample)


## More-manual annotations: ====
    clusterRefTab.hpc$cellType <- sce.hpc$cellType[match(clusterRefTab.hpc$merged, sce.hpc$collapsedCluster)]
    clusterRefTab.hpc$cellType <- as.character(clusterRefTab.hpc$cellType)
    # Rename 'Ambig.glial' -> 'Tcell' from what we've seen with downstream PB-modeled markers
    clusterRefTab.hpc$cellType[clusterRefTab.hpc$cellType=="Ambig.glial"] <- "Tcell"
    
    # Make new column for subclusers
    clusterRefTab.hpc$manual <- clusterRefTab.hpc$cellType
    # Inhib gets all split up
    clusterRefTab.hpc$manual <- ifelse(clusterRefTab.hpc$cellType == "Inhib",
                                         paste0(clusterRefTab.hpc$cellType, ".", c(1:5)),
                                         as.character(clusterRefTab.hpc$manual))
    # Excit 
    clusterRefTab.hpc$manual <- ifelse(clusterRefTab.hpc$origClust %in% c(9,25),
                                         paste0(clusterRefTab.hpc$cellType, ".1"),
                                         as.character(clusterRefTab.hpc$manual))
    
    clusterRefTab.hpc$manual <- ifelse(clusterRefTab.hpc$origClust %in% c(12,24),
                                       paste0(clusterRefTab.hpc$cellType, ".2"),
                                       as.character(clusterRefTab.hpc$manual))
    
    clusterRefTab.hpc$manual <- ifelse(clusterRefTab.hpc$origClust %in% c(15,30),
                                       paste0(clusterRefTab.hpc$cellType, ".3"),
                                       as.character(clusterRefTab.hpc$manual))
    
    clusterRefTab.hpc$manual <- ifelse(clusterRefTab.hpc$origClust %in% c(27),
                                       paste0(clusterRefTab.hpc$cellType, ".4"),
                                       as.character(clusterRefTab.hpc$manual))
    
    clusterRefTab.hpc$manual <- ifelse(clusterRefTab.hpc$origClust %in% c(14),
                                       paste0(clusterRefTab.hpc$cellType, ".5"),
                                       as.character(clusterRefTab.hpc$manual))

        ## All other glial types will be kept the same
    
    ## Post-hoc: Looks like Excit.5 truly inhibitory, and Inhib.1 truly excitatory
     #           (the latter is the 33 high-VCAN nuclei ID'd in this sample at pan-brain) 
    
        clusterRefTab.hpc$manual <- ifelse(clusterRefTab.hpc$manual == "Excit.5",
                                           "Inhib.1",
                                           as.character(clusterRefTab.hpc$manual))
        clusterRefTab.hpc$manual <- ifelse(clusterRefTab.hpc$origClust == 28,
                                           "Excit.5",
                                           as.character(clusterRefTab.hpc$manual))
        
        # --> THEN re-run the below
    
    # end subcluster annotation ====

    
## Add new annotations
sce.hpc$cellType.split <- clusterRefTab.hpc$manual[match(sce.hpc$prelimCluster,
                                                         clusterRefTab.hpc$origClust)]
sce.hpc$cellType.split <- factor(sce.hpc$cellType.split)

table(sce.hpc$cellType.split, sce.hpc$cellType)
    # good

table(sce.hpc$cellType.split) # (printing post-hoc-corrected annotations)
    #Ambig.lowNtrxts           Astro         Excit.1         Excit.2         Excit.3
    #            101            1343             116             117             310
    #        Excit.4         Excit.5         Inhib.1         Inhib.2         Inhib.3
    #             26              33              30              90             139
    #        Inhib.4         Inhib.5           Micro           Oligo             OPC
    #             55              56            1253            5885             864
    #          Tcell
    #             26

## Save these
save(sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda")



## Also print expression at this level of partitioning


# First remove "Ambig.lowNtrxts":
sce.hpc <- sce.hpc[ ,sce.hpc$cellType.split != "Ambig.lowNtrxts"]
sce.hpc$cellType.split <- droplevels(sce.hpc$cellType.split)

pdf("pdfs/regionSpecific_HPC-n3_marker-logExprs_cellTypesSplit_Apr2020.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.hpc, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:15], length(markers.mathys.custom[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(markers.mathys.custom)[i], " markers"))
  )
}
dev.off()


## Let's also re-plot reducedDims with new [broad & split] cell type annotations
#        (and rename old file with prefix 'zold_')
pdf("pdfs/regionSpecific_HPC-n3_reducedDims-with-collapsedClusters_Apr2020.pdf")
plotReducedDim(sce.hpc, dimred="PCA", ncomponents=5, colour_by="cellType", point_alpha=0.5)
plotTSNE(sce.hpc, colour_by="sample", point_size=3.5, point_alpha=0.5)
plotTSNE(sce.hpc, colour_by="prelimCluster", point_size=3.5, point_alpha=0.5)
plotTSNE(sce.hpc, colour_by="cellType", point_size=3.5, point_alpha=0.5)
plotTSNE(sce.hpc, colour_by="cellType.split", point_size=3.5, point_alpha=0.5)
plotTSNE(sce.hpc, colour_by="sum", point_size=3.5, point_alpha=0.5)
plotUMAP(sce.hpc, colour_by="cellType", point_size=3.5, point_alpha=0.5)
plotUMAP(sce.hpc, colour_by="cellType.split", point_size=3.5, point_alpha=0.5)
dev.off()


## And finally, for reference:
table(sce.hpc$cellType.split, sce.hpc$sample)
    #         hpc.5161 hpc.5212 hpc.5287
    # Astro        584      582      177
    # Excit.1        6        9      101
    # Excit.2      117        0        0
    # Excit.3        5      289       16
    # Excit.4        0       21        5
    # Excit.5       33        0        0
    # Inhib.1       30        0        0
    # Inhib.2       34       26       30
    # Inhib.3       72       37       30
    # Inhib.4       39       11        5
    # Inhib.5       28       18       10
    # Micro        520      536      197
    # Oligo       2594     2198     1093
    # OPC          395      263      206
    # Tcell         13       10        3


### Session info for 23Apr2021 ==============================
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
# [7] scran_1.18.5                scater_1.18.6              
# [9] ggplot2_3.3.3               EnsDb.Hsapiens.v86_2.99.0  
# [11] ensembldb_2.14.1            AnnotationFilter_1.14.0    
# [13] GenomicFeatures_1.42.3      AnnotationDbi_1.52.0       
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
# [9] rstudioapi_0.13           bit64_4.0.5              
# [11] fansi_0.4.2               xml2_1.3.2               
# [13] splines_4.0.4             R.methodsS3_1.8.1        
# [15] sparseMatrixStats_1.2.1   cachem_1.0.4             
# [17] Rsamtools_2.6.0           ResidualMatrix_1.0.0     
# [19] dbplyr_2.1.1              R.oo_1.24.0              
# [21] HDF5Array_1.18.1          compiler_4.0.4           
# [23] httr_1.4.2                dqrng_0.2.1              
# [25] assertthat_0.2.1          Matrix_1.3-2             
# [27] fastmap_1.1.0             lazyeval_0.2.2           
# [29] limma_3.46.0              BiocSingular_1.6.0       
# [31] prettyunits_1.1.1         tools_4.0.4              
# [33] rsvd_1.0.3                igraph_1.2.6             
# [35] gtable_0.3.0              glue_1.4.2               
# [37] GenomeInfoDbData_1.2.4    dplyr_1.0.5              
# [39] rappdirs_0.3.3            Rcpp_1.0.6               
# [41] vctrs_0.3.6               Biostrings_2.58.0        
# [43] rhdf5filters_1.2.0        rtracklayer_1.50.0       
# [45] DelayedMatrixStats_1.12.3 stringr_1.4.0            
# [47] beachmat_2.6.4            lifecycle_1.0.0          
# [49] irlba_2.3.3               statmod_1.4.35           
# [51] XML_3.99-0.6              edgeR_3.32.1             
# [53] zlibbioc_1.36.0           scales_1.1.1             
# [55] hms_1.0.0                 ProtGenerics_1.22.0      
# [57] rhdf5_2.34.0              RColorBrewer_1.1-2       
# [59] curl_4.3                  memoise_2.0.0            
# [61] gridExtra_2.3             segmented_1.3-3          
# [63] biomaRt_2.46.3            stringi_1.5.3            
# [65] RSQLite_2.2.7             BiocParallel_1.24.1      
# [67] rlang_0.4.10              pkgconfig_2.0.3          
# [69] bitops_1.0-6              lattice_0.20-41          
# [71] purrr_0.3.4               Rhdf5lib_1.12.1          
# [73] GenomicAlignments_1.26.0  bit_4.0.4                
# [75] tidyselect_1.1.0          magrittr_2.0.1           
# [77] R6_2.5.0                  generics_0.1.0           
# [79] DelayedArray_0.16.3       DBI_1.1.1                
# [81] pillar_1.6.0              withr_2.4.2              
# [83] RCurl_1.98-1.3            tibble_3.1.1             
# [85] crayon_1.4.1              utf8_1.2.1               
# [87] BiocFileCache_1.14.0      viridis_0.6.0            
# [89] progress_1.2.2            locfit_1.5-9.4           
# [91] grid_4.0.4                blob_1.2.1               
# [93] R.utils_2.10.1            openssl_1.4.3            
# [95] munsell_0.5.0             beeswarm_0.3.1           
# [97] viridisLite_0.4.0         vipor_0.4.5              
# [99] askpass_1.1 


