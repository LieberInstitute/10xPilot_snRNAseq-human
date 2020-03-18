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
    


### Picking up with optimally-defined PC space ===
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/NeuN-sortedNAc_n2_cleaned-combined_MNTMar2020.rda",
     verbose=TRUE)
    # sce.nac.all, chosen.hvgs.nac.all, pc.choice.n12


# How many PCs is optimal?:
metadata(pc.choice.n12)$chosen
    ## 204

## Assign this chosen ( PCs) to 'PCA_opt'
reducedDim(sce.nac.all, "PCA_opt") <- reducedDim(sce.nac.all, "PCA")[ ,1:(metadata(pc.choice.n12)$chosen)]


## t-SNE
set.seed(109)
sce.nac.all <- runTSNE(sce.nac.all, dimred="PCA_opt")


## UMAP
set.seed(109)
sce.nac.all <- runUMAP(sce.nac.all, dimred="PCA_opt")


## Add some phenodata
sce.nac.all$region <- ss(sce.nac.all$sample,".5",1)
sce.nac.all$donor <- paste0("Br",ss(sce.nac.all$sample,"y.",2))

df.processDate <- data.frame(sample=unique(sce.nac.all$sample))
df.processDate$date <- ifelse(df.processDate$sample =="hpc.5161","08May19", "23Jul19")
df.processDate$date

sce.nac.all$processDate <- ifelse(sce.nac.all$sample=="amy.5161", "04Sep19", "25Sep19")


## Make some pd to add to colData
ref.sampleInfo <- data.frame(sampleID = unique(sce.nac.all$sample))
ref.sampleInfo$realBatch <- ifelse(ref.sampleInfo$sampleID %in% c("dlpfc.5212", "hpc.5287", "nac.5161", "nac.5212"),
                                   "R2.Jul23", "R1.May08")
ref.sampleInfo$realBatch <- ifelse(ref.sampleInfo$sampleID %in% c("dlpfc.5161", "hpc.5212", "nac.5287", "amy.5161"),
                                   "R3.Sep04", ref.sampleInfo$realBatch)
ref.sampleInfo$realBatch <- ifelse(ref.sampleInfo$sampleID %in% c("amy.5212", "sacc.5161", "sacc.5212", "nac.5182", "nac.5207"),
                                   "R4.Sep25", ref.sampleInfo$realBatch)

ref.sampleInfo$protocol <- "Frank"
ref.sampleInfo$protocol <- ifelse(ref.sampleInfo$realBatch=="R2.Jul23", "pseudoSort", ref.sampleInfo$protocol)

# Add to sce.nac.all colData
sce.nac.all$region <- ss(sce.nac.all$sample,".5",1)
sce.nac.all$donor <- paste0("Br",substr(sce.nac.all$sample, start=nchar(sce.nac.all$sample)-3, stop=nchar(sce.nac.all$sample)))
sce.nac.all$processDate <- ref.sampleInfo$realBatch[match(sce.nac.all$sample, ref.sampleInfo$sampleID)]
sce.nac.all$protocol <- ref.sampleInfo$protocol[match(sce.nac.all$processDate, ref.sampleInfo$realBatch)]


# Save for now
save(sce.nac.all, chosen.hvgs.nac.all, pc.choice.n12, ref.sampleInfo,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/NeuN-sortedNAc_n2_cleaned-combined_MNTMar2020.rda")


### Clustering: Two-step ======================================================
### Step 1: Perform graph-based clustering in this optimal PC space
#         - take k=20 NN to build graph
snn.gr <- buildSNNGraph(sce.nac.all, k=20, use.dimred="PCA_opt")
clusters.k20 <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters.k20)
    ## 

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
clusIndexes.single <- unlist(clusIndexes[which(lengths(clusIndexes)==1)])
clusIndexes.mult <- clusIndexes[which(lengths(clusIndexes)>1)]
    

prelimCluster.PBcounts.mult <- sapply(clusIndexes.mult, function(ii){
  rowSums(assays(sce.nac.all)$counts[ ,ii])
  }
)
    ## of dim [1] 33538   242

prelimCluster.PBcounts <- cbind(prelimCluster.PBcounts.mult,
                                as.matrix(assays(sce.nac.all)$counts[ ,clusIndexes.single]))
colnames(prelimCluster.PBcounts)[243:247] <- as.character(243:247)

    # And btw...
    table(rowSums(prelimCluster.PBcounts)==0)  # 3500 TRUE


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
myplclust(tree.clusCollapsed, main="Pan-region (n=12) prelim-kNN-cluster relationships",
          cex.main=1, cex.lab=0.8, cex=0.6)


clust.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                               minClusterSize=2, deepSplit=1, cutHeight=900)


table(clust.treeCut)
unname(clust.treeCut[order.dendrogram(dend)])
    ## cutHeight at 900 looks best; proceed with this

## Add new labels to those prelimClusters cut off   - unneededfor pan-brain analysis
#clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)]==0)] <- max(clust.treeCut)+c(1)

labels_colors(dend) <- tableau20[clust.treeCut[order.dendrogram(dend)]]

# Print for future reference
pdf("pdfs/panBrain-n12_HC-prelimCluster-relationships_Feb2020.pdf", width=13, height=6)
par(cex=0.45, font=2, mar=c(5,7,4,0))
#plot(dend, main="Pan-region (n=12) prelim-kNN-cluster relationships", cex.main=2.5)
myplclust(tree.clusCollapsed, lab.col=tableau20[clust.treeCut],
          main="Pan-region (n=12) prelim-kNN-cluster relationships", cex.main=2.5)
dev.off()


# Make reference for new cluster assignment
clusterRefTab.nac.all <- data.frame(origClust=order.dendrogram(dend),
                                merged=clust.treeCut[order.dendrogram(dend)])


# Assign as 'collapsedCluster'
sce.nac.all$collapsedCluster <- factor(clusterRefTab.nac.all$merged[match(sce.nac.all$prelimCluster,
                                                                          clusterRefTab.nac.all$origClust)])

# Print some visualizations:
pdf("pdfs/panBrain-n12_reducedDims-with-collapsedClusters_Feb2020.pdf")
plotReducedDim(sce.nac.all, dimred="PCA", ncomponents=5, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.nac.all, colour_by="region", point_alpha=0.5)
plotTSNE(sce.nac.all, colour_by="processDate", point_alpha=0.5)
plotTSNE(sce.nac.all, colour_by="sample", point_alpha=0.5)
plotTSNE(sce.nac.all, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.nac.all, colour_by="sum", point_alpha=0.5)
plotUMAP(sce.nac.all, colour_by="collapsedCluster", point_alpha=0.5)
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

pdf("pdfs/panBrain-n12_marker-logExprs_collapsedClusters_Feb2020.pdf", height=6, width=12)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.nac.all, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="collapsedCluster", colour_by="collapsedCluster", point_alpha=0.5, point_size=.7,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:19], length(markers.mathys.custom[[i]])))
  )
}
dev.off()

# Any driven by low n transcripts, overall?
newClusIndex <- splitit(sce.nac.all$collapsedCluster)
sapply(newClusIndex, function(x) {quantile(sce.nac.all[,x]$sum)})
    #

## Add annotations, looking at marker gene expression
annotationTab.pan <- data.frame(cluster=c(1, 2, 3, 4, 5,
                                          6, 7, 8, 9, 10,
                                          11, 12, 13, 14, 15,
                                          16, 17, 18, 19),
                                cellType=c("Excit.1", "Inhib.1", "Excit.2", "Astro.1", "Inhib.2",
                                           "Excit.3", "Inhib.3", "Inhib.4", "Excit.4", "OPC",
                                           "Oligo","Micro","Inhib.5","Excit.5","Excit.6",
                                           "Ambig.hiVCAN","Excit.7","Astro.2","Excit.8")
                                )

sce.nac.all$cellType <- annotationTab.pan$cellType[match(sce.nac.all$collapsedCluster,
                                                         annotationTab.pan$cluster)]

## Save
save(sce.nac.all, chosen.hvgs.nac.all, pc.choice.n12, ref.sampleInfo, clusterRefTab.nac.all,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/NeuN-sortedNAc_n2_cleaned-combined_MNTMar2020.rda")




## Misc exploration ======================================================================

# Is 16 maybe pericytes? Looking at genes from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5057134/
plotExpression(sce.nac.all, exprs_values = "logcounts", features=c("VTN", "IFITM1"),
               x="collapsedCluster", colour_by="collapsedCluster", point_alpha=0.5, point_size=.7,
               add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                            geom = "crossbar", width = 0.3,
                                            colour=rep(tableau20[1:19], 2))
    ## (no expression of either of these)

## How well correlated is the expression of those Ambig.hiVCAN nuclei anyway?
 #    Let's say compared to "Astro.2", which also has a low distribution of n transcripts
 #    (collapsedCluster 18)
library(pheatmap)

table(sce.nac.all$cellType, sce.nac.all$region)

clustIdx.18 <- which(sce.nac.all$collapsedCluster==18)  # "Astro.2"
clustIdx.16 <- which(sce.nac.all$collapsedCluster==16)


quantile(cor(as.matrix(assay(sce.nac.all, "logcounts")[ ,clustIdx.18])), probs=seq(0.05,1,by=0.05))
    #         5%       10%       15%       20%       25%       30%       35%       40%
    #  0.1435597 0.1550650 0.1635447 0.1710807 0.1772498 0.1833153 0.1881462 0.1927772
    #        45%       50%       55%       60%       65%       70%       75%       80%
    #  0.1974163 0.2016515 0.2070635 0.2124137 0.2177567 0.2237208 0.2296916 0.2373880
    #        85%       90%       95%      100%
    #  0.2479811 0.2656165 0.3195865 1.0000000      - so actually pretty poor...

pheatmap(cor(as.matrix(assay(sce.nac.all, "logcounts")[ ,clustIdx.18])),
         labels_row=F, labels_col=F)


quantile(cor(as.matrix(assay(sce.nac.all, "logcounts")[ ,clustIdx.16])), probs=seq(0.05,1,by=0.05))
    #         5%       10%       15%       20%       25%       30%       35%       40%
    #  0.1139671 0.1440696 0.1624547 0.1762752 0.1880453 0.1962381 0.2028947 0.2085750
    #        45%       50%       55%       60%       65%       70%       75%       80%
    #  0.2163537 0.2201661 0.2272381 0.2393277 0.2517573 0.2621352 0.2746192 0.2888875
    #        85%       90%       95%      100%
    #  0.3063891 0.3572942 0.4383301 1.0000000      - a little better with this 'Ambig.hiVCAN'

pheatmap(cor(as.matrix(assay(sce.nac.all, "logcounts")[ ,clustIdx.inhib2])),
         labels_row=F, labels_col=F)


# Showing pairwise correlation is good with more robust clusters...
clustIdx.inhib2 <- which(sce.nac.all$cellType=="Inhib.2")
pheatmap(cor(as.matrix(assay(sce.nac.all, "logcounts")[ ,clustIdx.inhib2])),
         labels_row=F, labels_col=F)
    ## much better
quantile(cor(as.matrix(assay(sce.nac.all, "logcounts")[ ,clustIdx.inhib2])), probs=seq(0.05,1,by=0.05))
    #         5%       10%       15%       20%       25%       30%       35%       40%
    #  0.6184637 0.6865806 0.7100607 0.7239689 0.7342448 0.7425027 0.7498117 0.7564875
    #        45%       50%       55%       60%       65%       70%       75%       80%
    #  0.7629176 0.7689357 0.7750985 0.7810500 0.7871078 0.7932005 0.7995850 0.8064660
    #        85%       90%       95%      100%
    #  0.8138758 0.8225428 0.8345357 1.0000000

    # Is this because it's a bigger cluster?
    clustIdx.excit4 <- which(sce.nac.all$cellType=="Excit.4")
    pheatmap(cor(as.matrix(assay(sce.nac.all, "logcounts")[ ,clustIdx.excit4])),
         labels_row=F, labels_col=F)
    quantile(cor(as.matrix(assay(sce.nac.all, "logcounts")[ ,clustIdx.excit4])), probs=seq(0.05,1,by=0.05))
        ## less good as "Inhib.2", but still much better than collapsedCluster 16 & 18

    ## What about collapsedCluster 19?  On second glance, it expresses many markers...
     #      (notably VCAN; some PLP1; and more SLC17A6 than SLC17A7)
    clustIdx.excit8 <- which(sce.nac.all$cellType=="Excit.8")
    pheatmap(cor(as.matrix(assay(sce.nac.all, "logcounts")[ ,clustIdx.excit8])),
             labels_row=F, labels_col=F)
    quantile(cor(as.matrix(assay(sce.nac.all, "logcounts")[ ,clustIdx.excit8])), probs=seq(0.05,1,by=0.05))
        ## Overall still pretty good, but less than "Excit.4"
    
## For reference === === === === ===
table(sce.nac.all$cellType, sce.nac.all$region)
    
    #





