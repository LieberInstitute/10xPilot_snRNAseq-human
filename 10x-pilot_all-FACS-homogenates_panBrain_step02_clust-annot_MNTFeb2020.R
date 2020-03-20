### MNT 10x snRNA-seq workflow: step 02
###   **PAN-BRAIN analyses**
###     - (n=12) all regions from: Br5161 & Br5212 & Br5287
###     - Amyg, DLPFC, HPC, NAc, and sACC
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
  # Newest iterations for normalization: cbind, THEN take scaled LSFs computed on all nuclei
  # (i.e. no more MBN, because batch is so confounded with sample)

# Add $sample identity
for(i in 1:length(pilot.data)){
  pilot.data[[i]]$sample <- names(pilot.data)[i]
}

sce.all.n12 <- cbind(pilot.data[[1]], pilot.data[[2]], pilot.data[[3]], pilot.data[[4]],
                     pilot.data[[5]], pilot.data[[6]], pilot.data[[7]], pilot.data[[8]],
                     pilot.data[[9]], pilot.data[[10]], pilot.data[[11]], pilot.data[[12]])

# Remove $logcounts
assay(sce.all.n12, "logcounts") <- NULL
# Re-generate log-normalized counts
sce.all.n12 <- logNormCounts(sce.all.n12)

geneVar.all.n12 <- modelGeneVar(sce.all.n12)
chosen.hvgs.all.n12 <- geneVar.all.n12$bio > 0
sum(chosen.hvgs.all.n12)
    # [1] 8748    - kinda surprising...


### Dimensionality reduction ================================================================

# Run PCA, taking top 250 (instead of default 50 PCs)
set.seed(109)
sce.all.n12 <- runPCA(sce.all.n12, subset_row=chosen.hvgs.all.n12, ncomponents=250,
                  BSPARAM=BiocSingular::RandomParam())

# Save into a new data file, which will dedicate for pan-brain-analyses
save(sce.all.n12, chosen.hvgs.all.n12, file="rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda")



    # === === === === === === === ===
    ## 'getClusteredPCs()' evaluated in qsub mode (with 'R-batchJob_panBrain_optimalPCselxn_MNTFeb2020.R')
    #    --> saved into same .rda





### Picking up with optimally-defined PC space ===
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda",
     verbose=TRUE)
    # sce.all.n12, chosen.hvgs.all.n12, pc.choice.n12


# How many PCs is optimal?:
metadata(pc.choice.n12)$chosen
    ## 204

## Assign this chosen ( PCs) to 'PCA_opt'
reducedDim(sce.all.n12, "PCA_opt") <- reducedDim(sce.all.n12, "PCA")[ ,1:(metadata(pc.choice.n12)$chosen)]


## t-SNE
set.seed(109)
sce.all.n12 <- runTSNE(sce.all.n12, dimred="PCA_opt")


## UMAP
set.seed(109)
sce.all.n12 <- runUMAP(sce.all.n12, dimred="PCA_opt")


## Add some phenodata
sce.all.n12$region <- ss(sce.all.n12$sample,".5",1)
sce.all.n12$donor <- paste0("Br",ss(sce.all.n12$sample,"y.",2))

df.processDate <- data.frame(sample=unique(sce.all.n12$sample))
df.processDate$date <- ifelse(df.processDate$sample =="hpc.5161","08May19", "23Jul19")
df.processDate$date

sce.all.n12$processDate <- ifelse(sce.all.n12$sample=="amy.5161", "04Sep19", "25Sep19")


## Make some pd to add to colData
ref.sampleInfo <- data.frame(sampleID = unique(sce.all.n12$sample))
ref.sampleInfo$realBatch <- ifelse(ref.sampleInfo$sampleID %in% c("dlpfc.5212", "hpc.5287", "nac.5161", "nac.5212"),
                                   "R2.Jul23", "R1.May08")
ref.sampleInfo$realBatch <- ifelse(ref.sampleInfo$sampleID %in% c("dlpfc.5161", "hpc.5212", "nac.5287", "amy.5161"),
                                   "R3.Sep04", ref.sampleInfo$realBatch)
ref.sampleInfo$realBatch <- ifelse(ref.sampleInfo$sampleID %in% c("amy.5212", "sacc.5161", "sacc.5212", "nac.5182", "nac.5207"),
                                   "R4.Sep25", ref.sampleInfo$realBatch)

ref.sampleInfo$protocol <- "Frank"
ref.sampleInfo$protocol <- ifelse(ref.sampleInfo$realBatch=="R2.Jul23", "pseudoSort", ref.sampleInfo$protocol)

# Add to sce.all.n12 colData
sce.all.n12$region <- ss(sce.all.n12$sample,".5",1)
sce.all.n12$donor <- paste0("Br",substr(sce.all.n12$sample, start=nchar(sce.all.n12$sample)-3, stop=nchar(sce.all.n12$sample)))
sce.all.n12$processDate <- ref.sampleInfo$realBatch[match(sce.all.n12$sample, ref.sampleInfo$sampleID)]
sce.all.n12$protocol <- ref.sampleInfo$protocol[match(sce.all.n12$processDate, ref.sampleInfo$realBatch)]


# Save for now
save(sce.all.n12, chosen.hvgs.all.n12, pc.choice.n12, ref.sampleInfo,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda")


### Clustering: Two-step ======================================================
### Step 1: Perform graph-based clustering in this optimal PC space
#         - take k=20 NN to build graph
snn.gr <- buildSNNGraph(sce.all.n12, k=20, use.dimred="PCA_opt")
clusters.k20 <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters.k20)
    ## 

# Assign as 'prelimCluster'
sce.all.n12$prelimCluster <- factor(clusters.k20)


### Step 2: Hierarchical clustering of pseudo-bulked ("PB'd") counts with most robust normalization
#         (as determined in: 'side-Rscript_testingStep2_HC-normalizn-approaches_wAmygData_MNTJan2020.R')
#           ** That is, to pseudo-bulk (aka 'cluster-bulk') on raw counts, on all [non-zero] genes,
#              normalize with `librarySizeFactors()`, log2-transform, then perform HC'ing


## Preliminary cluster index for pseudo-bulking
# Will need to split out by single-vs-multiple-nuclei-containing prelimClusters bc
#`rowSums()` doesn't like the former
clusIndexes = splitit(sce.all.n12$prelimCluster)
clusIndexes.single <- unlist(clusIndexes[which(lengths(clusIndexes)==1)])
clusIndexes.mult <- clusIndexes[which(lengths(clusIndexes)>1)]
    

prelimCluster.PBcounts.mult <- sapply(clusIndexes.mult, function(ii){
  rowSums(assays(sce.all.n12)$counts[ ,ii])
  }
)
    ## of dim [1] 33538   242

prelimCluster.PBcounts <- cbind(prelimCluster.PBcounts.mult,
                                as.matrix(assays(sce.all.n12)$counts[ ,clusIndexes.single]))
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
clusterRefTab.all.n12 <- data.frame(origClust=order.dendrogram(dend),
                                merged=clust.treeCut[order.dendrogram(dend)])


# Assign as 'collapsedCluster'
sce.all.n12$collapsedCluster <- factor(clusterRefTab.all.n12$merged[match(sce.all.n12$prelimCluster,
                                                                          clusterRefTab.all.n12$origClust)])

# Print some visualizations:
pdf("pdfs/panBrain-n12_reducedDims-with-collapsedClusters_Feb2020.pdf")
plotReducedDim(sce.all.n12, dimred="PCA", ncomponents=5, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.all.n12, colour_by="region", point_alpha=0.5)
plotTSNE(sce.all.n12, colour_by="processDate", point_alpha=0.5)
plotTSNE(sce.all.n12, colour_by="sample", point_alpha=0.5)
plotTSNE(sce.all.n12, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.all.n12, colour_by="sum", point_alpha=0.5)
plotUMAP(sce.all.n12, colour_by="collapsedCluster", point_alpha=0.5)
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
    plotExpression(sce.all.n12, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="collapsedCluster", colour_by="collapsedCluster", point_alpha=0.5, point_size=.7,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:19], length(markers.mathys.custom[[i]])))
  )
}
dev.off()

# Any driven by low n transcripts, overall?
newClusIndex <- splitit(sce.all.n12$collapsedCluster)
sapply(newClusIndex, function(x) {quantile(sce.all.n12[,x]$sum)})
    #          1         2        3        4      5        6       7        8     9
    #0%     2919    101.00   4331.0   640.00    781    981.0  2841.0  1569.00  6934
    #25%   32795    213.25  31245.5  5109.25  29238  28993.0 26476.0 14009.25 12159
    #50%   49699  17730.50  41311.0  7567.00  38642  40330.5 35681.0 22541.50 18112
    #75%   67486  30153.50  56835.0 11068.50  50639  58606.5 44494.5 33442.75 33203
    #100% 196431 110271.00 127499.0 34702.00 121477 165583.0 77108.0 93619.00 68229
    #        10    11    12       13        14        15       16     17       18
    #0%    1139   328   126   822.00   4625.00    787.00   400.00   4507   299.00
    #25%   7511  4388  3194 24697.25  31875.00   9955.75   813.75  32680   661.75
    #50%  10415  6266  4655 33093.50  46351.00  19465.50  1063.00  43737   859.00
    #75%  13839  8531  6182 43756.00  61969.75  46144.50  2864.25  53954  1150.50
    #100% 34023 31312 25323 88046.00 118374.00 150461.00 25711.00 100139 18940.00
    #           19
    #0%    3641.00
    #25%   7541.00
    #50%  10025.00
    #75%  12912.25
    #100% 23657.00

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

sce.all.n12$cellType <- annotationTab.pan$cellType[match(sce.all.n12$collapsedCluster,
                                                         annotationTab.pan$cluster)]

## Save
save(sce.all.n12, chosen.hvgs.all.n12, pc.choice.n12, ref.sampleInfo, clusterRefTab.all.n12,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda")



### MNT 20Mar2020 === === ===
  # Add region-specific annotations

load("rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda", verbose=T)
    # sce.all.n12, chosen.hvgs.all.n12, pc.choice.n12, ref.sampleInfo, clusterRefTab.all.n12

## Amyg
load("rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    # sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo
    rm(chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy)

## DLPFC
load("rdas/regionSpecific_DLPFC-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    # sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc, ref.sampleInfo
    rm(chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc)

## HPC
load("rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    # sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo
    rm(chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc)
# Change an annotation
sce.hpc$cellType <- factor(gsub(pattern="Ambig.glial", "Tcell", sce.hpc$cellType))

## NAc
load("rdas/regionSpecific_NAc-n3_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    # sce.nac, chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac, ref.sampleInfo
    rm(chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac)

## sACC
load("rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    # sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc, ref.sampleInfo
    rm(chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc)


regions <- list(sce.amy, sce.dlpfc, sce.hpc, sce.nac, sce.sacc)

BC.sample.ii <- paste0(colnames(sce.all.n12),".",sce.all.n12$sample)

sce.all.n12$cellType.RS <- NA

for(r in regions){
  matchTo <- paste0(colnames(r),".",r$sample)
  BC.sample.sub <- BC.sample.ii %in% matchTo
  sce.all.n12$cellType.RS[BC.sample.sub] <- as.character(r$cellType[match(BC.sample.ii[BC.sample.sub], matchTo)])
}

table(sce.all.n12$cellType.RS)
sce.all.n12$cellType.RS <- factor(sce.all.n12$cellType.RS)
table(sce.all.n12$cellType.RS)


# Getting rid of sub-cell types in sACC sample and the pan-brain assignments to compare
table(ss(as.character(sce.all.n12$cellType.RS),"\\.",1) ==
        ss(as.character(sce.all.n12$cellType),"\\.",1))
    # FALSE  TRUE
    #   836 33234       - 33234/(836+ 33234) = 97.5% congruence

# Region specific
table(ss(as.character(sce.all.n12$cellType.RS),"\\.",1))
    #Ambig Astro Excit Inhib Micro Oligo   OPC Tcell
    #  414  3863  2927  2662  2987 18664  2527    26

table( ss(as.character(sce.all.n12$cellType),"\\.",1))
    #Ambig Astro Excit Inhib Micro Oligo   OPC
    #   32  3828  2848  3110  3077 18614  2561


## Pretty good - let's save
save(sce.all.n12, chosen.hvgs.all.n12, pc.choice.n12, ref.sampleInfo, clusterRefTab.all.n12,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda")




## A bit on "Ambig.hiVCAN"
CTids <- splitit(sce.all.n12$cellType)
sapply(CTids, function(x){quantile(sce.all.n12[ ,x]$sum)})
    # Ambig.hiVCAN  Astro.1  Astro.2 Excit.1  Excit.2  Excit.3 Excit.4   Excit.5
    # 0%         400.00   640.00   299.00    2919   4331.0    981.0    6934   4625.00
    # 25%        813.75  5109.25   661.75   32795  31245.5  28993.0   12159  31875.00
    # 50%       1063.00  7567.00   859.00   49699  41311.0  40330.5   18112  46351.00
    # 75%       2864.25 11068.50  1150.50   67486  56835.0  58606.5   33203  61969.75
    # 100%     25711.00 34702.00 18940.00  196431 127499.0 165583.0   68229 118374.00
    # Excit.6 Excit.7  Excit.8   Inhib.1 Inhib.2 Inhib.3  Inhib.4  Inhib.5
    # 0%      787.00    4507  3641.00    101.00     781  2841.0  1569.00   822.00
    # 25%    9955.75   32680  7541.00    213.25   29238 26476.0 14009.25 24697.25
    # 50%   19465.50   43737 10025.00  17730.50   38642 35681.0 22541.50 33093.50
    # 75%   46144.50   53954 12912.25  30153.50   50639 44494.5 33442.75 43756.00
    # 100% 150461.00  100139 23657.00 110271.00  121477 77108.0 93619.00 88046.00
    # Micro Oligo   OPC
    # 0%     126   328  1139
    # 25%   3194  4388  7511
    # 50%   4655  6266 10415
    # 75%   6182  8531 13839
    # 100% 25323 31312 34023      

sce.all.n12$cellType.RS[sce.all.n12$cellType=="Ambig.hiVCAN"]
    # [1] Inhib   Excit   OPC     OPC     OPC     OPC     OPC     OPC     OPC
    # [10] Oligo   Micro   OPC     OPC     OPC     OPC     Micro   OPC     OPC
    # [19] OPC     OPC     OPC     Micro   OPC     OPC     OPC     OPC     OPC
    # [28] OPC     Excit.2 Oligo   Excit.2 Excit.4

        # so looks like it's just mostly OPC nuclei that had less N transcripts captured


## Misc exploration ======================================================================

# Is 16 maybe pericytes? Looking at genes from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5057134/
plotExpression(sce.all.n12, exprs_values = "logcounts", features=c("VTN", "IFITM1"),
               x="collapsedCluster", colour_by="collapsedCluster", point_alpha=0.5, point_size=.7,
               add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                            geom = "crossbar", width = 0.3,
                                            colour=rep(tableau20[1:19], 2))
    ## (no expression of either of these)

## How well correlated is the expression of those Ambig.hiVCAN nuclei anyway?
 #    Let's say compared to "Astro.2", which also has a low distribution of n transcripts
 #    (collapsedCluster 18)
library(pheatmap)

table(sce.all.n12$cellType, sce.all.n12$region)

clustIdx.18 <- which(sce.all.n12$collapsedCluster==18)  # "Astro.2"
clustIdx.16 <- which(sce.all.n12$collapsedCluster==16)


quantile(cor(as.matrix(assay(sce.all.n12, "logcounts")[ ,clustIdx.18])), probs=seq(0.05,1,by=0.05))
    #         5%       10%       15%       20%       25%       30%       35%       40%
    #  0.1435597 0.1550650 0.1635447 0.1710807 0.1772498 0.1833153 0.1881462 0.1927772
    #        45%       50%       55%       60%       65%       70%       75%       80%
    #  0.1974163 0.2016515 0.2070635 0.2124137 0.2177567 0.2237208 0.2296916 0.2373880
    #        85%       90%       95%      100%
    #  0.2479811 0.2656165 0.3195865 1.0000000      - so actually pretty poor...

pheatmap(cor(as.matrix(assay(sce.all.n12, "logcounts")[ ,clustIdx.18])),
         labels_row=F, labels_col=F)


quantile(cor(as.matrix(assay(sce.all.n12, "logcounts")[ ,clustIdx.16])), probs=seq(0.05,1,by=0.05))
    #         5%       10%       15%       20%       25%       30%       35%       40%
    #  0.1139671 0.1440696 0.1624547 0.1762752 0.1880453 0.1962381 0.2028947 0.2085750
    #        45%       50%       55%       60%       65%       70%       75%       80%
    #  0.2163537 0.2201661 0.2272381 0.2393277 0.2517573 0.2621352 0.2746192 0.2888875
    #        85%       90%       95%      100%
    #  0.3063891 0.3572942 0.4383301 1.0000000      - a little better with this 'Ambig.hiVCAN'

pheatmap(cor(as.matrix(assay(sce.all.n12, "logcounts")[ ,clustIdx.inhib2])),
         labels_row=F, labels_col=F)


# Showing pairwise correlation is good with more robust clusters...
clustIdx.inhib2 <- which(sce.all.n12$cellType=="Inhib.2")
pheatmap(cor(as.matrix(assay(sce.all.n12, "logcounts")[ ,clustIdx.inhib2])),
         labels_row=F, labels_col=F)
    ## much better
quantile(cor(as.matrix(assay(sce.all.n12, "logcounts")[ ,clustIdx.inhib2])), probs=seq(0.05,1,by=0.05))
    #         5%       10%       15%       20%       25%       30%       35%       40%
    #  0.6184637 0.6865806 0.7100607 0.7239689 0.7342448 0.7425027 0.7498117 0.7564875
    #        45%       50%       55%       60%       65%       70%       75%       80%
    #  0.7629176 0.7689357 0.7750985 0.7810500 0.7871078 0.7932005 0.7995850 0.8064660
    #        85%       90%       95%      100%
    #  0.8138758 0.8225428 0.8345357 1.0000000

    # Is this because it's a bigger cluster?
    clustIdx.excit4 <- which(sce.all.n12$cellType=="Excit.4")
    pheatmap(cor(as.matrix(assay(sce.all.n12, "logcounts")[ ,clustIdx.excit4])),
         labels_row=F, labels_col=F)
    quantile(cor(as.matrix(assay(sce.all.n12, "logcounts")[ ,clustIdx.excit4])), probs=seq(0.05,1,by=0.05))
        ## less good as "Inhib.2", but still much better than collapsedCluster 16 & 18

    ## What about collapsedCluster 19?  On second glance, it expresses many markers...
     #      (notably VCAN; some PLP1; and more SLC17A6 than SLC17A7)
    clustIdx.excit8 <- which(sce.all.n12$cellType=="Excit.8")
    pheatmap(cor(as.matrix(assay(sce.all.n12, "logcounts")[ ,clustIdx.excit8])),
             labels_row=F, labels_col=F)
    quantile(cor(as.matrix(assay(sce.all.n12, "logcounts")[ ,clustIdx.excit8])), probs=seq(0.05,1,by=0.05))
        ## Overall still pretty good, but less than "Excit.4"
    
## For reference === === === === ===
table(sce.all.n12$cellType, sce.all.n12$region)
    
    #              amy dlpfc  hpc  nac sacc
    #Ambig.hiVCAN    1     3   23    0    5
    #Astro.1       842   486 1223  541  606
    #Astro.2         9    27   74    5   15
    #Excit.1        11   258    4    0  648
    #Excit.2         0   193    1    0  329
    #Excit.3       343     0  383    1    1
    #Excit.4         1     0    1   10   21
    #Excit.5         0    83    3    0  182
    #Excit.6         1     7  176    0    0
    #Excit.7         0    33    3    0   83
    #Excit.8        39     0   33    0    0
    #Inhib.1       183   228  250   71  218
    #Inhib.2        69   133   57   24  266
    #Inhib.3       113    45  110    0  115
    #Inhib.4        16   133   29    0  246
    #Inhib.5       130     1   31  638    4
    #Micro         790   272 1271  215  529
    #Oligo        3450  3228 5887 2804 3245
    #OPC           634   269  885  239  534






