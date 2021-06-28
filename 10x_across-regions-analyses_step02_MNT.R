### MNT 10x snRNA-seq workflow: step 02
###   **PAN-BRAIN analyses**
###     - (n=12) all regions from: Br5161 & Br5212 & Br5287
###     - Amyg, DLPFC, HPC, NAc, and sACC
### Initiated MNT 07Feb2020
### test.edit
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

# Print some visualizations:      (updated, below)
#pdf("pdfs/panBrain-n12_reducedDims-with-collapsedClusters_Feb2020.pdf")
plotReducedDim(sce.all.n12, dimred="PCA", ncomponents=5, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.all.n12, colour_by="region", point_alpha=0.5)
plotTSNE(sce.all.n12, colour_by="processDate", point_alpha=0.5)
plotTSNE(sce.all.n12, colour_by="sample", point_alpha=0.5)
plotTSNE(sce.all.n12, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.all.n12, colour_by="sum", point_alpha=0.5)
plotUMAP(sce.all.n12, colour_by="collapsedCluster", point_alpha=0.5)
#dev.off()



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
  'MSNs.pan' = c("PPP1R1B","BCL11B")
)

#pdf("pdfs/zold_panBrain-n12_marker-logExprs_collapsedClusters_Feb2020.pdf", height=6, width=12)
pdf("pdfs/panBrain-n12_marker-logExprs_collapsedClusters_Apr2020.pdf", height=6, width=12)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.all.n12, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:19], length(markers.mathys.custom[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(markers.mathys.custom)[i], " markers"))
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
# load("rdas/regionSpecific_DLPFC-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
#     # sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc, ref.sampleInfo
#     rm(chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc)
    #or
    load("rdas/regionSpecific_DLPFC-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
        # sce.dlpfc.st
    
## HPC
load("rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    # sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo
    rm(chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc)
# Change an annotation
sce.hpc$cellType <- factor(gsub(pattern="Ambig.glial", "Tcell", sce.hpc$cellType))

## NAc
#load("rdas/regionSpecific_NAc-n3_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
#    # sce.nac, chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac, ref.sampleInfo
#    rm(chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac)

    ## Actually use the annotations for these respective nuclei from the all-NAc analysis
    load("rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda", verbose=T)
        # sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all, ref.sampleInfo
        rm(chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all)
    
    # TEMPorarily set $cellType.final to $cellType and change casing of 'ambig.lowNtrxts'
    sce.nac.all$cellType <- sce.nac.all$cellType.final
    sce.nac.all$cellType <- ifelse(sce.nac.all$cellType=="ambig.lowNtrxts",
                                   "Ambig.lowNtrxts",
                                   as.character(sce.nac.all$cellType))
    # Then collapse D1 and D2 subclusters & otherwise-inhibitory
    sce.nac.all$cellType[grep("MSN.D1", sce.nac.all$cellType)] <- "MSN.D1"
    sce.nac.all$cellType[grep("MSN.D2", sce.nac.all$cellType)] <- "MSN.D2"
    sce.nac.all$cellType[grep("Inhib", sce.nac.all$cellType)] <- "Inhib"
    
    sce.nac.all$cellType <- factor(sce.nac.all$cellType)
    
## sACC
load("rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    # sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc, ref.sampleInfo
    rm(chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc)

    # Collapse as above
    sce.sacc$cellType <- as.character(sce.sacc$cellType)
    sce.sacc$cellType[grep("Inhib", sce.sacc$cellType)] <- "Inhib"
    sce.sacc$cellType[grep("Excit", sce.sacc$cellType)] <- "Excit"
    sce.sacc$cellType <- factor(sce.sacc$cellType)
    

## Assign to new colData column
regions <- list(sce.amy, sce.dlpfc, sce.hpc, sce.nac.all, sce.sacc)

BC.sample.ii <- paste0(colnames(sce.all.n12),".",sce.all.n12$sample)

#sce.all.n12$cellType.RS <- NA
    # or
    sce.all.n12$cellType.RS.sub <- NA

for(r in regions){
  matchTo <- paste0(colnames(r),".",r$sample)
  BC.sample.sub <- BC.sample.ii %in% matchTo
  sce.all.n12$cellType.RS[BC.sample.sub] <- as.character(r$cellType[match(BC.sample.ii[BC.sample.sub], matchTo)])
}

table(sce.all.n12$cellType.RS)
sce.all.n12$cellType.RS <- factor(sce.all.n12$cellType.RS)
table(sce.all.n12$cellType.RS)


## Added chunk 11May2020: add in NON-collapsed region-specific annotations ===
    # First, just load all objects but don't collapse, above

    # Temp - call everything $cellType.split (if not already)
    sce.nac.all$cellType.split <- sce.nac.all$cellType.final
    sce.sacc$cellType.split <- sce.sacc$cellType
    
    regions <- list(sce.amy, sce.dlpfc.st, sce.hpc, sce.nac.all, sce.sacc)
    names(regions) <- c("amy", "dlpfc", "hpc", "nac", "sacc")
    BC.sample.ii <- paste0(colnames(sce.all.n12),".",sce.all.n12$sample)
    
    # Add that annotation
    for(r in names(regions)){
      matchTo <- paste0(colnames(regions[[r]]),".",regions[[r]]$sample)
      BC.sample.sub <- BC.sample.ii %in% matchTo
      # Shorten
      regions[[r]]$cellType.split <- gsub("Excit", "Ex", regions[[r]]$cellType.split)
      regions[[r]]$cellType.split <- gsub("Inhib", "In", regions[[r]]$cellType.split)
      # Add region
      regions[[r]]$cellType.split <- paste0(regions[[r]]$cellType.split, "_", r)
      sce.all.n12$cellType.RS.sub[BC.sample.sub] <- as.character(regions[[r]]$cellType.split[match(BC.sample.ii[BC.sample.sub], matchTo)])
    }
    sce.all.n12$cellType.RS.sub <- factor(sce.all.n12$cellType.RS.sub)
    unique(sce.all.n12$cellType.RS.sub)    
        # 73 == 68 + 5 'Ambig.lowNtrxts' for each region.  good.



# Getting rid of sub-cell types in sACC sample and the pan-brain assignments to compare
table(gsub("MSN","Inhib",ss(as.character(sce.all.n12$cellType.RS),"\\.",1)) ==
        ss(as.character(sce.all.n12$cellType),"\\.",1))
    # FALSE  TRUE
    #   866 33204    - 33204/(866+ 33204) = 97.5% congruence

# Region specific
table(ss(as.character(sce.all.n12$cellType.RS),"\\.",1))
    # Ambig Astro Excit Inhib Micro   MSN Oligo   OPC Tcell
    #   445  3864  2927  2019  2956   642 18664  2527    26

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

        ## Further evidence it's best to drop 'Excit.4' from downstream analyses === ===
         # (given $cellType.RS.sub annotations)
        sce.all.n12$cellType.RS.sub[which(sce.all.n12$cellType=="Excit.4")]
            # [1] In.5_amy     In.1_hpc     MSN.D2.2_nac In.2_nac     MSN.D1.3_nac
            # [6] MSN.D2.2_nac MSN.D2.1_nac MSN.D2.1_nac MSN.D2.2_nac MSN.D1.4_nac
            # [11] MSN.D1.4_nac MSN.D1.4_nac Ex.1_sacc    Ex.1_sacc    Ex.1_sacc
            # [16] Ex.1_sacc    Ex.1_sacc    Ex.1_sacc    Ex.1_sacc    Ex.1_sacc
            # [21] Ex.3_sacc    Ex.1_sacc    Ex.3_sacc    Ex.1_sacc    In.1_sacc
            # [26] Ex.3_sacc    Ex.1_sacc    Ex.1_sacc    Ex.1_sacc    Ex.1_sacc
            # [31] Ex.4_sacc    Ex.1_sacc    Ex.1_sacc
    
    
    ## What about collapsedCluster 19?  On second glance, it expresses many markers...
     #      (notably VCAN; some PLP1; and more SLC17A6 than SLC17A7)
    clustIdx.excit8 <- which(sce.all.n12$cellType=="Excit.8")
    pheatmap(cor(as.matrix(assay(sce.all.n12, "logcounts")[ ,clustIdx.excit8])),
             labels_row=F, labels_col=F)
    quantile(cor(as.matrix(assay(sce.all.n12, "logcounts")[ ,clustIdx.excit8])), probs=seq(0.05,1,by=0.05))
        ## Overall still pretty good, but less than "Excit.4"
    
## For reference ========================================
table(sce.all.n12$cellType, sce.all.n12$region)
    
    #              amy dlpfc  hpc  nac sacc           & MNT notes === ===
    #Ambig.hiVCAN    1     3   23    0    5           < 50 nuclei - maybe drop, or ignore at least
    #Astro.1       842   486 1223  541  606
    #Astro.2         9    27   74    5   15
    #Excit.1        11   258    4    0  648           cortical
    #Excit.2         0   193    1    0  329           cortical
    #Excit.3       343     0  383    1    1           non-cort.
    #Excit.4         1     0    1   10   21           < 50 nuclei - maybe drop, or ignore at least
    #Excit.5         0    83    3    0  182           cortical
    #Excit.6         1     7  176    0    0           HPC
    #Excit.7         0    33    3    0   83           cortical
    #Excit.8        39     0   33    0    0           non-cort.
    #Inhib.1       183   228  250   71  218           all
    #Inhib.2        69   133   57   24  266           all
    #Inhib.3       113    45  110    0  115           non-NAc
    #Inhib.4        16   133   29    0  246           non-NAc
    #Inhib.5       130     1   31  638    4           NAc-heavy (but a good amount amyg; some hpc)
    #Micro         790   272 1271  215  529
    #Oligo        3450  3228 5887 2804 3245
    #OPC           634   269  885  239  534


# With new added [collapsed] region-specific annotations
table(sce.all.n12$cellType, sce.all.n12$cellType.RS)
    #              Ambig.lowNtrxts Astro Excit Inhib Micro MSN.D1 MSN.D2 Oligo   OPC Tcell
    # Ambig.hiVCAN               0     0     4     1     3      0      0     2    22     0
    # Astro.1                    0  3682     0     0     0      0      0    16     0     0
    # Astro.2                    2   128     0     0     0      0      0     0     0     0
    # Excit.1                    0     0   921     0     0      0      0     0     0     0
    # Excit.2                    0     0   523     0     0      0      0     0     0     0
    # Excit.3                    0     0   726     2     0      0      0     0     0     0
    # Excit.4                    0     0    22     2     0      4      5     0     0     0
    # Excit.5                    0     0   268     0     0      0      0     0     0     0
    # Excit.6                    0     0   184     0     0      0      0     0     0     0
    # Excit.7                    0     0   119     0     0      0      0     0     0     0
    # Excit.8                    0     0     0    72     0      0      0     0     0     0
    # Inhib.1                  368     4     0   574     2      0      0     0     2     0
    # Inhib.2                    2     0     0   547     0      0      0     0     0     0
    # Inhib.3                    0     0     0   383     0      0      0     0     0     0
    # Inhib.4                    0     0     1   423     0      0      0     0     0     0
    # Inhib.5                    0     0   157    14     0    461    172     0     0     0
    # Micro                     66     2     0     0  2948      0      0    35     0    26
    # Oligo                      1     1     2     1     3      0      0 18606     0     0
    # OPC                        6    47     0     0     0      0      0     5  2503     0



### Re-print some visualizations - MNT 10Apr2020 ====================
# With original coordinates and annotation
pdf("pdfs/panBrain-n12_reducedDims-with-collapsedClusters_Apr2020.pdf")
plotReducedDim(sce.all.n12, dimred="PCA", ncomponents=5, colour_by="cellType", point_alpha=0.5)
plotTSNE(sce.all.n12, colour_by="region", point_alpha=0.5, point_size=2.5) + ggtitle("t-SNE on opt PCs (d=204)")
plotTSNE(sce.all.n12, colour_by="processDate", point_alpha=0.5, point_size=2.5) + ggtitle("t-SNE on opt PCs (d=204)")
plotTSNE(sce.all.n12, colour_by="sample", point_alpha=0.5, point_size=2.5) + ggtitle("t-SNE on opt PCs (d=204)")
plotTSNE(sce.all.n12, colour_by="cellType", point_alpha=0.5, point_size=2.5) + ggtitle("t-SNE on opt PCs (d=204)")
# Region-specific annotation
plotTSNE(sce.all.n12, colour_by="cellType.RS", point_alpha=0.5, point_size=2.5) + ggtitle("t-SNE on opt PCs (d=204): region-specific annot")
plotTSNE(sce.all.n12, colour_by="sum", point_alpha=0.5, point_size=2.5) + ggtitle("t-SNE on opt PCs (d=204)")
plotUMAP(sce.all.n12, colour_by="cellType", point_alpha=0.5, point_size=2.5) + ggtitle("UMAP on opt PCs (d=204)")
dev.off()

# How many PCs?
head(attr(reducedDim(sce.all.n12, "PCA"), "percentVar"), n=75)
    # [1] 19.93610257  8.56551975  5.48499030  2.69104457  2.29910336  0.86823624  0.71376460
    # [8]  0.48377800  0.32700635  0.31751426  0.29901816  0.25911566  0.22686401  0.21400937
    # [15]  0.20303113  0.18768708  0.18136268  0.17008113  0.15380352  0.14366776  0.13666573
    # [22]  0.12527399  0.11727016  0.11141297  0.10514416  0.10317094  0.09743560  0.09296292
    # [29]  0.08796466  0.08415014  0.08363044  0.08132186  0.08054339  0.07687419  0.07353210
    # [36]  0.07313253  0.07035361  0.06841083  0.06810964  0.06737097  0.06464408  0.06338695
    # [43]  0.06261076  0.06170788  0.06062411  0.05978423  0.05923917  0.05845309  0.05738356
    # [50]  0.05623525  0.05588523  0.05546358  0.05420676  0.05329303  0.05254012  0.05199590
    # [57]  0.05148879  0.05073573  0.05023337  0.04904533  0.04896795  0.04883193  0.04804923
    # [64]  0.04767100  0.04691248  0.04680456  0.04633970  0.04582236  0.04560561  0.04480624
    # [71]  0.04469499  0.04421348  0.04352501  0.04339655  0.04307707

# 0.05% var or greater
reducedDim(sce.all.n12, "PCA_59") <- reducedDim(sce.all.n12, "PCA")[ ,c(1:59)]
# 0.1% var or greater
reducedDim(sce.all.n12, "PCA_26") <- reducedDim(sce.all.n12, "PCA")[ ,c(1:26)]

# First remove this reducedDim bc this has caused trouble previously
reducedDim(sce.all.n12, "TSNE") <- NULL

## 59 PCs tsNE === (this one looks better than 26 PCs actually)
set.seed(109)
sce.all.tsne.59pcs <- runTSNE(sce.all.n12, dimred="PCA_59")

save(sce.all.tsne.59pcs, file="rdas/ztemp_panBrain-n12_SCE-with-tSNEon59PCs_MNT.rda")
rm(sce.all.tsne.59pcs)


# MNT 16Apr: Deciding to remove the clusters won't focus on for plotting:
    #'Ambig.hiVCAN' & 'Excit.4' & those .RS-annot'd as 'Ambig.lowNtrxts'
sce.all.tsne.59pcs <- sce.all.tsne.59pcs[ ,sce.all.tsne.59pcs$cellType.RS != "Ambig.lowNtrxts"] # 445
sce.all.tsne.59pcs$cellType.RS <- droplevels(sce.all.tsne.59pcs$cellType.RS)

sce.all.tsne.59pcs <- sce.all.tsne.59pcs[ ,sce.all.tsne.59pcs$cellType != "Ambig.hiVCAN"] # 32 nuclei
sce.all.tsne.59pcs <- sce.all.tsne.59pcs[ ,sce.all.tsne.59pcs$cellType != "Excit.4"]  # 33 nuclei
sce.all.tsne.59pcs$cellType <- droplevels(sce.all.tsne.59pcs$cellType)

# Add broad cell type taken from pan-brain annotation
sce.all.tsne.59pcs$cellType.broad <- ss(as.character(sce.all.tsne.59pcs$cellType), "\\.", 1)

#pdf("pdfs/exploration/ztemp_panBrain-n12_TSNEon59PCs_MNT.pdf")
pdf("pdfs/pubFigures/panBrain-n12_tSNEon59PCs_3x3PCA_smallClustersDropped_MNTApr2020.pdf", width=9)
plotTSNE(sce.all.tsne.59pcs, colour_by="cellType", point_alpha=0.5, point_size=4.0,
         text_by="cellType.broad", text_size=8, theme_size=18) +
  ggtitle("t-SNE on top 59 PCs (pan-brain annot.)") + theme(plot.title = element_text(size=19))

plotTSNE(sce.all.tsne.59pcs, colour_by="cellType.broad", point_alpha=0.5, point_size=4.0,
         text_by="cellType.broad", text_size=7, theme_size=22) +
  ggtitle("t-SNE on top 59 PCs (broad pan-brain annot.)") + theme(plot.title = element_text(size=18))

plotTSNE(sce.all.tsne.59pcs, colour_by="cellType.RS", point_alpha=0.5, point_size=4.0,
         text_by="cellType.RS", text_size=5.5, theme_size=17) +
  ggtitle("t-SNE on top 59 PCs (region-specific annot.)") + theme(plot.title = element_text(size=18))

# Top 3 PCs
plotReducedDim(sce.all.tsne.59pcs, dimred="PCA", ncomponents=3, colour_by="cellType.broad",
               point_alpha=0.5, theme_size=15, add_legend=FALSE) + 
  ggtitle("Top three PCs (broad pan-brain annot.)") + theme(plot.title = element_text(size=18))
dev.off()


    ### MNT update 31Aug2020 === ===
      # Facet some different iterations of this 'best' tSNE: maybe by region
    
    pdf("pdfs/pubFigures/panBrain-n12_tSNEon59PCs_faceted_Aug2020.pdf", width=9)
    plotTSNE(sce.all.tsne.59pcs, colour_by="region", point_alpha=0.5, point_size=4.0, theme_size=22) +
      facet_wrap(~ sce.all.tsne.59pcs$region)
      ggtitle("t-SNE on top 59 PCs (broad pan-brain annot.)") + theme(plot.title = element_text(size=18))
    dev.off()
    
    ## More manually to have shadow of those for each region ======
    custom.cols <- c("DLPFC"=tableau20[1],
                  "sACC"=tableau20[3],
                  "HPC"=tableau20[5],
                  "AMY"=tableau20[7],
                  "NAc"=tableau20[9])
    
    sce.temp <- sce.all.tsne.59pcs
    
    ## DLPFC
    # Reorder to plot the region nuclei last
    sce.temp <- cbind(sce.temp[ ,sce.temp$region!="dlpfc"],
                      sce.temp[ ,sce.temp$region=="dlpfc"])
    # Grey out other regions
    sce.temp$reg.temp <- ifelse(sce.temp$region=="dlpfc", "DLPFC", NA)
    
    p.dlpfc <- plotTSNE(sce.temp, colour_by="reg.temp", point_alpha=0.9, point_size=3.5, theme_size=15,
             add_legend=FALSE) +
      scale_fill_manual(values=custom.cols["DLPFC"]) + ggtitle("DLPFC") +
      theme(plot.title = element_text(size=30),
            axis.title = element_text(size=0),
            axis.text = element_text(size=20))
    
    ## HPC
    # Reorder to plot the region nuclei last
    sce.temp <- cbind(sce.temp[ ,sce.temp$region!="hpc"],
                      sce.temp[ ,sce.temp$region=="hpc"])
    # Grey out other regions
    sce.temp$reg.temp <- ifelse(sce.temp$region=="hpc", "HPC", NA)
    
    p.hpc <- plotTSNE(sce.temp, colour_by="reg.temp", point_alpha=0.9, point_size=3.5, theme_size=15,
                        add_legend=FALSE) +
      scale_fill_manual(values=custom.cols["HPC"]) + ggtitle("HIPPO") +
      theme(plot.title = element_text(size=30),
            axis.title = element_text(size=0),
            axis.text = element_text(size=20))
    
    ## sACC
    # Reorder to plot the region nuclei last
    sce.temp <- cbind(sce.temp[ ,sce.temp$region!="sacc"],
                      sce.temp[ ,sce.temp$region=="sacc"])
    # Grey out other regions
    sce.temp$reg.temp <- ifelse(sce.temp$region=="sacc", "sACC", NA)
    
    p.sacc <- plotTSNE(sce.temp, colour_by="reg.temp", point_alpha=0.9, point_size=3.5, theme_size=15,
                      add_legend=FALSE) +
      scale_fill_manual(values=custom.cols["sACC"]) + ggtitle("sACC") +
      theme(plot.title = element_text(size=30),
            axis.title = element_text(size=0),
            axis.text = element_text(size=20))
    
    ## AMY
    # Reorder to plot the region nuclei last
    sce.temp <- cbind(sce.temp[ ,sce.temp$region!="amy"],
                      sce.temp[ ,sce.temp$region=="amy"])
    # Grey out other regions
    sce.temp$reg.temp <- ifelse(sce.temp$region=="amy", "AMY", NA)
    
    p.amy <- plotTSNE(sce.temp, colour_by="reg.temp", point_alpha=0.9, point_size=3.5, theme_size=15,
                       add_legend=FALSE) +
      scale_fill_manual(values=custom.cols["AMY"]) + ggtitle("AMY") +
      theme(plot.title = element_text(size=30),
            axis.title = element_text(size=0),
            axis.text = element_text(size=20))
    
    ## NAc
    # Reorder to plot the region nuclei last
    sce.temp <- cbind(sce.temp[ ,sce.temp$region!="nac"],
                      sce.temp[ ,sce.temp$region=="nac"])
    # Grey out other regions
    sce.temp$reg.temp <- ifelse(sce.temp$region=="nac", "NAc", NA)
    
    p.nac <- plotTSNE(sce.temp, colour_by="reg.temp", point_alpha=0.9, point_size=3.5, theme_size=15,
                      add_legend=FALSE) +
      scale_fill_manual(values=custom.cols["NAc"]) + ggtitle("NAc") +
      theme(plot.title = element_text(size=30),
            axis.title = element_text(size=0),
            axis.text = element_text(size=20))
    
            ## end region-colored t-SNEs ========
    
    
    ## All nuclei (pan-brain annotation, from above) ===
    p.full <- plotTSNE(sce.all.tsne.59pcs, colour_by="cellType", point_alpha=0.5, point_size=4.0,
             text_by="cellType.broad", text_size=8, theme_size=24) +
      ggtitle("t-SNE on top 59 PCs (pan-brain annot.)") + theme(plot.title = element_text(size=28))
    
    lay <- rbind(c(1,1,2),
                 c(1,1,3),
                 c(6,5,4))
    
    pdf("pdfs/pubFigures/panBrain-n12_tSNEon59PCs_faceted_v2_Aug2020.pdf", width=13.5, height=12.5)
    grid.arrange(grobs=list(p.full,
                         p.nac,
                         p.amy,
                         p.hpc,
                         p.dlpfc,
                         p.sacc),
                 layout_matrix=lay)
    dev.off()


        
## Print broad marker heatmap of pan-brain-defined clusters === === ===
load("rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda", verbose=T)
    #sce.all.n12, chosen.hvgs.all.n12, pc.choice.n12, ref.sampleInfo, clusterRefTab.all.n12

# As decided for marker detection, remove the clusters that won't focus on:
#     'Ambig.hiVCAN' & 'Excit.4' & those .RS-annot'd as 'Ambig.lowNtrxts'
sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType.RS != "Ambig.lowNtrxts"] # 445
sce.all.n12$cellType.RS <- droplevels(sce.all.n12$cellType.RS)

sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType != "Ambig.hiVCAN"] # 32 nuclei
sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType != "Excit.4"]  # 33 nuclei
sce.all.n12$cellType <- droplevels(sce.all.n12$cellType)

cell.idx <- splitit(sce.all.n12$cellType)
dat <- as.matrix(assay(sce.all.n12, "logcounts"))


genes <- c('SNAP25','SLC17A6','SLC17A7','GAD1','GAD2','AQP4','GFAP','C3','CD74','MBP','PDGFRA','VCAN','CLDN5','FLT1','SKAP1','TRAC')
current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMeans(dat[genes, ii])))

current_dat <- current_dat[ ,c(3:14, 1,2, 15:17)]

pdf("pdfs/pubFigures/heatmap-geneExprs_panBrain-annot_mean-broadMarkers_MNT.pdf")
pheatmap(current_dat, cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(0.02, 4, length.out = 101),
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "OrRd"))(100),
         fontsize_row = 17.5, fontsize_col = 17.5,
         main="Broad cell type marker expression (pan-brain annot.)", fontsize=11)
dev.off()




## 26 PCs tSNE (ignore) ================
set.seed(109)
sce.all.tsne.26pcs <- runTSNE(sce.all.n12, dimred="PCA_26")

# Print these tests
pdf("pdfs/exploration/ztemp_panBrain-n12_TSNEon26PCs_MNT.pdf")
plotTSNE(sce.all.tsne.26pcs, colour_by="cellType", point_alpha=0.5, point_size=2.5) + ggtitle("t-SNE on top 26 PCs (pan-brain annot)")
plotTSNE(sce.all.tsne.26pcs, colour_by="cellType.RS", point_alpha=0.5, point_size=2.5) + ggtitle("t-SNE on top 26 PCs (region-specific annot)")
dev.off()

save(sce.all.tsne.26pcs, file="rdas/ztemp_panBrain-n12_SCE-with-tSNEon26PCs_MNT.rda")
rm(sce.all.tsne.26pcs)






