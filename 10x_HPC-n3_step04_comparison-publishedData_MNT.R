### MNT 10x snRNA-seq workflow: step 04
###   **Region-specific analyses**
###     - (3x) HPC samples from: Br5161 & Br5212 & Br5287
###     - Setup and comparison to Habib, et al (DroNc-seq paper)
### test.edit MNT 21Jul2020 new MacBook Pro
#####################################################################

library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(batchelor)
library(DropletUtils)
library(jaffelab)
library(readr)
library(readxl)
library(lattice)
library(RColorBrewer)
library(pheatmap)

### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")
# ===


## Taken from AnJa's 'process_habib_counts.R' in the 'path':
# First
path <- '/dcl01/ajaffe/data/lab/singleCell/habib/'
    ## read counts
    geneCounts = read_delim(paste0(path,"GTEx_droncseq_hip_pcf/GTEx_droncseq_hip_pcf.umi_counts.txt.gz"),
                            delim = "\t")
    geneCounts = as.data.frame(geneCounts)
    rownames(geneCounts) = geneCounts$X1
    geneCounts$X1= NULL
    geneCounts = as.matrix(geneCounts)
    
    ## more info from supp table
    cellPheno = read_excel(paste0(path,"nmeth.4407-S10.xlsx"), skip=20)
    colnames(cellPheno) = c("CellID", "NumGenes","NumTx", "ClusterID", "ClusterName")
    cellPheno = cellPheno[match(colnames(geneCounts), cellPheno$CellID),]
    cellPheno = as.data.frame(cellPheno)
    cellPheno$SampleID = ss(cellPheno$CellID, "_")
    
    ## tsneInfo
    tsneInfo = read_delim(paste0(path,"GTEx_droncseq_hip_pcf/GTEx_droncseq_hip_pcf.tsne.txt.gz"),
                          delim = "\t")
    colnames(tsneInfo)[1] = "CellID"
    identical(tsneInfo$CellID, cellPheno$CellID) #  TRUE
    # cellPheno = cbind(cellPheno, tsneInfo[,2:3])        ## Will import into reducedDims
tsneInfo <- tsneInfo[ ,-1]
    
    ## read in sample info
    subjPheno = read_excel(paste0(path,"Human_Samples_To_Tissue.xlsx"))
    colnames(subjPheno) = c("SampleID", "TissueID")
    subjPheno = as.data.frame(subjPheno)
    
    subjPheno$SampleID = gsub("_", "-", subjPheno$SampleID)
    subjPheno$SampleID[subjPheno$SampleID == "hHP3"] = "hHP3b"
    subjPheno$TissueID[subjPheno$SampleID == "PFC-CD"] = "SM-4W9H6"
    
    ## check
    all(subjPheno$SampleID %in% cellPheno$SampleID)
    cellPheno$TissueID = subjPheno$TissueID[match(cellPheno$SampleID, subjPheno$SampleID)]
    length(unique(subjPheno$TissueID)) # 7
    
    ## and donor data
    donorPheno = read_excel(paste0(path,"nmeth.4407-S9.xlsx"),skip=2)
    donorPheno = donorPheno[2:8, ]
    colnames(donorPheno) = c("TissueID", "DonorID", "Material",
                             "RIN_Bulk", "Age", "Gender", "Race", "TissueSite", "PMI","Hardy")
    donorPheno$Region = ifelse(	donorPheno$TissueSite == "Brain - Hippocampus", "HIPPO", "DLPFC")
    donorPheno = as.data.frame(donorPheno)
    
    mm = match(cellPheno$TissueID, donorPheno$TissueID)
    cellPheno$Region = donorPheno$Region[mm]
    cellPheno$Age = donorPheno$Age[mm]
    cellPheno$RIN_Bulk = donorPheno$RIN_Bulk[mm]
    cellPheno$PMI = donorPheno$PMI[mm]
    cellPheno$DonorID = donorPheno$DonorID[mm]
    rownames(cellPheno) = cellPheno$CellID

    
## Now turn into SCE
sce.habib <- SingleCellExperiment(assays = list(counts=geneCounts),
                                  colData = cellPheno)
reducedDim(sce.habib, "TSNE.given") <- as.matrix(tsneInfo)

LSFvec <- librarySizeFactors(sce.habib)
sce.habib <- logNormCounts(sce.habib, size_factors=LSFvec)

# Plot some fun genes
pdf("pdfs/exploration/Habib_DroNc-seq/HIPPO-broadMarkers-by-annotated-cellType_MNT.pdf", height=6, width=8)
plotExpression(sce.habib, exprs_values = "logcounts", features=c("SNAP25","SYT1","MBP","VCAN",
                                                                  "CSF1R","AQP4","VTN","CLDN5"),
               x="ClusterName", colour_by="ClusterName", point_alpha=0.5, point_size=.7, ncol=4,
               add_legend=F, theme_size=8) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                          geom = "crossbar", width = 0.3,
                                                          colour=rep(tableau20[1:15], 8)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


# Save this in project dir
save(sce.habib, file="rdas/zSCE_habib-dlpfc-hippo_MNT.rda")

# Plot the provided tSNE with these colors - have to trick 
reducedDim(sce.habib, "TSNE") <- reducedDim(sce.habib, "TSNE.given")

pdf("pdfs/exploration/Habib_DroNc-seq/HIPPO-given-tSNE-by-annotated-cellType_MNT.pdf")
plotTSNE(sce.habib, colour_by="Region", point_alpha=0.5, point_size=2.5)
plotTSNE(sce.habib, colour_by="ClusterName", point_alpha=0.5, point_size=2.5)
dev.off()


# Distribution of annotated clusters across regions
table(sce.habib$ClusterName, sce.habib$Region)
    #              DLPFC HIPPO
    # ASC1           534   670
    # ASC2           523   182
    # END            116   138
    # exCA1           71   350
    # exCA3           82   663
    # exDG          1037   416    *interesting...
    # exPFC1        2964   140
    # exPFC2         258    39
    # GABA1          704   188
    # GABA2          620   203
    # MG             118   271
    # NSC            131    70
    # ODC1          1268  1697
    # OPC            380   304
    # Unclassified   230   596


## Also btw:
table(sce.habib$SampleID, sce.habib$DonorID)
    ##                1    2    3    4    5
    # hCc           780    0    0    0    0
    # hCd           812    0    0    0    0
    # hCe           733    0    0    0    0
    # hCf           779    0    0    0    0
    # hHP1          486    0    0    0    0
    # hHP2          320    0    0    0    0
    # hHP2a           0  144    0    0    0
    # hHP2b           0  166    0    0    0
    # hHP2c           0  115    0    0    0
    # hHP3b           0   47    0    0    0
    # HP2-A           0    0  717    0    0
    # HP2-B           0    0  674    0    0
    # HP3-A           0    0 1012    0    0
    # HP3-B           0    0    0 1014    0
    # HP3-B-united    0    0    0 1232    0
    # humanPFCa     752    0    0    0    0
    # humanPFCb     659    0    0    0    0
    # PFC-CD          0  918    0    0    0
    # PFC2-A1         0    0    0    0 1092
    # PFC2-A2         0    0    0    0  953
    # PFC2-A3         0    0    0    0 1085
    # PFC2-A5         0    0    0    0  473

        # Looks like 'DonorID' ~= 'TissueID'
        table(sce.habib$TissueID, sce.habib$DonorID)
        #             1    2    3    4    5
        # SM-4RGJU  806    0    0    0    0
        # SM-4RGKC 4515    0    0    0    0
        # SM-4W9GN    0  472    0    0    0
        # SM-4W9H6    0  918    0    0    0
        # SM-DG7EP    0    0 2403    0    0
        # SM-DG7EQ    0    0    0 2246    0
        # SM-DIQCU    0    0    0    0 3603


### -> Split DLPFC & HPC nuclei and analyze, separately, as if LIBD dataset ============
load("rdas/zSCE_habib-dlpfc-hippo_MNT.rda", verbose=T)
    # sce.habib

table(sce.habib$Region)
    #DLPFC HIPPO
    # 9036  5927

sce.habib.dlpfc <- sce.habib[ ,sce.habib$Region == "DLPFC"]
sce.habib.hpc <- sce.habib[ ,sce.habib$Region == "HIPPO"]

# Re-calculate logcounts
LSF.dlpfc <- librarySizeFactors(sce.habib.dlpfc)
assay(sce.habib.dlpfc, "logcounts") <- NULL
sce.habib.dlpfc <- logNormCounts(sce.habib.dlpfc, size_factors=LSF.dlpfc)

LSF.hpc <- librarySizeFactors(sce.habib.hpc)
assay(sce.habib.hpc, "logcounts") <- NULL
sce.habib.hpc <- logNormCounts(sce.habib.hpc, size_factors=LSF.hpc)

## Save these
save(sce.habib, sce.habib.dlpfc, sce.habib.hpc, file="rdas/zSCE_habib-dlpfc-hippo_MNT.rda")

rm(sce.habib)



### Work with HPC for now (have other, better DLPFC datasets) ==============================

# Determine HVGs
geneVar.hpc <- modelGeneVar(sce.habib.hpc)
chosen.hvgs.hpc <- geneVar.hpc$bio > 0
sum(chosen.hvgs.hpc)

# Run PCA, taking top 100 (instead of default 50 PCs)
set.seed(109)
sce.habib.hpc <- runPCA(sce.habib.hpc, subset_row=chosen.hvgs.hpc, ncomponents=200,
                  BSPARAM=BiocSingular::RandomParam())

sum(attr(reducedDim(sce.habib.hpc, "PCA"), "percentVar"))
    # [1] 31.26207 with 200 PCs - let's just go with this;; (22.37983 with 100 PCs)


## t-SNE
set.seed(109)
sce.habib.hpc <- runTSNE(sce.habib.hpc, dimred="PCA")
  
    # How does this look?
    plotTSNE(sce.habib.hpc, colour_by="TissueID") # looks like fireworks...
    plotReducedDim(sce.habib.hpc, dimred="PCA", ncomponents=3, colour_by="TissueID", point_alpha=0.5)
        # Some veeeeery strong donor effects, or at least bias...


### Clustering: Two-step === === === ===
### Step 1: Perform graph-based clustering in this optimal PC space
#         - take k=20 NN to build graph
snn.gr <- buildSNNGraph(sce.habib.hpc, k=20, use.dimred="PCA")
clusters.k20 <- igraph::cluster_walktrap(snn.gr)$membership
table(clusters.k20)
    ##   1    2    3    4    5    6    7    8    9   10
     # 837 2368   77   37   42    6  608 1754   66  132

# Assign as 'prelimCluster'
sce.habib.hpc$prelimCluster <- factor(clusters.k20)

# Is sample driving this 'high-res' clustering at this level?
table(sce.habib.hpc$prelimCluster, sce.habib.hpc$TissueID)  
    # not bad

# Visualize this
plotTSNE(sce.habib.hpc, colour_by="prelimCluster")  # Actually pretty good.  Let's keep this


    ### Step 2: Hierarchical clustering of pseudo-bulked ("PB'd") counts with most robust normalization
    #         (as determined in: 'side-Rscript_testingStep2_HC-normalizn-approaches_wAmygData_MNTJan2020.R')
        ## *** UNNEEDED - proceed with 'prelimCluster'

# Save this for now
save(sce.habib.hpc, file="rdas/zSCE_habib_HPC-only_MNT.rda")


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

pdf("pdfs/exploration/Habib_DroNc-seq/HIPPO-MNTclusters_marker-logExprs_Jul2020.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.habib.hpc, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
#                   x="prelimCluster", colour_by="prelimCluster", point_alpha=0.5, point_size=.7,
                   x="cellType.mnt", colour_by="cellType.mnt", point_alpha=0.5, point_size=.7,
                   add_legend=F) + stat_summary(fun = median, fun.min = median, fun.max = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium, length(markers.mathys.custom[[i]]))) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
        ggtitle(label=paste0(names(markers.mathys.custom)[i], " markers"))
  )
}
dev.off()


## Add quick annotations - some of them don't seem to be any broad cell type...
annotationTab.hpc <- data.frame(cluster=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                                cellType=c("Glia.mixed", "Neuron.mixed", "Endo.1", "Unknown.small", "Endo.2",
                                           "Endo.small", "Astro.1", "Oligo", "Astro.2", "Unknown")
                                )

sce.habib.hpc$cellType.mnt <- annotationTab.hpc$cellType[match(sce.habib.hpc$prelimCluster,
                                                     annotationTab.hpc$cluster)]
sce.habib.hpc$cellType.mnt <- factor(sce.habib.hpc$cellType.mnt)

# Save this
save(sce.habib.hpc, chosen.hvgs.hpc, annotationTab.hpc, file="rdas/zSCE_habib_HPC-only_MNT.rda")

sce.habib.hpc$ClusterName.habib <- sce.habib.hpc$ClusterName

# Plot some reducedDims
pdf("pdfs/exploration/Habib_DroNc-seq/HIPPO-reducedDims-with-MNTclusters_Jul2020.pdf")
plotReducedDim(sce.habib.hpc, dimred="PCA", ncomponents=4, colour_by="TissueID", point_alpha=0.5)
plotReducedDim(sce.habib.hpc, dimred="PCA", ncomponents=4, colour_by="cellType.mnt", point_alpha=0.5)
plotTSNE(sce.habib.hpc, colour_by="TissueID", point_size=3.5, point_alpha=0.5)
plotTSNE(sce.habib.hpc, colour_by="cellType.mnt", point_size=3.5, point_alpha=0.5)
plotTSNE(sce.habib.hpc, colour_by="ClusterName.habib", point_size=3.5, point_alpha=0.5)
plotTSNE(sce.habib.hpc, colour_by="NumTx", point_size=3.5, point_alpha=0.5)
dev.off()


# Btw
    table(sce.habib.hpc$cellType.mnt, sce.habib.hpc$TissueID)
        #                SM-4RGJU SM-4W9GN SM-DG7EP SM-DG7EQ
        # Astro.1             62       40      244      262
        # Astro.2              0        0       53       13
        # Endo.1              12        3       17       45
        # Endo.2               9        3       10       20
        # Endo.small           2        1        3        0
        # Glia.mixed         125       45      351      316
        # Neuron.mixed       323      235     1047      763
        # Oligo              190       59      678      827
        # Unknown             62       70        0        0
        # Unknown.small       21       16        0        0
    
    round(prop.table(table(sce.habib.hpc$cellType.mnt)), 3)
        #    Astro.1       Astro.2        Endo.1        Endo.2    Endo.small
        #      0.103         0.011         0.013         0.007         0.001
        # Glia.mixed  Neuron.mixed         Oligo       Unknown Unknown.small
        #      0.141         0.400         0.296         0.022         0.006



## Define markers for these clusters and t's to compare LIBD HPC to ========================
    
    ## Btw for HPC, $DonorID == $TissueID

### Run ANOVA real quick, just on HVGs ===
library(edgeR)
library(doMC)
registerDoMC(cores=4)

mat = assays(sce.habib.hpc)$logcounts[chosen.hvgs.hpc, ]

## do regression
varCompAnalysis.splitDonor = foreach(i = 1:nrow(mat)) %dopar% {
  if(i %% 100 == 0) cat(".")
  #fit = lm(as.numeric(mat[i,]) ~ cellType.mnt + TissueID + PMI,
  fit = lm(as.numeric(mat[i,]) ~ cellType.mnt + SampleID + PMI,    # or with more split 'SampleID'
           data=colData(sce.habib.hpc))
  full = anova(fit)
  fullSS = full$"Sum Sq"
  signif(cbind(full, PctExp = fullSS/sum(fullSS)*100), 3)
}
    ## Looks like couldn't estimate effect of PMI, even with modeling on just 'TissueID' 

names(varCompAnalysis) = rownames(mat)
names(varCompAnalysis.splitDonor) = rownames(mat)

## make boxplot
varExpl = t(sapply(varCompAnalysis, function(x) x[,"PctExp"]))
colnames(varExpl) = rownames(varCompAnalysis[[1]])

varExpl.splitDonor = t(sapply(varCompAnalysis.splitDonor, function(x) x[,"PctExp"]))
colnames(varExpl.splitDonor) = rownames(varCompAnalysis.splitDonor[[1]])


pdf("pdfs/exploration/HIPPO-anova_MNTclusters_Jul2020.pdf")
boxplot(varExpl, main="ANOVA on human HIPPO DroNc-seq \n (Habib, et al. Nat. Methods 2017)",
        ylab="Percent Var explained (%))")
boxplot(varExpl.splitDonor, main="ANOVA on human HIPPO DroNc-seq \n (Habib, et al. Nat. Methods 2017)",
        ylab="Percent Var explained (%))")
dev.off()
    # SampleID only slightly better...

save(varCompAnalysis, varCompAnalysis.splitDonor, file="./zs-habib-HIPPO_anova_output_MNT.rda")

apply(varExpl, 2, quantile)
    #      cellType.mnt TissueID Residuals
    # 0%        0.00922 9.68e-06      16.0
    # 25%       0.10300 2.74e-02      99.2
    # 50%       0.24800 5.21e-02      99.7
    # 75%       0.75100 9.71e-02      99.8
    # 100%     84.00000 8.84e+00     100.0

apply(varExpl.splitDonor, 2, quantile)
    #      cellType.mnt SampleID Residuals        -> use this model for `findMarkers()` (no 'PMI')
    # 0%        0.00922  0.00614      15.9
    # 25%       0.10300  0.09620      99.0
    # 50%       0.24800  0.14600      99.6
    # 75%       0.75100  0.22200      99.8
    # 100%     84.00000  9.47000      99.9



### Modeling for subcluster-specific genes - cluster-vs-all test, only === === ===
  # (only doing this iteration because this yields t-statistics to compare to human)

table(rowSums(assay(sce.habib.hpc, "counts"))==0)  # 3962 TRUE

# Drop genes with all 0's
sce.habib.hpc <- sce.habib.hpc[!rowSums(assay(sce.habib.hpc, "counts"))==0, ]

# Model unwanted effects
mod <- with(colData(sce.habib.hpc), model.matrix(~ SampleID))
mod <- mod[ ,-1]

markers.habibHPC.t.1vAll <- list()
for(i in levels(sce.habib.hpc$cellType.mnt)){
  # Make temporary contrast
  sce.habib.hpc$contrast <- ifelse(sce.habib.hpc$cellType.mnt==i, 1, 0)
  # Test cluster vs. all
  markers.habibHPC.t.1vAll[[i]] <- findMarkers(sce.habib.hpc, groups=sce.habib.hpc$contrast,
                                            assay.type="logcounts", design=mod, test="t",
                                            direction="up", pval.type="all", full.stats=T)
}

## Then, temp set of stats to get the standardized logFC
    temp.1vAll <- list()
    for(i in levels(sce.habib.hpc$cellType.mnt)){
      # Make temporary contrast
      sce.habib.hpc$contrast <- ifelse(sce.habib.hpc$cellType.mnt==i, 1, 0)
      # Test cluster vs. all
      temp.1vAll[[i]] <- findMarkers(sce.habib.hpc, groups=sce.habib.hpc$contrast,
                                     assay.type="logcounts", design=mod, test="t",
                                     std.lfc=TRUE,
                                     direction="up", pval.type="all", full.stats=T)
    }

## For some reason all the results are in the second List entry (first is always empty)

# Replace that empty slot with the entry with the actul stats
markers.habibHPC.t.1vAll <- lapply(markers.habibHPC.t.1vAll, function(x){ x[[2]] })
# Same for that with std.lfc
temp.1vAll <- lapply(temp.1vAll, function(x){ x[[2]] })

# Now just pull from the 'stats.0' DataFrame column
markers.habibHPC.t.1vAll <- lapply(markers.habibHPC.t.1vAll, function(x){ x$stats.0 })
temp.1vAll <- lapply(temp.1vAll, function(x){ x$stats.0 })

# Re-name std.lfc column and add to the first result
for(i in names(temp.1vAll)){
  colnames(temp.1vAll[[i]])[1] <- "std.logFC"
  markers.habibHPC.t.1vAll[[i]] <- cbind(markers.habibHPC.t.1vAll[[i]], temp.1vAll[[i]]$std.logFC)
  # Oh the colname is kept weird
  colnames(markers.habibHPC.t.1vAll[[i]])[4] <- "std.logFC"
  # Then re-organize
  markers.habibHPC.t.1vAll[[i]] <- markers.habibHPC.t.1vAll[[i]][ ,c("logFC","std.logFC","log.p.value","log.FDR")]
}


## Let's save this along with the previous pairwise results
save(markers.habibHPC.t.1vAll, file="rdas/zs-habib_markers-stats_MNTclusters_findMarkers-SN-LEVEL_Jul2020.rda")



sapply(markers.habibHPC.t.1vAll, function(x){table(x$log.FDR < log10(0.000001))})
    #      Astro.1 Astro.2 Endo.1 Endo.2 Endo.small Glia.mixed Neuron.mixed Oligo
    # FALSE   26844   26957  26485  27242      27699      27325        23849 27006
    # TRUE     1305    1192   1664    907        450        824         4300  1143
    #       Unknown Unknown.small
    # FALSE   26413         27070
    # TRUE     1736          1079


## Print those top 40 into table so can compare to human markers
markerList.t.habibHPC.1vAll <- lapply(markers.habibHPC.t.1vAll, function(x){
  rownames(x)[x[ ,"log.FDR"] < log10(0.000001)]
  }
)

top40genes <- sapply(markerList.t.habibHPC.1vAll, function(x) head(x, n=40))

write.csv(top40genes,"tables/forRef_top40genesLists_habib-HIPPO_MNTclusters_SN-LEVEL-tests_Jul2020.csv")


## Print top 20 marker genes to better help ID these ============
# Print these to pngs
markerList.t.habibHPC.1vAll <- lapply(markers.habibHPC.t.1vAll, function(x){
  rownames(x)[x[ ,"log.FDR"] < log10(0.000001)]
  }
)
top20genes <- lapply(markerList.t.habibHPC.1vAll, function(x) head(x, n=20))

# Print
for(i in names(top20genes)){
  png(paste0("pdfs/exploration/Habib_DroNc-seq/HIPPO-MNTcluster-markers_t-sn-level_1vALL_top20markers-",
             i,"_logExprs_Jul2020.png"), height=900, width=1200)
  print(
    plotExpression(sce.habib.hpc, exprs_values = "logcounts", features=top20genes[[i]],
                   x="cellType.mnt", colour_by="cellType.mnt", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium, 20)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " (MNT cluster ID) top 20 markers"))
  )
  dev.off()
}


cellType.idx <- splitit(sce.habib.hpc$cellType.mnt)
sapply(cellType.idx, function(x){quantile(sce.habib.hpc[ ,x]$NumTx)})
    #      Astro.1 Astro.2 Endo.1  Endo.2 Endo.small Glia.mixed Neuron.mixed  Oligo
    # 0%     236.0  249.00    520  269.00     324.00        225       226.00  233.0
    # 25%    343.0  388.00    784  355.25     355.50        299       654.00  308.0
    # 50%    443.5  617.00   1105  569.50     402.50        369       887.50  380.5
    # 75%    646.0  948.75   1604  949.00     526.75        512      1312.25  523.0
    # 100%  5025.0 3359.00   9441 2019.00    1012.00       8092      8902.00 4283.0
    #      Unknown Unknown.small
    # 0%     556.0           592
    # 25%    774.0           855
    # 50%   1108.0          1065
    # 75%   1534.5          1766
    # 100%  4004.0          4515    - so relatively, these two are quite trxnally active...
    #                                 However, ~all their top markers are ribosomal pt's and,
    #                                 randomly, XIST ('Unknown')
    #                               - HOWEVER, they are all male...:
    
    # From Supplementary Table 6:
        # SM-4RGJU	1	Tissue:Fresh Frozen Tissue	7.3	60	Male	White	Brain - Hippocampus
        # SM-4W9GN	2	Tissue:Fresh Frozen Tissue	6.8	53	Male	White	Brain - Hippocampus
        # SM-DG7EP	3	Tissue:Fresh Frozen Tissue	7.1	59	Male	White	Brain - Hippocampus
        # SM-DG7EQ	4	Tissue:Fresh Frozen Tissue	6.9	65	Male	White	Brain - Hippocampus
    
    sapply(markers.habibHPC.t.1vAll, function(x) which(rownames(x)=="XIST"))
        # Astro.1       Astro.2        Endo.1        Endo.2    Endo.small
        #   26436          5629         26169         18614         27868
        # Glia.mixed  Neuron.mixed         Oligo       Unknown Unknown.small
        #      26538         27914         26231            14            57
    
    sapply(markers.habibHPC.t.1vAll, function(x) which(rownames(x)=="TSIX"))
        # Astro.1       Astro.2        Endo.1        Endo.2    Endo.small
        #   25146         22792         24966         25231         24531
        # Glia.mixed  Neuron.mixed         Oligo       Unknown Unknown.small
        #       7655         24818         14326           337          3573
    

### Comparison to LIBD HPC ==========
# Bring in human stats; create t's
load("rdas/markers-stats_HPC-n3_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.hpc.t.1vAll, markers.hpc.t.design, markers.hpc.wilcox.block
    rm(markers.hpc.t.design, markers.hpc.wilcox.block)

# Need to add t's with N nuclei used in constrasts
load("rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    # sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo
    rm(chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo)

# First drop "Ambig.lowNtrxts" (50 nuclei)
sce.hpc <- sce.hpc[ ,sce.hpc$cellType.split != "Ambig.lowNtrxts"]
sce.hpc$cellType.split <- droplevels(sce.hpc$cellType.split)

## As above, calculate and add t-statistic (= std.logFC * sqrt(N)) from contrasts
#      and fix row order to the first entry "Astro"
fixTo <- rownames(markers.hpc.t.1vAll[["Astro"]])

for(s in names(markers.hpc.t.1vAll)){
  markers.hpc.t.1vAll[[s]]$t.stat <- markers.hpc.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.hpc))
  markers.hpc.t.1vAll[[s]] <- markers.hpc.t.1vAll[[s]][fixTo, ]
}

# Pull out the t's
ts.hpc <- sapply(markers.hpc.t.1vAll, function(x){x$t.stat})
rownames(ts.hpc) <- fixTo



## Then for Habib et al. - fix row order to the first entry "Astro.1"
fixTo <- rownames(markers.habibHPC.t.1vAll[["Astro.1"]])

for(s in names(markers.habibHPC.t.1vAll)){
  markers.habibHPC.t.1vAll[[s]]$t.stat <- markers.habibHPC.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.habib.hpc))
  markers.habibHPC.t.1vAll[[s]] <- markers.habibHPC.t.1vAll[[s]][fixTo, ]
}

# Pull out the t's
ts.habib.hpc <- sapply(markers.habibHPC.t.1vAll, function(x){x$t.stat})
rownames(ts.habib.hpc) <- fixTo


## Take intersecting between two and subset/reorder
sharedGenes <- intersect(rownames(ts.habib.hpc), rownames(ts.hpc))
length(sharedGenes) #16,850

ts.habib.hpc <- ts.habib.hpc[sharedGenes, ]
ts.hpc <- ts.hpc[sharedGenes, ]


cor_t_hpc <- cor(ts.hpc, ts.habib.hpc)
rownames(cor_t_hpc) = paste0(rownames(cor_t_hpc),"_","libd")
colnames(cor_t_hpc) = paste0(colnames(cor_t_hpc),"_","habib")
range(cor_t_hpc)


## Heatmap
theSeq.all = seq(-.80, .80, by = 0.01)
my.col.all <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.all)-1)

pdf("pdfs/exploration/Habib_DroNc-seq/overlap_MNTclusters-to-LIBD-10x-pilot-splitClusters_Jul2020.pdf")
pheatmap(cor_t_hpc,
         color=my.col.all,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.all,
         fontsize=11, fontsize_row=14, fontsize_col=14,
         display_numbers=T, number_format="%.2f", fontsize_number=6.5,
         main="Correlation of cluster-specific t's to MNT-redefined \n clusters of (Habib et al. Nat. Methods 2017)")
dev.off()



head(rownames(markers.habibHPC.t.1vAll[["Unknown"]]),n=100)
# "HIST1H4C" "TXN"      "RPL39"    "MIF"      "SLIRP"    "RPS15A"
# "RPL27"    "RPS29"    "H2AFZ"    "SNRPE"    "RPS13"    "RPL23A"
# "RPS18"    "XIST"     "SNRPG"    "RPL34"    "RPL26"    "PTTG1"
# "COX6C"    "UQCR10"   "RPS7"     "COX7B"    "RPS27A"   "USMG5"
# "TOP2A"    "RPS3"     "HSPE1"    "ROMO1"    "NUSAP1"   "RPL32"
# "RPS14"    "ATP5J2"   "RPS3A"    "CKS2"     "RPL38"    "RPL35"
# "RPA3"     "HIST1H1E" "SNRPD2"   "SNRPD1"   "RPL35A"   "RPL26L1"
# "RPS19"    "RPS6"     "RPL27A"   "RPS21"    "RPL41"    "RPL37A"
# "RPS15"    "RPS4X"    "COX7A2"   "NDUFA1"   "ATP5E"    "COX7C"
#  "RPL31"    "RPL12"    "RPS23"    "ATPIF1"   "RPS12"    "NDUFA4"
# "ATP5I"    "RPL10A"   "UQCRQ"    "RPL24"    "RPLP1"    "C1QBP"
# "RPL10"    "SF3B14"   "SEC61G"   "RPL18A"   "SNHG5"    "RPS28"
# "MID1"     "RPS25"    "NPM1"     "NDUFC2"   "RPL15"    "RPL37"
# "DLGAP5"   "TMEM258"  "NDUFB6"   "RPL36AL"  "RPL23"    "ATP5J"
# "RPS16"    "RPS20"    "SHFM1"    "C7orf55"  "ATP5EP2"  "RPL29"
# "RPL36"    "SLC25A5"  "SEC61B"   "COX17"    "SNRPF"    "UBL5"
# "CENPW"    "POLR2L"   "RPL19"    "RPL30"
    ## Gene Ontology gives no direction for these, either...
     # (as per http://bioinformatics.sdstate.edu/go/)


### Compare MNT annotations to Habib annotations ============
table(sce.habib.hpc$ClusterName, sce.habib.hpc$cellType.mnt)
    #              Astro.1 Astro.2 Endo.1 Endo.2 Endo.small Glia.mixed Neuron.mixed
    # ASC1             501       0      0      0          0        142           16
    # ASC2              91       0      0      0          0         73           10
    # END                0       0     77     42          6         10            3
    # exCA1              0       0      0      0          0         17          326
    # exCA3              0       0      0      0          0          6          656
    # exDG               0       0      0      0          0          7          409
    # exPFC1             9       1      0      0          0         32           77
    # exPFC2             1       0      0      0          0         14            5
    # GABA1              0       0      0      0          0          6          181
    # GABA2              0       0      0      0          0          2          201
    # MG                 2       0      0      0          0        203           50
    # NSC                0      63      0      0          0          7            0
    # ODC1               0       0      0      0          0         82            7
    # OPC                1       0      0      0          0        189           90
    # Unclassified       3       2      0      0          0         47          337
    # 
    #              Oligo Unknown Unknown.small
    # ASC1            11       0             0
    # ASC2             8       0             0
    # END              0       0             0
    # exCA1            7       0             0
    # exCA3            1       0             0
    # exDG             0       0             0
    # exPFC1          21       0             0
    # exPFC2          19       0             0
    # GABA1            1       0             0
    # GABA2            0       0             0
    # MG              16       0             0
    # NSC              0       0             0
    # ODC1          1608       0             0
    # OPC             24       0             0
    # Unclassified    38     132            37    * good that my 'Unknown's -> 'Unclassified'

    # From their supplementary text/Methods, it looks like they saw this in mouse (Fig. S3)
    #    at leaset and "thus removed from subsequent analyses." - probably same thing happened
    #    for human dataset too, as there are no 'Unclass.' or anything on human tSNEs

    #    -->  Might be best to run through this again, AFTER dropping those 'Unclassified'
    #         cells...




### Revised cell type annotations for comparisons =======================================
  # MNT 26Jul2020    -> Use the provided annotations and drop 'Unclassified'
  #                     (just wanna see how this looks, even though 'exDG' had waaay more DLPFC...)

## As before:
# Drop genes with all 0's
sce.habib.hpc <- sce.habib.hpc[!rowSums(assay(sce.habib.hpc, "counts"))==0, ]

# Drop 'Unclassified'
sce.habib.hpc$ClusterName <- factor(sce.habib.hpc$ClusterName)
sce.habib.hpc <- sce.habib.hpc[ ,!sce.habib.hpc$ClusterName == "Unclassified"]
sce.habib.hpc$ClusterName <- droplevels(sce.habib.hpc$ClusterName)

# Model unwanted effects
mod <- with(colData(sce.habib.hpc), model.matrix(~ SampleID))
mod <- mod[ ,-1]

markers.habibHPC.t.1vAll.pub <- list()
for(i in levels(sce.habib.hpc$ClusterName)){
  # Make temporary contrast
  sce.habib.hpc$contrast <- ifelse(sce.habib.hpc$ClusterName==i, 1, 0)
  # Test cluster vs. all
  markers.habibHPC.t.1vAll.pub[[i]] <- findMarkers(sce.habib.hpc, groups=sce.habib.hpc$contrast,
                                               assay.type="logcounts", design=mod, test="t",
                                               direction="up", pval.type="all", full.stats=T)
}

    ## Then, temp set of stats to get the standardized logFC
    temp.1vAll <- list()
    for(i in levels(sce.habib.hpc$ClusterName)){
      # Make temporary contrast
      sce.habib.hpc$contrast <- ifelse(sce.habib.hpc$ClusterName==i, 1, 0)
      # Test cluster vs. all
      temp.1vAll[[i]] <- findMarkers(sce.habib.hpc, groups=sce.habib.hpc$contrast,
                                     assay.type="logcounts", design=mod, test="t",
                                     std.lfc=TRUE,
                                     direction="up", pval.type="all", full.stats=T)
    }

## For some reason all the results are in the second List entry (first is always empty)

# Replace that empty slot with the entry with the actul stats
markers.habibHPC.t.1vAll.pub <- lapply(markers.habibHPC.t.1vAll.pub, function(x){ x[[2]] })
# Same for that with std.lfc
temp.1vAll <- lapply(temp.1vAll, function(x){ x[[2]] })

# Now just pull from the 'stats.0' DataFrame column
markers.habibHPC.t.1vAll.pub <- lapply(markers.habibHPC.t.1vAll.pub, function(x){ x$stats.0 })
temp.1vAll <- lapply(temp.1vAll, function(x){ x$stats.0 })

# Re-name std.lfc column and add to the first result
for(i in names(temp.1vAll)){
  colnames(temp.1vAll[[i]])[1] <- "std.logFC"
  markers.habibHPC.t.1vAll.pub[[i]] <- cbind(markers.habibHPC.t.1vAll.pub[[i]], temp.1vAll[[i]]$std.logFC)
  # Oh the colname is kept weird
  colnames(markers.habibHPC.t.1vAll.pub[[i]])[4] <- "std.logFC"
  # Then re-organize
  markers.habibHPC.t.1vAll.pub[[i]] <- markers.habibHPC.t.1vAll.pub[[i]][ ,c("logFC","std.logFC","log.p.value","log.FDR")]
}


## Let's save this along with the previous pairwise results
save(markers.habibHPC.t.1vAll.pub, file="rdas/zs-habib_markers-stats_reportedClusters_findMarkers-SN-LEVEL_Jul2020.rda")


sapply(markers.habibHPC.t.1vAll.pub, function(x){table(x$log.FDR < log10(0.000001))})
    #        ASC1  ASC2   END exCA1 exCA3  exDG exPFC1 exPFC2 GABA1 GABA2    MG   NSC
    # FALSE 26948 27406 26224 25994 25923 26635  27820  27868 26926 26809 26575 26974
    # TRUE   1201   743  1925  2155  2226  1514    329    281  1223  1340  1574  1175
    #        ODC1   OPC
    # FALSE 27106 27083
    # TRUE   1043  1066


## Print those top 40 into table so can compare to human markers
markerList.t.habibHPC.1vAll <- lapply(markers.habibHPC.t.1vAll.pub, function(x){
  rownames(x)[x[ ,"log.FDR"] < log10(0.000001)]
  }
)

top40genes <- sapply(markerList.t.habibHPC.1vAll, function(x) head(x, n=40))
write.csv(top40genes,"tables/forRef_top40genesLists_habib-HIPPO_reportedClusters_SN-LEVEL-tests_Jul2020.csv")


### Comparison to LIBD HPC (again) ==========
# Bring in human stats; create t's
load("rdas/markers-stats_HPC-n3_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.hpc.t.1vAll, markers.hpc.t.design, markers.hpc.wilcox.block
    rm(markers.hpc.t.design, markers.hpc.wilcox.block)

# Need to add t's with N nuclei used in constrasts
load("rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    # sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo
    rm(chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo)

# First drop "Ambig.lowNtrxts" (50 nuclei)
sce.hpc <- sce.hpc[ ,sce.hpc$cellType.split != "Ambig.lowNtrxts"]
sce.hpc$cellType.split <- droplevels(sce.hpc$cellType.split)

## As above, calculate and add t-statistic (= std.logFC * sqrt(N)) from contrasts
#      and fix row order to the first entry "Astro"
fixTo <- rownames(markers.hpc.t.1vAll[["Astro"]])

for(s in names(markers.hpc.t.1vAll)){
  markers.hpc.t.1vAll[[s]]$t.stat <- markers.hpc.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.hpc))
  markers.hpc.t.1vAll[[s]] <- markers.hpc.t.1vAll[[s]][fixTo, ]
}

# Pull out the t's
ts.hpc <- sapply(markers.hpc.t.1vAll, function(x){x$t.stat})
rownames(ts.hpc) <- fixTo



## Then for Habib et al. - fix row order to the first entry "Astro.1"
fixTo <- rownames(markers.habibHPC.t.1vAll.pub[["ASC1"]])

for(s in names(markers.habibHPC.t.1vAll.pub)){
  markers.habibHPC.t.1vAll.pub[[s]]$t.stat <- markers.habibHPC.t.1vAll.pub[[s]]$std.logFC * sqrt(ncol(sce.habib.hpc))
  markers.habibHPC.t.1vAll.pub[[s]] <- markers.habibHPC.t.1vAll.pub[[s]][fixTo, ]
}

# Pull out the t's
ts.habib.hpc <- sapply(markers.habibHPC.t.1vAll.pub, function(x){x$t.stat})
rownames(ts.habib.hpc) <- fixTo


## Take intersecting between two and subset/reorder
sharedGenes <- intersect(rownames(ts.habib.hpc), rownames(ts.hpc))
length(sharedGenes) #16,850

ts.habib.hpc <- ts.habib.hpc[sharedGenes, ]
ts.hpc <- ts.hpc[sharedGenes, ]


cor_t_hpc <- cor(ts.hpc, ts.habib.hpc)
rownames(cor_t_hpc) = paste0(rownames(cor_t_hpc),"_","libd")
colnames(cor_t_hpc) = paste0(colnames(cor_t_hpc),"_","habib")
range(cor_t_hpc)
    # [1] -0.5068500  0.7739983

## Heatmap
theSeq.all = seq(-.80, .80, by = 0.01)
my.col.all <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.all)-1)

pdf("pdfs/exploration/Habib_DroNc-seq/overlap_reportedClusters-to-LIBD-10x-pilot-splitClusters_Jul2020.pdf")
pheatmap(cor_t_hpc,
         color=my.col.all,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.all,
         fontsize=11, fontsize_row=14, fontsize_col=14,
         display_numbers=T, number_format="%.2f", fontsize_number=6.5,
         legend_breaks=c(seq(-0.8,0.8,by=0.4)),
         main="Correlation of cluster-specific t's to reported clusters \n in (Habib et al. Nat. Methods 2017)")
dev.off()





