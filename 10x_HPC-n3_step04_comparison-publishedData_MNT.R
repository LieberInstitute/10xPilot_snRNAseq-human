### MNT 10x snRNA-seq workflow: step 04
###   **Region-specific analyses**
###     - (3x) HPC samples
###     - Setup and comparison to Habib, et al (DroNc-seq paper)
### Updated for revision MNT 2021
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



### Comparison to LIBD HPC ==============================
# Bring in human stats; create t's
load("rdas/revision/markers-stats_HPC-n3_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.hpc.t.pw, markers.hpc.t.1vAll, medianNon0.hpc
    rm(markers.hpc.t.pw, medianNon0.hpc)

# Need to add t's with N nuclei used in constrasts
load("rdas/revision/regionSpecific_HPC-n3_cleaned-combined_SCE_MNT2021.rda", verbose=T)
    # sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, ref.sampleInfo, clusterRefTab.hpc, annotationTab.hpc, cell_colors.hpc
    rm(chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo)

sce.hpc
    # class: SingleCellExperiment 
    # dim: 33538 10268 
    # metadata(2): merge.info pca.info
    # assays(2): counts logcounts
    # rownames(33538): MIR1302-2HG FAM138A ... AC213203.1 FAM231C
    # rowData names(6): gene_id gene_version ... gene_biotype Symbol.uniq
    # colnames(10268): AAACCCATCTGTCAGA-1 AAACCCATCTGTCGCT-1 ...
    #   TTTGGTTGTGGTCCGT-1 TTTGTTGCAGAAACCG-1
    # colData names(20): Sample Barcode ... collapsedCluster cellType
    # reducedDimNames(4): PCA_corrected PCA_opt TSNE UMAP
    # altExpNames(0):

table(sce.hpc$cellType)
    # Astro_A       Astro_B  drop.doublet drop.lowNTx_A drop.lowNTx_B 
    #     936           234             5           105            19 
    # Excit_A       Excit_B       Excit_C       Excit_D       Excit_E 
    #      87           421             6            35             6 
    # Excit_F       Excit_G       Excit_H       Inhib_A       Inhib_B 
    #      29             6            33           300            30 
    # Inhib_C       Inhib_D         Micro         Mural         Oligo 
    #       5            31          1161            43          5912 
    #     OPC       OPC_COP         Tcell 
    #     823            15            26 

# First drop flagged "drop." clusters
sce.hpc <- sce.hpc[ ,-grep("drop.", sce.hpc$cellType)]
sce.hpc$cellType <- droplevels(sce.hpc$cellType)

## As above, calculate and add t-statistic (= std.logFC * sqrt(N)) from contrasts
#      and fix row order to the first entry "Astro"
markers.hpc.enriched <- lapply(markers.hpc.t.1vAll, function(x){x[[2]]})

fixTo <- rownames(markers.hpc.enriched[["Astro_A"]])

for(s in names(markers.hpc.enriched)){
  markers.hpc.enriched[[s]]$t.stat <- markers.hpc.enriched[[s]]$std.logFC * sqrt(ncol(sce.hpc))
  markers.hpc.enriched[[s]] <- markers.hpc.enriched[[s]][fixTo, ]
}

# Pull out the t's
ts.hpc <- sapply(markers.hpc.enriched, function(x){x$t.stat})
rownames(ts.hpc) <- fixTo



## Then for Habib et al. - fix row order to the first entry "Astro.1" ======
load("rdas/zs-habib_markers-stats_reportedClusters_findMarkers-SN-LEVEL_Jul2020.rda", verbose=T)
    # markers.habibHPC.t.1vAll.pub
names(markers.habibHPC.t.1vAll.pub)

load("rdas/zSCE_habib-dlpfc-hippo_MNT.rda", verbose=T)
    # sce.habib, sce.habib.hpc, sce.habib.dlpfc
    rm(sce.habib, sce.habib.dlpfc)

sce.habib.hpc
    # class: SingleCellExperiment 
    # dim: 32111 5927 
    # metadata(0):
    # assays(2): counts logcounts
    # rownames(32111): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
    # rowData names(0):
    # colnames(5927): hHP1_AACACTATCTAC hHP1_CTACGCATCCAT ...
    #   HP3-B-united_GGTAATAAGTTG HP3-B-united_GCCCCAAAGGAT
    # colData names(13): CellID NumGenes ... DonorID sizeFactor
    # reducedDimNames(1): TSNE.given
    # altExpNames(0):
    
# Drop 'Unclassified', as above
sce.habib.hpc$ClusterName <- factor(sce.habib.hpc$ClusterName)
sce.habib.hpc <- sce.habib.hpc[ ,!sce.habib.hpc$ClusterName == "Unclassified"]
sce.habib.hpc$ClusterName <- droplevels(sce.habib.hpc$ClusterName)

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
length(sharedGenes) #16,854

ts.habib.hpc <- ts.habib.hpc[sharedGenes, ]
ts.hpc <- ts.hpc[sharedGenes, ]


cor_t_hpc <- cor(ts.hpc, ts.habib.hpc)
rownames(cor_t_hpc) = paste0(rownames(cor_t_hpc),"_","libd")
colnames(cor_t_hpc) = paste0(colnames(cor_t_hpc),"_","habib")
range(cor_t_hpc)
    # [1] -0.5096868  0.7770631

## Heatmap
theSeq.all = seq(-.80, .80, by = 0.01)
my.col.all <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.all)-1)

pdf("pdfs/revision/pubFigures/overlap_reportedClusters-to-LIBD-HPC_allExpressedGenes_MNT2021.pdf")
pheatmap(cor_t_hpc,
         color=my.col.all,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.all,
         fontsize=11, fontsize_row=14, fontsize_col=14,
         display_numbers=T, number_format="%.2f", fontsize_number=6.5,
         legend_breaks=c(seq(-0.8,0.8,by=0.4)),
         main="Correlation of cluster-specific t's to reported clusters \n in (Habib et al. Nat. Methods 2017)")
dev.off()


### Session info for 23Jul2021 =======================================
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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices datasets  utils     methods  
# [9] base     
# 
# other attached packages:
#   [1] pheatmap_1.0.12             RColorBrewer_1.1-2          lattice_0.20-41            
# [4] readxl_1.3.1                readr_1.4.0                 jaffelab_0.99.30           
# [7] rafalib_1.0.0               DropletUtils_1.10.3         batchelor_1.6.3            
# [10] scran_1.18.7                scater_1.18.6               ggplot2_3.3.3              
# [13] SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0 Biobase_2.50.0             
# [16] GenomicRanges_1.42.0        GenomeInfoDb_1.26.7         IRanges_2.24.1             
# [19] S4Vectors_0.28.1            BiocGenerics_0.36.1         MatrixGenerics_1.2.1       
# [22] matrixStats_0.58.0         
# 
# loaded via a namespace (and not attached):
#   [1] viridis_0.6.0             edgeR_3.32.1              BiocSingular_1.6.0       
# [4] splines_4.0.4             viridisLite_0.4.0         DelayedMatrixStats_1.12.3
# [7] scuttle_1.0.4             R.utils_2.10.1            assertthat_0.2.1         
# [10] statmod_1.4.35            dqrng_0.3.0               cellranger_1.1.0         
# [13] GenomeInfoDbData_1.2.4    vipor_0.4.5               pillar_1.6.0             
# [16] glue_1.4.2                limma_3.46.0              beachmat_2.6.4           
# [19] XVector_0.30.0            colorspace_2.0-0          Matrix_1.3-4             
# [22] R.oo_1.24.0               pkgconfig_2.0.3           zlibbioc_1.36.0          
# [25] purrr_0.3.4               scales_1.1.1              HDF5Array_1.18.1         
# [28] ResidualMatrix_1.0.0      BiocParallel_1.24.1       googledrive_1.0.1        
# [31] tibble_3.1.1              generics_0.1.0            ellipsis_0.3.2           
# [34] withr_2.4.2               magrittr_2.0.1            crayon_1.4.1             
# [37] R.methodsS3_1.8.1         fansi_0.4.2               segmented_1.3-4          
# [40] bluster_1.0.0             beeswarm_0.4.0            tools_4.0.4              
# [43] hms_1.0.0                 lifecycle_1.0.0           Rhdf5lib_1.12.1          
# [46] munsell_0.5.0             locfit_1.5-9.4            DelayedArray_0.16.3      
# [49] irlba_2.3.3               compiler_4.0.4            rsvd_1.0.5               
# [52] rlang_0.4.11              rhdf5_2.34.0              grid_4.0.4               
# [55] RCurl_1.98-1.3            rhdf5filters_1.2.0        BiocNeighbors_1.8.2      
# [58] igraph_1.2.6              bitops_1.0-7              gtable_0.3.0             
# [61] DBI_1.1.1                 R6_2.5.0                  gridExtra_2.3            
# [64] dplyr_1.0.5               utf8_1.2.1                ggbeeswarm_0.6.0         
# [67] Rcpp_1.0.6                vctrs_0.3.8               tidyselect_1.1.1         
# [70] sparseMatrixStats_1.2.1


