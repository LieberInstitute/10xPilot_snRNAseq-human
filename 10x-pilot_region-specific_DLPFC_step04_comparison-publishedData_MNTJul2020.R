### MNT 10x snRNA-seq workflow: step 04
###   **Region-specific analyses**
###     - (2x) DLPFC samples from: Br5161 & Br5212
###     - Setup and comparison to Mathys, et al (AZD snRNA-seq paper)
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
library(Matrix)
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

# Set the path
path <- '/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/'
# Read in data
pd = read.csv(paste0(path,"mathys/snRNAseqPFC_BA10_biospecimen_metadata.csv"), as.is=TRUE)
    # of dim 48 x 9
pheno = read.delim(paste0(path,"mathys/filtered_column_metadata.txt"), row.names = 1)
    # of dim 70634 x 6

## MNT explore:
colnames(pheno)
# [1] "projid"          "tsne1"           "tsne2"           "pre.cluster"
# [5] "broad.cell.type" "Subcluster"
#     - 48 unique 'projid's - these are synonymous with pd$specimenID--i.e. donors

table(pheno$broad.cell.type)
      # Ast   End    Ex    In   Mic   Oli   Opc   Per
      #3392   121 34976  9196  1920 18235  2627   167
table(pheno$pre.cluster)
    #     1     2     3     4     5     6     7     8     9    10    11    12    13
    # 15113 10612  5315  4938  4142  3615  3392  3122  3197  2450  2780  2627  1920
    #    14    15    16    17    19    20    21
    #  1789  1782  1698  1018   288   424   412   - I think they annotated, THEN
    #                                               looked at clustering w/IN each of these

## If you look at 'Subcluster' x donor, you see clear donor bias in these 'subclusters'
    # (this is acknowledged in paper)
rownames(as.matrix(table(pheno$Subcluster, pheno$projid)))
    # [1] "Ast0" "Ast1" "Ast2" "Ast3" "End1" "End2" "Ex0"  "Ex1"  "Ex11" "Ex12"
    # [11] "Ex14" "Ex2"  "Ex3"  "Ex4"  "Ex5"  "Ex6"  "Ex7"  "Ex8"  "Ex9"  "In0"
    # [21] "In1"  "In10" "In11" "In2"  "In3"  "In4"  "In5"  "In6"  "In7"  "In8"
    # [31] "In9"  "Mic0" "Mic1" "Mic2" "Mic3" "Oli0" "Oli1" "Oli3" "Oli4" "Oli5"
    # [41] "Opc0" "Opc1" "Opc2" "Per"
unname(as.matrix(table(pheno$Subcluster, pheno$projid)))  # only prints top 31
tail(unname(as.matrix(table(pheno$Subcluster, pheno$projid))), n=13)


dat = readMM(paste0(path,"mathys/filtered_count_matrix.mtx"))
    # is of 17926 x 70634 dim
genes = read.delim(paste0(path,"mathys/filtered_gene_row_names.txt"),header=FALSE,as.is=TRUE)

## add names
rownames(dat) = genes$V1
colnames(dat) = rownames(pheno)


## add more pheno
pheno$individualID = pd$individualID[match(pheno$projid, pd$projid)]

pd2 = read.csv(paste0(path,"mathys/ROSMAP_Clinical_2019-05_v3.csv"),as.is=TRUE)
pd2$age_death[pd2$age_death == "90+"] = 90
pd2$age_death = as.numeric(pd2$age_death)
pd2$Dx = factor(ifelse(pd2$age_first_ad_dx == "", "Control", "AD"),
                levels = c("Control", "AD"))


# Expands to dims of per-nucleus 'pheno'
pd2 = pd2[match(pheno$individualID, pd2$individualID),]
pheno$Dx = pd2$Dx
pheno$age_death = pd2$age_death
pheno$msex = pd2$msex
pheno$race = pd2$race

## Create SCE
sce.mathys <- SingleCellExperiment(assays = list(counts = dat),
                                   colData = pheno,
                                   rowData = genes)

save(sce.mathys, file="rdas/referenceDatasets/SCE_mathys-PFC-BA10_MNT.rda")



## Process, take provided annotations and make markers/stats
LSFvec <- librarySizeFactors(sce.mathys)
sce.mathys <- logNormCounts(sce.mathys, size_factors=LSFvec)

# Pull tSNE coordinates out of provided pd/pheno
reducedDim(sce.mathys, "TSNE.given") <- as.matrix(colData(sce.mathys)[ ,c("tsne1", "tsne2")])
# and take those out of the colData
colData(sce.mathys) <- colData(sce.mathys)[ ,-c(2:3)]

# Any all-0 genes? shouldn't be...
table(rowSums(assay(sce.mathys, "counts"))==0)  # 3 TRUE...
# Drop these
sce.mathys <- sce.mathys[!rowSums(assay(sce.mathys, "counts"))==0, ]


    # Re-save with normalized counts and assigned reducedDim slot
    save(sce.mathys, file="rdas/referenceDatasets/SCE_mathys-PFC-BA10_MNT.rda")

    
    

### Run ANOVA real quick (oh, just do in pseudo-bulked space with prev-processed data) ====
library(edgeR)
library(doMC)
registerDoMC(cores=4)

#mat = assays(sce.mathys)$logcounts
load("rdas/referenceDatasets/mathys_broad-pseudobulked.rda", verbose=T)
    # sce_pseudobulk.broad, corfit

mat = assays(sce_pseudobulk.broad)$logcounts

## do regression (** takes too long currently bc a lot of nuclei... might want to run in PB space)
varCompAnalysis.pb.ind = foreach(i = 1:nrow(mat)) %dopar% {
  if(i %% 100 == 0) cat(".")
  #fit = lm(as.numeric(mat[i,]) ~ broad.cell.type + Dx + age_death + msex + race,
  fit = lm(as.numeric(mat[i,]) ~ broad.cell.type + Dx + age_death + msex + race + individualID,
           data=colData(sce_pseudobulk.broad))
  full = anova(fit)
  fullSS = full$"Sum Sq"
  signif(cbind(full, PctExp = fullSS/sum(fullSS)*100), 3)
}

names(varCompAnalysis.pb) = rownames(mat)
names(varCompAnalysis.pb.ind) = rownames(mat)

## make boxplot
varExpl.pb = t(sapply(varCompAnalysis.pb, function(x) x[,"PctExp"]))
colnames(varExpl.pb) = rownames(varCompAnalysis.pb[[1]])

varExpl.pb.ind = t(sapply(varCompAnalysis.pb.ind, function(x) x[,"PctExp"]))
colnames(varExpl.pb.ind) = rownames(varCompAnalysis.pb.ind[[1]])


pdf("pdfs/exploration/anova_mathys-PFC-BA-10_pseudo-bulked_Jul2020.pdf")
boxplot(varExpl.pb, main="ANOVA on human PFC 10x-snRNA-seq \n (Mathys, et al. Nature 2019)",
        ylab="Percent Var explained (%))")
boxplot(varExpl.pb.ind, main="ANOVA on human PFC 10x-snRNA-seq \n (Mathys, et al. Nature 2019)",
        ylab="Percent Var explained (%))")
dev.off()
    # ok so DEF keep 'individualID'

save(varCompAnalysis.pb, varCompAnalysis.pb.ind, file="rdas/referenceDatasets/mathys_anova-output_pseudo-bulked_MNTJul2020.rda")

apply(varExpl.pb, 2, function(x){quantile(x, na.rm=T)})
    #      broad.cell.type       Dx age_death     msex     race Residuals
    # 0%             0.272 1.33e-09  3.79e-09 7.22e-11 2.91e-09      4.56
    # 25%           20.700 2.32e-02  2.64e-02 2.19e-02 1.41e-02     50.05
    # 50%           35.600 1.03e-01  1.19e-01 1.02e-01 6.45e-02     63.40
    # 75%           49.200 3.00e-01  3.30e-01 3.08e-01 1.87e-01     78.20
    # 100%          95.400 3.54e+00  4.32e+00 7.15e+01 9.21e+00     99.40

apply(varExpl.pb.ind, 2, function(x){quantile(x, na.rm=T)})
    #      broad.cell.type       Dx age_death     msex     race individualID Residuals
    # 0%             0.272 1.33e-09  3.79e-09 7.22e-11 2.91e-09        0.584      3.84
    # 25%           20.700 2.32e-02  2.64e-02 2.19e-02 1.41e-02        7.230     41.90
    # 50%           35.600 1.03e-01  1.19e-01 1.02e-01 6.45e-02        9.570     53.30
    # 75%           49.200 3.00e-01  3.30e-01 3.08e-01 1.87e-01       11.700     65.90
    # 100%          95.400 3.54e+00  4.32e+00 7.15e+01 9.21e+00       61.200     88.30



### Modeling for subcluster-specific genes - cluster-vs-all test, only === === ===
# (only doing this iteration because this yields t-statistics to compare to human)

  # *** Pseudo-bulk values bc this otherwise takes too long with 70k nuclei:
    # Make the pseudo-bulked SCE
    sce.mathys.PB <- aggregateAcrossCells(sce.mathys, ids=paste0(sce.mathys$individualID,":",sce.mathys$broad.cell.type),
                                         use_exprs_values="counts")
    sce.mathys.PB
    
    # Remove stored `sizeFactors()` because this will mess you up
    #     * Also, to be safe, can always provide manually-computed SFs:
    sizeFactors(sce.mathys.PB) <- NULL
    LSFvec <- librarySizeFactors(sce.mathys.PB)
    sce.mathys.PB <- logNormCounts(sce.mathys.PB, size_factors=LSFvec)

    # Save this along with the sn-level SCE
    save(sce.mathys, sce.mathys.PB, file="rdas/referenceDatasets/SCE_mathys-PFC-BA10_MNT.rda")
    
    
    # Drop genes with all 0's
    sce.mathys.PB <- sce.mathys.PB[!rowSums(assay(sce.mathys.PB, "counts"))==0, ]
        ## keeps 28128 genes


# Model unwanted effects
#mod <- with(colData(sce.mathys.PB), model.matrix(~ individualID + Dx + age_death + msex))
#mod <- with(colData(sce.mathys.PB), model.matrix(~ individualID + Dx + msex))
#mod <- with(colData(sce.mathys.PB), model.matrix(~ individualID + msex))
mod <- with(colData(sce.mathys.PB), model.matrix(~ individualID))
mod <- mod[ ,-1]
    # Error in .ranksafe_qr(full.design) : design matrix is not of full rank
        # get this when remove 'age_death' && 'Dx'... howEVER, looks like it runs fine only
        #     modeling on 'individualID'

# Factorize
sce.mathys.PB$broad.cell.type <- factor(sce.mathys.PB$broad.cell.type)

markers.mathysPFC.t.1vAll <- list()
for(i in levels(sce.mathys.PB$broad.cell.type)){
  # Make temporary contrast
  sce.mathys.PB$contrast <- ifelse(sce.mathys.PB$broad.cell.type==i, 1, 0)
  # Test cluster vs. all
  markers.mathysPFC.t.1vAll[[i]] <- findMarkers(sce.mathys.PB, groups=sce.mathys.PB$contrast,
                                               assay.type="logcounts", design=mod, test="t",
                                               direction="up", pval.type="all", full.stats=T)
}

## Then, temp set of stats to get the standardized logFC
    temp.1vAll <- list()
    for(i in levels(sce.mathys.PB$broad.cell.type)){
      # Make temporary contrast
      sce.mathys.PB$contrast <- ifelse(sce.mathys.PB$broad.cell.type==i, 1, 0)
      # Test cluster vs. all
      temp.1vAll[[i]] <- findMarkers(sce.mathys.PB, groups=sce.mathys.PB$contrast,
                                     assay.type="logcounts", design=mod, test="t",
                                     std.lfc=TRUE,
                                     direction="up", pval.type="all", full.stats=T)
    }

    ## For some reason all the results are in the second List entry (first is always empty)

# Replace that empty slot with the entry with the actul stats
markers.mathysPFC.t.1vAll <- lapply(markers.mathysPFC.t.1vAll, function(x){ x[[2]] })
# Same for that with std.lfc
temp.1vAll <- lapply(temp.1vAll, function(x){ x[[2]] })

# Now just pull from the 'stats.0' DataFrame column
markers.mathysPFC.t.1vAll <- lapply(markers.mathysPFC.t.1vAll, function(x){ x$stats.0 })
temp.1vAll <- lapply(temp.1vAll, function(x){ x$stats.0 })

# Re-name std.lfc column and add to the first result
for(i in names(temp.1vAll)){
  colnames(temp.1vAll[[i]])[1] <- "std.logFC"
  markers.mathysPFC.t.1vAll[[i]] <- cbind(markers.mathysPFC.t.1vAll[[i]], temp.1vAll[[i]]$std.logFC)
  # Oh the colname is kept weird
  colnames(markers.mathysPFC.t.1vAll[[i]])[4] <- "std.logFC"
  # Then re-organize
  markers.mathysPFC.t.1vAll[[i]] <- markers.mathysPFC.t.1vAll[[i]][ ,c("logFC","std.logFC","log.p.value","log.FDR")]
}

sapply(markers.mathysPFC.t.1vAll, function(x){table(x$log.FDR < log10(0.000001))})
    #        Ast   End    Ex    In   Mic   Oli   Opc   Per
    # FALSE 14752 17760  7872 10544 16217 12414 15770 17732
    # TRUE   3171   163 10051  7379  1706  5509  2153   191

sapply(markers.mathysPFC.t.1vAll, function(x){head(rownames(x), n=20)})
    # ok looks pretty good!

## Let's save this along with the previous pairwise results
save(markers.mathysPFC.t.1vAll, file="rdas/referenceDatasets/zs-mathys_markers-stats_given-clusters_PB-findMarkers-SN-LEVEL_Aug2020.rda")






### Comparison to Velmeshev, et al (PFC & ACC) ========
## Load within-PFC statistics
load("/dcl01/ajaffe/data/lab/singleCell/velmeshev2019/analysis_MNT/markers-stats_velmeshev-et-al_ASD-cortex-withinRegion_findMarkers-SN-LEVEL_MNTAug2020.rda",
     verbose=T)
    # markers.asdVelm.t.pfc, markers.asdVelm.t.acc
    #rm(markers.asdVelm.t.acc)

# And corresponding SCE (to generate t.stat's)
load("/dcl01/ajaffe/data/lab/singleCell/velmeshev2019/analysis_MNT/SCE_asd-velmeshev-etal_MNT.rda", verbose=T)
# sce.asd

sce.asd.pfc <- sce.asd[ ,sce.asd$region=="PFC"]
sce.asd.acc <- sce.asd[ ,sce.asd$region=="ACC"]

# Need to convert Symbol in sce.dlpfc > EnsemblID, and also use n nuclei for t.stat
load("rdas/regionSpecific_DLPFC-n2_SCE_cellTypesSplit-fromST_Apr2020.rda", verbose=T)
    # sce.dlpfc.st, clusterRefTab.dlpfc, chosen.hvgs.dlpfc, ref.sampleInfo

# Drop "Ambig.lowNtrxts" (168 nuclei)
sce.dlpfc.st <- sce.dlpfc.st[ ,sce.dlpfc.st$cellType.split != "Ambig.lowNtrxts"]
sce.dlpfc.st$cellType.split <- droplevels(sce.dlpfc.st$cellType.split)

## Load LIBD DLPFC stats (don't need the pw result)
load("rdas/markers-stats_DLPFC_n2_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.dlpfc.t.1vAll, markers.dlpfc.t.design
    rm(markers.dlpfc.t.design)

for(i in names(markers.dlpfc.t.1vAll)){
  rownames(markers.dlpfc.t.1vAll[[i]]) <- rowData(sce.dlpfc.st)$ID[match(rownames(markers.dlpfc.t.1vAll[[i]]),
                                                                         rownames(sce.dlpfc.st))]
}

    # nrow(markers.asdVelm.t.pfc[[1]])  # 36501
    # nrow(markers.dlpfc.t.1vAll[[1]])  # 28111
    # length(intersect(rownames(markers.asdVelm.t.pfc[[1]]), rownames(markers.dlpfc.t.1vAll[[1]]))) # 26970
    # 
    # fixTo <- intersect(rownames(markers.asdVelm.t.pfc[[1]]), rownames(markers.dlpfc.t.1vAll[[1]]))


## Calculate and add t-statistic (= std.logFC * sqrt(N)) from contrasts
#      and fix row order to the first entry "Astro"
fixTo <- rownames(markers.dlpfc.t.1vAll[["Astro"]])
for(s in names(markers.dlpfc.t.1vAll)){
  markers.dlpfc.t.1vAll[[s]]$t.stat <- markers.dlpfc.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.dlpfc.st))
  markers.dlpfc.t.1vAll[[s]] <- markers.dlpfc.t.1vAll[[s]][fixTo, ]
}

# Pull out the t's
ts.dlpfc <- sapply(markers.dlpfc.t.1vAll, function(x){x$t.stat})
rownames(ts.dlpfc) <- fixTo



## Then for Velmeshev et al. - fix row order to the first entry "AST-FB"
fixTo <- rownames(markers.asdVelm.t.pfc[["AST-FB"]])

for(s in names(markers.asdVelm.t.pfc)){
  markers.asdVelm.t.pfc[[s]]$t.stat <- markers.asdVelm.t.pfc[[s]]$std.logFC * sqrt(ncol(sce.asd.pfc))
  markers.asdVelm.t.pfc[[s]] <- markers.asdVelm.t.pfc[[s]][fixTo, ]
}

# Pull out the t's
ts.velmeshev.pfc <- sapply(markers.asdVelm.t.pfc, function(x){x$t.stat})
rownames(ts.velmeshev.pfc) <- fixTo


## Take intersecting between two and subset/reorder
sharedGenes <- intersect(rownames(ts.velmeshev.pfc), rownames(ts.dlpfc))
length(sharedGenes) # 26,970

ts.velmeshev.pfc <- ts.velmeshev.pfc[sharedGenes, ]
ts.dlpfc <- ts.dlpfc[sharedGenes, ]


cor_t_dlpfc <- cor(ts.dlpfc, ts.velmeshev.pfc)
rownames(cor_t_dlpfc) = paste0(rownames(cor_t_dlpfc),"_","libd")
colnames(cor_t_dlpfc) = paste0(colnames(cor_t_dlpfc),"_","asd.pfc")
range(cor_t_dlpfc)

# For some reason not in alphabetical order like the others...
cor_t_dlpfc <- cor_t_dlpfc[sort(rownames(cor_t_dlpfc)), ]

## Heatmap
theSeq.all = seq(-.95, .95, by = 0.01)
my.col.all <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.all)-1)

pdf("pdfs/exploration/Velmeshev-ASD_pfc-acc/overlap-velmeshev-ASD-pfc_with_LIBD-10x-DLPFC_Aug2020.pdf")
pheatmap(cor_t_dlpfc,
         color=my.col.all,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.all,
         fontsize=9.5, fontsize_row=11, fontsize_col=11,
         display_numbers=T, number_format="%.2f", fontsize_number=5.0,
         legend_breaks=c(seq(-0.95,0.95,by=0.475)),
         main="Correlation of cluster-specific t's between LIBD DLPFC to \n PFC from (Velmeshev et al. Science 2019)")
dev.off()




### What if compared between both the .acc set of stats vs the .pfc?? =============

## Set up ACC t's
fixTo <- rownames(markers.asdVelm.t.acc[["AST-FB"]])

for(s in names(markers.asdVelm.t.acc)){
  markers.asdVelm.t.acc[[s]]$t.stat <- markers.asdVelm.t.acc[[s]]$std.logFC * sqrt(ncol(sce.asd.acc))
  markers.asdVelm.t.acc[[s]] <- markers.asdVelm.t.acc[[s]][fixTo, ]
}

# Pull out the t's
ts.velmeshev.acc <- sapply(markers.asdVelm.t.acc, function(x){x$t.stat})
rownames(ts.velmeshev.acc) <- fixTo

sharedGenes.all <- intersect(rownames(ts.velmeshev.pfc), rownames(ts.dlpfc))
sharedGenes.all <- intersect(sharedGenes.all, rownames(ts.velmeshev.acc))
    # of length 26,970

# Subset/order
ts.dlpfc <- ts.dlpfc[sharedGenes.all, ]
ts.velmeshev.pfc <- ts.velmeshev.pfc[sharedGenes.all, ]
ts.velmeshev.acc <- ts.velmeshev.acc[sharedGenes.all, ]

colnames(ts.velmeshev.pfc) <- paste0(colnames(ts.velmeshev.pfc),"_pfc")
colnames(ts.velmeshev.acc) <- paste0(colnames(ts.velmeshev.acc),"_acc")

ts.velmeshev.full <- cbind(ts.velmeshev.pfc, ts.velmeshev.acc)

cor_t_dlpfc.asd <- cor(ts.dlpfc, ts.velmeshev.full)
range(cor_t_dlpfc.asd)


## Heatmap
# Add some cluster info for add'l heatmap annotations
regionInfo <- data.frame(region=ss(colnames(ts.velmeshev.full), "_",2))
rownames(regionInfo) <- colnames(ts.velmeshev.full)


# Print
theSeq.all = seq(-.95, .95, by = 0.01)
my.col.all <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.all)-1)

pdf("pdfs/exploration/Velmeshev-ASD_pfc-acc/overlap-velmeshev-ASD-bothRegions_with_LIBD-10x-DLPFC_Aug2020.pdf", width=10)
pheatmap(cor_t_dlpfc.asd,
         color=my.col.all,
         annotation_col=regionInfo,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.all,
         fontsize=10, fontsize_row=11, fontsize_col=10,
         display_numbers=T, number_format="%.2f", fontsize_number=4.5,
         legend_breaks=c(seq(-0.95,0.95,by=0.475)),
         main="Correlation of cluster-specific t's between LIBD DLPFC to \n ACC & PFC from (Velmeshev et al. Science 2019)")
dev.off()





    
    
