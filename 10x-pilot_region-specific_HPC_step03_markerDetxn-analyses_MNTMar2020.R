### MNT 10x snRNA-seq workflow: step 03 - marker detection
###   **Region-specific analyses**
###     - (3x) HPC samples from: Br5161 & Br5212 & Br5287
### Initiated MNT 13Mar2020
#####################################################################

library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(batchelor)
library(DropletUtils)
library(jaffelab)
library(limma)

### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

# ===


### Cell type marker gene detection =======================================
#   ** Approach - pseudo-bulk on sample:cellType stratification, then treat as SCE, so that
#                 can use 'findMarkers()' function that utilizes different tests

load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    ## sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo

    # sce.hpc has 10444 nuclei

# First drop "Ambig.lowNtrxts" (101 nuclei)
sce.hpc <- sce.hpc[ ,sce.hpc$cellType != "Ambig.lowNtrxts"]
sce.hpc$cellType <- droplevels(sce.hpc$cellType)

# Then make the pseudo-bulked SCE
sce.hpc.PB <- aggregateAcrossCells(sce.hpc, ids=paste0(sce.hpc$sample,":",sce.hpc$cellType),
                                   use_exprs_values="counts")

# Drop genes with all 0's
sce.hpc.PB <- sce.hpc.PB[!rowSums(assay(sce.hpc.PB, "counts"))==0, ]
    ## keeps 28757 genes

# Remove stored `sizeFactors()` because this will mess you up
#     * Also, to be safe, can always provide manually-computed SFs:
sizeFactors(sce.hpc.PB) <- NULL
LSFvec <- librarySizeFactors(sce.hpc.PB)
sce.hpc.PB <- logNormCounts(sce.hpc.PB, size_factors=LSFvec)



## Find markers using stringent [max-p-value-of-all-pw-comparisons] test ('pval.type="all"') ===
# No 'design'
markers.hpc.t <- findMarkers(sce.hpc.PB, groups=sce.hpc.PB$cellType,
                                     assay.type="logcounts",
                                     direction="up", pval.type="all", full.stats=T)

sapply(markers.hpc.t, function(x){table(x$FDR<0.05)})
    # none


# With 'design='? (and go ahead and use normalized counts--"countsNormd")
design.PB <- model.matrix(~sce.hpc.PB$processDate)
design.PB <- design.PB[ , -1, drop=F] # 'drop=F' to keep as matrix - otherwise turns into numeric vec

# "logcounts"
markers.hpc.t.design.log <- findMarkers(sce.hpc.PB, groups=sce.hpc.PB$cellType,
                                        assay.type="logcounts", design=design.PB,
                                        direction="up", pval.type="all", full.stats=T)

sapply(markers.hpc.t.design.log, function(x){table(x$FDR<0.05)})
    ##      Ambig.glial Astro Excit Inhib Micro Oligo   OPC
    # FALSE       28281 28275 28605 28631 28014 28429 28628
    # TRUE          476   482   152   126   743   328   129


# "countsNormd"
assay(sce.hpc.PB, "countsNormd") <- t(apply(assay(sce.hpc.PB, "counts"), 1,
                                                   function(x) {x/LSFvec}))


markers.hpc.t.design.countsN <- findMarkers(sce.hpc.PB, groups=sce.hpc.PB$cellType,
                                            assay.type="countsNormd", design=design.PB,
                                            direction="up", pval.type="all", full.stats=T)

sapply(markers.hpc.t.design.countsN, function(x){table(x$FDR<0.05)})
    ##      Ambig.glial Astro Excit Inhib Micro Oligo   OPC
    # FALSE       28257 27588 28388 28284 27323 27815 28324
    # TRUE          500  1169   369   473  1434   942   433


markerList.PB.tDesign.hpc.countsN <- lapply(markers.hpc.t.design.countsN, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)
sum(lengths(markerList.PB.tDesign.hpc.countsN))
    ## 5320

markerList.PB.tDesign.hpc.log <- lapply(markers.hpc.t.design.log, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)

sum(lengths(markerList.PB.tDesign.hpc.log))
    ## 2436

length(intersect(unlist(markerList.PB.tDesign.hpc.log), unlist(markerList.PB.tDesign.hpc.countsN)))
    ## 2284


save(markers.hpc.t.design.log, markers.hpc.t.design.countsN,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_HPC-n3_findMarkers_MNTMar2020.rda")




###### Direct limma approach ####
#################################

## Load in sce.dlpfc and cluster:sample bulk
 # (already done up top)

# Clean up colData
colData(sce.hpc.PB) <- colData(sce.hpc.PB)[ ,c(13:17,19:20)]


## Extract the count data
mat <- assays(sce.hpc.PB)$logcounts

## Build a group model
mod <- with(colData(sce.hpc.PB), model.matrix(~ 0 + cellType))
colnames(mod) <- gsub('cellType', '', colnames(mod))

corfit <- duplicateCorrelation(mat, mod, block = sce.hpc.PB$donor)
corfit$consensus.correlation
    # [1] 0.03102348

fit <-
  lmFit(
    mat,
    design = mod,
    block = sce.hpc.PB$donor,
    correlation = corfit$consensus.correlation
  )
eb <- eBayes(fit)


## Contrasts for pairwise comparison
cellType_combs <- combn(colnames(mod), 2) # will be `choose(ncol(x), 2)` columns long, of course
cellType_contrasts <- apply(cellType_combs, 2, function(x) {
  z <- paste(x, collapse = '-')
  makeContrasts(contrasts = z, levels = mod)
})
rownames(cellType_contrasts) <- colnames(mod)
colnames(cellType_contrasts) <- apply(cellType_combs, 2, paste, collapse = '-')

eb_contrasts <- eBayes(contrasts.fit(fit, cellType_contrasts))

## Tabulating significant hits
pvals_contrasts <- eb_contrasts$p.value

data.frame(
  'FDRsig' = colSums(apply(pvals_contrasts, 2, p.adjust, 'fdr') < 0.05),
  'Pval10-6sig' = colSums(pvals_contrasts < 1e-6),
  'Pval10-8sig' = colSums(pvals_contrasts < 1e-8)
)
    #                   FDRsig Pval10.6sig Pval10.8sig
    # Ambig.glial-Astro  14189        5765        3949
    # Ambig.glial-Excit  15080        5712        3859
    # Ambig.glial-Inhib  14623        5890        3974
    # Ambig.glial-Micro  11800        4866        3279
    # Ambig.glial-Oligo  12518        4813        3268
    # Ambig.glial-OPC    13347        5529        3769
    # Astro-Excit         5860        1155         424
    # Astro-Inhib         5357        1048         399
    # Astro-Micro         7116        1755         788
    # Astro-Oligo         5275        1255         488
    # Astro-OPC           3899         720         243
    # Excit-Inhib         1417         193          57
    # Excit-Micro         8489        1995         919
    # Excit-Oligo         6448        1422         599
    # Excit-OPC           4654         785         271
    # Inhib-Micro         8105        1987         849
    # Inhib-Oligo         6127        1417         611
    # Inhib-OPC           3731         606         202
    # Micro-Oligo         5709        1353         606
    # Micro-OPC           6702        1679         766
    # Oligo-OPC           4269        1034         407



## Then each cellType vs the rest
cellType_idx <- splitit(sce.hpc.PB$cellType)

eb0_list <- lapply(cellType_idx, function(x) {
  res <- rep(0, ncol(sce.hpc.PB))
  res[x] <- 1
  m <- model.matrix(~ res)
  eBayes(
    lmFit(
      mat,
      design = m,
      block = sce.hpc.PB$donor,
      correlation = corfit$consensus.correlation
    )
  )
})

## Extract the p-values
pvals0_contrasts <- sapply(eb0_list, function(x) {
  x$p.value[, 2, drop = FALSE]
})

## Extract the tstats
t0_contrasts_cell <- sapply(eb0_list, function(x) {
  x$t[, 2, drop = FALSE]
})

# Any signif
data.frame(
  'FDRsig' = colSums(apply(pvals0_contrasts, 2, p.adjust, 'fdr') < 0.05),
  'Pval10-6sig' = colSums(pvals0_contrasts < 1e-6),
  'Pval10-8sig' = colSums(pvals0_contrasts < 1e-8)
)
    ##             FDRsig Pval10.6sig Pval10.8sig
    #  Ambig.glial  13699        3556        2416
    #  Astro          984         203          92
    #  Excit          638          91          39
    #  Inhib          429          93          48
    #  Micro         1380         383         229
    #  Oligo          491         108          34
    #  OPC            243          70          34

# For only (+) t-stats
data.frame(
  'FDRsig' = colSums(apply(pvals0_contrasts, 2, p.adjust, 'fdr') < 0.05 &
                       t0_contrasts_cell > 0),
  'Pval10-6sig' = colSums(pvals0_contrasts < 1e-6 &
                            t0_contrasts_cell > 0),
  'Pval10-8sig' = colSums(pvals0_contrasts < 1e-8 &
                            t0_contrasts_cell > 0)
)

    #            FDRsig Pval10.6sig Pval10.8sig
    #Ambig.glial   1598         174          78
    #Astro          955         202          92
    #Excit          597          89          39
    #Inhib          424          93          48
    #Micro         1318         377         225
    #Oligo          395          81          22
    #OPC            242          70          34




## Save for later
eb_contrasts.hpc.broad <- eb_contrasts
eb_list.hpc.broad <- eb0_list

save(eb_contrasts.hpc.broad, eb_list.hpc.broad, sce.hpc.PB,
     file = 'rdas/markers-stats_HPC-n3_manualContrasts_MNTMar2020.rda')




### Further explore - what are some of these "Ambig.glial" markers? === === === ===
  # -> print expression - do any of them look real?
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo
load("rdas/markers-stats_HPC-n3_manualContrasts_MNTMar2020.rda", verbose=T)
    # eb_contrasts.hpc.broad, eb_list.hpc.broad, sce.hpc.PB,

## Extract the p-values
pvals0_contrasts <- sapply(eb_list.hpc.broad, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts) <- rownames(sce.hpc.PB)

## Extract the tstats
t0_contrasts_cell <- sapply(eb_list.hpc.broad, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts_cell) <- rownames(sce.hpc.PB)

table(pvals0_contrasts[ ,"Ambig.glial"] <= 1e-8 & t0_contrasts_cell[ ,"Ambig.glial"] > 0)

potentialGenes <- (pvals0_contrasts[ ,"Ambig.glial"] <= 1e-8 & t0_contrasts_cell[ ,"Ambig.glial"] > 0)

quantile(t0_contrasts_cell[ ,"Ambig.glial"][potentialGenes])
    #0%       25%       50%       75%      100%
    #9.146918 10.183891 11.362352 16.448760 35.878306

upperQrtls <- potentialGenes & t0_contrasts_cell[ ,"Ambig.glial"] >= 11.362352

genesToPrint <- rownames(t0_contrasts_cell)[upperQrtls]
    ##  [1] "LCK"       "CD2"       "PYHIN1"    "FCRL6"     "SLAMF6"    "SLAMF1"
    #[7] "CD48"      "SLAMF7"    "LINC01871" "CD8B"      "OXNAD1"    "TRAT1"
    #[13] "NAA50"     "TIGIT"     "IL7R"      "GZMK"      "GZMA"      "ITK"
    #[19] "DOK2"      "PDE7A"     "LINC00861" "CD3E"      "CD3D"      "CD3G"
    #[25] "TRAC"      "LAT"       "ACAP1"     "SLFN12L"   "CCL5"      "CCL4"
    #[31] "IKZF3"     "TBX21"     "SKAP1"     "SKAP1-AS1" "SLA2"      "CNN2"
    #[37] "IL2RB"     "GRAP2"     "UBASH3A"


# First drop "Ambig.lowNtrxts" (101 nuclei)
sce.hpc <- sce.hpc[ ,sce.hpc$cellType != "Ambig.lowNtrxts"]
sce.hpc$cellType <- droplevels(sce.hpc$cellType)


pdf("pdfs/zExplore_HPC_ambig-glialCluster_markerExpression_MNTMar2020.pdf", height=3, width=4)
for(i in 1:length(genesToPrint)){  
  print(
    plotExpression(sce.hpc, exprs_values = "logcounts", features=genesToPrint[i],
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=tableau10medium[1:7])
        )
}
dev.off()

    ## Looks pretty real.  And these are very obviously T-cells...

pvals0_contrasts[ ,"Ambig.glial"][upperQrtls]
max(pvals0_contrasts[ ,"Ambig.glial"][upperQrtls]) # 2.12e-10
sum(pvals0_contrasts[ ,"Ambig.glial"] <= 2.2e-10 & t0_contrasts_cell[ ,"Ambig.glial"] > 0)
    #[1] 40

    ## So these are some of the best pvals too, even though pvals & t's are not necessarily super
     #      correlated
     # sapply(1:7, function(x){cor(pvals0_contrasts[ ,x], t0_contrasts_cell[ ,x])})


## For reference === === ===
table(sce.hpc$cellType, sce.hpc$sample)
    #                 hpc.5161 hpc.5212 hpc.5287
    # Ambig.glial           13       10        3
    # Ambig.lowNtrxts       48       39       14
    # Astro                584      582      177
    # Excit                158      319      122
    # Inhib                206       92       75
    # Micro                520      536      197
    # Oligo               2594     2198     1093
    # OPC                  395      263      206


### How does this compare to results of `findMarkers()`? =======

## Extract the p-values and compute fdrs
pvals0_contrasts <- sapply(eb_list.hpc.broad, function(x) {
  x$p.value[, 2, drop = FALSE]
})

fdrs0_contrasts = apply(pvals0_contrasts, 2, p.adjust, "fdr")

## Extract the tstats
t0_contrasts <- sapply(eb_list.hpc.broad, function(x) {
  x$t[, 2, drop = FALSE]
})

rownames(fdrs0_contrasts) <- rownames(sce.hpc.PB)
rownames(t0_contrasts) <- rownames(sce.hpc.PB)

markerList.PB.manual <- lapply(colnames(fdrs0_contrasts), function(x){
  rownames(fdrs0_contrasts)[fdrs0_contrasts[ ,x] < 0.05 & t0_contrasts[ ,x] > 0]
  # what if more stringent?
  #rownames(fdrs0_contrasts)[fdrs0_contrasts[ ,x] < 0.01 & t0_contrasts[ ,x] > 0]
})
names(markerList.PB.manual) <- colnames(fdrs0_contrasts)


## findMarkers() results - test was already just for up-regulated genes
sapply(markers.hpc.t.design.log, function(x){table(x$FDR<0.05)})

markerList.PB.hpc.tDesign.log <- lapply(markers.hpc.t.design.log, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)


lengths(markerList.PB.manual)



lengths(markerList.PB.hpc.tDesign.log)


sapply(names(markerList.PB.manual), function(x){
  length(intersect(markerList.PB.manual[[x]],
                   markerList.PB.hpc.tDesign.log[[x]]))}
)





