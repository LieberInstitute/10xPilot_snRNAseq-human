### MNT 10x snRNA-seq workflow: step 03 - marker detection
###   **Region-specific analyses**
###     - (2x) NAc samples from: Br5161 & Br5212 & Br5287
### Initiated MNT 12Feb2020
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


### Cell type marker gene detection ===
#   ** Approach - pseudo-bulk on sample:cellType stratification, then treat as SCE, so that
#                 can use 'findMarkers()' function that utilizes different tests

load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_NAc-n3_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    ## sce.nac, chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac, ref.sampleInfo

# First drop "Ambig.lowNtrxts" (52 nuclei)
sce.nac <- sce.nac[ ,sce.nac$cellType != "Ambig.lowNtrxts"] # 4496 nuclei afterwards
sce.nac$cellType <- droplevels(sce.nac$cellType)

# Then make the pseudo-bulked SCE
sce.nac.PB <- aggregateAcrossCells(sce.nac, ids=paste0(sce.nac$sample,":",sce.nac$cellType),
                                     use_exprs_values="counts")

# Drop genes with all 0's
sce.nac.PB <- sce.nac.PB[!rowSums(assay(sce.nac.PB, "counts"))==0, ]
## keeps 27618 genes

# Clean up colData
colData(sce.nac.PB) <- colData(sce.nac.PB)[ ,c(13:17,19:20)]

# Remove stored `sizeFactors()` because this will mess you up
#     * Also, to be safe, can always provide manually-computed SFs:
sizeFactors(sce.nac.PB) <- NULL
LSFvec <- librarySizeFactors(sce.nac.PB)
sce.nac.PB <- logNormCounts(sce.nac.PB, size_factors=LSFvec)


## Find markers using stringent [max-p-value-of-all-pw-comparisons] test ('pval.type="all"')
 #    - first try no 'design=' argument
markers.nac.t <- findMarkers(sce.nac.PB, groups=sce.nac.PB$cellType,
                               assay.type="logcounts",
                               direction="up", pval.type="all", full.stats=T)

sapply(markers.nac.t, function(x){table(x$FDR<0.05)})
    # none


## With 'design=' ? ===

# For reference:
table(sce.nac.PB$processDate, sce.nac.PB$sample)
    #         nac.5161 nac.5212 nac.5287
    #R2.Jul23        6        6        0
    #R3.Sep04        0        0        6

#design.PB <- model.matrix(~sce.nac.PB$sample + sce.nac.PB$processDate)
    ## This doesn't work with error:
     #    "Error in .ranksafe_qr(full.design) : design matrix is not of full rank"

## Add new combined design factor: sample:processDate
sce.nac.PB$sampleBatch <- paste0(sce.nac.PB$sample, ":", sce.nac.PB$processDate)
design.PB <- model.matrix(~sce.nac.PB$sampleBatch)
    ## Actually the processDate effect ends up being sucked up by sample...

design.PB <- design.PB[ , -1, drop=F] # 'drop=F' to keep as matrix (if one factor)
                                      #         - otherwise turns into numeric vec

markers.nac.t.design.counts <- findMarkers(sce.nac.PB, groups=sce.nac.PB$cellType,
                                             assay.type="counts", design=design.PB,
                                             direction="up", pval.type="all", full.stats=T)

sapply(markers.nac.t.design.counts, function(x){table(x$FDR<0.05)})
    ##      Astro Inhib.NRGNneg Inhib.NRGNpos Micro Oligo   OPC
    # FALSE 27617         27611         20691 27425 26715 27571
    # TRUE      1             7          6927   193   903    47



## normalized and transformed "logcounts"
markers.nac.t.design.log <- findMarkers(sce.nac.PB, groups=sce.nac.PB$cellType,
                                          assay.type="logcounts", design=design.PB,
                                          direction="up", pval.type="all", full.stats=T)


sapply(markers.nac.t.design.log, function(x){table(x$FDR<0.05)})
    ##      Astro Inhib.NRGNneg Inhib.NRGNpos Micro Oligo   OPC
    # FALSE 27395         27523         27068 26468 27292 27522
    # TRUE    223            95           550  1150   326    96

markerList.PB.tDesign.nac.log <- lapply(markers.nac.t.design.log, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)


## Normalized $counts, but not log-transformed
assay(sce.nac.PB, "countsNormd") <- t(apply(assay(sce.nac.PB, "counts"), 1, function(x) {x/LSFvec}))


markers.nac.t.design.countsN <- findMarkers(sce.nac.PB, groups=sce.nac.PB$cellType,
                                              assay.type="countsNormd", design=design.PB,
                                              direction="up", pval.type="all", full.stats=T)

sapply(markers.nac.t.design.countsN, function(x){table(x$FDR<0.05)})
    ##      Astro Inhib.NRGNneg Inhib.NRGNpos Micro Oligo   OPC
    # FALSE 26932         27237         26486 25658 26749 27236
    # TRUE    686           381          1132  1960   869   382

markerList.PB.tDesign.nac.countsN <- lapply(markers.nac.t.design.countsN, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)

# Between "logcounts" & "countsNormd"    
intersectingMarkers <- list()
for(i in names(markerList.PB.tDesign.nac.log)){
  intersectingMarkers[[i]] <- intersect(markerList.PB.tDesign.nac.countsN[[i]],
                                        markerList.PB.tDesign.nac.log[[i]])
}
lengths(intersectingMarkers)
    ## so most of those detected in "logcounts" are included in "countsNormd" (not surprising)


        ## sample:processDate in this still just breaks into three groups:
            #design.PB <- model.matrix(~sce.nac.PB$sample)
            #design.PB <- design.PB[ , -1, drop=F]            - this yields same result

        ## When try to block on 'processDate':
            # "Error in FUN(x, groups, ..., log.p = TRUE) : [\n] cannot specify both 'block' and 'design'"
        
        

save(markers.nac.t.design.log, markers.nac.t.design.countsN,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_NAc-n3-homs-only_findMarkers_MNTMar2020.rda")




###########################
### direct limma method ###
library(limma)

## Load in sce.dlpfc and cluster:sample bulk
# (already done up top)

## Extract the count data
mat <- assays(sce.nac.PB)$logcounts

## Build a group model
#mod <- with(colData(sce.nac.PB), model.matrix(~ 0 + cellType))
#colnames(mod) <- gsub('cellType', '', colnames(mod))

        #corfit.processDate <- duplicateCorrelation(mat, mod, block = sce.nac.PB$processDate)
        #corfit.processDate$consensus.correlation
        #    # [1] 0.005682906
        #
        #    ## With model adding 'sample'
        #    mod <- with(colData(sce.nac.PB), model.matrix(~ 0 + cellType + sample))
        #    corfit.processDate <- duplicateCorrelation(mat, mod, block = sce.nac.PB$processDate)
        #    corfit.processDate$consensus.correlation
        #        ## [1] 1    - yeaahh this isn't right - all $atanh.correlations == Inf......

#corfit <- duplicateCorrelation(mat, mod, block = sce.nac.PB$donor)
#corfit$consensus.correlation
    # [1] 0.05230507        - should always be on $sample as per AnJa (doesn't affect DLPFC/Amyg)


## With model adding 'processDate'
mod <- with(colData(sce.nac.PB), model.matrix(~ 0 + cellType + processDate))
colnames(mod) <- gsub('cellType', '', colnames(mod))
corfit <- duplicateCorrelation(mat, mod, block = sce.nac.PB$donor)
corfit$consensus.correlation
    ## 0.003670811



fit <-
  lmFit(
    mat,
    design = mod,
    block = sce.nac.PB$donor,
    correlation = corfit$consensus.correlation
  )
eb <- eBayes(fit)


## Contrasts for pairwise comparison
cellType_combs <- combn(colnames(mod)[1:6], 2)
cellType_contrasts <- apply(cellType_combs, 2, function(x) {
  z <- paste(x, collapse = '-')
  makeContrasts(contrasts = z, levels = mod[ ,1:6])
})
rownames(cellType_contrasts) <- colnames(mod)[1:6]
colnames(cellType_contrasts) <- apply(cellType_combs, 2, paste, collapse = '-')

fit.sub <- fit[ ,1:6]

eb_contrasts <- eBayes(contrasts.fit(fit.sub, cellType_contrasts))

## Tabulating significant hits
pvals_contrasts <- eb_contrasts$p.value

data.frame(
  'FDRsig' = colSums(apply(pvals_contrasts, 2, p.adjust, 'fdr') < 0.05),
  'Pval10-6sig' = colSums(pvals_contrasts < 1e-6),
  'Pval10-8sig' = colSums(pvals_contrasts < 1e-8)
)




## Then each cellType vs the rest
cellType_idx <- splitit(sce.nac.PB$cellType)

eb0_list <- lapply(cellType_idx, function(x) {
  res <- rep(0, ncol(sce.nac.PB))
  res[x] <- 1
  m <- model.matrix(~ res + processDate, data=colData(sce.nac.PB))
  eBayes(
    lmFit(
      mat,
      design = m,
      block = sce.nac.PB$donor,
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


data.frame(
  'FDRsig' = colSums(apply(pvals0_contrasts, 2, p.adjust, 'fdr') < 0.05),
  'Pval10-6sig' = colSums(pvals0_contrasts < 1e-6),
  'Pval10-8sig' = colSums(pvals0_contrasts < 1e-8)
)

# For only (+) t-stats
data.frame(
  'FDRsig' = colSums(apply(pvals0_contrasts, 2, p.adjust, 'fdr') < 0.05 &
                       t0_contrasts_cell > 0),
  'Pval10-6sig' = colSums(pvals0_contrasts < 1e-6 &
                            t0_contrasts_cell > 0),
  'Pval10-8sig' = colSums(pvals0_contrasts < 1e-8 &
                            t0_contrasts_cell > 0)
)

    ## Without t > 0 subset:
        ##              FDRsig Pval10.6sig Pval10.8sig
        # Astro            771         120          48
        # Inhib.NRGNneg    428          70          34
        # Inhib.NRGNpos   1953         276         111
        # Micro           4388         687         308
        # Oligo           1104         183          59
        # OPC              275          63          31
    
    ## With t > 0
        ##              FDRsig Pval10.6sig Pval10.8sig
        # Astro            609         112          46
        # Inhib.NRGNneg    293          61          33
        # Inhib.NRGNpos   1590         252         106
        # Micro           2406         471         237
        # Oligo            753         126          44
        # OPC              243          62          31




## Save for later
eb_contrasts.nac.broad <- eb_contrasts
eb_list.nac.broad <- eb0_list

save(eb_contrasts.nac.broad, eb_list.nac.broad, sce.nac.PB,
     file = 'rdas/markers-stats_NAc-n3-homs-only_manualContrasts_MNTMar2020.rda')



### MNT 11Mar2020: How does this compare to results of `findMarkers()`? === === === === ===

## Extract the p-values and compute fdrs
pvals0_contrasts <- sapply(eb_list.nac.broad, function(x) {
  x$p.value[, 2, drop = FALSE]
})

fdrs0_contrasts = apply(pvals0_contrasts, 2, p.adjust, "fdr")

## Extract the tstats
t0_contrasts <- sapply(eb_list.nac.broad, function(x) {
  x$t[, 2, drop = FALSE]
})

rownames(fdrs0_contrasts) <- rownames(sce.nac.PB)
rownames(t0_contrasts) <- rownames(sce.nac.PB)

markerList.PB.manual <- lapply(colnames(fdrs0_contrasts), function(x){
  rownames(fdrs0_contrasts)[fdrs0_contrasts[ ,x] < 0.05 & t0_contrasts[ ,x] > 0]
  # what if more stringent?
  #rownames(fdrs0_contrasts)[fdrs0_contrasts[ ,x] < 0.01 & t0_contrasts[ ,x] > 0]
})
names(markerList.PB.manual) <- colnames(fdrs0_contrasts)


## findMarkers() results - test was already just for up-regulated genes
sapply(markers.nac.t.design.log, function(x){table(x$FDR<0.05)})

markerList.PB.nac.tDesign.log <- lapply(markers.nac.t.design.log, function(x){
  rownames(x)[x$FDR < 0.05]
  })


lengths(markerList.PB.manual)
#   Astro Inhib.NRGNneg Inhib.NRGNpos      Micro     Oligo    OPC
#     609           293          1590       2406       753    243


lengths(markerList.PB.nac.tDesign.log)
#     Astro Inhib.NRGNneg Inhib.NRGNpos       Micro    Oligo    OPC
#       223            95           550        1150      326     96



sapply(names(markerList.PB.manual), function(x){
  length(intersect(markerList.PB.manual[[x]],
                   markerList.PB.nac.tDesign.log[[x]]))}
)
##    Astro Inhib.NRGNneg Inhib.NRGNpos       Micro     Oligo     OPC
#       207            91           537        1134       301      87




