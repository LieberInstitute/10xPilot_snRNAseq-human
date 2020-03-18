### MNT 10x snRNA-seq workflow: step 03 - marker detection
###   **Region-specific analyses**
###     - (2x) amygdala samples from: Br5161 & Br5212
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


### Cell type marker gene detection =======================================
#   ** Approach - pseudo-bulk on sample:cellType stratification, then treat as SCE, so that
#                 can use 'findMarkers()' function that utilizes different tests

load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    ## sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo

# First drop "Ambig.lowNtrxts" (50 nuclei)
sce.amy <- sce.amy[ ,sce.amy$cellType != "Ambig.lowNtrxts"]
sce.amy$cellType <- droplevels(sce.amy$cellType)

# Then make the pseudo-bulked SCE
sce.amy.PB <- aggregateAcrossCells(sce.amy, ids=paste0(sce.amy$sample,":",sce.amy$cellType),
                                   use_exprs_values="counts")

# Drop genes with all 0's
sce.amy.PB <- sce.amy.PB[!rowSums(assay(sce.amy.PB, "counts"))==0, ]
    ## keeps 28470 genes; 28464 if drop "Ambig.lowNtrxts" nuclei first

# Remove stored `sizeFactors()` because this will mess you up
#     * Also, to be safe, can always provide manually-computed SFs:
sizeFactors(sce.amy.PB) <- NULL
LSFvec <- librarySizeFactors(sce.amy.PB)
sce.amy.PB <- logNormCounts(sce.amy.PB, size_factors=LSFvec)

## Find markers using stringent [max-p-value-of-all-pw-comparisons] test ('pval.type="all"')
## SKIP THIS - MNT 05Mar2020 ================================
 #    - first try no 'design=' argument
markers.amy.t <- findMarkers(sce.amy.PB, groups=sce.amy.PB$cellType,
                             assay.type="logcounts",
                             direction="up", pval.type="all", full.stats=T)

sapply(markers.amy.t, function(x){table(x$FDR<0.05)})
    # none


## With 'design=' ? *** MUST DO ON $counts, NOT NORMALIZED COUNTS
design.PB <- model.matrix(~sce.amy.PB$sample)
design.PB <- design.PB[ , -1, drop=F] # 'drop=F' to keep as matrix - otherwise turns into numeric vec

markers.amy.t.design.counts <- findMarkers(sce.amy.PB, groups=sce.amy.PB$cellType,
                                           assay.type="counts", design=design.PB,
                                           direction="up", pval.type="all", full.stats=T)

sapply(markers.amy.t.design.counts, function(x){table(x$FDR<0.05)[2]})
    ##Ambig.lowNtrxts.NA         Astro.TRUE         Excit.TRUE         Inhib.TRUE
    #                NA                470               1231                 10
    #        Micro.TRUE         Oligo.TRUE           OPC.TRUE
    #               528               1903                 78



## normalized and transformed "logcounts"
markers.amy.t.design.log <- findMarkers(sce.amy.PB, groups=sce.amy.PB$cellType,
                                        assay.type="logcounts", design=design.PB,
                                        direction="up", pval.type="all", full.stats=T)

sapply(markers.amy.t.design.log, function(x){table(x$FDR<0.05)})
    ##      Ambig.lowNtrxts Astro Excit Inhib Micro Oligo   OPC
    #FALSE           28381 28241 28223 28427 27375 28229 28403
    #TRUE               89   229   247    43  1095   241    67

markerList.PB.tDesign.amy.log <- lapply(markers.amy.t.design.log, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)


## Normalized $counts, but not log-transformed
assay(sce.amy.PB, "countsNormd") <- t(apply(assay(sce.amy.PB, "counts"), 1, function(x) {x/LSFvec}))


markers.amy.t.design.countsN <- findMarkers(sce.amy.PB, groups=sce.amy.PB$cellType,
                                            assay.type="countsNormd", design=design.PB,
                                            direction="up", pval.type="all", full.stats=T)

sapply(markers.amy.t.design.countsN, function(x){table(x$FDR<0.05)})
    ##      Ambig.lowNtrxts Astro Excit Inhib Micro Oligo   OPC
    #FALSE           28393 27685 27803 28066 26531 27726 28138
    #TRUE               77   785   667   404  1939   744   332

markerList.PB.tDesign.amy.countsN <- lapply(markers.amy.t.design.countsN, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)

# Between "logcounts" & "countsNormd"    
intersectingMarkers <- list()
for(i in names(markerList.PB.tDesign.amy.log)){
  intersectingMarkers[[i]] <- intersect(markerList.PB.tDesign.amy.countsN[[i]],
                                        markerList.PB.tDesign.amy.log[[i]])
}
lengths(intersectingMarkers)
## so most of those detected in "logcounts" are included in "countsNormd" (not surprising)


      ## Using 'block=' is N/A in this PB'd setting:
      #    -> warnings such as "no within-block comparison between ____ and ____"



### end skip - with no Ambig cluster ===============================
#sce.amy.PBnoAmbig <- sce.amy.PB[ ,!sce.amy.PB$cellType=="Ambig.lowNtrxts"]

    ## MNT comment - replacing the below 'sce.dlpfc.PBnoAmbig', with just 'sce.dlpfc.PB',
    #     to make more follow-able

#assay(sce.amy.PBnoAmbig, "logcounts") <- NULL
#LSFvec.noAmbig <- librarySizeFactors(sce.amy.PBnoAmbig)
#sce.amy.PBnoAmbig <- logNormCounts(sce.amy.PBnoAmbig,
#                                   size_factors=LSFvec.noAmbig)

## Remove that level too
#sce.amy.PBnoAmbig$cellType <- factor(sce.amy.PBnoAmbig$cellType,
#                                     levels=unique(sce.amy.PBnoAmbig$cellType))

markers.amy.t.noAmbig <- findMarkers(sce.amy.PB, groups=sce.amy.PB$cellType,
                                     assay.type="logcounts",
                                     direction="up", pval.type="all", full.stats=T)

sapply(markers.amy.t.noAmbig, function(x){table(x$FDR<0.05)})
    # still none


# With 'design='? (and go ahead and use normalized counts--"countsNormd")
design.PB.noAmbig <- model.matrix(~sce.amy.PB$processDate)
design.PB.noAmbig <- design.PB.noAmbig[ , -1, drop=F] # 'drop=F' to keep as matrix - otherwise turns into numeric vec

# "logcounts"
#markers.amy.t.design.noAmbig.log <- findMarkers(sce.amy.PB, groups=sce.amy.PB$cellType,
markers.amy.t.design.log <- findMarkers(sce.amy.PB, groups=sce.amy.PB$cellType,
                                        assay.type="logcounts", design=design.PB.noAmbig,
                                        direction="up", pval.type="all", full.stats=T)

sapply(markers.amy.t.design.noAmbig.log, function(x){table(x$FDR<0.05)})
    ##      Astro Excit Inhib Micro Oligo   OPC
    # FALSE 28198 28261 28450 26407 28009 28425
    # TRUE    272   209    20  2063   461    45   - when first had rm'd all-0 genes (old)

    ##      Astro Excit Inhib Micro Oligo   OPC
    # FALSE 28192 28255 28444 26401 28001 28419
    # TRUE    272   209    20  2063   463    45   - when first dropping "Ambig.lowNtrxts" (almost same result)


# "countsNormd" - need to re-normalize after dropping 'Ambig.lowNtrxts'
#assay(sce.amy.PB, "countsNormd") <- NULL
assay(sce.amy.PB, "countsNormd") <- t(apply(assay(sce.amy.PB, "counts"), 1,
                                                   #function(x) {x/LSFvec.noAmbig}))
                                                   function(x) {x/LSFvec})) # from up top - MNT 05Mar2020


#markers.amy.t.design.noAmbig.countsN <- findMarkers(sce.amy.PB, groups=sce.amy.PB$cellType,
markers.amy.t.design.countsN <- findMarkers(sce.amy.PB, groups=sce.amy.PB$cellType,
                                            assay.type="countsNormd", design=design.PB.noAmbig,
                                            direction="up", pval.type="all", full.stats=T)

sapply(markers.amy.t.design.noAmbig.countsN, function(x){table(x$FDR<0.05)})
    ##      Astro Excit Inhib Micro Oligo   OPC
    # FALSE 27406 27792 28026 25476 27414 28116
    # TRUE   1064   678   444  2994  1056   354   - when first had rm'd all-0 genes (old)

    ##      Astro Excit Inhib Micro Oligo   OPC
    # FALSE 27400 27786 28020 25470 27408 28110
    # TRUE   1064   678   444  2994  1056   354   - when first dropping "Ambig.lowNtrxts" (same result)


markerList.PB.tDesign.amy.noAmbig.countsN <- lapply(markers.amy.t.design.noAmbig.countsN, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)


sum(lengths(markerList.PB.tDesign.amy.countsN))
    ## 4948
sum(lengths(markerList.PB.tDesign.amy.noAmbig.countsN))
    ## 6590
length(intersect(unlist(markerList.PB.tDesign.amy.countsN), unlist(markerList.PB.tDesign.amy.noAmbig.countsN)))
    ## 4440


save(markers.amy.t.design.log, markers.amy.t.design.countsN,
#     markers.amy.t.design.noAmbig.log, markers.amy.t.design.noAmbig.countsN,
#     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_Amyg_n2_MNTFeb2020.rda")
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_Amyg-n2_findMarkers_MNTMar2020.rda")




###### Direct limma approach ####
#################################

## Load in sce.dlpfc and cluster:sample bulk
 # (already done up top)

# Clean up colData
colData(sce.amy.PB) <- colData(sce.amy.PB)[ ,c(13:17,19:20)]


## Extract the count data
mat <- assays(sce.amy.PB)$logcounts

## Build a group model
mod <- with(colData(sce.amy.PB), model.matrix(~ 0 + cellType))
colnames(mod) <- gsub('cellType', '', colnames(mod))

corfit <- duplicateCorrelation(mat, mod, block = sce.amy.PB$donor)
corfit$consensus.correlation
    # [1] 0.01202754

fit <-
  lmFit(
    mat,
    design = mod,
    block = sce.amy.PB$donor,
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
#             FDRsig Pval10.6sig Pval10.8sig
#Astro-Excit   8182         498          38
#Astro-Inhib   7443         441          32
#Astro-Micro  10503         895         112
#Astro-Oligo   6991         562          58
#Astro-OPC     5197         299          25
#Excit-Inhib   1808          51           1
#Excit-Micro  12477        1160         126
#Excit-Oligo   9305         739          59
#Excit-OPC     6736         264           9
#Inhib-Micro  11356        1078         108
#Inhib-Oligo   8215         668          45
#Inhib-OPC     5038         196           7
#Micro-Oligo   8878         980         128
#Micro-OPC    10207         951         115
#Oligo-OPC     6159         446          4



## Then each cellType vs the rest
cellType_idx <- splitit(sce.amy.PB$cellType)

eb0_list <- lapply(cellType_idx, function(x) {
  res <- rep(0, ncol(sce.amy.PB))
  res[x] <- 1
  m <- model.matrix(~ res)
  eBayes(
    lmFit(
      mat,
      design = m,
      block = sce.amy.PB$donor,
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
#FDRsig Pval10.6sig Pval10.8sig
#Astro    781          78          18
#Excit    631          67          19
#Inhib    155          28           4
#Micro   4860         562         127
#Oligo   1024         100           8
#OPC      154          27           5

## With t > 0
#      FDRsig Pval10.6sig Pval10.8sig
#Astro    697          76          18
#Excit    556          66          19
#Inhib    152          28           4
#Micro   2941         431         107
#Oligo    629          71           5
#OPC      139          26           5




## Save for later
eb_contrasts.amy.broad <- eb_contrasts
eb_list.amy.broad <- eb0_list

save(eb_contrasts.amy.broad, eb_list.amy.broad, sce.amy.PB,
     file = 'rdas/markers-stats_Amyg-n2_manualContrasts_MNTMar2020.rda')



### MNT 11Mar2020: How does this compare to results of `findMarkers()`? === === === === ===

## Extract the p-values and compute fdrs
pvals0_contrasts <- sapply(eb_list.amy.broad, function(x) {
  x$p.value[, 2, drop = FALSE]
})

fdrs0_contrasts = apply(pvals0_contrasts, 2, p.adjust, "fdr")

## Extract the tstats
t0_contrasts <- sapply(eb_list.amy.broad, function(x) {
  x$t[, 2, drop = FALSE]
})

rownames(fdrs0_contrasts) <- rownames(sce.amy.PB)
rownames(t0_contrasts) <- rownames(sce.amy.PB)

markerList.PB.manual <- lapply(colnames(fdrs0_contrasts), function(x){
  #rownames(fdrs0_contrasts)[fdrs0_contrasts[ ,x] < 0.05 & t0_contrasts[ ,x] > 0]
  # what if more stringent?
  rownames(fdrs0_contrasts)[fdrs0_contrasts[ ,x] < 0.01 & t0_contrasts[ ,x] > 0]
})
names(markerList.PB.manual) <- colnames(fdrs0_contrasts)


## findMarkers() results - test was already just for up-regulated genes
sapply(markers.amy.t.design.log, function(x){table(x$FDR<0.05)})

markerList.PB.amy.tDesign.log <- lapply(markers.amy.t.design.log, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)


lengths(markerList.PB.manual)
    # Astro Excit Inhib Micro Oligo   OPC
    #   697   556   152  2941   629   139
    
    ## With FDR<0.01:
    # Astro Excit Inhib Micro Oligo   OPC
    #   332   244    64  1761   332    79


lengths(markerList.PB.amy.tDesign.log)
    #Astro Excit Inhib Micro Oligo   OPC
    #  272   209    20  2063   463    45

sapply(names(markerList.PB.manual), function(x){
  length(intersect(markerList.PB.manual[[x]],
                   markerList.PB.amy.tDesign.log[[x]]))}
)
    ## Astro Excit Inhib Micro Oligo   OPC
    #    206   137    20  1857   336    35

    ## overlap with manual method's FDR<0.01 cutoff
    # Astro Excit Inhib Micro Oligo   OPC
    #   164   107    18  1430   233    29


### MNT aside, 05Mar2020 ===============
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)


# First drop "Ambig.lowNtrxts" (168 nuclei)
sce.amy.noAmbig <- sce.amy[ ,sce.amy$cellType != "Ambig.lowNtrxts"]
sce.amy.noAmbig$cellType <- droplevels(sce.amy.noAmbig$cellType)
# Then make the pseudo-bulked SCE
sce.amy.PBafterDrop <- aggregateAcrossCells(sce.amy.noAmbig, ids=paste0(sce.amy.noAmbig$sample,":",sce.amy.noAmbig$cellType),
                                              use_exprs_values="counts")
# Drop genes with all 0's
sce.amy.PBafterDrop <- sce.amy.PBafterDrop[!rowSums(assay(sce.amy.PBafterDrop, "counts"))==0, ]
    ## 28464 remaining genes

## OR

# Just make the pseudo-bulked SCE, without dropping that cluster
sce.amy.PB <- aggregateAcrossCells(sce.amy, ids=paste0(sce.amy$sample,":",sce.amy$cellType),
                                     use_exprs_values="counts")
# Drop genes with all 0's
sce.amy.PB <- sce.amy.PB[!rowSums(assay(sce.amy.PB, "counts"))==0, ]
    ## 28470

## Then
genesOfInterest <- setdiff(rownames(sce.amy.PB), rownames(sce.amy.PBafterDrop))
#  [1] "FOXC1"      "CSAG1"      "AP000879.1" "TBX3"       "FOXL1"        "FOXS1"

rowSums(assay(sce.amy.PB, "counts")[genesOfInterest, ])
## all 1's# 166 zeros

    ##    - so this is definitely just a poor, lowly-captured cluster; explains the poor
    ##      intra-cluster correlation, compared to other clusters



## And old...
# For BoG abstract ===
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_Amyg_n2_MNTFeb2020.rda", verbose=T)
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_DLPFC_n2_MNTFeb2020.rda", verbose=T)

intersectingMarkers <- list()
for(i in names(markerList.PB.tDesign.amy)){
  intersectingMarkers[[i]] <- intersect(markerList.PB.tDesign.amy[[i]],
                                        markerList.PB.tDesign.dlpfc[[i]])
}
# CEACAM21 a microglial marker - turns out to have been implicated in SCZ
#     (Jewish-Israeli familial study - I don't think in large GWAS though)

specificMarkers.amy <- list()
for(i in names(markerList.PB.tDesign.amy)){
  specificMarkers.amy[[i]] <- setdiff(markerList.PB.tDesign.amy[[i]],
                                      markerList.PB.tDesign.dlpfc[[i]])  
}

specificMarkers.dlpfc <- list()
for(i in names(markerList.PB.tDesign.amy)){
  specificMarkers.dlpfc[[i]] <- setdiff(markerList.PB.tDesign.dlpfc[[i]],
                                        markerList.PB.tDesign.amy[[i]])
}