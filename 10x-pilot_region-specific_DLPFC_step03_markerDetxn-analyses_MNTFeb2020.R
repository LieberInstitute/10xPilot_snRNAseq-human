### MNT 10x snRNA-seq workflow: step 03 - marker detection
###   **Region-specific analyses**
###     - (2x) DLPFC samples from: Br5161 & Br5212
### Initiated MNT 12Feb2020 - modified MNT 05Mar2020
### Modification notes: First drop "Ambig.lowNtrxts" cluster
###         (see bottom chunk for deets)
#####################################################################

library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(batchelor)
library(DropletUtils)
library(jaffelab)
library(limma)
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


### Cell type marker gene detection =======================================
#   ** Approach - pseudo-bulk on sample:cellType stratification, then treat as SCE, so that
#                 can use 'findMarkers()' function that utilizes different tests

load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_DLPFC-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    ## sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc, ref.sampleInfo

# First drop "Ambig.lowNtrxts" (168 nuclei)
sce.dlpfc <- sce.dlpfc[ ,sce.dlpfc$cellType != "Ambig.lowNtrxts"]
sce.dlpfc$cellType <- droplevels(sce.dlpfc$cellType)

# Then make the pseudo-bulked SCE
sce.dlpfc.PB <- aggregateAcrossCells(sce.dlpfc, ids=paste0(sce.dlpfc$sample,":",sce.dlpfc$cellType),
                                     use_exprs_values="counts")

# Drop genes with all 0's
sce.dlpfc.PB <- sce.dlpfc.PB[!rowSums(assay(sce.dlpfc.PB, "counts"))==0, ]
    ## keeps 28128 genes; 28111 if drop Ambig.lowNtrxts first

# Remove stored `sizeFactors()` because this will mess you up
#     * Also, to be safe, can always provide manually-computed SFs:
sizeFactors(sce.dlpfc.PB) <- NULL
LSFvec <- librarySizeFactors(sce.dlpfc.PB)
sce.dlpfc.PB <- logNormCounts(sce.dlpfc.PB, size_factors=LSFvec)

### Find markers using stringent [max-p-value-of-all-pw-comparisons] test ('pval.type="all"') === === ===
## SKIP THIS - MNT 05Mar2020 ================================
 #    - first try no 'design=' argument
markers.dlpfc.t <- findMarkers(sce.dlpfc.PB, groups=sce.dlpfc.PB$cellType,
                               assay.type="logcounts",
                               direction="up", pval.type="all", full.stats=T)

sapply(markers.dlpfc.t, function(x){table(x$FDR<0.05)})
    # none


## With 'design=' ? *** MUST DO ON $counts, NOT NORMALIZED COUNTS
design.PB <- model.matrix(~sce.dlpfc.PB$sample)
design.PB <- design.PB[ , -1, drop=F] # 'drop=F' to keep as matrix - otherwise turns into numeric vec

markers.dlpfc.t.design.counts <- findMarkers(sce.dlpfc.PB, groups=sce.dlpfc.PB$cellType,
                                             assay.type="counts", design=design.PB,
                                             direction="up", pval.type="all", full.stats=T)

sapply(markers.dlpfc.t.design.counts, function(x){table(x$FDR<0.05)})
    ##      Ambig.lowNtrxts Astro Excit Inhib Micro Oligo   OPC
    #FALSE           28127 28121 27994 28114 28113 28104 28126
    #TRUE                1     7   134    14    15    24     2



## normalized and transformed "logcounts"
markers.dlpfc.t.design.log <- findMarkers(sce.dlpfc.PB, groups=sce.dlpfc.PB$cellType,
                                          assay.type="logcounts", design=design.PB,
                                          direction="up", pval.type="all", full.stats=T)

sapply(markers.dlpfc.t.design.log, function(x){table(x$FDR<0.05)})
    ##       Ambig.lowNtrxts Astro Excit Inhib Micro Oligo   OPC
    # FALSE           28077 28050 27843 28100 27735 28053 28097
    # TRUE               51    78   285    28   393    75    31

markerList.PB.tDesign.dlpfc.log <- lapply(markers.dlpfc.t.design.log, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)


## Normalized $counts, but not log-transformed
assay(sce.dlpfc.PB, "countsNormd") <- t(apply(assay(sce.dlpfc.PB, "counts"), 1, function(x) {x/LSFvec}))


markers.dlpfc.t.design.countsN <- findMarkers(sce.dlpfc.PB, groups=sce.dlpfc.PB$cellType,
                                              assay.type="countsNormd", design=design.PB,
                                              direction="up", pval.type="all", full.stats=T)

sapply(markers.dlpfc.t.design.countsN, function(x){table(x$FDR<0.05)})
    ##  Ambig.lowNtrxts Astro Excit Inhib Micro Oligo   OPC
    #FALSE       27600 27576 27576 27990 27028 27594 27899
    #TRUE          528   552   552   138  1100   534   229

markerList.PB.tDesign.dlpfc.countsN <- lapply(markers.dlpfc.t.design.countsN, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)

# Between "logcounts" & "countsNormd"    
intersectingMarkers <- list()
for(i in names(markerList.PB.tDesign.dlpfc.log)){
  intersectingMarkers[[i]] <- intersect(markerList.PB.tDesign.dlpfc.countsN[[i]],
                                        markerList.PB.tDesign.dlpfc.log[[i]])
}
lengths(intersectingMarkers)
    ## so most of those detected in "logcounts" are included in "countsNormd" (not surprising)


        ## Using 'block=' is N/A in this PB'd setting:
        #    -> warnings such as "no within-block comparison between ____ and ____"
        


### end skip - with no Ambig cluster ===============================
#sce.dlpfc.PBnoAmbig <- sce.dlpfc.PB[ ,!sce.dlpfc.PB$cellType=="Ambig.lowNtrxts"]

    ## MNT comment - replacing the below 'sce.dlpfc.PBnoAmbig', with just 'sce.dlpfc.PB',
     #     to make more follow-able


#sce.dlpfc.PB$logcounts <- NULL
#LSFvec.noAmbig <- librarySizeFactors(sce.dlpfc.PB)
#sce.dlpfc.PB <- logNormCounts(sce.dlpfc.PB,
#                              size_factors=LSFvec.noAmbig)

## Remove that level too
#sce.dlpfc.PB$cellType <- factor(sce.dlpfc.PB$cellType,
#                                       levels=unique(sce.dlpfc.PB$cellType))

markers.dlpfc.t.noAmbig <- findMarkers(sce.dlpfc.PB, groups=sce.dlpfc.PB$cellType,
                                       assay.type="logcounts",
                                       direction="up", pval.type="all", full.stats=T)

sapply(markers.dlpfc.t.noAmbig, function(x){table(x$FDR<0.05)})
    # still none


# With 'design='? (and go ahead and use normalized counts--"countsNormd")
design.PB.noAmbig <- model.matrix(~sce.dlpfc.PB$processDate)
    ## same result if you used $sample or $donor (for DLPFC, where it's confounded)
design.PB.noAmbig <- design.PB.noAmbig[ , -1, drop=F] # 'drop=F' to keep as matrix - otherwise turns into numeric vec

# "logcounts"
#markers.dlpfc.t.design.noAmbig.log <- findMarkers(sce.dlpfc.PB, groups=sce.dlpfc.PB$cellType,
markers.dlpfc.t.design.log <- findMarkers(sce.dlpfc.PB, groups=sce.dlpfc.PB$cellType,
                                          assay.type="logcounts", design=design.PB.noAmbig,
                                          direction="up", pval.type="all", full.stats=T)

sapply(markers.dlpfc.t.design.noAmbig.log, function(x){table(x$FDR<0.05)})
    ##       Astro Excit Inhib Micro Oligo   OPC
    #  FALSE 28074 27912 28112 27219 27997 28118
    #  TRUE     54   216    16   909   131    10   - when first had rm'd all-0 genes (old)

    ##       Astro Excit Inhib Micro Oligo   OPC
    #  FALSE 28057 27895 28095 27202 27980 28101
    #  TRUE     54   216    16   909   131    10   - when first dropping "Ambig.lowNtrxts" (same result)


# "countsNormd" - need to re-normalize after dropping 'Ambig.lowNtrxts'
#assay(sce.dlpfc.PB, "countsNormd") <- NULL
assay(sce.dlpfc.PB, "countsNormd") <- t(apply(assay(sce.dlpfc.PB, "counts"), 1,
                                                     #function(x) {x/LSFvec.noAmbig}))
                                                     function(x) {x/LSFvec})) # from up top - MNT 05Mar2020

#markers.dlpfc.t.design.noAmbig.countsN <- findMarkers(sce.dlpfc.PB, groups=sce.dlpfc.PB$cellType,
markers.dlpfc.t.design.countsN <- findMarkers(sce.dlpfc.PB, groups=sce.dlpfc.PB$cellType,
                                              assay.type="countsNormd", design=design.PB.noAmbig,
                                              direction="up", pval.type="all", full.stats=T)

sapply(markers.dlpfc.t.design.noAmbig.countsN, function(x){table(x$FDR<0.05)})
    ##      Astro Excit Inhib Micro Oligo   OPC
    #FALSE 27470 27467 27980 26384 27428 27888
    #TRUE    658   661   148  1744   700   240   - when first had rm'd all-0 genes (old)

    ##      Astro Excit Inhib Micro Oligo   OPC
    # FALSE 27453 27450 27957 26367 27411 27871
    # TRUE    658   661   154  1744   700   240  - when first dropping "Ambig.lowNtrxts" (slightly diff)


markerList.PB.tDesign.dlpfc.noAmbig.countsN <- lapply(markers.dlpfc.t.design.noAmbig.countsN, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)


sum(lengths(markerList.PB.tDesign.dlpfc.countsN))
    ## 3633
sum(lengths(markerList.PB.tDesign.dlpfc.noAmbig.countsN))
    ## 4151
length(intersect(unlist(markerList.PB.tDesign.dlpfc.countsN), unlist(markerList.PB.tDesign.dlpfc.noAmbig.countsN)))
    ## 2794


save(markers.dlpfc.t.design.log, markers.dlpfc.t.design.countsN,
#     markers.dlpfc.t.design.noAmbig.log, markers.dlpfc.t.design.noAmbig.countsN,
#     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_DLPFC_n2_MNTFeb2020.rda")
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_DLPFC_n2_findMarkers_MNTMar2020.rda")


# For reference
#write.csv(markerList.PB.tDesign.dlpfc.log[["Ambig.lowNtrxts"]], row.names = F,
#          file="temp_ambigDLPFC_markers.csv")


    ## For reference: can get rowRanges/full rowData from spatial SCE?
    #load("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Human_DLPFC_Visium_processedData_rseList.rda",
    #     verbose=TRUE)
    ## rseList, geom_spatial
    
    names(rseList)  # numeric IDs of length 12
    rse.eg.ST <- rseList[["151507"]]
    head(rowData(rse.eg.ST))  # has 'gene_id', 'gene_version', 'gene_name','gene_source', 'gene_biotype'
    
    table(rowData(rse.eg.ST)$gene_id == rowData(sce.dlpfc)$ID)
        ## all TRUE
    


## Gene list enrichment analyses - Nvm ====
    load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_DLPFC_n2_MNTFeb2020.rda",
         verbose=T)
    rm(markers.dlpfc.t.design.countsN, markers.dlpfc.t.design.log)
        # Since poor correlation of the ambig.lowNtrxts cluster
        # (this cluster seems to be just driven by low captured libraries)
    rm(markers.dlpfc.t.design.noAmbig.countsN)  # let's just stay conservative with stats
    head(markers.dlpfc.t.design.noAmbig.log[["Excit"]])
        ## There's no single test statistic so let's just do manual modeling
         #    for AnJa's ST enrichment approach
    rm(markers.dlpfc.t.design.noAmbig.log)
    ## ====

## Load in sce.dlpfc and cluster:sample bulk
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_DLPFC-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc, ref.sampleInfo
        
# First drop "Ambig.lowNtrxts" (168 nuclei)
sce.dlpfc <- sce.dlpfc[ ,sce.dlpfc$cellType != "Ambig.lowNtrxts"]
sce.dlpfc$cellType <- droplevels(sce.dlpfc$cellType)

# Then make the pseudo-bulked SCE
sce.dlpfc.PB <- aggregateAcrossCells(sce.dlpfc, ids=paste0(sce.dlpfc$sample,":",sce.dlpfc$cellType),
                                     use_exprs_values="counts")

# Drop genes with all 0's
sce.dlpfc.PB <- sce.dlpfc.PB[!rowSums(assay(sce.dlpfc.PB, "counts"))==0, ]
    ## keeps 28111 genes

# Remove stored `sizeFactors()` because this will mess you up
#     * Also, to be safe, can always provide manually-computed SFs:
sizeFactors(sce.dlpfc.PB) <- NULL
LSFvec <- librarySizeFactors(sce.dlpfc.PB)
sce.dlpfc.PB <- logNormCounts(sce.dlpfc.PB, size_factors=LSFvec)

# Clean up colData
colData(sce.dlpfc.PB) <- colData(sce.dlpfc.PB)[ ,c(13:17,19:20)]


        ## Actually what is this?:
        load("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/rda/tstats_Human_DLPFC_snRNAseq_Nguyen.Rdata",
             verbose=T)
            # tstats_Human_DLPFC_snRNAseq_Nguyen_topLayer
        
        dim(tstats_Human_DLPFC_snRNAseq_Nguyen_topLayer)  # 692 x 31


## Load gene lists from ST work === === === ===
load("rdas/geneLists-fromSTpaper_forEnrichTests_MNT.rda", verbose=T)
    # geneLists.fromST
        ## 37 lists



### Adapted from `/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/layer_specificity.R`
  #     ** will eventually want to put into side-script; clean up; then add into this script


###### Direct limma approach ####
#################################

## Extract the count data
mat <- assays(sce.dlpfc.PB)$logcounts

## Build a group model
mod <- with(colData(sce.dlpfc.PB), model.matrix(~ 0 + cellType))
colnames(mod) <- gsub('cellType', '', colnames(mod))

corfit <- duplicateCorrelation(mat, mod, block = sce.dlpfc.PB$donor)
corfit$consensus.correlation
    # [1] 0.03835962

fit <-
  lmFit(
    mat,
    design = mod,
    block = sce.dlpfc.PB$donor,
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
colnames(cellType_contrasts) <-
  apply(cellType_combs, 2, paste, collapse = '-')

eb_contrasts <- eBayes(contrasts.fit(fit, cellType_contrasts))

## Tabulating significant hits
pvals_contrasts <- eb_contrasts$p.value

data.frame(
  'FDRsig' = colSums(apply(pvals_contrasts, 2, p.adjust, 'fdr') < 0.05),
  'Pval10-6sig' = colSums(pvals_contrasts < 1e-6),
  'Pval10-8sig' = colSums(pvals_contrasts < 1e-8)
)
    #             FDRsig Pval10.6sig Pval10.8sig
    #Astro-Excit   5636         247           7
    #Astro-Inhib   4573         183           4
    #Astro-Micro   6878         499          40
    #Astro-Oligo   4487         233           6
    #Astro-OPC     3211         131           6
    #Excit-Inhib    907          28           1
    #Excit-Micro   9029         586          43
    #Excit-Oligo   5984         306          11
    #Excit-OPC     4311         162           4
    #Inhib-Micro   8129         534          45
    #Inhib-Oligo   5061         261          11
    #Inhib-OPC     2971         113           3
    #Micro-Oligo   6063         464          37
    #Micro-OPC     6558         501          49
    #Oligo-OPC     3617         205           9



## Then each cellType vs the rest
cellType_idx <- splitit(sce.dlpfc.PB$cellType)

eb0_list <- lapply(cellType_idx, function(x) {
  res <- rep(0, ncol(sce.dlpfc.PB))
  res[x] <- 1
  m <- model.matrix(~ res)
  eBayes(
    lmFit(
      mat,
      design = m,
      block = sce.dlpfc.PB$donor,
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
    #Astro    537          49          16
    #Excit    690          90          26
    #Inhib    123          28           5
    #Micro   3589         359         116
    #Oligo    770          71          13
    #OPC      113          24           4

## With t > 0
    #      FDRsig Pval10.6sig Pval10.8sig
    #Astro    461          47          16
    #Excit    635          89          26
    #Inhib    121          28           5
    #Micro   1909         264         101
    #Oligo    480          42          10
    #OPC       97          24           4




## Save for later
eb_contrasts.dlpfc.broad <- eb_contrasts
eb_list.dlpfc.broad <- eb0_list

save(eb_contrasts.dlpfc.broad, eb_list.dlpfc.broad, sce.dlpfc.PB,
     file = 'rdas/markers-stats_DLPFC_n2_manualContrasts_MNTMar2020.rda')


### How does this compare to results of `findMarkers()`? === === === === ===
load("rdas/markers-stats_DLPFC_n2_manualContrasts_MNTMar2020.rda", verbose=T)
    # eb_contrasts.dlpfc.broad, eb_list.dlpfc.broad, sce.dlpfc.PB

## Extract the p-values and compute fdrs
pvals0_contrasts <- sapply(eb_list.dlpfc.broad, function(x) {
  x$p.value[, 2, drop = FALSE]
})

fdrs0_contrasts = apply(pvals0_contrasts, 2, p.adjust, "fdr")

## Extract the tstats
t0_contrasts <- sapply(eb_list.dlpfc.broad, function(x) {
  x$t[, 2, drop = FALSE]
})

rownames(fdrs0_contrasts) <- rownames(sce.dlpfc.PB)
rownames(t0_contrasts) <- rownames(sce.dlpfc.PB)

markerList.PB.manual <- lapply(colnames(fdrs0_contrasts), function(x){
  rownames(fdrs0_contrasts)[fdrs0_contrasts[ ,x] < 0.05 & t0_contrasts[ ,x] > 0]
})
names(markerList.PB.manual) <- colnames(fdrs0_contrasts)


## findMarkers() results - test was already just for up-regulated genes
sapply(markers.dlpfc.t.design.log, function(x){table(x$FDR<0.05)})

markerList.PB.dlpfc.tDesign.log <- lapply(markers.dlpfc.t.design.log, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)


lengths(markerList.PB.manual)
    # Astro Excit Inhib Micro Oligo   OPC
    #   461   635   121  1909   480    97

lengths(markerList.PB.dlpfc.tDesign.log)
    # Astro Excit Inhib Micro Oligo   OPC
    #    54   216    16   909   131    10

sapply(names(markerList.PB.manual), function(x){
  length(intersect(markerList.PB.manual[[x]],
                   markerList.PB.dlpfc.tDesign.log[[x]]))}
  )
    ## Astro Excit Inhib Micro Oligo   OPC
     #    50   203    14   853   114     7      - pretty good overlap between the latter and the manual-limma method



### MNT aside, 05Mar2020 ===========
  # Noticed there are 17 genes which are 0 expression AFTER first dropping "Ambig.lowNtrxts" cluster

load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_DLPFC-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    ## sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc, ref.sampleInfo

# First drop "Ambig.lowNtrxts" (168 nuclei)
sce.dlpfc.noAmbig <- sce.dlpfc[ ,sce.dlpfc$cellType != "Ambig.lowNtrxts"]
sce.dlpfc.noAmbig$cellType <- droplevels(sce.dlpfc.noAmbig$cellType)
# Then make the pseudo-bulked SCE
sce.dlpfc.PBafterDrop <- aggregateAcrossCells(sce.dlpfc.noAmbig, ids=paste0(sce.dlpfc.noAmbig$sample,":",sce.dlpfc.noAmbig$cellType),
                                     use_exprs_values="counts")
# Drop genes with all 0's
sce.dlpfc.PBafterDrop <- sce.dlpfc.PBafterDrop[!rowSums(assay(sce.dlpfc.PBafterDrop, "counts"))==0, ]

## OR

# Just make the pseudo-bulked SCE, without dropping that cluster
sce.dlpfc.PB <- aggregateAcrossCells(sce.dlpfc, ids=paste0(sce.dlpfc$sample,":",sce.dlpfc$cellType),
                                              use_exprs_values="counts")
# Drop genes with all 0's
sce.dlpfc.PB <- sce.dlpfc.PB[!rowSums(assay(sce.dlpfc.PB, "counts"))==0, ]

## Then
genesOfInterest <- setdiff(rownames(sce.dlpfc.PB), rownames(sce.dlpfc.PBafterDrop))
#  [1] "FOXD2"      "FASLG"      "FGFBP2"     "TRGV5"      "TRBC1"
#  [6] "GPR174"     "PRF1"       "CLEC2B"     "KRT72"      "KRT73"
#  [11] "IFNG"       "TBX3"       "GZMB"       "TMEM30B"    "AC106782.6"
#  [16] "S1PR4"      "NKG7

rowSums(assay(sce.dlpfc.PB, "counts")[genesOfInterest, ])
    ## mostly 1's; highest is 11 for PRF11
table(assay(sce.dlpfc, "counts")["PRF1",sce.dlpfc$cellType=="Ambig.lowNtrxts"] ==0) # 166 zeros

    ##    - so this is definitely just a poor, lowly-captured cluster; explains the poor
    ##      intra-cluster correlation, compared to other clusters





