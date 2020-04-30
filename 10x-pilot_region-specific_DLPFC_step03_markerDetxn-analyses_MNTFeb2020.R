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
    







###### Direct limma approach ####
#################################
    
    ### Adapted from `/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/
      #               Layer_Guesses/layer_specificity.R`

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




### Top markers to print / potentially test with RNA-scope === === === ===
load('/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_DLPFC_n2_manualContrasts_MNTMar2020.rda',
     verbose=T)
    # eb_contrasts.dlpfc.broad, eb_list.dlpfc.broad, sce.dlpfc.PB

# follow chunk 'How does this compare to results of `findMarkers()`?' for fdr & t mats


markerList.PB.manual <- lapply(colnames(fdrs0_contrasts), function(x){
  rownames(fdrs0_contrasts)[fdrs0_contrasts[ ,x] < 0.001 & t0_contrasts[ ,x] > 0]
})
names(markerList.PB.manual) <- colnames(fdrs0_contrasts)
lengths(markerList.PB.manual)
    # Astro Excit Inhib Micro Oligo   OPC
    #    76   130    28   494    74    22

markerTs.fdr.001 <- lapply(colnames(fdrs0_contrasts), function(x){
  as.matrix(t0_contrasts[fdrs0_contrasts[ ,x] < 0.001 & t0_contrasts[ ,x] > 0, x])
})

names(markerTs.fdr.001) <- colnames(fdrs0_contrasts)

markerList.sorted <- lapply(markerTs.fdr.001, function(x){
  x[,1][order(x, decreasing=TRUE)]
})

genes2plot <- lapply(markerList.sorted, function(x){head(x, n=20)})


## Let's plot some expression of these to see how much are 'real' (not driven by outliers)
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_DLPFC-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc, ref.sampleInfo
    rm(chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc, ref.sampleInfo)

# As before, first drop "Ambig.lowNtrxts" (168 nuclei)
sce.dlpfc <- sce.dlpfc[ ,sce.dlpfc$cellType != "Ambig.lowNtrxts"]
sce.dlpfc$cellType <- droplevels(sce.dlpfc$cellType)


pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_DLPFC-n2_top20markers_logExprs_Mar2020.pdf", height=7.5, width=9.5)
for(i in 1:length(genes2plot)){
  print(
    plotExpression(sce.dlpfc, exprs_values = "logcounts", features=c(names(genes2plot[[i]])),
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:6], length(genes2plot[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(genes2plot)[i], " top 20 markers"))
  )
}
dev.off()


## What if just subset on protein-coding first?
library(rtracklayer)
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/zref_genes-GTF-fromGRCh38-3.0.0_33538.rda", verbose=T)
    # gtf

table(gtf$gene_biotype)
    #        antisense       IG_C_gene IG_C_pseudogene       IG_D_gene       IG_J_gene
    #             5497              14               9              37              18
    #  IG_J_pseudogene       IG_V_gene IG_V_pseudogene         lincRNA  protein_coding
    #                3             144             188            7484           19912
    #        TR_C_gene       TR_D_gene       TR_J_gene TR_J_pseudogene       TR_V_gene
    #                6               4              79               4             106
    #  TR_V_pseudogene
    #               33

table(rownames(sce.dlpfc) %in% gtf$gene_name)
    # FALSE  TRUE
    #    48 33490    - probably because of the `uniquify`
table(rowData(sce.dlpfc)$Symbol %in% gtf$gene_name)
    #  TRUE
    # 33538

# Are they the same order?
table(rowData(sce.dlpfc)$ID == gtf$gene_id) # all TRUE

table(!rowSums(assay(sce.dlpfc, "counts"))==0)  # 28111     - good
keepVec <- !rowSums(assay(sce.dlpfc, "counts"))==0

gtf <- gtf[keepVec, ]
# Then
table(gtf$gene_id == rowData(sce.dlpfc.PB)$ID)  # all 28111 TRUE      - good

## Make pt-coding list
markerList.sorted.pt <- lapply(markerList.sorted, function(x){
  x[names(x) %in% gtf$gene_name[gtf$gene_biotype=="protein_coding"]]
})

lengths(markerList.sorted)
    # Astro Excit Inhib Micro Oligo   OPC
    #    76   130    28   494    74    22

lengths(markerList.sorted.pt)
    #Astro Excit Inhib Micro Oligo   OPC
    #   51    36    19   397    40     9



genes2plot.pt <- lapply(markerList.sorted.pt, function(x){head(x, n=20)})

# Plot these
pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_DLPFC-n2_top20markers_logExprs_pt-coding_Mar2020.pdf", height=7.5, width=9.5)
for(i in 1:length(genes2plot.pt)){
  print(
    plotExpression(sce.dlpfc, exprs_values = "logcounts", features=c(names(genes2plot.pt[[i]])),
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:6], length(genes2plot.pt[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(genes2plot.pt)[i], " top 20 protein-coding markers"))
  )
}
dev.off()


## How much they intersect with the top protein-coding-agnostic set?
sapply(names(genes2plot), function(x){intersect(names(genes2plot[[x]]), names(genes2plot.pt[[x]]))})
    # $Astro
    #   [1] "SIX5"     "OTOS"     "ZIC5"     "CRACR2B"  "OTX1"     "KCNE5"
    #   [7] "C6orf223" "P2RY2"    "DMRTA1"   "TMPRSS3"  "RASL12"   "IGFN1"
    #   [13] "TFAP2C"   "GDPD2"    "C22orf31" "LAMA1"
    # 
    # $Excit
    #   [1] "OR14I1"   "CD200R1L" "TNNT2"    "ROS1"     "C4orf54"
    # 
    # $Inhib
    #   [1] "DLX5"    "DLX2"    "PRLHR"   "DLX6"    "DLX1"    "SLC32A1" "SP9"
    #   [8] "LHX6"    "NKX6-3"  "CHRNA2"  "PLSCR5"  "DEPDC1"  "CARD10"
    # 
    # $Micro
    #   [1] "RNASE6"   "IL1B"     "VSIG4"    "MPEG1"    "RGS18"    "RGS1"
    #   [7] "CCL3"     "GIMAP6"   "TREM2"    "NCF4"     "TREML1"   "IL1A"
    #   [13] "SERPINA1"
    # 
    # $Oligo
    #   [1] "FFAR1"   "HOXD1"   "GPIHBP1" "SMIM6"   "NGFR"    "LYRM9"   "SLC5A11"
    #   [8] "RNASE1"  "TMEM88B" "NIPAL4"
    # 
    # $OPC
    #   [1] "NR0B1"   "CPXM1"   "KCNG4"   "DCAF4L2" "B3GNT7"  "COL20A1" "KIF18B"
    #   [8] "CCKAR"




### Modeling of spatially-registered neuronal subtypes ================================
  # Added MNT late-Mar/Apr2020

  ## Extract the count data
  mat <- assays(sce.dlpfc.PB)$logcounts
  


## Load SCE with new info
load("rdas/regionSpecific_DLPFC-n2_SCE_cellTypesSplit-fromST_Apr2020.rda", verbose=T)
    # sce.dlpfc.st, clusterRefTab.dlpfc, chosen.hvgs.dlpfc, ref.sampleInfo

table(sce.dlpfc.st$cellType.split)

# First drop "Ambig.lowNtrxts" (168 nuclei)
sce.dlpfc.st <- sce.dlpfc.st[ ,sce.dlpfc.st$cellType.split != "Ambig.lowNtrxts"]
sce.dlpfc.st$cellType.split <- droplevels(sce.dlpfc.st$cellType.split)

# Then make the pseudo-bulked SCE
sce.dlpfc.st.PB <- aggregateAcrossCells(sce.dlpfc.st, ids=paste0(sce.dlpfc.st$sample,":",sce.dlpfc.st$cellType.split),
                                    use_exprs_values="counts")

# Clean up colData
colData(sce.dlpfc.st.PB) <- colData(sce.dlpfc.st.PB)[ ,c(13:17,19:21)]

# Drop genes with all 0's
sce.dlpfc.st.PB <- sce.dlpfc.st.PB[!rowSums(assay(sce.dlpfc.st.PB, "counts"))==0, ]
    ## keeps 28111 genes

# Remove stored `sizeFactors()` because this will mess you up
#     * Also, to be safe, can always provide manually-computed SFs:
sizeFactors(sce.dlpfc.st.PB) <- NULL
LSFvec <- librarySizeFactors(sce.dlpfc.st.PB)
sce.dlpfc.st.PB <- logNormCounts(sce.dlpfc.st.PB, size_factors=LSFvec)

## Extract the count data
mat <- assays(sce.dlpfc.st.PB)$logcounts

## Build a group model
mod <- with(colData(sce.dlpfc.st.PB), model.matrix(~ 0 + cellType.split))
colnames(mod) <- gsub('cellType.split', '', colnames(mod))

corfit <- duplicateCorrelation(mat, mod, block = sce.dlpfc.st.PB$donor)
corfit$consensus.correlation
    # [1] 0.03655833
    # (0.03835962 at the broad cell type level, so not much difference..!)


## Test each cellType vs the rest - only for neuronal, bc this wouldn't affect glial stats
cellType_idx <- splitit(sce.dlpfc.st.PB$cellType.split)

## Will have to do this by excitatory or inhib. subtypes, separately
eb0_list_neurons <- list()
# Excitatory
for(k in names(cellType_idx)[ss(names(cellType_idx),"\\.",1) %in% c("Excit")]){
  # Subtype of interest
  subtype <- rep(0, ncol(sce.dlpfc.st.PB))
  subtype[ cellType_idx[[k]] ] <- 1
  # Add broad excitatory coef to model              # or don't..!
  #broadtype <- rep(0, ncol(sce.dlpfc.st.PB))
  #broadtype[grep("Excit", colnames(mat))] <- "Excit"
  
  m <- model.matrix(~ subtype)# + broadtype)
  
  eb0_list_neurons[[k]] <- eBayes(
    lmFit(
      mat,
      design = m,
      block = sce.dlpfc.st.PB$donor,
      correlation = corfit$consensus.correlation
    )
  )
}

# Inhibitory
for(k in names(cellType_idx)[ss(names(cellType_idx),"\\.",1) %in% c("Inhib")]){
  # Subtype of interest
  subtype <- rep(0, ncol(sce.dlpfc.st.PB))
  subtype[ cellType_idx[[k]] ] <- 1
  # Add broad excitatory coef to model              # or don't..!
  #broadtype <- rep(0, ncol(sce.dlpfc.st.PB))
  #broadtype[grep("Inhib", colnames(mat))] <- "Inhib"
  
  m <- model.matrix(~ subtype)# + broadtype)
  
  eb0_list_neurons[[k]] <- eBayes(
    lmFit(
      mat,
      design = m,
      block = sce.dlpfc.st.PB$donor,
      correlation = corfit$consensus.correlation
    )
  )
}

    ## with broad term included -> 'eb0_list_neurons.broad'


## Extract the p-values
pvals0_contrasts.st <- sapply(eb0_list_neurons, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts.st) <- rownames(sce.dlpfc.st.PB)

## Extract the tstats
t0_contrasts.st <- sapply(eb0_list_neurons, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts.st) <- rownames(sce.dlpfc.st.PB)

# Any signif with (+) t's?
data.frame(
  'FDRsig.05' = colSums(apply(pvals0_contrasts.st, 2, p.adjust, 'fdr') < 0.05 & 
                          t0_contrasts.st > 0),
  'FDRsig.01' = colSums(apply(pvals0_contrasts.st, 2, p.adjust, 'fdr') < 0.01 & 
                          t0_contrasts.st > 0),
  'Pval10-6sig' = colSums(pvals0_contrasts.st < 1e-6 & 
                            t0_contrasts.st > 0),
  'Pval10-8sig' = colSums(pvals0_contrasts.st < 1e-8 & 
                            t0_contrasts.st > 0)
)
    #                FDRsig.05 FDRsig.01 Pval10.6sig Pval10.8sig
    # Excit.ambig          240       176          24          14
    # Excit.L2:3           140        60          33          13
    # Excit.L3:4           117        28          18           5
    # Excit.L4:5            96        44          29           6
    # Excit.L5             210        97          48          26
    # Excit.L5:6           117        49          29          13
    # Excit.L6.broad        71        17           8           5
    # Inhib.1              217       136          72          36
    # Inhib.2              602       289           9           0
    # Inhib.3              113        36          11           7
    # Inhib.4               91        41          19           6
    # Inhib.5              101        43          30          16
    # Inhib.6              182        94          58          30


    ## Note that before, at the broad-cell-type-level:

            ## With t > 0
            #      FDRsig Pval10.6sig Pval10.8sig
            #Astro    461          47          16
            #Excit    635          89          26
            #Inhib    121          28           5     - very poor numbers here.
            #Micro   1909         264         101
            #Oligo    480          42          10
            #OPC       97          24           4



    ## ** If DO include broad cell type in model:
    #                FDRsig.05 FDRsig.01 Pval10.6sig Pval10.8sig
    # Excit.ambig          212       158          24          13
    # Excit.L2:3            64        38          23           9
    # Excit.L3:4            39        20           9           3
    # Excit.L4:5            47        34          19           6
    # Excit.L5             172        77          43          21
    # Excit.L5:6            67        34          23          11
    # Excit.L6.broad        30         7           7           5
    # Inhib.1              221       127          66          35
    # Inhib.2              498       194           4           0
    # Inhib.3               72        21          10           5
    # Inhib.4               92        23          14           4
    # Inhib.5               80        39          28          14
    # Inhib.6              163        87          56          26


          ## Exploration of effects of including broad neuronal type term ====
          t0_contrasts.st.broad <- sapply(eb0_list_neurons.broad, function(x) {
              x$t[, 2, drop = FALSE]
            })
          rownames(t0_contrasts.st.broad) <- rownames(sce.dlpfc.st.PB)
          
          col.pal = brewer.pal(10,"RdBu")
          
          ## Within-approach t's
          pdf("pdfs/exploration/zTemp_sub-clusterTs_DLPFC_Mar2020.pdf")
          pheatmap(round(cor(t0_contrasts.st.broad), 3), main="With broad neuronal term", display_numbers=T,
                   cluster_cols=F, cluster_rows=F, color=col.pal, breaks=seq(-1,1,by=.2))
          pheatmap(round(cor(t0_contrasts.st), 3), main="No broad neuronal term", display_numbers=T,
                   cluster_cols=F, cluster_rows=F, color=col.pal, breaks=seq(-1,1,by=.2))
          dev.off()
          
          
          broadTypeEffect <- sapply(eb0_list_neurons.broad, function(x) {x$coefficients[ ,2]})
          
          apply(broadTypeEffect, 2, quantile)
          
          apply(broadTypeEffect, 2, mean)
          
          
          # What about between the two approaches??
          subtypes <-colnames(t0_contrasts.st)
          sapply(subtypes, function(x){cor(t0_contrasts.st[ ,x], t0_contrasts.st.broad[ ,x])})
              # Excit.ambig     Excit.L2:3     Excit.L3:4     Excit.L4:5       Excit.L5
              # 0.9292742      0.8743395      0.8710827      0.8393692      0.8927272
              # Excit.L5:6 Excit.L6.broad        Inhib.1        Inhib.2        Inhib.3
              # 0.8643328      0.8698543      0.9664003      0.9620088      0.9203056
              # Inhib.4        Inhib.5        Inhib.6
              # 0.8663071      0.8517854      0.8657143
          
          # end explore ====


## (as before)  -   let's just go with fdr < 0.05 for now
fdrs0_contrasts.st = apply(pvals0_contrasts.st, 2, p.adjust, "fdr")
rownames(fdrs0_contrasts.st) <- rownames(sce.dlpfc.st.PB)


markerList.PB.manual <- lapply(colnames(fdrs0_contrasts.st), function(x){
  rownames(fdrs0_contrasts.st)[fdrs0_contrasts.st[ ,x] < 0.05 & t0_contrasts.st[ ,x] > 0]
})
names(markerList.PB.manual) <- colnames(fdrs0_contrasts.st)

# Get t's
markerTs.fdr.05 <- lapply(colnames(fdrs0_contrasts.st), function(x){
  as.matrix(t0_contrasts.st[fdrs0_contrasts.st[ ,x] < 0.05 & t0_contrasts.st[ ,x] > 0, x])
})

names(markerTs.fdr.05) <- colnames(fdrs0_contrasts.st)

# Sort for largest difference
markerList.sorted <- lapply(markerTs.fdr.05, function(x){
  x[,1][order(x, decreasing=TRUE)]
})

genes2plot <- lapply(markerList.sorted, function(x){head(x, n=20)})



#pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_DLPFC-n2_top20markers-SUBtypes_logExprs_Apr2020.pdf", height=7.5, width=9.5)
pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_DLPFC-n2_top20markers-SUBtypes_noBroadTerm_Apr2020.pdf", height=7.5, width=9.5)
for(i in 1:length(genes2plot)){
  print(
    plotExpression(sce.dlpfc.st, exprs_values = "logcounts", features=c(names(genes2plot[[i]])),
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F, theme_size=8) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:17], length(genes2plot[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(genes2plot)[i], " top 20 markers"))
  )
}
dev.off()



## What if just subset on protein-coding first? ===
library(rtracklayer)
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/zref_genes-GTF-fromGRCh38-3.0.0_33538.rda", verbose=T)
# gtf

table(rownames(sce.dlpfc.st) %in% gtf$gene_name)
# FALSE  TRUE
#    48 33490    - probably because of the `uniquify`
table(rowData(sce.dlpfc.st)$Symbol %in% gtf$gene_name)
#  TRUE
# 33538

# Are they the same order?
table(rowData(sce.dlpfc.st)$ID == gtf$gene_id) # all TRUE

table(!rowSums(assay(sce.dlpfc.st, "counts"))==0)  # 28111     - good
keepVec <- !rowSums(assay(sce.dlpfc.st, "counts"))==0

gtf <- gtf[keepVec, ]
# Then
table(gtf$gene_id == rowData(sce.dlpfc.st.PB)$ID)  # all 28111 TRUE      - good

## Make pt-coding list
markerList.sorted.pt <- lapply(markerList.sorted, function(x){
  x[names(x) %in% gtf$gene_name[gtf$gene_biotype=="protein_coding"]]
})

lengths(markerList.sorted)

lengths(markerList.sorted.pt)
    #   Excit.ambig     Excit.L2:3     Excit.L3:4     Excit.L4:5       Excit.L5
    #            72             39             45             34             85
    #    Excit.L5:6 Excit.L6.broad        Inhib.1        Inhib.2        Inhib.3
    #            39             24             93            488             45
    #       Inhib.4        Inhib.5        Inhib.6
    #            27             48             81
    

    ## ** If DO include broad cell type in model:
    #


genes2plot.pt <- lapply(markerList.sorted.pt, function(x){head(x, n=20)})

# Plot these
#pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_DLPFC-n2_top20markers-SUBtypes_logExprs_pt-coding_Apr2020.pdf", height=7.5, width=9.5)
pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_DLPFC-n2_top20markers-SUBtypes_noBroadTerm_pt-coding_Apr2020.pdf", height=7.5, width=9.5)
for(i in 1:length(genes2plot.pt)){
  print(
    plotExpression(sce.dlpfc.st, exprs_values = "logcounts", features=c(names(genes2plot.pt[[i]])),
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F, theme_size=8) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:17], length(genes2plot.pt[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(genes2plot.pt)[i], " top 20 protein-coding markers"))
  )
}
dev.off()



# Save stats
eb_list.dlpfc.neuronalSubs.simple <- eb0_list_neurons
save(eb_list.dlpfc.neuronalSubs.simple, sce.dlpfc.st.PB,
     file="rdas/markers-stats_DLPFC_n2_manualContrasts_neuronalSubs_noBroadTerm_MNTApr2020.rda")


# And results from same WITH modeling broad neuronal type
eb_list.dlpfc.neuronalSubs <- eb0_list_neurons.broad
save(eb_list.dlpfc.neuronalSubs, sce.dlpfc.st.PB,
     file="rdas/markers-stats_DLPFC_n2_manualContrasts_neuronalSubs_MNTApr2020.rda")


# Write 'genes2plot's to a csv
names(genes2plot.pt) <- paste0(names(genes2plot.pt),"_pt")
top20genes <- cbind(sapply(genes2plot, names), sapply(genes2plot.pt, names))
top20genes <- top20genes[ ,sort(colnames(top20genes))]

write.csv(top20genes, file="tables/top20genesLists_DLPFC-n2_cellTypesSplit.csv")






### Single-nucleus-level tests for cell-type-specific genes ================================
  # MNT 23Apr2020 - added after seeing much better results in pan-brain analysis

## Load SCE with new info
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_DLPFC-n2_SCE_cellTypesSplit-fromST_Apr2020.rda", verbose=T)
    # sce.dlpfc.st, clusterRefTab.dlpfc, chosen.hvgs.dlpfc, ref.sampleInfo

table(sce.dlpfc.st$cellType.split)

# First drop "Ambig.lowNtrxts" (168 nuclei)
sce.dlpfc.st <- sce.dlpfc.st[ ,sce.dlpfc.st$cellType.split != "Ambig.lowNtrxts"]
sce.dlpfc.st$cellType.split <- droplevels(sce.dlpfc.st$cellType.split)

# Remove 0 genes across all nuclei
sce.dlpfc.st <- sce.dlpfc.st[!rowSums(assay(sce.dlpfc.st, "counts"))==0, ]  # keeps same 28111 genes


## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.dlpfc.st), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # 'drop=F' to keep as matrix - otherwise turns into numeric vec
    ## Get: "Error in .ranksafe_qr(full.design) : design matrix is not of full rank"
     #   if try to put 0 intercept

    ## MNT comment: actually doesn't matter

# Run pairwise t-tests
markers.dlpfc.t.design <- findMarkers(sce.dlpfc.st, groups=sce.dlpfc.st$cellType.split,
                                    assay.type="logcounts", design=mod, test="t",
                                    direction="up", pval.type="all", full.stats=T)

sapply(markers.dlpfc.t.design, function(x){table(x$FDR<0.05)})
    #      Astro Excit.ambig Excit.L2:3 Excit.L3:4 Excit.L4:5 Excit.L5 Excit.L5:6
    # FALSE 33308       33463      33536      33443      33513    33315      33490
    # TRUE    230          75          2         95         25      223         48
    #       Excit.L6.broad Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5 Inhib.6 Micro Oligo   OPC
    # FALSE          33501   33302   33347   33405   33513   33509   33512 33006 33390 33384
    # TRUE              37     236     191     133      25      29      26   532   148   154


## WMW: Blocking on donor (this test doesn't take 'design=' argument) ===
markers.dlpfc.wilcox.block <- findMarkers(sce.dlpfc.st, groups=sce.dlpfc.st$cellType.split,
                                        assay.type="logcounts", block=sce.dlpfc.st$donor, test="wilcox",
                                        direction="up", pval.type="all", full.stats=T)
    
    # no warnings as in pan-brain analyses, but NO results of FDR<0.05...:
sapply(markers.dlpfc.wilcox.block, function(x){table(x$FDR<0.05)})
    
    sapply(names(markers.dlpfc.wilcox.block), function(x){quantile(markers.dlpfc.wilcox.block[[x]]$FDR)})
        # Astro Excit.ambig Excit.L2:3 Excit.L3:4 Excit.L4:5  Excit.L5 Excit.L5:6 Excit.L6.broad
        # 0%   0.3332772           1          1          1          1 0.6095081          1              1
        # 25%  1.0000000           1          1          1          1 1.0000000          1              1
        # 50%  1.0000000           1          1          1          1 1.0000000          1              1
        # 75%  1.0000000           1          1          1          1 1.0000000          1              1
        # 100% 1.0000000           1          1          1          1 1.0000000          1              1
        # Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5 Inhib.6     Micro     Oligo       OPC
        # 0%         1       1       1       1       1       1 0.1850964 0.1071413 0.1395478
        # 25%        1       1       1       1       1       1 1.0000000 1.0000000 1.0000000
        # 50%        1       1       1       1       1       1 1.0000000 1.0000000 1.0000000
        # 75%        1       1       1       1       1       1 1.0000000 1.0000000 1.0000000
        # 100%       1       1       1       1       1       1 1.0000000 1.0000000 1.0000000


## Binomial ===
markers.dlpfc.binom.block <- findMarkers(sce.dlpfc.st, groups=sce.dlpfc.st$cellType.split,
                                       assay.type="logcounts", block=sce.dlpfc.st$donor, test="binom",
                                       direction="up", pval.type="all", full.stats=T)

sapply(markers.dlpfc.binom.block, function(x){table(x$FDR<0.05)})
# none - and these are ALL 1's across the board....

## Save all these for future reference
save(markers.dlpfc.t.design, #markers.dlpfc.wilcox.block, markers.dlpfc.binom.block,
     file="rdas/markers-stats_DLPFC_n2_findMarkers-SN-LEVEL_MNTApr2020.rda")


# Print these to pngs
markerList.t <- lapply(markers.dlpfc.t.design, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)

genes.top40.t <- lapply(markerList.t, function(x){head(x, n=40)})
smallerSets <- c("Excit.L4:5", "Inhib.4", "Inhib.5", "Inhib.6")

#dir.create("pdfs/exploration/DLPFC/")
## ~40+ marker genes
for(i in setdiff(names(genes.top40.t), c(smallerSets, "Excit.L2:3"))){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/DLPFC/DLPFC_t-sn-level-top40markers-",gsub(":",".",i),"_logExprs_Apr2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.dlpfc.st, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:17], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers: single-nucleus-level p.w. t-tests"))
  )
  dev.off()
}

## < 40 ('smallerSets')
for(i in smallerSets){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/DLPFC/DLPFC_t-sn-level-top40markers-",gsub(":",".",i),"_logExprs_Apr2020.png"), height=1900/2, width=1200)
  print(
    plotExpression(sce.dlpfc.st, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:17], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top markers with FDR < 0.05: single-nucleus-level p.w. t-tests"))
  )
  dev.off()
}

## 'Excit.L2:3'
png("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/DLPFC/DLPFC_t-sn-level-top40markers-Excit.L2.3_logExprs_Apr2020.png", height=1900/8, width=1200/3)
print(
  plotExpression(sce.dlpfc.st, exprs_values = "logcounts", features=genes.top40.t[["Excit.L2:3"]],
                 x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, #ncol=5,
                 add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                              geom = "crossbar", width = 0.3,
                                              colour=rep(tableau20[1:17], length(genes.top40.t[["Excit.L2:3"]]))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 10)) +  
    ggtitle(label="Excit.L2:3 top markers with FDR < 0.05: single-nucleus-level p.w. t-tests")
)
dev.off()






### Cluster-vs-all single-nucleus-level iteration ================================
  # MNT 30Apr2020

## Load SCE with new info
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_DLPFC-n2_SCE_cellTypesSplit-fromST_Apr2020.rda", verbose=T)
    # sce.dlpfc.st, clusterRefTab.dlpfc, chosen.hvgs.dlpfc, ref.sampleInfo

table(sce.dlpfc.st$cellType.split)

# First drop "Ambig.lowNtrxts" (168 nuclei)
sce.dlpfc.st <- sce.dlpfc.st[ ,sce.dlpfc.st$cellType.split != "Ambig.lowNtrxts"]
sce.dlpfc.st$cellType.split <- droplevels(sce.dlpfc.st$cellType.split)

# Remove 0 genes across all nuclei
sce.dlpfc.st <- sce.dlpfc.st[!rowSums(assay(sce.dlpfc.st, "counts"))==0, ]  # keeps same 28111 genes


## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.dlpfc.st), model.matrix(~ donor))

markers.dlpfc.t.1vAll <- list()
for(i in unique(sce.dlpfc.st$cellType.split)){
  # Make temporary contrast
  sce.dlpfc.st$contrast <- ifelse(sce.dlpfc.st$cellType.split==i, 1, 0)
  # Test cluster vs. all
  markers.dlpfc.t.1vAll[[i]] <- findMarkers(sce.dlpfc.st, groups=sce.dlpfc.st$contrast,
                                assay.type="logcounts", design=mod, test="t",
                                direction="up", pval.type="all", full.stats=T)
}
    # There were 17 warnings (use warnings() to see them)
    #   [17x]:
    #   Warning messages:
    #   1: In .fit_lm_internal(x, subset.row, groups, design = design,  ... :
    #     automatically removed intercept column

class(markers.dlpfc.t.1vAll[["Astro"]]) # SimpleList
dim(markers.dlpfc.t.1vAll[["Astro"]]) # NULL
head(markers.dlpfc.t.1vAll[["Astro"]])
    ## For some reason all the results are in the second List entry (first is always empty)

    table(markers.dlpfc.t.1vAll[["Astro"]][[2]]$FDR<0.05) # 4892
    table(markers.dlpfc.t.1vAll[["Astro"]][[2]]$stats.0$log.FDR < log10(0.05))  # 6143
        ## Why are these different?  Though probably should take the former, because
         # there is only one test for each of these clusters and not sure
         # WHAT these are...

    head(markers.dlpfc.t.1vAll[["Oligo"]][[2]])
        ## Nice, MBP and PLP1 are in the top 6
    

sapply(markers.dlpfc.t.1vAll, function(x){
  table(x[[2]]$stats.0$log.FDR < log10(.001))
  })
    #       Oligo Astro Inhib.4 Excit.L4:5 Micro Inhib.6   OPC Excit.L2:3 Excit.ambig
    # FALSE 25366 23223   21034      19061 23999   20956 23802      22479       20869
    # TRUE   2745  4888    7077       9050  4112    7155  4309       5632        7242
    #       Excit.L3:4 Excit.L5:6 Excit.L6.broad Inhib.5 Excit.L5 Inhib.1 Inhib.2
    # FALSE      23916      22320          21501   23246    23878   26699   26540
    # TRUE        4195       5791           6610    4865     4233    1412    1571
    #       Inhib.3
    # FALSE   25766
    # TRUE     2345


## Let's save this along with the previous pairwise results
load("rdas/markers-stats_DLPFC_n2_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
save(markers.dlpfc.t.1vAll, markers.dlpfc.t.design,
     file="rdas/markers-stats_DLPFC_n2_findMarkers-SN-LEVEL_MNTApr2020.rda")


# Replace that empty slot with the entry with the actul stats
markers.dlpfc.t.1vAll <- lapply(markers.dlpfc.t.1vAll, function(x){x[[2]]})


## Print these to pngs
markerList.t.1vAll <- lapply(markers.dlpfc.t.1vAll, function(x){
  rownames(x)[x$stats.0[ ,"log.FDR"] < log10(0.000001)]
  }
)
genes.top40.t <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/DLPFC/DLPFC_t-sn-level_1vALL_top40markers-",gsub(":",".",i),"_logExprs_Apr2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.dlpfc.st, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:17], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers: single-nucleus-level p.w. t-tests, cluster-vs-all"))
  )
  dev.off()
}


## How do they intersect?
markerList.t.pw <- lapply(markers.dlpfc.t.design, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)

# From pairwise t-tests, FDR < 0.05
lengths(markerList.t.pw)

# From cluster-vs-all others, FDR < 1e6
lengths(markerList.t.1vAll)

# Intersection
sapply(names(markerList.t.pw), function(c){
  length(intersect(markerList.t.pw[[c]],
                   markerList.t.1vAll[[c]]))
})



## Checking Berger method for pval.type=="all" ====
    ps.chosen <- markers.dlpfc.t.design[["Astro"]]$p.value
    ps.chosen <- log10(ps.chosen)
    head(ps.chosen)
        # [1] -95.30511 -78.89036 -63.47255 -43.07298 -38.99470 -38.68466
    
    ps.pairwise <- cbind(sapply(c(3:18), function(x){markers.dlpfc.t.design[["Astro"]][ ,x]$log.p.value}))
    head(ps.pairwise)
    
    maxpick <- apply(ps.pairwise, 1, which.max)
    
    max.p <- vector()
    for(i in 1:length(maxpick)){
      max.p[i] <- ps.pairwise[i, maxpick[i]]
    }
    
    head(max.p, n=20)
    head(ps.chosen, n=20) # not the same as the above..
    cor(max.p, ps.chosen) # == 1
        # interesting. So I guess it's NOT just the max p. but this transformed somehow
    




