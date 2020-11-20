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
    ## same result if you used $sample or $donor (for HPC, where it's confounded)
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

## Load in sce.hpc and cluster:sample bulk
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




### Top markers to print / potentially test with RNA-scope === === === ===
load('/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_HPC-n3_manualContrasts_MNTMar2020.rda',
     verbose=T)
    # eb_contrasts.hpc.broad, eb_list.hpc.broad, sce.hpc.PB

# follow chunk 'How does this compare to results of `findMarkers()`?' for fdr & t mats
colnames(pvals0_contrasts)[1] <- "Tcell"
colnames(fdrs0_contrasts)[1] <- "Tcell"
colnames(t0_contrasts)[1] <- "Tcell"

markerList.PB.manual <- lapply(colnames(fdrs0_contrasts), function(x){
  rownames(fdrs0_contrasts)[fdrs0_contrasts[ ,x] < 0.001 & t0_contrasts[ ,x] > 0]
})
names(markerList.PB.manual) <- colnames(fdrs0_contrasts)
lengths(markerList.PB.manual)
    # Tcell Astro Excit Inhib Micro Oligo   OPC
    #   660   306   114   116   526   109    83

markerTs.fdr.001 <- lapply(colnames(fdrs0_contrasts), function(x){
  as.matrix(t0_contrasts[fdrs0_contrasts[ ,x] < 0.001 & t0_contrasts[ ,x] > 0, x])
})

names(markerTs.fdr.001) <- colnames(fdrs0_contrasts)

markerList.sorted <- lapply(markerTs.fdr.001, function(x){
  x[,1][order(x, decreasing=TRUE)]
})

genes2plot <- lapply(markerList.sorted, function(x){head(x, n=20)})




## Let's plot some expression of these to see how much are 'real' (not driven by outliers)
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo
    rm(chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo)

# As before, first drop "Ambig.lowNtrxts" (101 nuclei)
sce.hpc <- sce.hpc[ ,sce.hpc$cellType != "Ambig.lowNtrxts"]
sce.hpc$cellType <- droplevels(sce.hpc$cellType)
# Then rename "Ambig.glial" to "Tcell" (26 nuclei)
#     (A posteriori - from downstream marker exploration)
sce.hpc$cellType <- factor(gsub(pattern="Ambig.glial", "Tcell", sce.hpc$cellType))


pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_HPC-n3_top20markers_logExprs_Mar2020.pdf", height=7.5, width=9.5)
for(i in 1:length(genes2plot)){
  print(
    plotExpression(sce.hpc, exprs_values = "logcounts", features=c(names(genes2plot[[i]])),
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:7], length(genes2plot[[i]]))) +
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

table(rownames(sce.hpc) %in% gtf$gene_name)
# FALSE  TRUE
#    48 33490    - probably because of the `uniquify`
table(rowData(sce.hpc)$Symbol %in% gtf$gene_name)
#  TRUE
# 33538

# Are they the same order?
table(rowData(sce.hpc)$ID == gtf$gene_id) # all TRUE

table(!rowSums(assay(sce.hpc, "counts"))==0)  # 28757     - good
keepVec <- !rowSums(assay(sce.hpc, "counts"))==0

gtf <- gtf[keepVec, ]
# Then
table(gtf$gene_id == rowData(sce.hpc.PB)$ID)  # all 28757 TRUE      - good

## Make pt-coding list
markerList.sorted.pt <- lapply(markerList.sorted, function(x){
  x[names(x) %in% gtf$gene_name[gtf$gene_biotype=="protein_coding"]]
})

lengths(markerList.sorted)
# Tcell Astro Excit Inhib Micro Oligo   OPC
#   660   306   114   116   526   109    83

lengths(markerList.sorted.pt)
# Tcell Astro Excit Inhib Micro Oligo   OPC
#   621   187    49    70   341    44    41



genes2plot.pt <- lapply(markerList.sorted.pt, function(x){head(x, n=20)})

# Plot these
pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_HPC-n3_top20markers_logExprs_pt-coding_Mar2020.pdf", height=7.5, width=9.5)
for(i in 1:length(genes2plot.pt)){
  print(
    plotExpression(sce.hpc, exprs_values = "logcounts", features=c(names(genes2plot.pt[[i]])),
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:7], length(genes2plot.pt[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(genes2plot.pt)[i], " top 20 protein-coding markers"))
  )
}
dev.off()


## How much they intersect with the top protein-coding-agnostic set?
sapply(names(genes2plot), function(x){intersect(names(genes2plot[[x]]), names(genes2plot.pt[[x]]))})
    # $Tcell
    # [1] "CD2"     "FCRL6"   "GZMK"    "UBASH3A" "TBX21"   "CD3D"    "IL7R"
    # [8] "SLAMF7"  "SLAMF6"  "CD3E"    "GZMA"    "SLAMF1"  "IL2RB"   "CCL4"
    # [15] "SLFN12L" "CD3G"    "ACAP1"
    # 
    # $Astro
    # [1] "SLC2A4"   "SIX5"     "OTX1"     "NFATC4"   "CYP4F11"  "OTOS"
    # [7] "TFAP2C"   "TMEM200B" "GDPD2"    "RASL12"
    # 
    # $Excit
    # [1] "EMX1"    "KCNV1"   "MS4A8"   "OR14I1"  "DNAJC5G" "GHSR"    "GRM2"
    # [8] "KNCN"    "LRRC10B" "SPANXN4" "P2RX2"
    # 
    # $Inhib
    # [1] "NKX2-1" "PLSCR5" "CRH"    "SP8"    "DLX2"   "LHX6"   "CHRNA2" "DLX5"
    # [9] "NKX6-3" "HTR3A"  "PRLHR"  "EREG"   "TRH"    "DLX1"   "OR4D1"  "OPN5"
    # 
    # $Micro
    # [1] "TREML1"  "IL1A"    "LILRB2"  "ASCL4"   "FCGR1A"  "TLR7"    "IL1B"
    # [8] "VSIG4"   "ANKRD22" "TLR8"    "CD300C"
    # 
    # $Oligo
    # [1] "FXYD4"     "PIP"       "LYRM9"     "NGFR"      "TRPV6"     "SECISBP2L"
    # 
    # $OPC
    # [1] "DCAF4L2" "GDF6"    "KCNG4"   "NR0B1"   "CSPG4"   "TIMP4"   "COL20A1"
    # [8] "GPR17"


# Write 'genes2plot's to a csv
names(genes2plot.pt) <- paste0(names(genes2plot.pt),"_pt")
top20genes <- cbind(sapply(genes2plot, names), sapply(genes2plot.pt, names))
top20genes <- top20genes[ ,sort(colnames(top20genes))]

write.csv(top20genes, file="tables/top20genesLists_HPC-n3_cellTypes.csv")





    ### ========================== ###
    ### SINGLE-NUCLEUS-LEVEL TESTS ###
    ### ========================== ###


### Single-nucleus-level tests for cell-type-specific genes ================================
# MNT 23Apr2020 - added after seeing much better results in pan-brain analysis

## Load SCE with new info
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
# sce.hpc, clusterRefTab.hpc, chosen.hvgs.hpc, ref.sampleInfo

table(sce.hpc$cellType.split)

# First drop "Ambig.lowNtrxts" (101 nuclei)
sce.hpc <- sce.hpc[ ,sce.hpc$cellType.split != "Ambig.lowNtrxts"]
sce.hpc$cellType.split <- droplevels(sce.hpc$cellType.split)

# Remove 0 genes across all nuclei
sce.hpc <- sce.hpc[!rowSums(assay(sce.hpc, "counts"))==0, ]  # keeps same 28757 genes


## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.hpc), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`


# Run pairwise t-tests
markers.hpc.t.design <- findMarkers(sce.hpc, groups=sce.hpc$cellType.split,
                                      assay.type="logcounts", design=mod, test="t",
                                      direction="up", pval.type="all", full.stats=T)

sapply(markers.hpc.t.design, function(x){table(x$FDR<0.05)})
    #       Astro Excit.1 Excit.2 Excit.3 Excit.4 Excit.5 Inhib.1 Inhib.2 Inhib.3
    # FALSE 28244   28435   28540   28498   28140   28394   28552   28629   28593
    # TRUE    513     322     217     259     617     363     205     128     164
    #       Inhib.4 Inhib.5 Micro Oligo   OPC Tcell
    # FALSE   28601   28595 28093 28356 28506 28006
    # TRUE      156     162   664   401   251   751


## WMW: Blocking on donor (this test doesn't take 'design=' argument) ===
markers.hpc.wilcox.block <- findMarkers(sce.hpc, groups=sce.hpc$cellType.split,
                                          assay.type="logcounts", block=sce.hpc$donor, test="wilcox",
                                          direction="up", pval.type="all", full.stats=T)

# no warnings as in pan-brain analyses, but NO results of FDR<0.05...:
sapply(markers.hpc.wilcox.block, function(x){table(x$FDR<0.05)})
      # Actually some decent results but many subclusters with 0 hits


## Binomial ===
markers.hpc.binom.block <- findMarkers(sce.hpc, groups=sce.hpc$cellType.split,
                                         assay.type="logcounts", block=sce.hpc$donor, test="binom",
                                         direction="up", pval.type="all", full.stats=T)

sapply(markers.hpc.binom.block, function(x){table(x$FDR<0.05)})
    # only a couple dozen hits for glia, only - disregard these

## Save all these for future reference
save(markers.hpc.t.design, markers.hpc.wilcox.block, #markers.hpc.binom.block,
     file="rdas/markers-stats_HPC-n3_findMarkers-SN-LEVEL_MNTApr2020.rda")


# Print these to pngs
markerList.t <- lapply(markers.hpc.t.design, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)

genes.top40.t <- lapply(markerList.t, function(x){head(x, n=40)})


#dir.create("pdfs/exploration/HPC/")
for(i in names(genes.top40.t)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/HPC/HPC_t-sn-level_pairwise_top40markers-", i, "_logExprs_Apr2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.hpc, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:15], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers: single-nucleus-level p.w. t-tests"))
  )
  dev.off()
}






### Cluster-vs-all single-nucleus-level iteration ================================
# MNT 30Apr2020

## Load SCE with new info
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo

table(sce.hpc$cellType.split)

# First drop "Ambig.lowNtrxts" (101 nuclei)
sce.hpc <- sce.hpc[ ,sce.hpc$cellType.split != "Ambig.lowNtrxts"]
sce.hpc$cellType.split <- droplevels(sce.hpc$cellType.split)

# Remove 0 genes across all nuclei
sce.hpc <- sce.hpc[!rowSums(assay(sce.hpc, "counts"))==0, ]  # keeps same 28757 genes


## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.hpc), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`


markers.hpc.t.1vAll <- list()
for(i in levels(sce.hpc$cellType.split)){
  # Make temporary contrast
  sce.hpc$contrast <- ifelse(sce.hpc$cellType.split==i, 1, 0)
  # Test cluster vs. all
  markers.hpc.t.1vAll[[i]] <- findMarkers(sce.hpc, groups=sce.hpc$contrast,
                                            assay.type="logcounts", design=mod, test="t",
                                            direction="up", pval.type="all", full.stats=T)
}

    ## Then, temp set of stats to get the standardized logFC
    temp.1vAll <- list()
    for(i in levels(sce.hpc$cellType.split)){
      # Make temporary contrast
      sce.hpc$contrast <- ifelse(sce.hpc$cellType.split==i, 1, 0)
      # Test cluster vs. all
      temp.1vAll[[i]] <- findMarkers(sce.hpc, groups=sce.hpc$contrast,
                                     assay.type="logcounts", design=mod, test="t",
                                     std.lfc=TRUE,
                                     direction="up", pval.type="all", full.stats=T)
    }



    ## As with DLPFC, for some reason all the results are in the
     #    second List entry (first is always empty)

head(markers.hpc.t.1vAll[["Oligo"]][[2]])
    ## Nice, MBP and PLP1 are again in the top 6


sapply(markers.hpc.t.1vAll, function(x){
  table(x[[2]]$stats.0$log.FDR < log10(.001))
})
    #       Oligo Micro   OPC Inhib.5 Inhib.2 Astro Inhib.3 Excit.2 Inhib.4 Tcell
    # FALSE 24914 22821 23236   24401   23540 21612   21436   22608   24460 26858
    # TRUE   3843  5936  5521    4356    5217  7145    7321    6149    4297  1899
    #       Inhib.1 Excit.5 Excit.3 Excit.1 Excit.4
    # FALSE   25456   25913   19431   23726   24170
    # TRUE     3301    2844    9326    5031    4587



# Replace that empty slot with the entry with the actul stats
markers.hpc.t.1vAll <- lapply(markers.hpc.t.1vAll, function(x){ x[[2]] })
# Same for that with std.lfc
temp.1vAll <- lapply(temp.1vAll, function(x){ x[[2]] })

# Now just pull from the 'stats.0' DataFrame column
markers.hpc.t.1vAll <- lapply(markers.hpc.t.1vAll, function(x){ x$stats.0 })
temp.1vAll <- lapply(temp.1vAll, function(x){ x$stats.0 })

# Re-name std.lfc column and add to the first result
for(i in names(temp.1vAll)){
  colnames(temp.1vAll[[i]])[1] <- "std.logFC"
  markers.hpc.t.1vAll[[i]] <- cbind(markers.hpc.t.1vAll[[i]], temp.1vAll[[i]]$std.logFC)
  # Oh the colname is kept weird
  colnames(markers.hpc.t.1vAll[[i]])[4] <- "std.logFC"
  # Then re-organize
  markers.hpc.t.1vAll[[i]] <- markers.hpc.t.1vAll[[i]][ ,c("logFC","std.logFC","log.p.value","log.FDR")]
}





## Let's save this along with the previous pairwise results
save(markers.hpc.t.1vAll, markers.hpc.t.design, markers.hpc.wilcox.block,
#     file="rdas/markers-stats_HPC-n3_findMarkers-SN-LEVEL_MNTApr2020.rda")
#     (deleting this older version - doesn't have the std.lfc result)
     file="rdas/markers-stats_HPC-n3_findMarkers-SN-LEVEL_MNTMay2020.rda")


## Print these to pngs
markerList.t.1vAll <- lapply(markers.hpc.t.1vAll, function(x){
  rownames(x)[x$log.FDR < log10(0.000001)]
 }
)
genes.top40.t <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/HPC/HPC_t-sn-level_1vALL_top40markers-",gsub(":",".",i),"_logExprs_Apr2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.hpc, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:15], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers: single-nucleus-level p.w. t-tests, cluster-vs-all"))
  )
  dev.off()
}


## How do they intersect?
markerList.t.pw <- lapply(markers.hpc.t.design, function(x){
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

    # Of top 40's:
    sapply(names(markerList.t.pw), function(c){
      length(intersect(lapply(markerList.t.pw, function(l){head(l,n=40)})[[c]],
                       lapply(markerList.t.1vAll, function(l){head(l,n=40)})[[c]]
                       ))
    })
    #   Astro Excit.1 Excit.2 Excit.3 Excit.4 Excit.5 Inhib.1 Inhib.2 Inhib.3 Inhib.4
    #      26      24      18      15      28      33      22      10      22      23
    # Inhib.5   Micro   Oligo     OPC   Tcell
    #      19      30      26      20      39


    
## Write these top 40 lists to a csv
names(markerList.t.pw) <- paste0(names(markerList.t.pw),"_pw")
names(markerList.t.1vAll) <- paste0(names(markerList.t.1vAll),"_1vAll")

top40genes <- cbind(sapply(markerList.t.pw, function(x) head(x, n=40)),
                    sapply(markerList.t.1vAll, function(y) head(y, n=40)))
top40genes <- top40genes[ ,sort(colnames(top40genes))]

write.csv(top40genes, file="tables/top40genesLists_HPC-n3_cellType.split_SN-LEVEL-tests_May2020.csv",
          row.names=FALSE)




### MNT add 18Nov2020 =================================
  # -> What if add param/requirement that for any given subcluster, median expression has to > 0?
load("rdas/markers-stats_HPC-n3_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.hpc.t.1vAll, markers.hpc.t.design, markers.hpc.wilcox.block

## Load SCE 
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo

table(sce.hpc$cellType.split)

# First drop "Ambig.lowNtrxts" (101 nuclei)
sce.hpc <- sce.hpc[ ,sce.hpc$cellType.split != "Ambig.lowNtrxts"]
sce.hpc$cellType.split <- droplevels(sce.hpc$cellType.split)

# Remove 0 genes across all nuclei
sce.hpc <- sce.hpc[!rowSums(assay(sce.hpc, "counts"))==0, ]


## Make list of Boolean param / cell subtype ===
#medianNon0.idx <- list()
cellSubtype.idx <- splitit(sce.hpc$cellType.split)
medianNon0.idx <- lapply(cellSubtype.idx, function(x){
  apply(as.matrix(assay(sce.hpc, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

lengths(medianNon0.idx)
sapply(medianNon0.idx, head)

# Add these to the stats for each set of markers
for(i in names(markers.hpc.t.1vAll)){
  markers.hpc.t.1vAll[[i]] <- cbind(markers.hpc.t.1vAll[[i]],
                                    medianNon0.idx[[i]][match(rownames(markers.hpc.t.1vAll[[i]]),
                                                           names(medianNon0.idx[[i]]))])
  colnames(markers.hpc.t.1vAll[[i]])[5] <- "non0median"
}


## Use these restrictions to print (to png) a refined top 40, as before ===
markerList.t.1vAll <- lapply(markers.hpc.t.1vAll, function(x){
  rownames(x)[x$log.FDR < log10(0.000001) & x$non0median==TRUE]
  }
)
    # lengths(markerList.t.1vAll)     # ( **without $non0median==TRUE restriction )
        #   Astro Excit.1 Excit.2 Excit.3 Excit.4 Excit.5 Inhib.1 Inhib.2 Inhib.3 Inhib.4
        #    5668    3876    4581    7414    3246    2033    2314    3679    5184    2806
        # Inhib.5   Micro   Oligo     OPC   Tcell
        #    2962    4934    3323    4182    1406

lengths(markerList.t.1vAll)
    #   Astro Excit.1 Excit.2 Excit.3 Excit.4 Excit.5 Inhib.1 Inhib.2 Inhib.3 Inhib.4
    #     847    1958    2659    3412    2000     832    1594    2111    2487    1861
    # Inhib.5   Micro   Oligo     OPC   Tcell
    #    1993     802     953    1065     354

genes.top40.t <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/HPC/HPC_t-sn-level_1vALL_top40markers-REFINED-",gsub(":",".",i),"_logExprs_Nov2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.hpc, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:15], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers, refined: single-nucleus-level p.w. t-tests, cluster-vs-all"))
  )
  dev.off()
}



## Do the same with the pairwise result ('markers.hpc.t.design') === === ===
# Add these to the stats for each set of markers
for(i in names(markers.hpc.t.design)){
  markers.hpc.t.design[[i]] <- cbind(markers.hpc.t.design[[i]],
                                    medianNon0.idx[[i]][match(rownames(markers.hpc.t.design[[i]]),
                                                           names(medianNon0.idx[[i]]))])
  colnames(markers.hpc.t.design[[i]])[17] <- "non0median"
}

markerList.t <- lapply(markers.hpc.t.design, function(x){
  rownames(x)[x$FDR < 0.05 & x$non0median==TRUE]
  }
)
    # lengths(markerList.t)     # ( **without $non0median==TRUE restriction )
        #   Astro Excit.1 Excit.2 Excit.3 Excit.4 Excit.5 Inhib.1 Inhib.2 Inhib.3 Inhib.4
        #     513     322     217     259     617     363     205     128     164     156
        # Inhib.5   Micro   Oligo     OPC   Tcell
        #     162     664     401     251     751

lengths(markerList.t)
    #  Astro Excit.1 Excit.2 Excit.3 Excit.4 Excit.5 Inhib.1 Inhib.2 Inhib.3 Inhib.4
    #    249     167      56     146     296      97      76      27      62      66
    #Inhib.5   Micro   Oligo     OPC   Tcell
    #     69     282     338     157     178


genes.top40.t <- lapply(markerList.t, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/HPC/HPC_t-sn-level_pairwise_top40markers-REFINED-", i, "_logExprs_Nov2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.hpc, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:15], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers, refined: single-nucleus-level p.w. t-tests"))
  )
  dev.off()
}

## Then write a new CSV of these refined top 40 genes ===
names(markerList.t) <- paste0(names(markerList.t),"_pw")
names(markerList.t.1vAll) <- paste0(names(markerList.t.1vAll),"_1vAll")

# Many of the PW results don't have at least 40 markers:
extend.idx <- names(which(lengths(markerList.t) < 40))
for(i in extend.idx){
  markerList.t[[i]] <- c(markerList.t[[i]], rep("", 40-length(markerList.t[[i]])))
}

top40genes <- cbind(sapply(markerList.t, function(x) head(x, n=40)),
                    sapply(markerList.t.1vAll, function(y) head(y, n=40)))
top40genes <- top40genes[ ,sort(colnames(top40genes))]

write.csv(top40genes, file="tables/top40genesLists-REFINED_HPC-n3_cellType.split_Nov2020.csv",
          row.names=FALSE)






