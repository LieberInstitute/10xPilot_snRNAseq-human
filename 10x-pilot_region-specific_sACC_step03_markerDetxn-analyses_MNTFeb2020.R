### MNT 10x snRNA-seq workflow: step 03 - marker detection
###   **Region-specific analyses**
###     - (2x) sACC samples from: Br5161 & Br5212
### Initiated MNT 12Feb2020
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

load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    ## sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc, ref.sampleInfo

# sce.sacc has 7047 nuclei

# First drop "Ambig.lowNtrxts" (43 nuclei)
sce.sacc <- sce.sacc[ ,sce.sacc$cellType != "Ambig.lowNtrxts"]
sce.sacc$cellType <- droplevels(sce.sacc$cellType)

# **This region was unique in that its collapsed clusters had multiple sub-clusters
#   - because all the others' were collapsed on broad cell types, we'll do the same here
#     (further exploration at the non-broad-level can be done later)

sce.sacc$cellType.broad <- as.factor(ss(as.character(sce.sacc$cellType), "\\.", 1))

# Then make the pseudo-bulked SCE
sce.sacc.PB <- aggregateAcrossCells(sce.sacc, ids=paste0(sce.sacc$sample,":",sce.sacc$cellType.broad),
                                   use_exprs_values="counts")

# Clean up colData
colData(sce.sacc.PB) <- colData(sce.sacc.PB)[ ,c(13:17,19:21)]

# Drop genes with all 0's
sce.sacc.PB <- sce.sacc.PB[!rowSums(assay(sce.sacc.PB, "counts"))==0, ]
    ## keeps 28774 genes

# Remove stored `sizeFactors()` because this will mess you up
#     * Also, to be safe, can always provide manually-computed SFs:
sizeFactors(sce.sacc.PB) <- NULL
LSFvec <- librarySizeFactors(sce.sacc.PB)
sce.sacc.PB <- logNormCounts(sce.sacc.PB, size_factors=LSFvec)



## Find markers using stringent [max-p-value-of-all-pw-comparisons] test ('pval.type="all"') ===
# No 'design'
markers.sacc.t <- findMarkers(sce.sacc.PB, groups=sce.sacc.PB$cellType.broad,
                             assay.type="logcounts",
                             direction="up", pval.type="all", full.stats=T)

sapply(markers.sacc.t, function(x){table(x$FDR<0.05)})
    # none


# With 'design='? (and go ahead and use normalized counts--"countsNormd")

# Add new combined design factor: sample:processDate
    # (these samples were processed on the same day, whereas other regions we model on processDate,
    #  which is usually confounded with $sample/$donor anyhow)
sce.sacc.PB$sampleBatch <- paste0(sce.sacc.PB$sample, ":", sce.sacc.PB$processDate)
design.PB <- model.matrix(~sce.sacc.PB$sampleBatch)
design.PB <- design.PB[ , -1, drop=F] # 'drop=F' to keep as matrix - otherwise turns into numeric vec

# "logcounts"
markers.sacc.t.design.log <- findMarkers(sce.sacc.PB, groups=sce.sacc.PB$cellType.broad,
                                        assay.type="logcounts", design=design.PB,
                                        direction="up", pval.type="all", full.stats=T)

sapply(markers.sacc.t.design.log, function(x){table(x$FDR<0.05)})
    ##      Astro Excit Inhib Micro Oligo   OPC
    # FALSE 28468 28361 28707 26898 28281 28720
    # TRUE    306   413    67  1876   493    54


# "countsNormd"
assay(sce.sacc.PB, "countsNormd") <- t(apply(assay(sce.sacc.PB, "counts"), 1,
                                            function(x) {x/LSFvec}))


markers.sacc.t.design.countsN <- findMarkers(sce.sacc.PB, groups=sce.sacc.PB$cellType.broad,
                                            assay.type="countsNormd", design=design.PB,
                                            direction="up", pval.type="all", full.stats=T)

sapply(markers.sacc.t.design.countsN, function(x){table(x$FDR<0.05)})
    ##      Astro Excit Inhib Micro Oligo   OPC
    # FALSE 27909 27530 28230 25967 27715 28366
    # TRUE    865  1244   544  2807  1059   408


markerList.PB.tDesign.sacc.countsN <- lapply(markers.sacc.t.design.countsN, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)
sum(lengths(markerList.PB.tDesign.sacc.countsN))
    ## 6927

markerList.PB.tDesign.sacc.log <- lapply(markers.sacc.t.design.log, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)

sum(lengths(markerList.PB.tDesign.sacc.log))
    ## 3209

length(intersect(unlist(markerList.PB.tDesign.sacc.log), unlist(markerList.PB.tDesign.sacc.countsN)))
    ## 3040


save(markers.sacc.t.design.log, markers.sacc.t.design.countsN,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_sACC-n2_findMarkers_MNTMar2020.rda")




###### Direct limma approach ####
#################################

## Load in sce.dlpfc and cluster:sample bulk
# (already done up top)


## Extract the count data
mat <- assays(sce.sacc.PB)$logcounts

## Build a group model
mod <- with(colData(sce.sacc.PB), model.matrix(~ 0 + cellType.broad))
colnames(mod) <- gsub('cellType.broad', '', colnames(mod))

corfit <- duplicateCorrelation(mat, mod, block = sce.sacc.PB$donor)
corfit$consensus.correlation
# [1] 0.00735874

fit <-
  lmFit(
    mat,
    design = mod,
    block = sce.sacc.PB$donor,
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
    #            FDRsig Pval10.6sig Pval10.8sig
    #Astro-Excit   8838         563          42
    #Astro-Inhib   7766         460          23
    #Astro-Micro   9752         959         130
    #Astro-Oligo   7085         571          54
    #Astro-OPC     5178         260          10
    #Excit-Inhib   2032          71           1
    #Excit-Micro  12924        1159         128
    #Excit-Oligo   9529         795          57
    #Excit-OPC     7096         361          30
    #Inhib-Micro  11807        1073         103
    #Inhib-Oligo   8369         677          49
    #Inhib-OPC     5585         253          16
    #Micro-Oligo   8986        1039         123
    #Micro-OPC     9918         903         101
    #Oligo-OPC     6379         483          28



## Then each cellType vs the rest
cellType_idx <- splitit(sce.sacc.PB$cellType.broad)

eb0_list <- lapply(cellType_idx, function(x) {
  res <- rep(0, ncol(sce.sacc.PB))
  res[x] <- 1
  m <- model.matrix(~ res)
  eBayes(
    lmFit(
      mat,
      design = m,
      block = sce.sacc.PB$donor,
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
    ##      FDRsig Pval10.6sig Pval10.8sig
    # Astro    579          76          17
    # Excit    845          93          15
    # Inhib    150          24           5
    # Micro   4812         467         126
    # Oligo   1062          88           7
    # OPC      149          28           5

# For only (+) t-stats
data.frame(
  'FDRsig' = colSums(apply(pvals0_contrasts, 2, p.adjust, 'fdr') < 0.05 &
                       t0_contrasts_cell > 0),
  'Pval10-6sig' = colSums(pvals0_contrasts < 1e-6 &
                            t0_contrasts_cell > 0),
  'Pval10-8sig' = colSums(pvals0_contrasts < 1e-8 &
                            t0_contrasts_cell > 0)
)

    ##      FDRsig Pval10.6sig Pval10.8sig
    # Astro    491          75          17
    # Excit    758          91          15
    # Inhib    146          24           5
    # Micro   2681         330         103
    # Oligo    633          54           4
    # OPC      131          28           5




## Save for later
eb_contrasts.sacc.broad <- eb_contrasts
eb_list.sacc.broad <- eb0_list

save(eb_contrasts.sacc.broad, eb_list.sacc.broad, sce.sacc.PB,
     file = 'rdas/markers-stats_sACC-n2_manualContrasts_MNTMar2020.rda')



## For reference === === ===
table(sce.sacc$cellType.broad, sce.sacc$sample)
    #       sacc.5161 sacc.5212
    # Astro       188       435
    # Excit       444       822
    # Inhib       224       606
    # Micro       235       267
    # Oligo      1852      1400
    # OPC         227       304


table(sce.sacc$cellType, sce.sacc$sample)
    #         sacc.5161 sacc.5212
    # Astro         188       435
    # Excit.1       187       399
    # Excit.2       151       259
    # Excit.3        68       117
    # Excit.4        38        47
    # Inhib.1       142       373
    # Inhib.2        82       233
    # Micro         235       267
    # Oligo        1852      1400
    # OPC           227       304


### How does this compare to results of `findMarkers()`? =======

## Extract the p-values and compute fdrs
pvals0_contrasts <- sapply(eb_list.sacc.broad, function(x) {
  x$p.value[, 2, drop = FALSE]
})

fdrs0_contrasts = apply(pvals0_contrasts, 2, p.adjust, "fdr")

## Extract the tstats
t0_contrasts <- sapply(eb_list.sacc.broad, function(x) {
  x$t[, 2, drop = FALSE]
})

rownames(fdrs0_contrasts) <- rownames(sce.sacc.PB)
rownames(t0_contrasts) <- rownames(sce.sacc.PB)

markerList.PB.manual <- lapply(colnames(fdrs0_contrasts), function(x){
  rownames(fdrs0_contrasts)[fdrs0_contrasts[ ,x] < 0.05 & t0_contrasts[ ,x] > 0]
  # what if more stringent?
  #rownames(fdrs0_contrasts)[fdrs0_contrasts[ ,x] < 0.01 & t0_contrasts[ ,x] > 0]
})
names(markerList.PB.manual) <- colnames(fdrs0_contrasts)


## findMarkers() results - test was already just for up-regulated genes
sapply(markers.sacc.t.design.log, function(x){table(x$FDR<0.05)})

markerList.PB.sacc.tDesign.log <- lapply(markers.sacc.t.design.log, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)


lengths(markerList.PB.manual)
    # Astro Excit Inhib Micro Oligo   OPC
    #   491   758   146  2681   633   131


lengths(markerList.PB.sacc.tDesign.log)
    # Astro Excit Inhib Micro Oligo   OPC
    #   306   413    67  1876   493    54

sapply(names(markerList.PB.manual), function(x){
  length(intersect(markerList.PB.manual[[x]],
                   markerList.PB.sacc.tDesign.log[[x]]))}
)
    # Astro Excit Inhib Micro Oligo   OPC
    #   194   314    45  1691   349    38






### Top markers to print / potentially test with RNA-scope === === === ===
load('/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_sACC-n2_manualContrasts_MNTMar2020.rda',
     verbose=T)
    # eb_contrasts.sacc.broad, eb_list.sacc.broad, sce.sacc.PB

# follow chunk 'How does this compare to results of `findMarkers()`?' for fdr & t mats


markerList.PB.manual <- lapply(colnames(fdrs0_contrasts), function(x){
  rownames(fdrs0_contrasts)[fdrs0_contrasts[ ,x] < 0.001 & t0_contrasts[ ,x] > 0]
})
names(markerList.PB.manual) <- colnames(fdrs0_contrasts)
lengths(markerList.PB.manual)
# Astro Excit Inhib Micro Oligo   OPC
#   100   139    19   750   100    28

markerTs.fdr.001 <- lapply(colnames(fdrs0_contrasts), function(x){
  as.matrix(t0_contrasts[fdrs0_contrasts[ ,x] < 0.001 & t0_contrasts[ ,x] > 0, x])
})

names(markerTs.fdr.001) <- colnames(fdrs0_contrasts)

markerList.sorted <- lapply(markerTs.fdr.001, function(x){
  x[,1][order(x, decreasing=TRUE)]
})

genes2plot <- lapply(markerList.sorted, function(x){head(x, n=20)})


## Let's plot some expression of these to see how much are 'real' (not driven by outliers)
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc, ref.sampleInfo
    rm(chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc, ref.sampleInfo)

# As before, first drop "Ambig.lowNtrxts" (43 nuclei)
sce.sacc <- sce.sacc[ ,sce.sacc$cellType != "Ambig.lowNtrxts"]
sce.sacc$cellType <- droplevels(sce.sacc$cellType)


# Collapse to broad cell types, since markers were defined at that level
sce.sacc$cellType.broad <- as.factor(ss(as.character(sce.sacc$cellType), "\\.", 1))


pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_sACC-n2_top20markers_logExprs_Mar2020.pdf", height=7.5, width=9.5)
for(i in 1:length(genes2plot)){
  print(
    plotExpression(sce.sacc, exprs_values = "logcounts", features=c(names(genes2plot[[i]])),
                   x="cellType.broad", colour_by="cellType.broad", point_alpha=0.5, point_size=.7, ncol=5,
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

table(rownames(sce.sacc) %in% gtf$gene_name)
    # FALSE  TRUE
    #    48 33490    - probably because of the `uniquify`
table(rowData(sce.sacc)$Symbol %in% gtf$gene_name)
    #  TRUE
    # 33538

# Are they the same order?
table(rowData(sce.sacc)$ID == gtf$gene_id) # all TRUE

table(!rowSums(assay(sce.sacc, "counts"))==0)  # 28774     - good
keepVec <- !rowSums(assay(sce.sacc, "counts"))==0

gtf <- gtf[keepVec, ]
# Then
table(gtf$gene_id == rowData(sce.sacc.PB)$ID)  # all 28774 TRUE      - good

## Make pt-coding list
markerList.sorted.pt <- lapply(markerList.sorted, function(x){
  x[names(x) %in% gtf$gene_name[gtf$gene_biotype=="protein_coding"]]
})

lengths(markerList.sorted)
# Astro Excit Inhib Micro Oligo   OPC
#   100   139    19   750   100    28

lengths(markerList.sorted.pt)
# Astro Excit Inhib Micro Oligo   OPC
#    58    37    12   615    68    11



genes2plot.pt <- lapply(markerList.sorted.pt, function(x){head(x, n=20)})

# Plot these
pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_sACC-n2_top20markers_logExprs_pt-coding_Mar2020.pdf", height=7.5, width=9.5)
for(i in 1:length(genes2plot.pt)){
  print(
    plotExpression(sce.sacc, exprs_values = "logcounts", features=c(names(genes2plot.pt[[i]])),
                   x="cellType.broad", colour_by="cellType.broad", point_alpha=0.5, point_size=.7, ncol=5,
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
#  [1] "FCN2"   "SLC2A4" "FHL5"   "CERS1"  "POU4F1"
# 
# $Excit
#  [1] "IZUMO3" "PRTN3"
# 
# $Inhib
#  [1] "DLX2"   "NKX2-1" "PLSCR5" "CHRNA2" "GAD2"   "NKX6-3" "PRLHR"  "EXOC1L"
#  [9] "SLC6A2" "FCRL4"  "SP9"    "ANO1"
# 
# $Micro
#  [1] "CD7"     "TOR4A"   "CCR5"    "CCL3L1"  "LILRA1"  "GAPT"    "SRGN"
#  [8] "LILRA4"  "NCKAP1L" "MYO1F"   "TNFSF10" "MS4A14"
# 
# $Oligo
#  [1] "HSD17B3"    "AC034102.2" "CCP110"     "LPGAT1"     "AC026316.5"
#  [6] "CTNNA3"     "AC079594.2" "HHIP"       "CLDND1"     "LRP2"
#  [11] "IFNA6"
# 
# $OPC
#  [1] "KCNG4"  "FOXD2"  "TM4SF1" "SGCA"   "HTRA3"  "HAS2"   "ACTL10" "TSSK2"
#  [9] "NR0B1"








### Modeling of neuronal subtypes ============================================================
## First load previous sce.sacc.PB & Re-run duplicate correlation
 #            - want to use at the level of just broad cellType:
load('/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_sACC-n2_manualContrasts_MNTMar2020.rda',
     verbose=T)
    # eb_contrasts.sacc.broad, eb_list.sacc.broad, sce.sacc.PB

    ## Extract the count data
    mat <- assays(sce.sacc.PB)$logcounts
    
    ## Build a group model
    mod <- with(colData(sce.sacc.PB), model.matrix(~ 0 + cellType.broad))
    colnames(mod) <- gsub('cellType.broad', '', colnames(mod))
    
    corfit <- duplicateCorrelation(mat, mod, block = sce.sacc.PB$donor)
    corfit$consensus.correlation
        # [1] 0.00735874

    rm(mat, mod, sce.sacc.PB)

### Now set up the collapsed-cluster-level PB'd object:
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc, ref.sampleInfo
    rm(chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc, ref.sampleInfo)
    # sce.sacc has 7047 nuclei

# First drop "Ambig.lowNtrxts" (43 nuclei)
sce.sacc <- sce.sacc[ ,sce.sacc$cellType != "Ambig.lowNtrxts"]
sce.sacc$cellType <- droplevels(sce.sacc$cellType)

# Add 'broad' cellType factor so can perform additive modeling
sce.sacc$cellType.broad <- as.factor(ss(as.character(sce.sacc$cellType), "\\.", 1))


# Then make the pseudo-bulked SCE
sce.sacc.PB <- aggregateAcrossCells(sce.sacc, ids=paste0(sce.sacc$sample,":",sce.sacc$cellType),
                                    use_exprs_values="counts")

# Clean up colData
colData(sce.sacc.PB) <- colData(sce.sacc.PB)[ ,c(13:17,19:21)]

# Drop genes with all 0's
sce.sacc.PB <- sce.sacc.PB[!rowSums(assay(sce.sacc.PB, "counts"))==0, ]
## keeps 28774 genes

# Remove stored `sizeFactors()` because this will mess you up
#     * Also, to be safe, can always provide manually-computed SFs:
sizeFactors(sce.sacc.PB) <- NULL
LSFvec <- librarySizeFactors(sce.sacc.PB)
sce.sacc.PB <- logNormCounts(sce.sacc.PB, size_factors=LSFvec)



## Extract the count data
mat <- assays(sce.sacc.PB)$logcounts

## Each cellType vs the rest - only for neuronal, bc this wouldn't affect glial stats
cellType_idx <- splitit(sce.sacc.PB$cellType)


## Will have to do this by excitatory or inhib. subtypes, separately
eb0_list_neurons <- list()
# Excitatory
for(k in names(cellType_idx)[ss(names(cellType_idx),"\\.",1) %in% c("Excit")]){
  # Subtype of interest
  subtype <- rep(0, ncol(sce.sacc.PB))
  subtype[ cellType_idx[[k]] ] <- 1
  # Add broad excitatory coef to model
  broadtype <- rep(0, ncol(sce.sacc.PB))
  broadtype[grep("Excit", colnames(mat))] <- "Excit"
  
  m <- model.matrix(~ subtype + broadtype)
  
  eb0_list_neurons[[k]] <- eBayes(
                             lmFit(
                               mat,
                               design = m,
                               block = sce.sacc.PB$donor,
                               correlation = corfit$consensus.correlation
                               )
                             )
}

# Inhibitory
for(k in names(cellType_idx)[ss(names(cellType_idx),"\\.",1) %in% c("Inhib")]){
  # Subtype of interest
  subtype <- rep(0, ncol(sce.sacc.PB))
  subtype[ cellType_idx[[k]] ] <- 1
  # Add broad excitatory coef to model
  broadtype <- rep(0, ncol(sce.sacc.PB))
  broadtype[grep("Inhib", colnames(mat))] <- "Inhib"
  
  m <- model.matrix(~ subtype + broadtype)
  
  eb0_list_neurons[[k]] <- eBayes(
    lmFit(
      mat,
      design = m,
      block = sce.sacc.PB$donor,
      correlation = corfit$consensus.correlation
    )
  )
}


## Extract the p-values
pvals0_contrasts.st <- sapply(eb0_list_neurons, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts.st) <- rownames(sce.sacc.PB)

## Extract the tstats
t0_contrasts.st <- sapply(eb0_list_neurons, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts.st) <- rownames(sce.sacc.PB)

# Any signif
data.frame(
  'FDRsig.05' = colSums(apply(pvals0_contrasts.st, 2, p.adjust, 'fdr') < 0.05),
  'FDRsig.01' = colSums(apply(pvals0_contrasts.st, 2, p.adjust, 'fdr') < 0.01),
  'Pval10-6sig' = colSums(pvals0_contrasts.st < 1e-6),
  'Pval10-8sig' = colSums(pvals0_contrasts.st < 1e-8)
)
    #         FDRsig.05 FDRsig.01 Pval10.6sig Pval10.8sig
    # Excit.1        48        24          15           1
    # Excit.2        57        20          11           3
    # Excit.3        48        32          15           3
    # Excit.4       400       179          56          18
    # Inhib.1        68        45          14           8
    # Inhib.2        68        45          14           8

    #                         ^ go for printing top 20 of these


## (as before)
fdrs0_contrasts.st = apply(pvals0_contrasts.st, 2, p.adjust, "fdr")
rownames(fdrs0_contrasts.st) <- rownames(sce.sacc.PB)


markerList.PB.manual <- lapply(colnames(fdrs0_contrasts.st), function(x){
  rownames(fdrs0_contrasts.st)[fdrs0_contrasts.st[ ,x] < 0.01 & t0_contrasts.st[ ,x] > 0]
})
names(markerList.PB.manual) <- colnames(fdrs0_contrasts.st)
lengths(markerList.PB.manual)
# Excit.1 Excit.2 Excit.3 Excit.4 Inhib.1 Inhib.2
#      24      20      31      94      15      30

markerTs.fdr.01 <- lapply(colnames(fdrs0_contrasts.st), function(x){
  as.matrix(t0_contrasts.st[fdrs0_contrasts.st[ ,x] < 0.01 & t0_contrasts.st[ ,x] > 0, x])
})

names(markerTs.fdr.01) <- colnames(fdrs0_contrasts.st)

markerList.sorted <- lapply(markerTs.fdr.01, function(x){
  x[,1][order(x, decreasing=TRUE)]
})

genes2plot <- lapply(markerList.sorted, function(x){head(x, n=20)})


pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_sACC-n2_top20markers-SUBtypes_logExprs_Mar2020.pdf", height=7.5, width=9.5)
for(i in 1:length(genes2plot)){
  print(
    plotExpression(sce.sacc, exprs_values = "logcounts", features=c(names(genes2plot[[i]])),
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:10], length(genes2plot[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(genes2plot)[i], " top 20 markers"))
  )
}
dev.off()



