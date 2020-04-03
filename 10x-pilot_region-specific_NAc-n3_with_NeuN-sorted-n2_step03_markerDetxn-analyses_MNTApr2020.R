### MNT 10x snRNA-seq workflow: step 03 - marker detection
###   **Region-specific analyses**
###     - (3x) NAc samples from: Br5161 & Br5212 & Br5287
###     - (2x) NeuN-sorted samples from: Br5207 & Br5182
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


    ## MNT comment 02Apr2020 === === === === === ===
    ## - run `findMarkers()` later - go for manual limma modeling for now
    ## === === === === === === === === === === === =


load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda",
     verbose=T)
    # sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all, ref.sampleInfo

###### Direct limma approach ####
#################################
table(sce.nac.all$cellType.split)
table(sce.nac.all$cellType.split, sce.nac.all$sample) 
    # nac.5161 nac.5212 nac.5287 nac.neun.5182 nac.neun.5207
    # ambig.lowNtrxts       19       42       22             7             3
    # Astro                149      384       12             0             0
    # Inhib.1                1        3        0            16             5
    # Inhib.2                1        1        1            42            11
    # Inhib.3                7        7        9            86           167
    # Inhib.4                9        8        4           104            58
    # Micro                 72       72       37             0             0
    # MSN.broad              0      266        0             0             0
    # MSN.D1.1               2        0        0           117            13
    # MSN.D1.2              10        3        0           285             3
    # MSN.D1.3              17        8        6           369           319
    # MSN.D1.4             178        2       72          1505          1829
    # MSN.D2.1               9        6        3           134           148
    # MSN.D2.2              41       14        5          1602          1870
    # Oligo               1454      854      499             0             0
    # OPC                   98      104       37             0             0

    #           - modeling might end up being a little weird....

# First drop "ambig.lowNtrxts" (93 nuclei)
sce.nac.all <- sce.nac.all[ ,sce.nac.all$cellType.split != "ambig.lowNtrxts"]
sce.nac.all$cellType.split <- droplevels(sce.nac.all$cellType.split)

# Then make the pseudo-bulked SCE
sce.nac.all.PB <- aggregateAcrossCells(sce.nac.all, ids=paste0(sce.nac.all$sample,":",sce.nac.all$cellType.split),
                                        use_exprs_values="counts")
    ## of 33538 x 59 dims

# Clean up colData
colData(sce.nac.all.PB) <- colData(sce.nac.all.PB)[ ,c(13:17,19:22)]

# Drop genes with all 0's
sce.nac.all.PB <- sce.nac.all.PB[!rowSums(assay(sce.nac.all.PB, "counts"))==0, ]
    ## keeps 29236 genes

# Remove stored `sizeFactors()` because this will mess you up
#     * Also, to be safe, can always provide manually-computed SFs:
sizeFactors(sce.nac.all.PB) <- NULL
LSFvec <- librarySizeFactors(sce.nac.all.PB)
sce.nac.all.PB <- logNormCounts(sce.nac.all.PB, size_factors=LSFvec)

## Extract the count data
mat <- assays(sce.nac.all.PB)$logcounts

## Build a group model * FOR THIS REGION ADDING 'processDate'
mod <- with(colData(sce.nac.all.PB), model.matrix(~ 0 + cellType.split + processDate))
colnames(mod) <- gsub('cellType.split', '', colnames(mod))

corfit <- duplicateCorrelation(mat, mod, block = sce.nac.all.PB$donor)
corfit$consensus.correlation
    ## [1] 0.02344084

## (other scripts have pairwise tests -- too many levels to test that here)


## Then each cellType vs the rest
cellType_idx <- splitit(sce.nac.all.PB$cellType.split)

eb0_list <- lapply(cellType_idx, function(x) {
  res <- rep(0, ncol(sce.nac.all.PB))
  res[x] <- 1
  m <- model.matrix(~ res + sce.nac.all.PB$processDate)
  eBayes(
    lmFit(
      mat,
      design = m,
      block = sce.nac.all.PB$donor,
      correlation = corfit$consensus.correlation
    )
  )
})

## Extract the p-values
pvals0_contrasts <- sapply(eb0_list, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts) <- rownames(sce.nac.all.PB)

fdrs0_contrasts = apply(pvals0_contrasts, 2, p.adjust, "fdr")

## Extract the tstats
t0_contrasts_cell <- sapply(eb0_list, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts_cell) <- rownames(sce.nac.all.PB)


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
    #           FDRsig Pval10.6sig Pval10.8sig
    # Astro       1713         390         235
    # Inhib.1      351          28           9
    # Inhib.2     4827         596         237
    # Inhib.3      391          85          42
    # Inhib.4      443         104          49
    # Micro       3351         972         688
    # MSN.broad    152          81          13
    # MSN.D1.1     101          13           2
    # MSN.D1.2     190          17           2
    # MSN.D1.3      88          20          11
    # MSN.D1.4     493          53          12
    # MSN.D2.1     178          27          10
    # MSN.D2.2     198          19           5
    # Oligo       1732         420         251
    # OPC          860         181         113


## With t > 0
    #           FDRsig Pval10.6sig Pval10.8sig
    # Astro       1496         371         228
    # Inhib.1       97          13           6
    # Inhib.2      220           5           0
    # Inhib.3      351          78          39
    # Inhib.4      400         101          49
    # Micro       2201         692         514
    # MSN.broad    152          81          13
    # MSN.D1.1      88          12           2
    # MSN.D1.2      65          12           1
    # MSN.D1.3      88          20          11
    # MSN.D1.4     485          53          12
    # MSN.D2.1     132          24          10
    # MSN.D2.2     194          19           5
    # Oligo       1402         370         227
    # OPC          804         179         112




## Save for later
#eb_contrasts.dlpfc.broad <- eb_contrasts
eb_list.nac.all <- eb0_list
corfit.nac.all <- corfit

save(eb_list.nac.all, sce.nac.all.PB, corfit.nac.all,
     file = 'rdas/markers-stats_NAc_all-n5_manualContrasts_MNTApr2020.rda')




### Top markers to print / potentially test with RNA-scope === === === ===
load('/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_NAc_all-n5_manualContrasts_MNTApr2020.rda',
     verbose=T)
    # eb_list.nac.all, sce.nac.all.PB, corfit.nac.all

# (set up p0/fdrs0 & t0 lists, above)

# Let's take those with fdr < 0.01
markerList.PB.manual <- lapply(colnames(fdrs0_contrasts), function(x){
  rownames(fdrs0_contrasts)[fdrs0_contrasts[ ,x] < 0.01 & t0_contrasts_cell[ ,x] > 0]
})
names(markerList.PB.manual) <- colnames(fdrs0_contrasts)
lengths(markerList.PB.manual)
    #    Astro   Inhib.1   Inhib.2   Inhib.3   Inhib.4     Micro MSN.broad  MSN.D1.1
    #     1010        25        66       191       230      1583       128        30
    # MSN.D1.2  MSN.D1.3  MSN.D1.4  MSN.D2.1  MSN.D2.2     Oligo       OPC
    #       23        40       200        53        57       876       492

markerTs.fdr.01 <- lapply(colnames(fdrs0_contrasts), function(x){
  as.matrix(t0_contrasts_cell[fdrs0_contrasts[ ,x] < 0.01 & t0_contrasts_cell[ ,x] > 0, x])
})

names(markerTs.fdr.01) <- colnames(fdrs0_contrasts)

markerList.sorted <- lapply(markerTs.fdr.01, function(x){
  x[,1][order(x, decreasing=TRUE)]
})

genes2plot <- lapply(markerList.sorted, function(x){head(x, n=20)})


## Let's plot some expression of these to see how much are 'real' (not driven by outliers)
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda",
     verbose=T)
    # sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all, ref.sampleInfo
    rm(chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all, ref.sampleInfo)

# As before, first drop "Ambig.lowNtrxts" (168 nuclei)
sce.nac.all <- sce.nac.all[ ,sce.nac.all$cellType != "Ambig.lowNtrxts"]
sce.nac.all$cellType <- droplevels(sce.nac.all$cellType)


pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_NAc-ALL-n5_top20markers_logExprs_Apr2020.pdf", height=7.5, width=9.5)
for(i in 1:length(genes2plot)){
  print(
    plotExpression(sce.nac.all, exprs_values = "logcounts", features=c(names(genes2plot[[i]])),
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:15], length(genes2plot[[i]]))) +
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

table(rownames(sce.nac.all) %in% gtf$gene_name)
    # FALSE  TRUE
    #    48 33490    - probably because of the `uniquify`
table(rowData(sce.nac.all)$Symbol %in% gtf$gene_name)
    #  TRUE
    # 33538

# Are they the same order?
table(rowData(sce.nac.all)$ID == gtf$gene_id) # all TRUE

table(!rowSums(assay(sce.nac.all, "counts"))==0)  # 29236     - good
keepVec <- !rowSums(assay(sce.nac.all, "counts"))==0

gtf <- gtf[keepVec, ]
# Then
table(gtf$gene_id == rowData(sce.nac.all.PB)$ID)  # all 29236 TRUE      - good

## Make pt-coding list
markerList.sorted.pt <- lapply(markerList.sorted, function(x){
  x[names(x) %in% gtf$gene_name[gtf$gene_biotype=="protein_coding"]]
})

lengths(markerList.sorted)

lengths(markerList.sorted.pt)
    #    Astro   Inhib.1   Inhib.2   Inhib.3   Inhib.4     Micro MSN.broad  MSN.D1.1
    #      555        18        56       107       142      1138        48        16
    # MSN.D1.2  MSN.D1.3  MSN.D1.4  MSN.D2.1  MSN.D2.2     Oligo       OPC
    #       15        15        55        19        13       463       264



genes2plot.pt <- lapply(markerList.sorted.pt, function(x){head(x, n=20)})

# Plot these
pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_NAc-ALL-n5_top20markers_logExprs_pt-coding_Apr2020.pdf", height=7.5, width=9.5)
for(i in 1:length(genes2plot.pt)){
  print(
    plotExpression(sce.nac.all, exprs_values = "logcounts", features=c(names(genes2plot.pt[[i]])),
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:15], length(genes2plot.pt[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(genes2plot.pt)[i], " top 20 protein-coding markers"))
  )
}
dev.off()


## How much they intersect with the top protein-coding-agnostic set?
sapply(names(genes2plot), function(x){intersect(names(genes2plot[[x]]), names(genes2plot.pt[[x]]))})
    #



