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
    ## same result if you used $sample or $donor (for Amyg, where it's confounded)
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
  rownames(fdrs0_contrasts)[fdrs0_contrasts[ ,x] < 0.05 & t0_contrasts[ ,x] > 0]
  # what if more stringent?
  #rownames(fdrs0_contrasts)[fdrs0_contrasts[ ,x] < 0.01 & t0_contrasts[ ,x] > 0]
})
names(markerList.PB.manual) <- colnames(fdrs0_contrasts)
lengths(markerList.PB.manual)


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



### Top markers to print / potentially test with RNA-scope === === === ===
load('/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_Amyg-n2_manualContrasts_MNTMar2020.rda',
     verbose=T)
    # eb_contrasts.amy.broad, eb_list.amy.broad, sce.amy.PB

    # follow chunk 'How does this compare to results of `findMarkers()`?' for fdr & t mats

# Take FDR < 0.005 (instead of < 0.001 as in some other regions)
markerList.PB.manual <- lapply(colnames(fdrs0_contrasts), function(x){
  rownames(fdrs0_contrasts)[fdrs0_contrasts[ ,x] < 0.005 & t0_contrasts[ ,x] > 0]
})
names(markerList.PB.manual) <- colnames(fdrs0_contrasts)
lengths(markerList.PB.manual)

    # FDR < 0.005
        # Astro Excit Inhib Micro Oligo   OPC
        #   258   170    49  1406   266    53
    
    # (FDR < 0.001)
        # Astro Excit Inhib Micro Oligo   OPC
        #   113    92    28   894   149    26

markerTs.fdr.005 <- lapply(colnames(fdrs0_contrasts), function(x){
  as.matrix(t0_contrasts[fdrs0_contrasts[ ,x] < 0.005 & t0_contrasts[ ,x] > 0, x])
})

names(markerTs.fdr.005) <- colnames(fdrs0_contrasts)

markerList.sorted <- lapply(markerTs.fdr.005, function(x){
  x[,1][order(x, decreasing=TRUE)]
})

genes2plot <- lapply(markerList.sorted, function(x){head(x, n=20)})


## Let's plot some expression of these to see how much are 'real' (not driven by outliers)
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo
    rm(chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo)

# As before, first drop "Ambig.lowNtrxts" (50 nuclei)
sce.amy <- sce.amy[ ,sce.amy$cellType != "Ambig.lowNtrxts"]
sce.amy$cellType <- droplevels(sce.amy$cellType)


#pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_Amyg-n2_top20markers_logExprs_Mar2020.pdf", height=7.5, width=9.5)
    ## 'Mar' iteration renamed with 'zold_' prefix - limited with too strict FDR cutoff such that top 20 pt-coding markers was < 20...
pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_Amyg-n2_top20markers_logExprs_Apr2020.pdf", height=7.5, width=9.5)
for(i in 1:length(genes2plot)){
  print(
    plotExpression(sce.amy, exprs_values = "logcounts", features=c(names(genes2plot[[i]])),
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

table(rownames(sce.amy) %in% gtf$gene_name)
    # FALSE  TRUE
    #    48 33490    - probably because of the `uniquify`
table(rowData(sce.amy)$Symbol %in% gtf$gene_name)
    #  TRUE
    # 33538

# Are they the same order?
table(rowData(sce.amy)$ID == gtf$gene_id) # all TRUE

table(!rowSums(assay(sce.amy, "counts"))==0)  # 28464     - good
keepVec <- !rowSums(assay(sce.amy, "counts"))==0

gtf <- gtf[keepVec, ]
# Then
table(gtf$gene_id == rowData(sce.amy.PB)$ID)  # all 28464 TRUE      - good

## Make pt-coding list
markerList.sorted.pt <- lapply(markerList.sorted, function(x){
  x[names(x) %in% gtf$gene_name[gtf$gene_biotype=="protein_coding"]]
})

lengths(markerList.sorted)
    # (above)

lengths(markerList.sorted.pt)
    # Astro Excit Inhib Micro Oligo   OPC
    #   154    63    31  1132   191    26



genes2plot.pt <- lapply(markerList.sorted.pt, function(x){head(x, n=20)})

# Plot these
#pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_Amyg-n2_top20markers_logExprs_pt-coding_Mar2020.pdf", height=7.5, width=9.5)
    ## 'Mar' iteration renamed with 'zold_' prefix - limited with too strict FDR cutoff such that top 20 pt-coding markers was < 20...
pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_Amyg-n2_top20markers_logExprs_pt-coding_Apr2020.pdf", height=7.5, width=9.5)
for(i in 1:length(genes2plot.pt)){
  print(
    plotExpression(sce.amy, exprs_values = "logcounts", features=c(names(genes2plot.pt[[i]])),
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
    # [1] "IGFN1"   "SLC14A1" "BBOX1"   "TMPRSS3" "MRVI1"   "NEUROG2" "LAMA1"
    # [8] "SIX5"    "YAP1"
    # 
    # $Excit
    # [1] "EMX1"   "OR14I1"
    # 
    # $Inhib
    # [1] "NKX2-1" "CYP1A1" "NMUR1"  "TRPM5"  "IBSP"   "LYPD6B" "NXPH2"  "ANO1"
    # [9] "CHRNA2" "VIP"    "PLSCR5" "DRD5"   "SP8"
    # 
    # $Micro
    # [1] "MYO1F"    "BIN2"     "CCL20"    "MSLN"     "CD300C"   "PTGS1"
    # [7] "PIK3AP1"  "SIGLEC10" "HLA-DOA"  "NCF2"     "PRDM11"   "TFEC"
    # [13] "FCGR3A"   "CASP1"
    # 
    # $Oligo
    # [1] "IFNA2"      "CLDND1"     "CTNNA3"     "QDPR"       "OPALIN"
    # [6] "GREM1"      "CD22"       "HHIP"       "CCP110"     "MAG"
    # [11] "LDB3"       "FXYD4"      "POU2AF1"    "AC026316.5"
    # 
    # $OPC
    # [1] "CPXM1"  "KCNG4"  "GDF6"   "GCOM1"  "CCKAR"  "CSPG4"  "ATP2C2" "NEU4"


# Write 'genes2plot's to a csv
names(genes2plot.pt) <- paste0(names(genes2plot.pt),"_pt")
        # * Only relevant where < 20 genes at a given FDR (but was being too strict in the first place):
        # genes2plot.pt.names <- sapply(genes2plot.pt, function(x){
        #                           sapply(1:20, function(i){ifelse(!is.na(x[i]), names(x)[i], NA)})
        #                         })
        # 
        # top20genes <- cbind(sapply(genes2plot, names), genes2plot.pt.names)
        # rownames(top20genes) <- NULL
top20genes <- cbind(sapply(genes2plot, names), sapply(genes2plot.pt, names))
top20genes <- top20genes[ ,sort(colnames(top20genes))]

write.csv(top20genes, file="tables/top20genesLists_Amyg-n2_cellTypes.csv")



### ========================== ###
### SINGLE-NUCLEUS-LEVEL TESTS ###
### ========================== ###


### Single-nucleus-level tests for cell-type-specific genes ================================

## Load SCE with new info
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo

table(sce.amy$cellType.split)

# First drop "Ambig.lowNtrxts" (50 nuclei)
sce.amy <- sce.amy[ ,sce.amy$cellType.split != "Ambig.lowNtrxts"]
sce.amy$cellType.split <- droplevels(sce.amy$cellType.split)

# Remove 0 genes across all nuclei
sce.amy <- sce.amy[!rowSums(assay(sce.amy, "counts"))==0, ]  # keeps same 28464 genes


## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.amy), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`


# Run pairwise t-tests
markers.amy.t.design <- findMarkers(sce.amy, groups=sce.amy$cellType.split,
                                    assay.type="logcounts", design=mod, test="t",
                                    direction="up", pval.type="all", full.stats=T)

sapply(markers.amy.t.design, function(x){table(x$FDR<0.05)})
    #      Astro Excit.1 Excit.2 Excit.3 Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5
    # FALSE 27751   28366   27714   27811   28433   28241   28134   28321   28215
    # TRUE    713      98     750     653      31     223     330     143     249
    # Micro Oligo   OPC
    # FALSE 27475 28030 28170
    # TRUE    989   434   294


## WMW: Blocking on donor (this test doesn't take 'design=' argument) ===
markers.amy.wilcox.block <- findMarkers(sce.amy, groups=sce.amy$cellType.split,
                                        assay.type="logcounts", block=sce.amy$donor, test="wilcox",
                                        direction="up", pval.type="all", full.stats=T)

# no warnings as in pan-brain analyses, but NO results of FDR<0.05...:
sapply(markers.amy.wilcox.block, function(x){table(x$FDR<0.05)})
    # Actually some decent results but 'Inhib.1' has none


## Binomial ===
markers.amy.binom.block <- findMarkers(sce.amy, groups=sce.amy$cellType.split,
                                       assay.type="logcounts", block=sce.amy$donor, test="binom",
                                       direction="up", pval.type="all", full.stats=T)

sapply(markers.amy.binom.block, function(x){table(x$FDR<0.05)})
# only a couple dozen hits for glia, only - disregard these

## Save all these for future reference
save(markers.amy.t.design, markers.amy.wilcox.block, #markers.amy.binom.block,
     file="rdas/markers-stats_Amyg-n2_findMarkers-SN-LEVEL_MNTMay2020.rda")


        # Btw - some have 0 p.value's and FDR's
        head(markers.amy.t.design[["Excit.2"]][ ,1:2])
        head(markers.amy.t.design[["Excit.3"]][ ,1:2])
        


# Print these to pngs
markerList.t <- lapply(markers.amy.t.design, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)

genes.top40.t <- lapply(markerList.t, function(x){head(x, n=40)})


#dir.create("pdfs/exploration/Amyg/")
for(i in names(genes.top40.t)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/Amyg/Amyg_t-sn-level_pairwise_top40markers-", i, "_logExprs_May2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.amy, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:12], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers: single-nucleus-level p.w. t-tests"))
  )
  dev.off()
}






### Cluster-vs-all single-nucleus-level iteration ======

## Load SCE with new info
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo

table(sce.amy$cellType.split)

# First drop "Ambig.lowNtrxts" (50 nuclei)
sce.amy <- sce.amy[ ,sce.amy$cellType.split != "Ambig.lowNtrxts"]
sce.amy$cellType.split <- droplevels(sce.amy$cellType.split)

# Remove 0 genes across all nuclei
sce.amy <- sce.amy[!rowSums(assay(sce.amy, "counts"))==0, ]  # keeps same 28464 genes


## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.amy), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

markers.amy.t.1vAll <- list()
for(i in levels(sce.amy$cellType.split)){
  # Make temporary contrast
  sce.amy$contrast <- ifelse(sce.amy$cellType.split==i, 1, 0)
  # Test cluster vs. all
  markers.amy.t.1vAll[[i]] <- findMarkers(sce.amy, groups=sce.amy$contrast,
                                            assay.type="logcounts", design=mod, test="t",
                                            direction="up", pval.type="all", full.stats=T)
}

    ## Then, temp set of stats to get the standardized logFC
    temp.1vAll <- list()
    for(i in levels(sce.amy$cellType.split)){
      # Make temporary contrast
      sce.amy$contrast <- ifelse(sce.amy$cellType.split==i, 1, 0)
      # Test cluster vs. all
      temp.1vAll[[i]] <- findMarkers(sce.amy, groups=sce.amy$contrast,
                                     assay.type="logcounts", design=mod, test="t",
                                     std.lfc=TRUE,
                                     direction="up", pval.type="all", full.stats=T)
    }


## For some reason all the results are in the second List entry (first is always empty)

# Replace that empty slot with the entry with the actul stats
markers.amy.t.1vAll <- lapply(markers.amy.t.1vAll, function(x){ x[[2]] })
# Same for that with std.lfc
temp.1vAll <- lapply(temp.1vAll, function(x){ x[[2]] })

# Now just pull from the 'stats.0' DataFrame column
markers.amy.t.1vAll <- lapply(markers.amy.t.1vAll, function(x){ x$stats.0 })
temp.1vAll <- lapply(temp.1vAll, function(x){ x$stats.0 })

# Re-name std.lfc column and add to the first result
for(i in names(temp.1vAll)){
  colnames(temp.1vAll[[i]])[1] <- "std.logFC"
  markers.amy.t.1vAll[[i]] <- cbind(markers.amy.t.1vAll[[i]], temp.1vAll[[i]]$std.logFC)
  # Oh the colname is kept weird
  colnames(markers.amy.t.1vAll[[i]])[4] <- "std.logFC"
  # Then re-organize
  markers.amy.t.1vAll[[i]] <- markers.amy.t.1vAll[[i]][ ,c("logFC","std.logFC","log.p.value","log.FDR")]
}


## Let's save this along with the previous pairwise results
save(markers.amy.t.1vAll, markers.amy.t.design, markers.amy.wilcox.block,
     file="rdas/markers-stats_Amyg-n2_findMarkers-SN-LEVEL_MNTMay2020.rda")




## Print these to pngs
markerList.t.1vAll <- lapply(markers.amy.t.1vAll, function(x){
  rownames(x)[x[ ,"log.FDR"] < log10(0.001)]
  }
)
genes.top40.t.1vAll <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t.1vAll)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/Amyg/Amyg_t-sn-level_1vALL_top40markers-",gsub(":",".",i),"_logExprs_May2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.amy, exprs_values = "logcounts", features=genes.top40.t.1vAll[[i]],
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:12], length(genes.top40.t.1vAll[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers: single-nucleus-level p.w. t-tests, cluster-vs-all"))
  )
  dev.off()
}






## How do these top 40 intersect? ===
sapply(names(genes.top40.t), function(c){
  length(intersect(genes.top40.t[[c]],
                   genes.top40.t.1vAll[[c]]))
})
    #  Astro Excit.1 Excit.2 Excit.3 Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5   Micro
    #     35       5      37      32       3      21      32      15      24      37
    #  Oligo     OPC
    #     30      27



## Write these top 40 lists to a csv
names(markerList.t) <- paste0(names(markerList.t),"_pw")
names(markerList.t.1vAll) <- paste0(names(markerList.t.1vAll),"_1vAll")

# PW result for "Inhib.1" doesn't have 40 markers:
markerList.t[["Inhib.1_pw"]] <- c(markerList.t[["Inhib.1_pw"]], rep("",9))

top40genes <- cbind(sapply(markerList.t, function(x) head(x, n=40)),
                    sapply(markerList.t.1vAll, function(y) head(y, n=40)))
top40genes <- top40genes[ ,sort(colnames(top40genes))]

write.csv(top40genes, file="tables/top40genesLists_Amyg-n2_cellType.split_SN-LEVEL-tests_May2020.csv",
          row.names=FALSE)



## 15May2020 for AnJa - print t's and FDRs ===
# Bring in human stats; create t's
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_Amyg-n2_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.amy.t.1vAll, markers.amy.t.design, markers.amy.wilcox.block
    rm(markers.amy.t.design, markers.amy.wilcox.block)
    
# Need to add t's with N nuclei used in constrasts
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    #sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo
    rm(chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy,ref.sampleInfo)

# First drop "ambig.lowNtrxts" (93 nuclei)
sce.amy <- sce.amy[ ,sce.amy$cellType.split != "Ambig.lowNtrxts"]
sce.amy$cellType.split <- droplevels(sce.amy$cellType.split)

## As above, calculate and add t-statistic (= std.logFC * sqrt(N)) from contrasts
#      and fix row order to the first entry "Astro"
fixTo <- rownames(markers.amy.t.1vAll[["Astro"]])
for(s in names(markers.amy.t.1vAll)){
  markers.amy.t.1vAll[[s]]$t.stat <- markers.amy.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.amy))
  markers.amy.t.1vAll[[s]] <- markers.amy.t.1vAll[[s]][fixTo, ]
}
# Pull out the t's
ts.amy <- sapply(markers.amy.t.1vAll, function(x){x$t.stat})
rownames(ts.amy) <- fixTo
# Change to EnsemblID
rownames(ts.amy) <- rowData(sce.amy)$ID[match(rownames(ts.amy), rownames(sce.amy))]

# logFDRs
logFDRs.amy <- sapply(markers.amy.t.1vAll, function(x){x$log.FDR})
rownames(logFDRs.amy) <- fixTo
# EnsemblID
rownames(logFDRs.amy) <- rownames(ts.amy)

save(ts.amy, logFDRs.amy, file="rdas/zForAnJa_AMY-sn-level-markerStats_MNT.rda")

apply(logFDRs.amy, 2, function(x) {table(x<log10(1e-6))})
    # Astro Excit.1 Excit.2 Excit.3 Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5
    # FALSE 23194   19604   25841   25201   21979   23391   26181   26787   24329
    # TRUE   5270    8860    2623    3263    6485    5073    2283    1677    4135
    # Micro Oligo   OPC
    # FALSE 24520 25698 24702
    # TRUE   3944  2766  3762



### For fig: Plot some top markers in vlnplot array (12Jun2020) === === === ===
  #     * had been done on local machine, but realized chunk isn't in script
  #       - updating 31Aug2020

load("rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
load("rdas/markers-stats_Amyg-n2_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.amy.t.1vAll, markers.amy.t.design, markers.amy.wilcox.block
    # Focus on the pairwise result (".design") bc more specific
    rm(markers.amy.t.1vAll, markers.amy.wilcox.block)

# First drop "ambig.lowNtrxts" (50 nuclei)
sce.amy <- sce.amy[ ,sce.amy$cellType.split != "Ambig.lowNtrxts"]
sce.amy$cellType.split <- droplevels(sce.amy$cellType.split)


# Take top two for broad glia
topToPrint <- as.data.frame(sapply(markers.amy.t.design, function(x) {head(rownames(x),n=2)}))
topToPrint <- c(topToPrint$Astro, c("NPTX1", "SLC30A3", # Excit.1
                                    "SLC17A6", "SOX4", "SOX11", #Excit.2
                                    "MCHR2", "CDH22", # Excit.3
                                    "CCK", "CALB2", "KIT", # Inhib.1/2/4
                                    "CNTNAP3", "CNTNAP3B", "CALB1", # Inhib.3
                                    "NPFFR2", "TLL1"), # Inhib.5
                topToPrint$Micro, topToPrint$Oligo, topToPrint$OPC)

table(topToPrint %in% rownames(sce.amy)) # good

# With top 2 per glial
pdf("pdfs/pubFigures/Amyg_topMarkers-ARRAY_logExprs_Jun2020_v1.pdf", height=17, width=4)
print(
  plotExpression(sce.amy, exprs_values = "logcounts", features=topToPrint,
                 x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=1,
                 add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                              geom = "crossbar", width = 0.3,
                                              colour=rep(tableau20[1:12], length(topToPrint))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12), plot.title = element_text(size = 25)) +  
    ggtitle(label="AMY top marker array (with glial)")
)
dev.off()

# Neuronal markers only (highlighted in paper)
topToPrint <- c("NRN1", "NPTX1", "SLC30A3", # Excit.1
  "SLC17A6", "VCAN", #Excit.2
  "MCHR2", "RBFOX3", # Excit.3
  "CCK", "CALB2", "KIT", # Inhib.1/2/4
  "CRH", "CALB1", # Inhib.3
  "NPFFR2", "TLL1")

pdf("pdfs/pubFigures/Amyg_topMarkers-ARRAY_logExprs_Jun2020_v2.pdf", height=11, width=4.5)
print(
  plotExpression(sce.amy, exprs_values = "logcounts", features=topToPrint,
                 x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=1,
                 add_legend=F, scales="free_y") + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                              geom = "crossbar", width = 0.3,
                                              colour=rep(tableau20[1:12], length(topToPrint))) +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 13),
          axis.text.y = element_text(size = 6.5),
          axis.title.y = element_text(angle = 90, size = 16),
          plot.title = element_text(size = 20),
          panel.grid.major=element_line(colour="grey95", size=0.8),
          panel.grid.minor=element_line(colour="grey95", size=0.4))# +  
    #ggtitle(label="AMY top marker array (neuronal only)")
)
dev.off()





### MNT add 18Nov2020 =================================
# -> What if add param/requirement that for any given subcluster, median expression has to > 0?
load("rdas/markers-stats_Amyg-n2_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.amy.t.1vAll, markers.amy.t.design, markers.amy.wilcox.block

## Load SCE 
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo

table(sce.amy$cellType.split)

# First drop "Ambig.lowNtrxts" (50 nuclei)
sce.amy <- sce.amy[ ,sce.amy$cellType.split != "Ambig.lowNtrxts"]
sce.amy$cellType.split <- droplevels(sce.amy$cellType.split)

# Remove 0 genes across all nuclei
sce.amy <- sce.amy[!rowSums(assay(sce.amy, "counts"))==0, ]


## Make list of Boolean param / cell subtype ===
cellSubtype.idx <- splitit(sce.amy$cellType.split)
medianNon0.idx <- lapply(cellSubtype.idx, function(x){
  apply(as.matrix(assay(sce.amy, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

lengths(medianNon0.idx)
sapply(medianNon0.idx, head)

# Add these to the stats for each set of markers
for(i in names(markers.amy.t.1vAll)){
  markers.amy.t.1vAll[[i]] <- cbind(markers.amy.t.1vAll[[i]],
                                     medianNon0.idx[[i]][match(rownames(markers.amy.t.1vAll[[i]]),
                                                               names(medianNon0.idx[[i]]))])
  colnames(markers.amy.t.1vAll[[i]])[5] <- "non0median"
}


## Use these restrictions to print (to png) a refined top 40, as before ===
markerList.t.1vAll <- lapply(markers.amy.t.1vAll, function(x){
  rownames(x)[x$log.FDR < log10(0.000001) & x$non0median==TRUE]
  }
)
    # lengths(markerList.t.1vAll)     # ( **without $non0median==TRUE restriction )
        #  Astro Excit.1 Excit.2 Excit.3 Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5   Micro
        #   5270    8860    2623    3263    6485    5073    2283    1677    4135    3944
        #  Oligo     OPC
        #   2766    3762

lengths(markerList.t.1vAll)
    #  Astro Excit.1 Excit.2 Excit.3 Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5   Micro
    #   1326    4344     948    1795    3424    3031    1647    1121    2500     823
    #  Oligo     OPC
    #    976    1361

genes.top40.t <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/Amyg/Amyg_t-sn-level_1vALL_top40markers-REFINED-",gsub(":",".",i),"_logExprs_Nov2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.amy, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:12], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers, refined: single-nucleus-level p.w. t-tests, cluster-vs-all"))
  )
  dev.off()
}



## Do the same with the pairwise result ('markers.amy.t.design') === === ===
# Add these to the stats for each set of markers
for(i in names(markers.amy.t.design)){
  markers.amy.t.design[[i]] <- cbind(markers.amy.t.design[[i]],
                                      medianNon0.idx[[i]][match(rownames(markers.amy.t.design[[i]]),
                                                                names(medianNon0.idx[[i]]))])
  colnames(markers.amy.t.design[[i]])[14] <- "non0median"
}

markerList.t <- lapply(markers.amy.t.design, function(x){
  rownames(x)[x$FDR < 0.05 & x$non0median==TRUE]
  }
)
    # lengths(markerList.t)     # ( **without $non0median==TRUE restriction )
        #  Astro Excit.1 Excit.2 Excit.3 Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5   Micro
        #    713      98     750     653      31     223     330     143     249     989
        #  Oligo     OPC
        #    434     294

lengths(markerList.t)
    #  Astro Excit.1 Excit.2 Excit.3 Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5   Micro
    #    403      44     226     350       7     142     182      34      86     405
    #  Oligo     OPC
    #    371     186


genes.top40.t <- lapply(markerList.t, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/Amyg/Amyg_t-sn-level_pairwise_top40markers-REFINED-", i, "_logExprs_Nov2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.amy, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:12], length(genes.top40.t[[i]]))) +
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

write.csv(top40genes, file="tables/top40genesLists-REFINED_Amyg-n2_cellType.split_Nov2020.csv",
          row.names=FALSE)












### Aside - difference b/tw dropping 'Ambig.lowNtrxts' before or after PB'ing ===========
  # MNT 05Mar2020
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)


# First drop "Ambig.lowNtrxts" (50 nuclei)
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