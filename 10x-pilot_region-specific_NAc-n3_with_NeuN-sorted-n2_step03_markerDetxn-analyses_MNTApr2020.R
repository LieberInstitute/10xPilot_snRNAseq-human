### MNT 10x snRNA-seq workflow: step 03 - marker detection
###   **Region-specific analyses**
###     - (3x) NAc samples from: Br5161 & Br5212 & Br5287
###     - (2x) NeuN-sorted samples from: Br5207 & Br5182
###   ** This final iteration includes the 5212-specific 'MSN.broad' resolved
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
table(sce.nac.all$cellType.final)
table(sce.nac.all$cellType.final, sce.nac.all$sample) 
    #                 nac.5161 nac.5212 nac.5287 nac.neun.5182 nac.neun.5207
    # ambig.lowNtrxts       19       42       22             7             3
    # Astro                149      384       12             0             0
    # Inhib.1                1        3        0            16             5
    # Inhib.2                1        1        1            42            11
    # Inhib.3                7        7        9            86           167
    # Inhib.4                9        8        4           104            58
    # Micro                 72       72       37             0             0
    # MSN.D1.1               2        0        0           117            13
    # MSN.D1.2              10        3        0           285             3
    # MSN.D1.3              17        8        6           369           319
    # MSN.D1.4             178      169       72          1505          1829
    # MSN.D2.1               9        6        3           134           148
    # MSN.D2.2              41      113        5          1602          1870
    # Oligo               1454      854      499             0             0
    # OPC                   98      104       37             0             0

# First drop "ambig.lowNtrxts" (93 nuclei)
sce.nac.all <- sce.nac.all[ ,sce.nac.all$cellType.final != "ambig.lowNtrxts"]
sce.nac.all$cellType.final <- droplevels(sce.nac.all$cellType.final)

# Then make the pseudo-bulked SCE
sce.nac.all.PB <- aggregateAcrossCells(sce.nac.all, ids=paste0(sce.nac.all$sample,":",sce.nac.all$cellType.final),
                                        use_exprs_values="counts")
    ## of 33538 x 58 dims

# Clean up colData
colData(sce.nac.all.PB) <- colData(sce.nac.all.PB)[ ,c(13:17,19:25)]

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
mod <- with(colData(sce.nac.all.PB), model.matrix(~ 0 + cellType.final + processDate))
colnames(mod) <- gsub('cellType.final', '', colnames(mod))

corfit <- duplicateCorrelation(mat, mod, block = sce.nac.all.PB$donor)
corfit$consensus.correlation
    ## [1] 0.01500338

## (other scripts have pairwise tests -- too many levels to test that here)


## Then each cellType vs the rest
cellType_idx <- splitit(sce.nac.all.PB$cellType.final)

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
    #         FDRsig Pval10.6sig Pval10.8sig
    # Astro      1670         380         224
    # Inhib.1     387          27           9
    # Inhib.2    4983         604         243
    # Inhib.3     351          78          41
    # Inhib.4     430          95          48
    # Micro      3356         966         675
    # MSN.D1.1     93          13           2
    # MSN.D1.2    188          16           0
    # MSN.D1.3     79          18          11
    # MSN.D1.4   1029         128          47
    # MSN.D2.1    173          28          10
    # MSN.D2.2    307          32           9
    # Oligo      1707         412         246
    # OPC         822         172         112


## With t > 0
    #          FDRsig Pval10.6sig Pval10.8sig
    # Astro      1462         361         218
    # Inhib.1      99          13           6
    # Inhib.2     215           3           0
    # Inhib.3     314          71          38
    # Inhib.4     386          91          48
    # Micro      2171         684         506
    # MSN.D1.1     80          12           2
    # MSN.D1.2     58          10           0
    # MSN.D1.3     79          18          11
    # MSN.D1.4   1027         128          47
    # MSN.D2.1    122          23           9
    # MSN.D2.2    300          32           9
    # Oligo      1384         360         221
    # OPC         771         170         111




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

# As above:
    ## Extract the p-values
    pvals0_contrasts <- sapply(eb_list.nac.all, function(x) {
      x$p.value[, 2, drop = FALSE]
    })
    rownames(pvals0_contrasts) <- rownames(sce.nac.all.PB)
    
    fdrs0_contrasts = apply(pvals0_contrasts, 2, p.adjust, "fdr")
    
    ## Extract the tstats
    t0_contrasts_cell <- sapply(eb_list.nac.all, function(x) {
      x$t[, 2, drop = FALSE]
    })
    rownames(t0_contrasts_cell) <- rownames(sce.nac.all.PB)

    
# Let's take those with fdr < 0.05
markerList.PB.manual <- lapply(colnames(fdrs0_contrasts), function(x){
  rownames(fdrs0_contrasts)[fdrs0_contrasts[ ,x] < 0.05 & t0_contrasts_cell[ ,x] > 0]
})
names(markerList.PB.manual) <- colnames(fdrs0_contrasts)
lengths(markerList.PB.manual)
    #    Astro  Inhib.1  Inhib.2  Inhib.3  Inhib.4    Micro MSN.D1.1 MSN.D1.2
    #     1462       99      215      314      386     2171       80       58
    # MSN.D1.3 MSN.D1.4 MSN.D2.1 MSN.D2.2    Oligo      OPC
    #       79     1027      122      300     1384      771

markerTs.fdr.05 <- lapply(colnames(fdrs0_contrasts), function(x){
  as.matrix(t0_contrasts_cell[fdrs0_contrasts[ ,x] < 0.05 & t0_contrasts_cell[ ,x] > 0, x])
})

names(markerTs.fdr.05) <- colnames(fdrs0_contrasts)

markerList.sorted <- lapply(markerTs.fdr.05, function(x){
  x[,1][order(x, decreasing=TRUE)]
})

genes2plot <- lapply(markerList.sorted, function(x){head(x, n=20)})

    ## MNT addition: print more for c("Inhib.2","MSN.D1.3", "MSN.D1.4", "MSN.D2.2")
    genes2plot.more <- lapply(markerList.sorted, function(x){head(x, n=40)})
    genes2plot.more <- list(genes2plot.more[["Inhib.2"]], genes2plot.more[["MSN.D1.3"]],
                            genes2plot.more[["MSN.D1.4"]], genes2plot.more[["MSN.D2.2"]])
    names(genes2plot.more) <- c("Inhib.2","MSN.D1.3", "MSN.D1.4", "MSN.D2.2")

## Let's plot some expression of these to see how much are 'real' (not driven by outliers)
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda",
     verbose=T)
    # sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all, ref.sampleInfo
    rm(chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all, ref.sampleInfo)

# As before, first drop "ambig.lowNtrxts" (93 nuclei)
sce.nac.all <- sce.nac.all[ ,sce.nac.all$cellType.final != "ambig.lowNtrxts"]
sce.nac.all$cellType.final <- droplevels(sce.nac.all$cellType.final)


pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_NAc-ALL-n5_top20markers_logExprs_Apr2020.pdf", height=7.5, width=9.5)
for(i in 1:length(genes2plot)){
  print(
    plotExpression(sce.nac.all, exprs_values = "logcounts", features=c(names(genes2plot[[i]])),
                   x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:14], length(genes2plot[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(genes2plot)[i], " top 20 markers"))
  )
}
dev.off()


    ## 'genes2plot.more' (added 10Apr2020):
    library(grDevices)
    for(i in names(genes2plot.more)){
    png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_NAc-ALL-n5_top40markers-",i,"_logExprs_Apr2020.png"), height=1900, width=1200)
      print(
        plotExpression(sce.nac.all, exprs_values = "logcounts", features=c(names(genes2plot.more[[i]])),
                       x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7, ncol=5,
                       add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                    geom = "crossbar", width = 0.3,
                                                    colour=rep(tableau20[1:14], length(genes2plot.more[[i]]))) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
          ggtitle(label=paste0(i, " top 40 markers"))
      )
    dev.off()
    }
    
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
    #    Astro  Inhib.1  Inhib.2  Inhib.3  Inhib.4    Micro MSN.D1.1 MSN.D1.2
    #      822       59      156      178      243     1583       28       32
    # MSN.D1.3 MSN.D1.4 MSN.D2.1 MSN.D2.2    Oligo      OPC
    #       36      324       52      103      740      409



genes2plot.pt <- lapply(markerList.sorted.pt, function(x){head(x, n=20)})

    ## MNT addition: print more for c("Inhib.2","MSN.D1.3", "MSN.D1.4", "MSN.D2.2")
    genes2plot.pt.more <- lapply(markerList.sorted.pt, function(x){head(x, n=40)})
    genes2plot.pt.more <- list(genes2plot.pt.more[["Inhib.2"]], genes2plot.pt.more[["MSN.D1.3"]],
                            genes2plot.pt.more[["MSN.D1.4"]], genes2plot.pt.more[["MSN.D2.2"]])
    names(genes2plot.pt.more) <- c("Inhib.2","MSN.D1.3", "MSN.D1.4", "MSN.D2.2")

# Plot these
pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_NAc-ALL-n5_top20markers_logExprs_pt-coding_Apr2020.pdf", height=7.5, width=9.5)
for(i in 1:length(genes2plot.pt)){
  print(
    plotExpression(sce.nac.all, exprs_values = "logcounts", features=c(names(genes2plot.pt[[i]])),
                   x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:14], length(genes2plot.pt[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(names(genes2plot.pt)[i], " top 20 protein-coding markers"))
  )
}
dev.off()

    ## 'genes2plot.pt.more' (added 10Apr2020):
    for(i in names(genes2plot.pt.more)){
      png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_NAc-ALL-n5_top40markers-",i,"_logExprs_pt-coding_Apr2020.png"), height=1900, width=1200)
      print(
        plotExpression(sce.nac.all, exprs_values = "logcounts", features=c(names(genes2plot.pt.more[[i]])),
                       x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7, ncol=5,
                       add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                    geom = "crossbar", width = 0.3,
                                                    colour=rep(tableau20[1:14], length(genes2plot.pt.more[[i]]))) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
          ggtitle(label=paste0(i, " top 40 protein-coding markers"))
      )
      dev.off()
    }



## How much they intersect with the top protein-coding-agnostic set?
sapply(names(genes2plot), function(x){intersect(names(genes2plot[[x]]), names(genes2plot.pt[[x]]))})
    #$Astro
      # [1] "FEZF1"   "KCNJ8"   "RASL12"  "TFAP2C"  "GLI1"    "PRODH"   "S1PR1"
      # [8] "TMPRSS3" "SLC2A4"  "IGFN1"   "IL33"
      # 
    # $Inhib.1
      # [1] "NPY"        "PNOC"       "ST6GALNAC2" "SST"        "HDAC4"
      # [6] "FBN2"       "KRTAP5-10"  "DDC"        "NOS1"       "DRD5"
      # [11] "SLC27A6"    "GPR156"     "TLL2"       "PCDH18"     "GDPD2"
      # 
    # $Inhib.2
      # [1] "NDUFB6" "PPP3CB" "SLC5A7" "SP8"    "BMP3"   "ACBD6"  "ACAT1"  "FCRL4"
      # [9] "VIP"    "SPTLC1" "FRMD7"  "CDC5L"  "ECEL1"  "TTC37"  "CDK14"  "TNS4"
      # [17] "KCNH5"  "ZRANB2"
      # 
    # $Inhib.3
      # [1] "EXOC1L" "LMO7DN" "SP5"    "OR2I1P" "PTHLH"  "PLSCR5" "KMO"    "DOK7"
      # 
    # $Inhib.4
      # [1] "CER1"    "TRH"     "RORC"    "CHRNA3"  "NPR3"    "COL2A1"  "GLP1R"
      # [8] "CBLN1"   "TAC3"    "COL13A1" "RERGL"
      # 
    # $Micro
      # [1] "SUCNR1"   "FCGR2B"   "HHEX"     "MPEG1"    "LILRA4"   "SERPINA1"
      # [7] "TAL1"     "VSIG4"    "GAPT"     "GIMAP7"   "NCF4"     "CEBPE"
      # [13] "CD300C"   "CCR1"     "LILRB2"   "RGS18"
      # 
    # $MSN.D1.1
      # [1] "CAPSL"    "C15orf62" "HBD"      "CCIN"     "SHISA8"   "SCGB1D4"  "FOXD4L4"
      # [8] "NPFFR2"   "GABRQ"
      # 
    # $MSN.D1.2
      # [1] "TMEM131" "EBF1"    "UMODL1"  "PCDHB1"  "DIO3"    "DWORF"   "TAF1L"
      # [8] "PDCD1"   "MAT1A"   "NCOA1"   "ASPG"    "TTC34"   "GIMAP5"
      # 
    # $MSN.D1.3
      # [1] "XDH"     "GPR26"   "LRRC55"  "RXFP1"   "LECT2"   "CCDC172"
      # 
    # $MSN.D1.4
      # [1] "OR51E2"  "IL17C"   "DLX4"    "CAVIN3"  "MINDY4B" "MCIDAS"  "MSLNL"
      # 
    # $MSN.D2.1
      # [1] "TH"      "ATP10A"  "C3orf52" "TTN"     "MROH2A"
      # 
    # $MSN.D2.2
      # [1] "SDR16C5"  "TMPRSS13"
      # 
    # $Oligo
      # [1] "FFAR1"      "FCRLA"      "NGFR"       "TMEM88B"    "IFNA2"
      # [6] "GPIHBP1"    "AC034102.2" "HOXD1"
      # 
    # $OPC
      # [1] "CPXM1"   "DCAF4L2" "KCNG4"   "TM4SF1"  "KLF17"   "NR0B1"   "FMO3"
      # [8] "TIMP4"   "GPR17"   "CYP2A13" "EMILIN3" "CSPG4"   "BGN"     "GDF6"
      # [15] "HMX1"


# Write 'genes2plot's to a csv
names(genes2plot.pt) <- paste0(names(genes2plot.pt),"_pt")
top20genes <- cbind(sapply(genes2plot, names), sapply(genes2plot.pt, names))
top20genes <- top20genes[ ,sort(colnames(top20genes))]

write.csv(top20genes, file="tables/top20genesLists_NAc-n5_cellType.final.csv")


png("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/regionSpecific_NAc-ALL-n5_top40markers-selectSubtypes_logExprs_Apr2020.png", height=360, width=480)
plotExpression(sce.nac.all, exprs_values = "logcounts", features="SEC63",
              x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7,
              add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                           geom = "crossbar", width = 0.3,
                                           colour=tableau20[1:14]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

