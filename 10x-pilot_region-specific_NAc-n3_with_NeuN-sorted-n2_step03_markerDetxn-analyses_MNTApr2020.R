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


pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/zold_regionSpecific_NAc-ALL-n5_top20markers_logExprs_Apr2020.pdf", height=7.5, width=9.5)
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
    png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/zold_regionSpecific_NAc-ALL-n5_top40markers-",i,"_logExprs_Apr2020.png"), height=1900, width=1200)
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
pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/zold_regionSpecific_NAc-ALL-n5_top20markers_logExprs_pt-coding_Apr2020.pdf", height=7.5, width=9.5)
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
      png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/zold_regionSpecific_NAc-ALL-n5_top40markers-",i,"_logExprs_pt-coding_Apr2020.png"), height=1900, width=1200)
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




## Markers pulled from Gokce, et al (doi: 10.1016/j.celrep.2016.06.059) =========
markers.gokce <- list(
  "D1.MSN" = c("Tac1","Drd1","Asic4","Slc35d3","Pdyn","Sfxn1","Nrxn1"),
                      # ^ edited from 'Drd1a', Accn4
  "D2.MSN" = c("Penk","Adora2a","Drd2","Gpr6","Grik3","Gpr52","Gnas"),
                      # ^ edited from 'A2a'
  "D1.Pcdh8" = c("Pcdh8","Adarb2","Tacr1","Tac1" ,"Nrxn2","Sema3e","Sema4a","Sema5a","Sema5b",
                 "Sema6d","Pcdh7","Ptprg","Ptprm","Ptpro","Ptpru","TAC3","Elavl4",
                                                                  # ^ edited from 'Tac2'
                 "Khdrbs3","Rbm20","Aff2","Lrpprc","Celf4",
                 # Depleted set:
                 "Nlgn1", "Calb1"),
  "D1.Foxp1" = c("Foxp1","Camk4"),
  "D2.Htr7" = c("Htr7","AGTR1","Penk","Tac1","Ptprt","Ngfr","Grik3","Cacng5",
                        # ^ edited from 'Agtr1a'
                "Tmeff2","Sox9","Sp8","Runx1","Mafb","Litaf",
                # Depleted set:
                "Cacna2d3","Synpr"),
  "D2.Synpr" = c("Synpr"),
  "gradient" = c("Dner","Cxcl14","Tnnt1","Meis2","Cartpt","Kcnip1","Calb1",
                 "Crym","Cnr1","Nnat","Gfra1","Wfs1","Th")
)

markers.gokce <- lapply(markers.gokce, toupper)

# Load SCE
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda",
     verbose=T)
    # sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all, ref.sampleInfo


# First drop "ambig.lowNtrxts" (93 nuclei)
sce.nac.all <- sce.nac.all[ ,sce.nac.all$cellType.final != "ambig.lowNtrxts"]
sce.nac.all$cellType.final <- droplevels(sce.nac.all$cellType.final)

# Which are there?
sapply(markers.gokce, function(x){x %in% rownames(sce.nac.all)})  # most of them

lapply(sapply(markers.gokce, function(x){x %in% rownames(sce.nac.all)}),  # most of them
       function(n){which(n==FALSE)})
    # So 'Drd1a', 'Accn4', 'A2a',       'Tac2',       'Agtr1a'


    # Exploring/identifying homologous gene names ====
    load("/dcl01/ajaffe/data/lab/singleCell/day_rat_snRNAseq/SCE_rat-NAc-PBd_w_matchingHsap-NAc-PBd_HomoloGene.IDs_MNT.rda", verbose=T)
    # sce.rat.PBsub, sce.hsap.PBsub, Readme
    
        # 'sce.hsap.PBsub' has 'HomoloGene.ID'

    hom = read.delim("http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt",
                     as.is=TRUE)
    
    hom_mm <- hom[hom$Common.Organism.Name == "mouse, laboratory", ]
    
    c('Drd1a', 'Accn4', 'A2a','Tac2','Agtr1a') %in% hom_mm$Symbol
    # FALSE     FALSE   FALSE   TRUE    TRUE
    # Drd1      Asic4  Adora2a
    
    # Then for 'Tac2'
    hom_mm$HomoloGene.ID[which(hom_mm$Symbol=="Tac2")]  # 7560
    rowData(sce.hsap.PBsub)$Symbol[rowData(sce.hsap.PBsub)$HomoloGene.ID==7560] # none...
    
    hom_hs <- hom[hom$Common.Organism.Name == "human", ]
    hom_hs[hom_hs$HomoloGene.ID==7560, ]  # symbol is TAC3
    'TAC3' %in% rowData(sce.nac.all)$Symbol # TRUE
        # ahhh so this one just didn't have a shared homolog with rat, I guess
    
    
    hom_mm$HomoloGene.ID[which(hom_mm$Symbol=="Agtr1a")]
    rowData(sce.hsap.PBsub)$Symbol[rowData(sce.hsap.PBsub)$HomoloGene.ID==3556]
        # AGTR1
    # end find synonyms =======




## Let's make a new dir and files for these graphics, since these are of various length
#dir.create("pdfs/exploration/gokce-etal_markers/")

for(i in names(markers.gokce)){
  pdf(paste0("./pdfs/exploration/gokce-etal_markers/",i,"-mouseStriatum-markers_human-NAcExpression_Apr2020.pdf"), height=2.6, width=3)
  # Print each gene's expression in its own page of the pdf
  for(g in 1:length(markers.gokce[[i]])){
    print(
      plotExpression(sce.nac.all, exprs_values = "logcounts", features=c(markers.gokce[[i]][g]),
                     x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7,
                     add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                  geom = "crossbar", width = 0.3,
                                                  colour=tableau20[1:14]) +
        theme(axis.text.x = element_text(angle=90, hjust=1, size=5.5), axis.text.y = element_text(size=7.5),
              plot.title = element_text(size=7)) +  
        ggtitle(label=paste0(i, " markers in human NAc subclusters: ", markers.gokce[[i]][g]))
    )
  }
  dev.off()
}


    
    
    ### ========================== ###
    ### SINGLE-NUCLEUS-LEVEL TESTS ###
    ### ========================== ###


### Single-nucleus level marker detection? =========
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda",
     verbose=T)
    # sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all, ref.sampleInfo

# First drop "ambig.lowNtrxts" (93 nuclei)
sce.nac.all <- sce.nac.all[ ,sce.nac.all$cellType.final != "ambig.lowNtrxts"]
sce.nac.all$cellType.final <- droplevels(sce.nac.all$cellType.final)
    
# Drop genes with all 0's
sce.nac.all <- sce.nac.all[!rowSums(assay(sce.nac.all, "counts"))==0, ]
    ## keeps 29236 genes


## Traditional t-test with design as in PB'd/limma approach ===
#mod <- with(colData(sce.nac.all), model.matrix(~ processDate + donor))
    # Error in .ranksafe_qr(full.design) : design matrix is not of full rank
    # -> just try processDate bc it describes more var (at least at PB-pan-brain level)
mod <- with(colData(sce.nac.all), model.matrix(~ processDate))
mod <- mod[ ,-1]

markers.nac.t.design <- findMarkers(sce.nac.all, groups=sce.nac.all$cellType.final,
                                    assay.type="logcounts", design=mod, test="t",
                                    direction="up", pval.type="all", full.stats=T)

sapply(markers.nac.t.design, function(x){table(x$FDR<0.05)})
    #      Astro Inhib.1 Inhib.2 Inhib.3 Inhib.4 Micro MSN.D1.1 MSN.D1.2 MSN.D1.3
    # FALSE 28405   28896   29057   28980   28953 27819    29154    28941    29149
    # TRUE    831     340     179     256     283  1417       82      295       87
    #       MSN.D1.4 MSN.D2.1 MSN.D2.2 Oligo   OPC
    # FALSE    29175    29006    29170 28743 28864
    # TRUE        61      230       66   493   372


## WMW: Blocking on sample (this test doesn't take 'design=' argument) ===
markers.nac.wilcox.block <- findMarkers(sce.nac.all, groups=sce.nac.all$cellType.final,
                                        assay.type="logcounts", block=sce.nac.all$donor, test="wilcox",
                                        direction="up", pval.type="all", full.stats=T)


sapply(markers.nac.wilcox.block, function(x){table(x$FDR<0.05)})
    ## none for about ~1/3 of these... none for Micros or Oligos...


## Binomial ===
markers.nac.binom.block <- findMarkers(sce.nac.all, groups=sce.nac.all$cellType.final,
                                       assay.type="logcounts", block=sce.nac.all$donor, test="binom",
                                       direction="up", pval.type="all", full.stats=T)


sapply(markers.nac.binom.block, function(x){table(x$FDR<0.05)})
    ## even worse than WMW... basically none


## Save all these for future reference
save(markers.nac.t.design, #markers.nac.wilcox.block, markers.nac.binom.block,
     file="rdas/markers-stats_NAc-n5_findMarkers-SN-LEVEL_MNTApr2020.rda")


## Print these to PNGs
markerList.t <- lapply(markers.nac.t.design, function(x){
    rownames(x)[x$FDR < 0.05]
  }
)

# Take top 40
genes.top40.t <- lapply(markerList.t, function(x){head(x, n=40)})

#dir.create("pdfs/exploration/NAc-n5-markers/")
for(i in names(genes.top40.t)){
  png(paste0("pdfs/exploration/NAc-n5-markers/NAc-all-n5_t-sn-level-top40markers-",i,"_logExprs_Apr2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.nac.all, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:14], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers: single-nucleus-level p.w. t-tests"))
  )
  dev.off()
}



### Cluster-vs-all single-nucleus-level iteration ======

## Load SCE with new info
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda",
     verbose=T)
    # sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all, ref.sampleInfo

# First drop "ambig.lowNtrxts" (93 nuclei)
sce.nac.all <- sce.nac.all[ ,sce.nac.all$cellType.final != "ambig.lowNtrxts"]
sce.nac.all$cellType.final <- droplevels(sce.nac.all$cellType.final)

# Drop genes with all 0's
sce.nac.all <- sce.nac.all[!rowSums(assay(sce.nac.all, "counts"))==0, ]  ## keeps 29236 genes


## Traditional t-test with design as in PB'd/limma approach ===
# mod <- with(colData(sce.nac), model.matrix(~ donor))
# mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`


#mod <- with(colData(sce.nac.all), model.matrix(~ processDate + donor))
    # Error in .ranksafe_qr(full.design) : design matrix is not of full rank
    # -> just try processDate bc it describes more var (at least at PB-pan-brain level)
mod <- with(colData(sce.nac.all), model.matrix(~ processDate))
mod <- mod[ ,-1]


markers.nac.t.1vAll <- list()
for(i in levels(sce.nac.all$cellType.final)){
  # Make temporary contrast
  sce.nac.all$contrast <- ifelse(sce.nac.all$cellType.final==i, 1, 0)
  # Test cluster vs. all
  markers.nac.t.1vAll[[i]] <- findMarkers(sce.nac.all, groups=sce.nac.all$contrast,
                                          assay.type="logcounts", design=mod, test="t",
                                          direction="up", pval.type="all", full.stats=T)
}

    ## Then, temp set of stats to get the standardized logFC
    temp.1vAll <- list()
    for(i in levels(sce.nac.all$cellType.final)){
      # Make temporary contrast
      sce.nac.all$contrast <- ifelse(sce.nac.all$cellType.final==i, 1, 0)
      # Test cluster vs. all
      temp.1vAll[[i]] <- findMarkers(sce.nac.all, groups=sce.nac.all$contrast,
                                     assay.type="logcounts", design=mod, test="t",
                                     std.lfc=TRUE,
                                     direction="up", pval.type="all", full.stats=T)
    }


## For some reason all the results are in the second List entry (first is always empty)

# Replace that empty slot with the entry with the actul stats
markers.nac.t.1vAll <- lapply(markers.nac.t.1vAll, function(x){ x[[2]] })
# Same for that with std.lfc
temp.1vAll <- lapply(temp.1vAll, function(x){ x[[2]] })

# Now just pull from the 'stats.0' DataFrame column
markers.nac.t.1vAll <- lapply(markers.nac.t.1vAll, function(x){ x$stats.0 })
temp.1vAll <- lapply(temp.1vAll, function(x){ x$stats.0 })

# Re-name std.lfc column and add to the first result
for(i in names(temp.1vAll)){
  colnames(temp.1vAll[[i]])[1] <- "std.logFC"
  markers.nac.t.1vAll[[i]] <- cbind(markers.nac.t.1vAll[[i]], temp.1vAll[[i]]$std.logFC)
  # Oh the colname is kept weird
  colnames(markers.nac.t.1vAll[[i]])[4] <- "std.logFC"
  # Then re-organize
  markers.nac.t.1vAll[[i]] <- markers.nac.t.1vAll[[i]][ ,c("logFC","std.logFC","log.p.value","log.FDR")]
}


## Let's save this along with the previous pairwise results
save(markers.nac.t.design, markers.nac.t.1vAll,
     file="rdas/markers-stats_NAc-n5_findMarkers-SN-LEVEL_MNTApr2020.rda")



## Print these to pngs
markerList.t.1vAll <- lapply(markers.nac.t.1vAll, function(x){
    rownames(x)[x[ ,"log.FDR"] < log10(0.001)]
  }
)
genes.top40.t.1vAll <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t.1vAll)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/NAc-n5-markers/NAc-all-n5_t-sn-level_1vALL_top40markers-",i,"_logExprs_May2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.nac.all, exprs_values = "logcounts", features=genes.top40.t.1vAll[[i]],
                   x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:14], length(genes.top40.t.1vAll[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers: single-nucleus-level t-tests, cluster-vs-all-others"))
  )
  dev.off()
}



## How do these top 40 intersect? ===
sapply(names(genes.top40.t), function(c){
  length(intersect(genes.top40.t[[c]],
                   genes.top40.t.1vAll[[c]]))
})
    #     Astro  Inhib.1  Inhib.2  Inhib.3  Inhib.4    Micro MSN.D1.1 MSN.D1.2
    #        35       30       16       26       22       40       17       24
    #  MSN.D1.3 MSN.D1.4 MSN.D2.1 MSN.D2.2    Oligo      OPC
    #         7        9       23       16       32       29


## Write these top 40 lists to a csv
names(markerList.t) <- paste0(names(markerList.t),"_pw")
names(markerList.t.1vAll) <- paste0(names(markerList.t.1vAll),"_1vAll")

top40genes <- cbind(sapply(markerList.t, function(x) head(x, n=40)),
                    sapply(markerList.t.1vAll, function(y) head(y, n=40)))
top40genes <- top40genes[ ,sort(colnames(top40genes))]

write.csv(top40genes, file="tables/top40genesLists_NAc-n5_cellType.final_SN-LEVEL-tests_May2020.csv",
          row.names=FALSE)



## Make marker array for Supp figure (MNT suggested panel A) ===========
load("rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda", verbose=T)
    #sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all, ref.sampleInfo

load("rdas/markers-stats_NAc-n5_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.nac.t.design, markers.nac.t.1vAll

# First make interneuron subset
sce.nac.int <- sce.nac.all[ ,grep("Inhib.", sce.nac.all$cellType.final)]
sce.nac.int$cellType.final <- droplevels(sce.nac.int$cellType.final)


# Take top four for 4 inhib. interneuron pops
topToPrint <- as.data.frame(sapply(markers.nac.t.1vAll, function(x) {head(rownames(x),n=4)}))
topToPrint <- topToPrint[grep("Inhib.", names(topToPrint))]

# Manual assignment, bc 'Inhib.2' is mostly driven by noise (but 0 median)
topToPrint["Inhib.2"] <- c("KCNJ6", "SDK1", "LRFN2","NR2F1-AS1")

table(unlist(topToPrint) %in% rownames(sce.nac.int)) # good

# Print
pdf("pdfs/pubFigures/suppFig_NAc_interneuron-marker-array_MNTSep2020.pdf", height=8, width=5.5)
print(
  plotExpression(sce.nac.int, exprs_values = "logcounts", features=c(t(topToPrint)),
                 x="cellType.final", colour_by="cellType.final", point_alpha=0.6, point_size=1.5, ncol=4,
                 add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                              geom = "crossbar", width = 0.3,
                                              colour=rep(tableau10medium[1:4], length(unlist(topToPrint)))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12), plot.title = element_text(size = 25)) +  
    ggtitle(label="Inhib.1          Inhib.2           Inhib.3          Inhib.4") + xlab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 13),
          axis.title.y = element_text(angle = 90, size = 16),
          plot.title = element_text(size = 15),
          panel.grid.major=element_line(colour="grey95", size=0.8),
          panel.grid.minor=element_line(colour="grey95", size=0.4))
)
dev.off()


## For MNT version panel B
pdf("pdfs/pubFigures/suppFig_NAc_interneuron-experiment-panelB_MNTSep2020.pdf", height=3.2, width=4.5)
print(
  plotExpression(sce.nac.int, exprs_values = "logcounts", features=c("GAD1", "KIT", "PTHLH", "PVALB"),
                 x="cellType.final", colour_by="cellType.final", point_alpha=0.7, point_size=1.2, ncol=2,
                 add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                              geom = "crossbar", width = 0.3,
                                              colour=rep(tableau10medium[1:4], 4)) +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9.5),
          axis.title.y = element_text(angle = 90, size = 10),
          panel.grid.major=element_line(colour="grey95", size=0.8),
          panel.grid.minor=element_line(colour="grey95", size=0.4))
)
dev.off()


## Supp Fig 5: Other markers pointed out in text ===
genes2print <- c("DRD1", "DRD2", "CASZ1", "GPR6", "EBF1", "GRM8")

# First drop "ambig.lowNtrxts" (93 nuclei)
sce.nac.all <- sce.nac.all[ ,sce.nac.all$cellType.final != "ambig.lowNtrxts"]
sce.nac.all$cellType.final <- droplevels(sce.nac.all$cellType.final)

pdf("pdfs/pubFigures/suppFig_NAc_other-MSN-markers_MNTSep2020.pdf", height=4.5, width=6.5)
print(
  plotExpression(sce.nac.all, exprs_values = "logcounts", features=genes2print,
                 x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=1.0, ncol=2,
                 add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                              geom = "crossbar", width = 0.3,
                                              colour=rep(tableau20[1:14], length(genes2print))) +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9.5),
          axis.title.y = element_text(angle = 90, size = 10),
          panel.grid.major=element_line(colour="grey95", size=0.8),
          panel.grid.minor=element_line(colour="grey95", size=0.4))
)
dev.off()



## Top markers for D1.4 / D2.2 often co-expressed ===
load("rdas/markers-stats_NAc-n5_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.nac.t.design, markers.nac.t.1vAll

cell.idx <- splitit(sce.nac.all$cellType.final)
dat <- as.matrix(assay(sce.nac.all, "logcounts"))
genes <- head(rownames(markers.nac.t.1vAll[["MSN.D1.4"]]), n=40)

pdf('pdfs/pubFigures/suppFigure_heatmap-Exprs_NAc-n5_top40-D1.4markers-1vAlltest_MNTSep2020.pdf',
    useDingbats=TRUE, height=6, width=10)
current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMeans(dat[genes, ii])))
pheatmap(t(current_dat), cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(0.02, 5, length.out = 101),
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "OrRd"))(100),
         fontsize_row = 20, fontsize_col=15)
dev.off()






### MNT add 18Nov2020 =================================
# -> What if add param/requirement that for any given subcluster, median expression has to > 0?
load("rdas/markers-stats_HPC-n3_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.nac.t.design, markers.nac.t.1vAll

## Load SCE 
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all, ref.sampleInfo

table(sce.nac.all$cellType.final)

# First drop "ambig.lowNtrxts" (93 nuclei)
sce.nac.all <- sce.nac.all[ ,sce.nac.all$cellType.final != "ambig.lowNtrxts"]
sce.nac.all$cellType.final <- droplevels(sce.nac.all$cellType.final)

# Remove 0 genes across all nuclei
sce.nac.all <- sce.nac.all[!rowSums(assay(sce.nac.all, "counts"))==0, ]


## Make list of Boolean param / cell subtype ===
cellSubtype.idx <- splitit(sce.nac.all$cellType.final)
medianNon0.idx <- lapply(cellSubtype.idx, function(x){
  apply(as.matrix(assay(sce.nac.all, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

lengths(medianNon0.idx)
sapply(medianNon0.idx, head)

# Add these to the stats for each set of markers
for(i in names(markers.nac.t.1vAll)){
  markers.nac.t.1vAll[[i]] <- cbind(markers.nac.t.1vAll[[i]],
                                    medianNon0.idx[[i]][match(rownames(markers.nac.t.1vAll[[i]]),
                                                              names(medianNon0.idx[[i]]))])
  colnames(markers.nac.t.1vAll[[i]])[5] <- "non0median"
}


## Use these restrictions to print (to png) a refined top 40, as before ===
markerList.t.1vAll <- lapply(markers.nac.t.1vAll, function(x){
  rownames(x)[x$log.FDR < log10(0.000001) & x$non0median==TRUE]
  }
)
    # lengths(markerList.t.1vAll)     # ( **without $non0median==TRUE restriction )
        #    Astro  Inhib.1  Inhib.2  Inhib.3  Inhib.4    Micro MSN.D1.1 MSN.D1.2
        #     4956     1295     1589     4917     4719     3656     1474     2468
        # MSN.D1.3 MSN.D1.4 MSN.D2.1 MSN.D2.2    Oligo      OPC
        #     2434     5893     2027     3106     2248     3259

lengths(markerList.t.1vAll)
    #    Astro  Inhib.1  Inhib.2  Inhib.3  Inhib.4    Micro MSN.D1.1 MSN.D1.2
    #     1214      690      561     1953     1716      763      700      989
    # MSN.D1.3 MSN.D1.4 MSN.D2.1 MSN.D2.2    Oligo      OPC
    #     1173     3054     1069     1912      700     1086

genes.top40.t <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/NAc-n5-markers/NAc_t-sn-level_1vALL_top40markers-REFINED-",gsub(":",".",i),"_logExprs_Nov2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.nac.all, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:14], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers, refined: single-nucleus-level p.w. t-tests, cluster-vs-all"))
  )
  dev.off()
}



## Do the same with the pairwise result ('markers.nac.t.design') === === ===
# Add these to the stats for each set of markers
for(i in names(markers.nac.t.design)){
  markers.nac.t.design[[i]] <- cbind(markers.nac.t.design[[i]],
                                     medianNon0.idx[[i]][match(rownames(markers.nac.t.design[[i]]),
                                                               names(medianNon0.idx[[i]]))])
  colnames(markers.nac.t.design[[i]])[16] <- "non0median"
}

markerList.t <- lapply(markers.nac.t.design, function(x){
  rownames(x)[x$FDR < 0.05 & x$non0median==TRUE]
  }
)
    # lengths(markerList.t)     # ( **without $non0median==TRUE restriction )
        #   Astro  Inhib.1  Inhib.2  Inhib.3  Inhib.4    Micro MSN.D1.1 MSN.D1.2
        #     831      340      179      256      283     1417       82      295
        #MSN.D1.3 MSN.D1.4 MSN.D2.1 MSN.D2.2    Oligo      OPC
        #      87       61      230       66      493      372

lengths(markerList.t)
    #   Astro  Inhib.1  Inhib.2  Inhib.3  Inhib.4    Micro MSN.D1.1 MSN.D1.2
    #     378      110        9      145      105      402       45      201
    #MSN.D1.3 MSN.D1.4 MSN.D2.1 MSN.D2.2    Oligo      OPC
    #      64       54      143       57      333      187


genes.top40.t <- lapply(markerList.t, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/NAc-n5-markers/NAc_t-sn-level_pairwise_top40markers-REFINED-", i, "_logExprs_Nov2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.nac.all, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:14], length(genes.top40.t[[i]]))) +
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

write.csv(top40genes, file="tables/top40genesLists-REFINED_NAc-n5_cellType.final_Nov2020.csv",
          row.names=FALSE)









