### LAH 10x snRNA-seq workflow: step 03 - marker detection
###   **Region-specific analyses**
###     - (3x) DLPFC samples
### 25May2021
#####################################################################

library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(batchelor)
library(DropletUtils)
library(jaffelab)
library(limma)
library(DeconvoBuddies)
library(here)

source("plotExpressionCustom.R")

### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

# ===


## Load SCE with new info
load(here("rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda"), verbose=T)
# sce.dlpfc
# chosen.hvgs.dlpfc
# pc.choice.dlpfc
# clusterRefTab.dlpfc
# ref.sampleInfo
# annotationTab.dlpfc

table(sce.dlpfc$cellType)

# Remove 0 genes across all nuclei
sce.dlpfc <- sce.dlpfc[!rowSums(assay(sce.dlpfc, "counts"))==0, ]  # 29310

## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.dlpfc), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`


# Run pairwise t-tests
markers.t.design <- findMarkers(sce.dlpfc, groups=sce.dlpfc$cellType,
                                      assay.type="logcounts", design=mod, test="t",
                                      direction="up", pval.type="all", full.stats=T)

sapply(markers.t.design, function(x){table(x$FDR<0.05)})
#       ambig.glial_A ambig.glial_B Astro Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Inhib_A Inhib_B Inhib_C Inhib_D
# FALSE         28848         28793 28988   29277   29232   29271   29203   29281   29255   29253   29303   29299   29240
# TRUE            462           517   322      33      78      39     107      29      55      57       7      11      70
#       Inhib_E Inhib_F Micro Oligo   OPC
# FALSE   29103   29155 28854 29080 29136
# TRUE      207     155   456   230   174


## WMW: Blocking on donor (this test doesn't take 'design=' argument) ===
markers.wilcox.block <- findMarkers(sce.dlpfc, groups=sce.dlpfc$cellType,
                                          assay.type="logcounts", block=sce.dlpfc$donor, test="wilcox",
                                          direction="up", pval.type="all", full.stats=T)

# no warnings as in pan-brain analyses, but NO results of FDR<0.05...:
sapply(markers.wilcox.block, function(x){table(x$FDR<0.05)})
      # Actually some decent results but many subclusters with 0 hits
# $Micro
# FALSE  TRUE 
# 29268    42 
# $Oligo
# FALSE  TRUE 
# 29224    86 
# $OPC
# FALSE  TRUE 
# 29293    17

## Binomial ===
markers.binom.block <- findMarkers(sce.dlpfc, groups=sce.dlpfc$cellType,
                                   assay.type="logcounts", block=sce.dlpfc$donor, test="binom",
                                   direction="up", pval.type="all", full.stats=T)

sapply(markers.binom.block, function(x){table(x$FDR<0.05)})
    # All FALSE

## Save all these for future reference
save(markers.t.design, markers.wilcox.block, #markers.binom.block,
     file="rdas/revision/markers-stats_DLPFC-n3_findMarkers-SN-LEVEL_LAHMay2021.rda")


# Print these to pngs
markerList.t <- lapply(markers.t.design, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)

genes.top40.t <- lapply(markerList.t, function(x){head(x, n=40)})

# $ambig.glial_A - T cells
# "LILRB5" Leukocyte    
# "RUNX3" 
# "SIGLEC1" - expressed only by a subpopulation of macrophages   
# "CD163" - exclusively expressed in monocytes and macrophages
# "IQGAP2"  
# "CD2"  peripheral blood T-cells      "ITK"  play a role in T-cell proliferation and differentiation   
#  "NKG7"      
# [9] "PRF1"       "GPR174"     "F13A1"      "CLEC2B"     "DOK2"       "CST7"       "SLFN12L"    "LINC00243" 
# [17] "LINC00861"  "IL32"       "LINC01839"  "GZMB"       "TRGV5"      "FGFBP2"     "MS4A6A"     "CD52"      
# [25] "MYO1G"      "SKAP1"      "TBC1D10C"   "VNN1"       "CCL3"       "CCR7"       "IL7R"       "JAML"      
# [33] "ZAP70"      "SH2D2A"     "LILRB2"     "PTPN22"     "ANXA1"      "SLAMF6"     "AL109767.1" "MARCO"     
# 
# $ambig.glial_B - Mural
# [1] "CARMN"      "EBF1"       "ITIH5"      "TBX3"       *"COL1A2"     "FOXC2"      "FOXC1"      "AC092957.1"
# [9] "NOTCH3"     "TBX2"       "PGR"        "LRRC32"     "TBX18"      "STARD8"     "LINC02147"  "TINAGL1"   
# [17] "TFPI"       *"CFH"        "SLC6A13"    "BMP5"       "ARHGAP29"   "FHL5"       "LINC02172"  "NOSTRIN"   
# [25] "PEAR1"      "LEF1"       "CLDN5"      "A4GALT"     "FOXS1"      "GJC1"       "SLC12A7"    "COL6A2"    
# [33] "FOXF2"      "ADGRF5"     "LINC02188"  "FOXL1"      "MYLK2"      "SOX18"      "AL353780.1" "GGT5"   

# genes.top40.t$ambig.glial_B[genes.top40.t$ambig.glial_B %in% mural]
# [1] "COL1A2"     "CFH"        "SLC6A13"    "TBX18"      "ITIH5"      "AC092957.1" "ARHGAP29"   "EBF1"

#dir.create("pdfs/revision/DLPFC/")
for(i in names(genes.top40.t)){
  png(paste0("pdfs/revision/DLPFC/DLPFC_t-sn-level_pairwise_top40markers-", i, "_logExprs_LAH2021.png"), height=1900, width=1200)
  print(
    plotExpressionCustom(sce = sce.dlpfc,
                         features = genes.top40.t[[i]], 
                         features_name = i, 
                         anno_name = "cellType")
  )
  dev.off()
}




### Cluster-vs-all single-nucleus-level iteration ================================
# MNT 30Apr2020

## Load SCE with new info
load(here("rdas/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda"),
     verbose=T)
    # sce.dlpfc, chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo

table(sce.dlpfc$cellType)

# First drop "Ambig.lowNtrxts" (101 nuclei)
sce.dlpfc <- sce.dlpfc[ ,sce.dlpfc$cellType != "Ambig.lowNtrxts"]
sce.dlpfc$cellType <- droplevels(sce.dlpfc$cellType)

# Remove 0 genes across all nuclei
sce.dlpfc <- sce.dlpfc[!rowSums(assay(sce.dlpfc, "counts"))==0, ]  # keeps same 28757 genes


## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.dlpfc), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`


markers.t.1vAll <- list()
for(i in levels(sce.dlpfc$cellType)){
  # Make temporary contrast
  sce.dlpfc$contrast <- ifelse(sce.dlpfc$cellType==i, 1, 0)
  # Test cluster vs. all
  markers.t.1vAll[[i]] <- findMarkers(sce.dlpfc, groups=sce.dlpfc$contrast,
                                            assay.type="logcounts", design=mod, test="t",
                                            direction="up", pval.type="all", full.stats=T)
}

    ## Then, temp set of stats to get the standardized logFC
    temp.1vAll <- list()
    for(i in levels(sce.dlpfc$cellType)){
      # Make temporary contrast
      sce.dlpfc$contrast <- ifelse(sce.dlpfc$cellType==i, 1, 0)
      # Test cluster vs. all
      temp.1vAll[[i]] <- findMarkers(sce.dlpfc, groups=sce.dlpfc$contrast,
                                     assay.type="logcounts", design=mod, test="t",
                                     std.lfc=TRUE,
                                     direction="up", pval.type="all", full.stats=T)
    }



    ## As with DLPFC, for some reason all the results are in the
     #    second List entry (first is always empty)

head(markers.t.1vAll[["Oligo"]][[2]])
    ## Nice, MBP and PLP1 are again in the top 6


sapply(markers.t.1vAll, function(x){
  table(x[[2]]$stats.0$log.FDR < log10(.001))
})
    #       Oligo Micro   OPC Inhib.5 Inhib.2 Astro Inhib.3 Excit.2 Inhib.4 Tcell
    # FALSE 24914 22821 23236   24401   23540 21612   21436   22608   24460 26858
    # TRUE   3843  5936  5521    4356    5217  7145    7321    6149    4297  1899
    #       Inhib.1 Excit.5 Excit.3 Excit.1 Excit.4
    # FALSE   25456   25913   19431   23726   24170
    # TRUE     3301    2844    9326    5031    4587



# Replace that empty slot with the entry with the actul stats
markers.t.1vAll <- lapply(markers.t.1vAll, function(x){ x[[2]] })
# Same for that with std.lfc
temp.1vAll <- lapply(temp.1vAll, function(x){ x[[2]] })

# Now just pull from the 'stats.0' DataFrame column
markers.t.1vAll <- lapply(markers.t.1vAll, function(x){ x$stats.0 })
temp.1vAll <- lapply(temp.1vAll, function(x){ x$stats.0 })

# Re-name std.lfc column and add to the first result
for(i in names(temp.1vAll)){
  colnames(temp.1vAll[[i]])[1] <- "std.logFC"
  markers.t.1vAll[[i]] <- cbind(markers.t.1vAll[[i]], temp.1vAll[[i]]$std.logFC)
  # Oh the colname is kept weird
  colnames(markers.t.1vAll[[i]])[4] <- "std.logFC"
  # Then re-organize
  markers.t.1vAll[[i]] <- markers.t.1vAll[[i]][ ,c("logFC","std.logFC","log.p.value","log.FDR")]
}





## Let's save this along with the previous pairwise results
save(markers.t.1vAll, markers.t.design, markers.wilcox.block,
#     file="rdas/revision/markers-stats_DLPFC-n3_findMarkers-SN-LEVEL_LAHMay2021.rda")
#     (deleting this older version - doesn't have the std.lfc result)
     file="rdas/revision/markers-stats_DLPFC-n3_findMarkers-SN-LEVEL_LAHMay2021.rda")


## Print these to pngs
markerList.t.1vAll <- lapply(markers.t.1vAll, function(x){
  rownames(x)[x$log.FDR < log10(0.000001)]
 }
)
genes.top40.t <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0(here(pdfs/DLPFC/DLPFC_t-sn-level_1vALL_top40markers-",gsub(":",".",i),"_logExprs_Apr2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.dlpfc, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:15], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers: single-nucleus-level p.w. t-tests, cluster-vs-all"))
  )
  dev.off()
}


## How do they intersect?
markerList.t.pw <- lapply(markers.t.design, function(x){
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

write.csv(top40genes, file="tables/top40genesLists_DLPFC-n3_cellType_SN-LEVEL-tests_May2020.csv",
          row.names=FALSE)




### MNT add 18Nov2020 =================================
  # -> What if add param/requirement that for any given subcluster, median expression has to > 0?
load("rdas/revision/markers-stats_DLPFC-n3_findMarkers-SN-LEVEL_LAHMay2021.rda", verbose=T)
    # markers.t.1vAll, markers.t.design, markers.wilcox.block

## Load SCE 
load(here("rdas/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda",
     verbose=T)
    # sce.dlpfc, chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo

table(sce.dlpfc$cellType)

# First drop "Ambig.lowNtrxts" (101 nuclei)
sce.dlpfc <- sce.dlpfc[ ,sce.dlpfc$cellType != "Ambig.lowNtrxts"]
sce.dlpfc$cellType <- droplevels(sce.dlpfc$cellType)

# Remove 0 genes across all nuclei
sce.dlpfc <- sce.dlpfc[!rowSums(assay(sce.dlpfc, "counts"))==0, ]


## Make list of Boolean param / cell subtype ===
cellSubtype.idx <- splitit(sce.dlpfc$cellType)
medianNon0.idx <- lapply(cellSubtype.idx, function(x){
  apply(as.matrix(assay(sce.dlpfc, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

lengths(medianNon0.idx)
sapply(medianNon0.idx, head)

# Add these to the stats for each set of markers
for(i in names(markers.t.1vAll)){
  markers.t.1vAll[[i]] <- cbind(markers.t.1vAll[[i]],
                                    medianNon0.idx[[i]][match(rownames(markers.t.1vAll[[i]]),
                                                           names(medianNon0.idx[[i]]))])
  colnames(markers.t.1vAll[[i]])[5] <- "non0median"
}


## Use these restrictions to print (to png) a refined top 40, as before ===
markerList.t.1vAll <- lapply(markers.t.1vAll, function(x){
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
  png(paste0(here(pdfs/DLPFC/DLPFC_t-sn-level_1vALL_top40markers-REFINED-",gsub(":",".",i),"_logExprs_Nov2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.dlpfc, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:15], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers, refined: single-nucleus-level p.w. t-tests, cluster-vs-all"))
  )
  dev.off()
}



## Do the same with the pairwise result ('markers.t.design') === === ===
# Add these to the stats for each set of markers
for(i in names(markers.t.design)){
  markers.t.design[[i]] <- cbind(markers.t.design[[i]],
                                    medianNon0.idx[[i]][match(rownames(markers.t.design[[i]]),
                                                           names(medianNon0.idx[[i]]))])
  colnames(markers.t.design[[i]])[17] <- "non0median"
}

markerList.t <- lapply(markers.t.design, function(x){
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
  png(paste0(here(pdfs/DLPFC/DLPFC_t-sn-level_pairwise_top40markers-REFINED-", i, "_logExprs_Nov2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.dlpfc, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
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

write.csv(top40genes, file="tables/top40genesLists-REFINED_DLPFC-n3_cellType_Nov2020.csv",
          row.names=FALSE)




## Aside: add in 't.stat' as in 'step04' analyses to save for LoHu/LeCo ===
for(s in names(markers.t.1vAll)){
  markers.t.1vAll[[s]]$t.stat <- markers.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.dlpfc))
}

save(markers.t.1vAll, markers.t.design, sce.dlpfc,
     file="rdas/revision/markerstats-and-SCE_DLPFC-n3_sn-level_cleaned_MNTNov2020.rda")




