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


## Load SCE with new info
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_HPC-n3_cleaned-combined_SCE_MNT2021.rda")
    # sce.hpc, clusterRefTab.hpc, chosen.hvgs.hpc, ref.sampleInfo

table(sce.hpc$cellType)
    #        ambig ambig.glial_A ambig.glial_B       Astro_A       Astro_B 
    #           43            15            19           936           234 
    #      Astro_C  drop.doublet       Excit_A       Excit_B       Excit_C 
    #          105             5           486            87             6 
    #      Excit_D       Excit_E       Excit_F         Inhib         Micro 
    #           35             6            38           331          1161 
    #        Oligo           OPC         Tcell 
    #         5912           823            26

# First drop "drop.doublet" (5 nuclei)
sce.hpc <- sce.hpc[ ,sce.hpc$cellType != "drop.doublet"]
sce.hpc$cellType <- droplevels(sce.hpc$cellType)

# Remove 0 genes across all nuclei
sce.hpc <- sce.hpc[!rowSums(assay(sce.hpc, "counts"))==0, ]  # keeps same 28768 genes


## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.hpc), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`



### Make list of Boolean param / cell subtype ===
# Will use this to assess more 'valid', non-noise-driving markers
cellSubtype.idx <- splitit(sce.hpc$cellType)
medianNon0.hpc <- lapply(cellSubtype.idx, function(x){
  apply(as.matrix(assay(sce.hpc, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

lengths(medianNon0.hpc)
sapply(medianNon0.hpc, table)




# Run pairwise t-tests
markers.hpc.t.pw <- findMarkers(sce.hpc, groups=sce.hpc$cellType,
                                assay.type="logcounts", design=mod, test="t",
                                direction="up", pval.type="all", full.stats=T)

sapply(markers.hpc.t.pw, function(x){table(x$FDR<0.05)})
    #       Astro Excit.1 Excit.2 Excit.3 Excit.4 Excit.5 Inhib.1 Inhib.2 Inhib.3
    # FALSE 28244   28435   28540   28498   28140   28394   28552   28629   28593
    # TRUE    513     322     217     259     617     363     205     128     164

    #       Inhib.4 Inhib.5 Micro Oligo   OPC Tcell
    # FALSE   28601   28595 28093 28356 28506 28006
    # TRUE      156     162   664   401   251   751


## WMW: Blocking on donor (this test doesn't take 'design=' argument) ===
markers.hpc.wilcox.block <- findMarkers(sce.hpc, groups=sce.hpc$cellType,
                                          assay.type="logcounts", block=sce.hpc$donor, test="wilcox",
                                          direction="up", pval.type="all", full.stats=T)

# no warnings as in pan-brain analyses, but NO results of FDR<0.05...:
sapply(markers.hpc.wilcox.block, function(x){table(x$FDR<0.05)})
      # Actually some decent results but many subclusters with 0 hits


## Binomial ===
markers.hpc.binom.block <- findMarkers(sce.hpc, groups=sce.hpc$cellType,
                                         assay.type="logcounts", block=sce.hpc$donor, test="binom",
                                         direction="up", pval.type="all", full.stats=T)

sapply(markers.hpc.binom.block, function(x){table(x$FDR<0.05)})
    # only a couple dozen hits for glia, only - disregard these


# Add respective 'non0median' column to the stats for each set of markers
#   (just the pw t-test for now)
for(i in names(markers.hpc.t.pw)){
  markers.hpc.t.pw[[i]] <- cbind(markers.hpc.t.pw[[i]],
                                 medianNon0.hpc[[i]][match(rownames(markers.hpc.t.pw[[i]]),
                                                           names(medianNon0.hpc[[i]]))])
  colnames(markers.hpc.t.pw[[i]])[20] <- "non0median"
}

sapply(markers.hpc.t.pw, function(x){table(x$FDR<0.05 & x$non0median == TRUE)["TRUE"]})

markerList.t.pw.logcounts <- lapply(markers.hpc.t.pw, function(x){
  rownames(x)[x$FDR < 0.05 & x$non0median==TRUE]
})

lengths(markerList.t.pw.logcounts)
    #   ambig ambig.glial_A ambig.glial_B       Astro_A       Astro_B       Astro_C       Excit_A 
    #      69           118             0           139            93             1            43 
    # Excit_B       Excit_C       Excit_D       Excit_E       Excit_F         Inhib         Micro 
    #     102           128            69            98            64            38           223 
    #   Oligo           OPC         Tcell 
    #     116            65           129 

# Explore some of the prelim results
plotExpression(sce.hpc, exprs_values = "logcounts", features=head(markerList.t.pw.logcounts[["ambig.glial_A"]], 15),
               x="cellType", colour_by="cellType", point_alpha=0.2, point_size=.7, ncol=5,
               add_legend=F, show_median=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25))

    # total UMI distribution?
    cellIdx <- splitit(sce.hpc$cellType)
    sapply(cellIdx, function(x){quantile(sce.hpc$sum[x])})
    sapply(cellIdx, function(x){round(quantile(sce.hpc$sizeFactor[x]),3)})
    
        # So 'ambig.glial_B' should definitely be a 'drop.lowNTx';
        # 'Astro_C' should be, too, since, while median total UMIs ~1500, their SF median is 0.18,
        #   so these the logcounts are scaled UP by a median of ~5x, inducing the single 'marker'

    # What about 'ambig' & 'ambig.glial_A'?
    head(markerList.t.pw.logcounts[["ambig"]], 40)
        # [1] "COL1A2"     "CFH"        "PHLDB2"     "SLC6A13"    "TBX18"     
        # [6] "ITIH5"      "AC092957.1" "ARHGAP29"   "COBLL1"     "CEMIP"     
        # [11] "EBF1"       "SLC13A4"    "SLC16A12"   "ABCA9"      "RBPMS"     
        # [16] "ADAMTS12"   "SLC6A12"    "FBLN1"      "MYO1B"      "COLEC12"
        #[21] "SVIL"        "ADAMTS9-AS2" "NID1"        "ECM2"        "PDGFRB"     
        # [26] "ARHGAP10"    "ATP1A2"      "KANK2"       "BICC1"       "IGFBP7"     
        # [31] "PLCB4"       "UACA"        "PARVA"       "TPM1"        "PLA2R1"     
        # [36] "RBMS3"       "FN1"         "VIM"         "YAP1"        "EYA1" 
            # These look like mural cells (either pericytes/VSMCs; unclear which)
      
    head(markerList.t.pw.logcounts[["ambig.glial_A"]], 40)
        # [1] "AC008080.4" "GPR17"      "ADAM33"     "BCAN"       "TNS3"      
        # [6] "LIMS2"      "BCAS1"      "AC006058.1" "MDFI"       "PRICKLE1"  
        # [11] "SEMA5B"     "ADAMTS20"   "SLC16A10"   "CENPJ"      "KCNS3"     
        # [16] "BMPER"      "ASIC4"      "NBAT1"      "AC044781.1" "SOX4"      
        # [21] "HRASLS"     "SMOC1"      "EPHB1"      "B3GNT7"     "DOCK6"     
        # [26] "RHBDL3"     "TMEM163"    "SIRT2"      "C9orf129"   "SCRG1"     
        # [31] "AC110285.1" "ANGPTL2"    "ELFN2"      "CRB1"       "RASGEF1B"  
        # [36] "REXO5"      "TNFRSF21"   "FRMD4A"     "IL17RB"     "KY"
            # These seem to be the 'differentiation-committed OPCs (COP)s' from
            #    Marques-Zeisel, et al. 2017 (doi: 10.1126/science.aaf6463)
    
        # MNT 21May2021: go back and re-annotate these accordingly


## Save all these for future reference ===
save(markers.hpc.t.pw.logcounts, #markers.hpc.wilcox.block, #markers.hpc.binom.block,
     medianNon0.hpc,
     file="rdas/revision/markers-stats_HPC-n3_findMarkers-SN-LEVEL_MNT2021.rda")



# Print these to pngs
markerList.t <- lapply(markers.hpc.t.pw, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)

genes.top40.t <- lapply(markerList.t, function(x){head(x, n=40)})


#dir.create("pdfs/revision/HPC/")
for(i in names(genes.top40.t)){
  png(paste0("pdfs/revision/HPC/HPC_t_pairwise_top40markers-", i, "_logExprs_MNT2021.png"), height=1900, width=1200)
  print(
    plotExpression(sce.hpc, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:15], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers: single-nucleus-level p.w. t-tests"))
  )
  dev.off()
}


# Save
save(markers.hpc.t.pw, markers.hpc.t.pw.logcounts,
     markers.hpc.wilcox.block, medianNon0.hpc,
     file="rdas/revision/markers-stats_HPC-n3_findMarkers-SN-LEVEL_MNT2021.rda")




### Cluster-vs-all single-nucleus-level iteration ================================
# MNT 30Apr2020

## Load SCE with new info
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo

table(sce.hpc$cellType)

# First drop "Ambig.lowNtrxts" (101 nuclei)
sce.hpc <- sce.hpc[ ,sce.hpc$cellType != "Ambig.lowNtrxts"]
sce.hpc$cellType <- droplevels(sce.hpc$cellType)

# Remove 0 genes across all nuclei
sce.hpc <- sce.hpc[!rowSums(assay(sce.hpc, "counts"))==0, ]  # keeps same 28757 genes


## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.hpc), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`


markers.hpc.t.1vAll <- list()
for(i in levels(sce.hpc$cellType)){
  # Make temporary contrast
  sce.hpc$contrast <- ifelse(sce.hpc$cellType==i, 1, 0)
  # Test cluster vs. all
  markers.hpc.t.1vAll[[i]] <- findMarkers(sce.hpc, groups=sce.hpc$contrast,
                                            assay.type="logcounts", design=mod, test="t",
                                            direction="up", pval.type="all", full.stats=T)
}

    ## Then, temp set of stats to get the standardized logFC
    temp.1vAll <- list()
    for(i in levels(sce.hpc$cellType)){
      # Make temporary contrast
      sce.hpc$contrast <- ifelse(sce.hpc$cellType==i, 1, 0)
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
save(markers.hpc.t.1vAll, markers.hpc.t.pw, markers.hpc.wilcox.block,
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
markerList.t.pw <- lapply(markers.hpc.t.pw, function(x){
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

write.csv(top40genes, file="tables/top40genesLists_HPC-n3_cellType_SN-LEVEL-tests_May2020.csv",
          row.names=FALSE)




### MNT add 18Nov2020 =================================
  # -> What if add param/requirement that for any given subcluster, median expression has to > 0?
load("rdas/markers-stats_HPC-n3_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.hpc.t.1vAll, markers.hpc.t.pw, markers.hpc.wilcox.block

## Load SCE 
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo

table(sce.hpc$cellType)

# First drop "Ambig.lowNtrxts" (101 nuclei)
sce.hpc <- sce.hpc[ ,sce.hpc$cellType != "Ambig.lowNtrxts"]
sce.hpc$cellType <- droplevels(sce.hpc$cellType)

# Remove 0 genes across all nuclei
sce.hpc <- sce.hpc[!rowSums(assay(sce.hpc, "counts"))==0, ]


## Make list of Boolean param / cell subtype ===
cellSubtype.idx <- splitit(sce.hpc$cellType)
medianNon0.hpc <- lapply(cellSubtype.idx, function(x){
  apply(as.matrix(assay(sce.hpc, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

lengths(medianNon0.hpc)
sapply(medianNon0.hpc, head)

# Add these to the stats for each set of markers
for(i in names(markers.hpc.t.1vAll)){
  markers.hpc.t.1vAll[[i]] <- cbind(markers.hpc.t.1vAll[[i]],
                                    medianNon0.hpc[[i]][match(rownames(markers.hpc.t.1vAll[[i]]),
                                                           names(medianNon0.hpc[[i]]))])
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
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:15], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers, refined: single-nucleus-level p.w. t-tests, cluster-vs-all"))
  )
  dev.off()
}



## Do the same with the pairwise result ('markers.hpc.t.pw') === === ===
# Add these to the stats for each set of markers
for(i in names(markers.hpc.t.pw)){
  markers.hpc.t.pw[[i]] <- cbind(markers.hpc.t.pw[[i]],
                                    medianNon0.hpc[[i]][match(rownames(markers.hpc.t.pw[[i]]),
                                                           names(medianNon0.hpc[[i]]))])
  colnames(markers.hpc.t.pw[[i]])[17] <- "non0median"
}

markerList.t <- lapply(markers.hpc.t.pw, function(x){
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

write.csv(top40genes, file="tables/top40genesLists-REFINED_HPC-n3_cellType_Nov2020.csv",
          row.names=FALSE)




## Aside: add in 't.stat' as in 'step04' analyses to save for LoHu/LeCo ===
for(s in names(markers.hpc.t.1vAll)){
  markers.hpc.t.1vAll[[s]]$t.stat <- markers.hpc.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.hpc))
}

save(markers.hpc.t.1vAll, markers.hpc.t.pw, sce.hpc,
     file="rdas/markerStats-and-SCE_HPC-n3_sn-level_cleaned_MNTNov2020.rda")



### Session info for 21May2021 ============
sessionInfo()
    # R version 4.0.4 RC (2021-02-08 r79975)
    # Platform: x86_64-pc-linux-gnu (64-bit)
    # Running under: CentOS Linux 7 (Core)
    # 
    # Matrix products: default
    # BLAS:   /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/lib/libRblas.so
    # LAPACK: /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/lib/libRlapack.so
    # 
    # locale:
    #   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    # [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    # [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    # [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    # [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    # [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    # 
    # attached base packages:
    #   [1] parallel  stats4    stats     graphics  grDevices datasets  utils    
    # [8] methods   base     
    # 
    # other attached packages:
    #   [1] limma_3.46.0                jaffelab_0.99.30           
    # [3] rafalib_1.0.0               DropletUtils_1.10.3        
    # [5] batchelor_1.6.2             scran_1.18.5               
    # [7] scater_1.18.6               ggplot2_3.3.3              
    # [9] EnsDb.Hsapiens.v86_2.99.0   ensembldb_2.14.1           
    # [11] AnnotationFilter_1.14.0     GenomicFeatures_1.42.3     
    # [13] AnnotationDbi_1.52.0        SingleCellExperiment_1.12.0
    # [15] SummarizedExperiment_1.20.0 Biobase_2.50.0             
    # [17] GenomicRanges_1.42.0        GenomeInfoDb_1.26.7        
    # [19] IRanges_2.24.1              S4Vectors_0.28.1           
    # [21] BiocGenerics_0.36.1         MatrixGenerics_1.2.1       
    # [23] matrixStats_0.58.0         
    # 
    # loaded via a namespace (and not attached):
    #   [1] googledrive_1.0.1         ggbeeswarm_0.6.0         
    # [3] colorspace_2.0-0          ellipsis_0.3.2           
    # [5] scuttle_1.0.4             bluster_1.0.0            
    # [7] XVector_0.30.0            BiocNeighbors_1.8.2      
    # [9] rstudioapi_0.13           farver_2.1.0             
    # [11] bit64_4.0.5               fansi_0.4.2              
    # [13] xml2_1.3.2                splines_4.0.4            
    # [15] R.methodsS3_1.8.1         sparseMatrixStats_1.2.1  
    # [17] cachem_1.0.4              Rsamtools_2.6.0          
    # [19] ResidualMatrix_1.0.0      dbplyr_2.1.1             
    # [21] R.oo_1.24.0               HDF5Array_1.18.1         
    # [23] compiler_4.0.4            httr_1.4.2               
    # [25] dqrng_0.2.1               assertthat_0.2.1         
    # [27] Matrix_1.3-2              fastmap_1.1.0            
    # [29] lazyeval_0.2.2            BiocSingular_1.6.0       
    # [31] prettyunits_1.1.1         tools_4.0.4              
    # [33] rsvd_1.0.3                igraph_1.2.6             
    # [35] gtable_0.3.0              glue_1.4.2               
    # [37] GenomeInfoDbData_1.2.4    dplyr_1.0.5              
    # [39] rappdirs_0.3.3            Rcpp_1.0.6               
    # [41] vctrs_0.3.6               Biostrings_2.58.0        
    # [43] rhdf5filters_1.2.0        rtracklayer_1.50.0       
    # [45] DelayedMatrixStats_1.12.3 stringr_1.4.0            
    # [47] beachmat_2.6.4            lifecycle_1.0.0          
    # [49] irlba_2.3.3               statmod_1.4.35           
    # [51] XML_3.99-0.6              edgeR_3.32.1             
    # [53] zlibbioc_1.36.0           scales_1.1.1             
    # [55] hms_1.0.0                 ProtGenerics_1.22.0      
    # [57] rhdf5_2.34.0              RColorBrewer_1.1-2       
    # [59] curl_4.3                  memoise_2.0.0            
    # [61] gridExtra_2.3             segmented_1.3-3          
    # [63] biomaRt_2.46.3            stringi_1.5.3            
    # [65] RSQLite_2.2.7             BiocParallel_1.24.1      
    # [67] rlang_0.4.10              pkgconfig_2.0.3          
    # [69] bitops_1.0-7              lattice_0.20-41          
    # [71] purrr_0.3.4               Rhdf5lib_1.12.1          
    # [73] labeling_0.4.2            GenomicAlignments_1.26.0 
    # [75] cowplot_1.1.1             bit_4.0.4                
    # [77] tidyselect_1.1.1          magrittr_2.0.1           
    # [79] R6_2.5.0                  generics_0.1.0           
    # [81] DelayedArray_0.16.3       DBI_1.1.1                
    # [83] pillar_1.6.0              withr_2.4.2              
    # [85] RCurl_1.98-1.3            tibble_3.1.1             
    # [87] crayon_1.4.1              utf8_1.2.1               
    # [89] BiocFileCache_1.14.0      viridis_0.6.0            
    # [91] progress_1.2.2            locfit_1.5-9.4           
    # [93] grid_4.0.4                blob_1.2.1               
    # [95] digest_0.6.27             R.utils_2.10.1           
    # [97] openssl_1.4.3             munsell_0.5.0            
    # [99] beeswarm_0.3.1            viridisLite_0.4.0        
    # [101] vipor_0.4.5               askpass_1.1
