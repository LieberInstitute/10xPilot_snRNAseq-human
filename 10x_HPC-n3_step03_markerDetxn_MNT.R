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




## Aside: add in 't.stat' as in 'step04' analyses to save for LoHu/LeCo ===
for(s in names(markers.hpc.t.1vAll)){
  markers.hpc.t.1vAll[[s]]$t.stat <- markers.hpc.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.hpc))
}

save(markers.hpc.t.1vAll, markers.hpc.t.design, sce.hpc,
     file="rdas/markerStats-and-SCE_HPC-n3_sn-level_cleaned_MNTNov2020.rda")




