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


## Load SCE with new info
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    ## sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc, ref.sampleInfo

table(sce.sacc$cellType)

# First drop "Ambig.lowNtrxts" (43 nuclei)
sce.sacc <- sce.sacc[ ,sce.sacc$cellType != "Ambig.lowNtrxts"]
sce.sacc$cellType <- droplevels(sce.sacc$cellType)

# Remove 0 genes across all nuclei
sce.sacc <- sce.sacc[!rowSums(assay(sce.sacc, "counts"))==0, ]  # keeps same 28774 genes


## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.sacc), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`


# Run pairwise t-tests
markers.sacc.t.design <- findMarkers(sce.sacc, groups=sce.sacc$cellType,
                                    assay.type="logcounts", design=mod, test="t",
                                    direction="up", pval.type="all", full.stats=T)

sapply(markers.sacc.t.design, function(x){table(x$FDR<0.05)})
    #       Astro Excit.1 Excit.2 Excit.3 Excit.4 Inhib.1 Inhib.2 Micro Oligo   OPC
    # FALSE 27821   28246   28493   28378   27967   28455   28282 27319 28059 28272
    # TRUE    953     528     281     396     807     319     492  1455   715   502


## WMW: Blocking on donor (this test doesn't take 'design=' argument) ===
markers.sacc.wilcox.block <- findMarkers(sce.sacc, groups=sce.sacc$cellType,
                                        assay.type="logcounts", block=sce.sacc$donor, test="wilcox",
                                        direction="up", pval.type="all", full.stats=T)

# no warnings as in pan-brain analyses, but NO results of FDR<0.05...:
sapply(markers.sacc.wilcox.block, function(x){table(x$FDR<0.05)})
    #       Astro Excit.1 Excit.2 Excit.3 Excit.4 Inhib.1 Inhib.2 Micro Oligo   OPC
    # FALSE 28312   28377   28615   28600   28337   28581   28560 28129 28247 28476
    # TRUE    462     397     159     174     437     193     214   645   527   298


## Binomial ===
markers.sacc.binom.block <- findMarkers(sce.sacc, groups=sce.sacc$cellType,
                                       assay.type="logcounts", block=sce.sacc$donor, test="binom",
                                       direction="up", pval.type="all", full.stats=T)

sapply(markers.sacc.binom.block, function(x){table(x$FDR<0.05)})
    #       Astro Excit.1 Excit.2 Excit.3 Excit.4 Inhib.1 Inhib.2 Micro Oligo   OPC
    # FALSE 28610   28623   28743   28742   28619   28729   28712 28496 28630 28674
    # TRUE    164     151      31      32     155      45      62   278   144   100

## Save all these for future reference
save(markers.sacc.t.design, markers.sacc.wilcox.block, markers.sacc.binom.block,
     file="rdas/markers-stats_sACC-n2_findMarkers-SN-LEVEL_MNTMay2020.rda")


    # Btw (as in other regions) - some have 0 p.value's and FDR's
    head(markers.sacc.t.design[["Excit.3"]][ ,1:2])
        # SULF1 top gene with "0" p.value & FDR, but this is clearly a great marker
        # gene for thsi cluster, so these are probably thresholded at some point


# Print these to pngs
markerList.t.pw <- lapply(markers.sacc.t.design, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)

genes.top40.t <- lapply(markerList.t.pw, function(x){head(x, n=40)})


#dir.create("pdfs/exploration/sACC/")
for(i in names(genes.top40.t)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/sACC/sACC_t-sn-level_pairwise_top40markers-", i, "_logExprs_May2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.sacc, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:10], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers: single-nucleus-level p.w. t-tests"))
  )
  dev.off()
}






### Cluster-vs-all single-nucleus-level iteration ======

## Load SCE with new info
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    ## sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc, ref.sampleInfo

table(sce.sacc$cellType)

# First drop "Ambig.lowNtrxts" (43 nuclei)
sce.sacc <- sce.sacc[ ,sce.sacc$cellType != "Ambig.lowNtrxts"]
sce.sacc$cellType <- droplevels(sce.sacc$cellType)

# Remove 0 genes across all nuclei
sce.sacc <- sce.sacc[!rowSums(assay(sce.sacc, "counts"))==0, ]  # keeps same 28774 genes


## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.sacc), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

markers.sacc.t.1vAll <- list()
for(i in levels(sce.sacc$cellType)){
  # Make temporary contrast
  sce.sacc$contrast <- ifelse(sce.sacc$cellType==i, 1, 0)
  # Test cluster vs. all
  markers.sacc.t.1vAll[[i]] <- findMarkers(sce.sacc, groups=sce.sacc$contrast,
                                          assay.type="logcounts", design=mod, test="t",
                                          direction="up", pval.type="all", full.stats=T)
}

    ## Then, temp set of stats to get the standardized logFC
    temp.1vAll <- list()
    for(i in levels(sce.sacc$cellType)){
      # Make temporary contrast
      sce.sacc$contrast <- ifelse(sce.sacc$cellType==i, 1, 0)
      # Test cluster vs. all
      temp.1vAll[[i]] <- findMarkers(sce.sacc, groups=sce.sacc$contrast,
                                     assay.type="logcounts", design=mod, test="t",
                                     std.lfc=TRUE,
                                     direction="up", pval.type="all", full.stats=T)
    }


## For some reason all the results are in the second List entry (first is always empty)

# Replace that empty slot with the entry with the actul stats
markers.sacc.t.1vAll <- lapply(markers.sacc.t.1vAll, function(x){ x[[2]] })
# Same for that with std.lfc
temp.1vAll <- lapply(temp.1vAll, function(x){ x[[2]] })

# Now just pull from the 'stats.0' DataFrame column
markers.sacc.t.1vAll <- lapply(markers.sacc.t.1vAll, function(x){ x$stats.0 })
temp.1vAll <- lapply(temp.1vAll, function(x){ x$stats.0 })

# Re-name std.lfc column and add to the first result
for(i in names(temp.1vAll)){
  colnames(temp.1vAll[[i]])[1] <- "std.logFC"
  markers.sacc.t.1vAll[[i]] <- cbind(markers.sacc.t.1vAll[[i]], temp.1vAll[[i]]$std.logFC)
  # Oh the colname is kept weird
  colnames(markers.sacc.t.1vAll[[i]])[4] <- "std.logFC"
  # Then re-organize
  markers.sacc.t.1vAll[[i]] <- markers.sacc.t.1vAll[[i]][ ,c("logFC","std.logFC","log.p.value","log.FDR")]
}


## Let's save this along with the previous pairwise results
save(markers.sacc.t.design, markers.sacc.wilcox.block, markers.sacc.binom.block,  # pairwise set
     markers.sacc.t.1vAll,
     file="rdas/markers-stats_sACC-n2_findMarkers-SN-LEVEL_MNTMay2020.rda")


## Print these to pngs
markerList.t.1vAll <- lapply(markers.sacc.t.1vAll, function(x){
  rownames(x)[x[ ,"log.FDR"] < log10(0.01)]
  }
)
genes.top40.t.1vAll <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t.1vAll)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/sACC/sACC_t-sn-level_1vALL_top40markers-",i,"_logExprs_May2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.sacc, exprs_values = "logcounts", features=genes.top40.t.1vAll[[i]],
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:10], length(genes.top40.t.1vAll[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers: single-nucleus-level t-tests, cluster-vs-all"))
  )
  dev.off()
}




## How do these top 40 intersect? ===
sapply(names(genes.top40.t), function(c){
  length(intersect(genes.top40.t[[c]],
                   genes.top40.t.1vAll[[c]]))
})
    #  Astro Excit.1 Excit.2 Excit.3 Excit.4 Inhib.1 Inhib.2   Micro   Oligo     OPC
    #     36      21      21      24      32      20      29      36      31      30



## Write these top 40 lists to a csv
names(markerList.t.pw) <- paste0(names(markerList.t.pw),"_pw")
names(markerList.t.1vAll) <- paste0(names(markerList.t.1vAll),"_1vAll")

# PW result for "Inhib.1" doesn't have 40 markers:
#markerList.t.pw[["Inhib.1_pw"]] <- c(markerList.t.pw[["Inhib.1_pw"]], rep("",9))

top40genes <- cbind(sapply(markerList.t.pw, function(x) head(x, n=40)),
                    sapply(markerList.t.1vAll, function(y) head(y, n=40)))
top40genes <- top40genes[ ,sort(colnames(top40genes))]

write.csv(top40genes, file="tables/top40genesLists_sACC-n2_cellType_SN-LEVEL-tests_May2020.csv",
          row.names=FALSE)




### MNT add 18Nov2020 =================================
# -> What if add param/requirement that for any given subcluster, median expression has to > 0?
load("rdas/markers-stats_sACC-n2_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.sacc.t.design, markers.sacc.wilcox.block, markers.sacc.binom.block, markers.sacc.t.1vAll   

## Load SCE 
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc, ref.sampleInfo

table(sce.sacc$cellType)

# First drop "Ambig.lowNtrxts" (43 nuclei)
sce.sacc <- sce.sacc[ ,sce.sacc$cellType != "Ambig.lowNtrxts"]
sce.sacc$cellType <- droplevels(sce.sacc$cellType)

# Remove 0 genes across all nuclei
sce.sacc <- sce.sacc[!rowSums(assay(sce.sacc, "counts"))==0, ]


## Make list of Boolean param / cell subtype ===
cellSubtype.idx <- splitit(sce.sacc$cellType)
medianNon0.idx <- lapply(cellSubtype.idx, function(x){
  apply(as.matrix(assay(sce.sacc, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

lengths(medianNon0.idx)
sapply(medianNon0.idx, head)

# Add these to the stats for each set of markers
for(i in names(markers.sacc.t.1vAll)){
  markers.sacc.t.1vAll[[i]] <- cbind(markers.sacc.t.1vAll[[i]],
                                    medianNon0.idx[[i]][match(rownames(markers.sacc.t.1vAll[[i]]),
                                                           names(medianNon0.idx[[i]]))])
  colnames(markers.sacc.t.1vAll[[i]])[5] <- "non0median"
}


## Use these restrictions to print (to png) a refined top 40, as before ===
markerList.t.1vAll <- lapply(markers.sacc.t.1vAll, function(x){
  rownames(x)[x$log.FDR < log10(0.000001) & x$non0median==TRUE]
  }
)
    # lengths(markerList.t.1vAll)     # ( **without $non0median==TRUE restriction )
        #  Astro Excit.1 Excit.2 Excit.3 Excit.4 Inhib.1 Inhib.2   Micro   Oligo     OPC
        #   3468    9706    8983    5462    4475    8084    7451    3341    2306    3077

lengths(markerList.t.1vAll)
    #  Astro Excit.1 Excit.2 Excit.3 Excit.4 Inhib.1 Inhib.2   Micro   Oligo     OPC
    #    835    4745    4621    3371    2894    3891    3558     714     843    1159

genes.top40.t <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/sACC/sACC_t-sn-level_1vALL_top40markers-REFINED-",gsub(":",".",i),"_logExprs_Nov2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.sacc, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:10], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers, refined: single-nucleus-level p.w. t-tests, cluster-vs-all"))
  )
  dev.off()
}



## Do the same with the pairwise result ('markers.sacc.t.design') === === ===
# Add these to the stats for each set of markers
for(i in names(markers.sacc.t.design)){
  markers.sacc.t.design[[i]] <- cbind(markers.sacc.t.design[[i]],
                                     medianNon0.idx[[i]][match(rownames(markers.sacc.t.design[[i]]),
                                                            names(medianNon0.idx[[i]]))])
  colnames(markers.sacc.t.design[[i]])[12] <- "non0median"
}

markerList.t <- lapply(markers.sacc.t.design, function(x){
  rownames(x)[x$FDR < 0.05 & x$non0median==TRUE]
  }
)
    # lengths(markerList.t)     # ( **without $non0median==TRUE restriction )
        #  Astro Excit.1 Excit.2 Excit.3 Excit.4 Inhib.1 Inhib.2   Micro   Oligo     OPC
        #    953     528     281     396     807     319     492    1455     715     502

lengths(markerList.t)
    #  Astro Excit.1 Excit.2 Excit.3 Excit.4 Inhib.1 Inhib.2   Micro   Oligo     OPC
    #    353     272     132     137     322     182     208     417     439     234


genes.top40.t <- lapply(markerList.t, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/sACC/sACC_t-sn-level_pairwise_top40markers-REFINED-", i, "_logExprs_Nov2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.sacc, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:10], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers, refined: single-nucleus-level p.w. t-tests"))
  )
  dev.off()
}

## Then write a new CSV of these refined top 40 genes ===
names(markerList.t) <- paste0(names(markerList.t),"_pw")
names(markerList.t.1vAll) <- paste0(names(markerList.t.1vAll),"_1vAll")

# # Many of the PW results don't have at least 40 markers:      - N/A for sACC
# extend.idx <- names(which(lengths(markerList.t) < 40))
# for(i in extend.idx){
#   markerList.t[[i]] <- c(markerList.t[[i]], rep("", 40-length(markerList.t[[i]])))
# }

top40genes <- cbind(sapply(markerList.t, function(x) head(x, n=40)),
                    sapply(markerList.t.1vAll, function(y) head(y, n=40)))
top40genes <- top40genes[ ,sort(colnames(top40genes))]

write.csv(top40genes, file="tables/top40genesLists-REFINED_sACC-n2_cellType.split_Nov2020.csv",
          row.names=FALSE)



## Aside: add in 't.stat' as in 'step04' analyses to save for LoHu/LeCo ===
for(s in names(markers.sacc.t.1vAll)){
  markers.sacc.t.1vAll[[s]]$t.stat <- markers.sacc.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.sacc))
}

save(markers.sacc.t.1vAll, markers.sacc.t.design, sce.sacc,
     file="rdas/markerStats-and-SCE_sACC-n2_sn-level_cleaned_MNTNov2020.rda")


