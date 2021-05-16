### MNT 10x snRNA-seq workflow: step 03 - marker detection
###   **Region-specific analyses**
###     - (2x) DLPFC samples from: Br5161 & Br5212
### Initiated MNT 12Feb2020 - modified MNT 05Mar2020
### Modification notes: First drop "Ambig.lowNtrxts" cluster
###         (see bottom chunk for deets)
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



### Single-nucleus-level tests for cell-type-specific genes ================================
  # MNT 23Apr2020 - added after seeing much better results in pan-brain analysis

## Load SCE with new info
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_DLPFC-n2_SCE_cellTypesSplit-fromST_Apr2020.rda", verbose=T)
    # sce.dlpfc.st, clusterRefTab.dlpfc, chosen.hvgs.dlpfc, ref.sampleInfo

table(sce.dlpfc.st$cellType.split)

# First drop "Ambig.lowNtrxts" (168 nuclei)
sce.dlpfc.st <- sce.dlpfc.st[ ,sce.dlpfc.st$cellType.split != "Ambig.lowNtrxts"]
sce.dlpfc.st$cellType.split <- droplevels(sce.dlpfc.st$cellType.split)

# Remove 0 genes across all nuclei
sce.dlpfc.st <- sce.dlpfc.st[!rowSums(assay(sce.dlpfc.st, "counts"))==0, ]  # keeps same 28111 genes


## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.dlpfc.st), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # 'drop=F' to keep as matrix - otherwise turns into numeric vec
    ## Get: "Error in .ranksafe_qr(full.design) : design matrix is not of full rank"
     #   if try to put 0 intercept

    ## MNT comment: actually doesn't matter

# Run pairwise t-tests
markers.dlpfc.t.design <- findMarkers(sce.dlpfc.st, groups=sce.dlpfc.st$cellType.split,
                                    assay.type="logcounts", design=mod, test="t",
                                    direction="up", pval.type="all", full.stats=T)

sapply(markers.dlpfc.t.design, function(x){table(x$FDR<0.05)})
    #       Astro Excit.ambig Excit.L2:3 Excit.L3:4 Excit.L4:5 Excit.L5 Excit.L5:6
    # FALSE 27877       28034      28109      28015      28086    27883      28063
    # TRUE    234          77          2         96         25      228         48
    #       Excit.L6.broad Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5 Inhib.6 Micro Oligo   OPC
    # FALSE          28073   27866   27911   27973   28084   28081   28083 27571 27958 27955
    # TRUE              38     245     200     138      27      30      28   540   153   156


## WMW: Blocking on donor (this test doesn't take 'design=' argument) ===
markers.dlpfc.wilcox.block <- findMarkers(sce.dlpfc.st, groups=sce.dlpfc.st$cellType.split,
                                        assay.type="logcounts", block=sce.dlpfc.st$donor, test="wilcox",
                                        direction="up", pval.type="all", full.stats=T)
    
    # no warnings as in pan-brain analyses, but NO results of FDR<0.05...:
sapply(markers.dlpfc.wilcox.block, function(x){table(x$FDR<0.05)})
    
    sapply(names(markers.dlpfc.wilcox.block), function(x){quantile(markers.dlpfc.wilcox.block[[x]]$FDR)})
        # Astro Excit.ambig Excit.L2:3 Excit.L3:4 Excit.L4:5  Excit.L5 Excit.L5:6 Excit.L6.broad
        # 0%   0.3332772           1          1          1          1 0.6095081          1              1
        # 25%  1.0000000           1          1          1          1 1.0000000          1              1
        # 50%  1.0000000           1          1          1          1 1.0000000          1              1
        # 75%  1.0000000           1          1          1          1 1.0000000          1              1
        # 100% 1.0000000           1          1          1          1 1.0000000          1              1
        # Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5 Inhib.6     Micro     Oligo       OPC
        # 0%         1       1       1       1       1       1 0.1850964 0.1071413 0.1395478
        # 25%        1       1       1       1       1       1 1.0000000 1.0000000 1.0000000
        # 50%        1       1       1       1       1       1 1.0000000 1.0000000 1.0000000
        # 75%        1       1       1       1       1       1 1.0000000 1.0000000 1.0000000
        # 100%       1       1       1       1       1       1 1.0000000 1.0000000 1.0000000


## Binomial ===
markers.dlpfc.binom.block <- findMarkers(sce.dlpfc.st, groups=sce.dlpfc.st$cellType.split,
                                       assay.type="logcounts", block=sce.dlpfc.st$donor, test="binom",
                                       direction="up", pval.type="all", full.stats=T)

sapply(markers.dlpfc.binom.block, function(x){table(x$FDR<0.05)})
    # none - and these are ALL 1's across the board....

## Save all these for future reference
save(markers.dlpfc.t.design, #markers.dlpfc.wilcox.block, markers.dlpfc.binom.block,
     file="rdas/markers-stats_DLPFC_n2_findMarkers-SN-LEVEL_MNTApr2020.rda")


# Print these to pngs
markerList.t <- lapply(markers.dlpfc.t.design, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)

genes.top40.t <- lapply(markerList.t, function(x){head(x, n=40)})
smallerSets <- c("Excit.L4:5", "Inhib.4", "Inhib.5", "Inhib.6")

#dir.create("pdfs/exploration/DLPFC/")
## ~40+ marker genes
for(i in setdiff(names(genes.top40.t), c(smallerSets, "Excit.L2:3"))){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/DLPFC/DLPFC_t-sn-level-top40markers-",gsub(":",".",i),"_logExprs_Apr2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.dlpfc.st, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:17], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers: single-nucleus-level p.w. t-tests"))
  )
  dev.off()
}

## < 40 ('smallerSets')
for(i in smallerSets){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/DLPFC/DLPFC_t-sn-level-top40markers-",gsub(":",".",i),"_logExprs_Apr2020.png"), height=1900/2, width=1200)
  print(
    plotExpression(sce.dlpfc.st, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:17], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top markers with FDR < 0.05: single-nucleus-level p.w. t-tests"))
  )
  dev.off()
}

## 'Excit.L2:3'
png("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/DLPFC/DLPFC_t-sn-level-top40markers-Excit.L2.3_logExprs_Apr2020.png", height=1900/8, width=1200/3)
print(
  plotExpression(sce.dlpfc.st, exprs_values = "logcounts", features=genes.top40.t[["Excit.L2:3"]],
                 x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, #ncol=5,
                 add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                              geom = "crossbar", width = 0.3,
                                              colour=rep(tableau20[1:17], length(genes.top40.t[["Excit.L2:3"]]))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 10)) +  
    ggtitle(label="Excit.L2:3 top markers with FDR < 0.05: single-nucleus-level p.w. t-tests")
)
dev.off()






### Cluster-vs-all single-nucleus-level iteration ================================
  # Update MNT 05May2020: add iteration to use and report the standardized logFC (Cohen's D)
  # -> see 'side-Rscript_markerDetxn-analyses_comparisons-to-sn-levelstats_Apr2020.R'
  #    which explores this

## Load SCE with new info
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_DLPFC-n2_SCE_cellTypesSplit-fromST_Apr2020.rda", verbose=T)
    # sce.dlpfc.st, clusterRefTab.dlpfc, chosen.hvgs.dlpfc, ref.sampleInfo

table(sce.dlpfc.st$cellType.split)

# First drop "Ambig.lowNtrxts" (168 nuclei)
sce.dlpfc.st <- sce.dlpfc.st[ ,sce.dlpfc.st$cellType.split != "Ambig.lowNtrxts"]
sce.dlpfc.st$cellType.split <- droplevels(sce.dlpfc.st$cellType.split)

# Remove 0 genes across all nuclei
sce.dlpfc.st <- sce.dlpfc.st[!rowSums(assay(sce.dlpfc.st, "counts"))==0, ]  # keeps same 28111 genes


## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.dlpfc.st), model.matrix(~ donor))

markers.dlpfc.t.1vAll <- list()
for(i in unique(sce.dlpfc.st$cellType.split)){
  # Make temporary contrast
  sce.dlpfc.st$contrast <- ifelse(sce.dlpfc.st$cellType.split==i, 1, 0)
  # Test cluster vs. all
  markers.dlpfc.t.1vAll[[i]] <- findMarkers(sce.dlpfc.st, groups=sce.dlpfc.st$contrast,
                                assay.type="logcounts", design=mod, test="t",
                                direction="up", pval.type="all", full.stats=T)
}
    # There were 17 warnings (use warnings() to see them)
    #   [17x]:
    #   Warning messages:
    #   1: In .fit_lm_internal(x, subset.row, groups, design = design,  ... :
    #     automatically removed intercept column

## Temp set of stats to get the standardized logFC
temp.1vAll <- list()
for(i in unique(sce.dlpfc.st$cellType.split)){
  # Make temporary contrast
  sce.dlpfc.st$contrast <- ifelse(sce.dlpfc.st$cellType.split==i, 1, 0)
  # Test cluster vs. all
  temp.1vAll[[i]] <- findMarkers(sce.dlpfc.st, groups=sce.dlpfc.st$contrast,
                                 assay.type="logcounts", design=mod, test="t",
                                 std.lfc=TRUE,
                                 direction="up", pval.type="all", full.stats=T)
}


class(markers.dlpfc.t.1vAll[["Astro"]]) # SimpleList
dim(markers.dlpfc.t.1vAll[["Astro"]]) # NULL
head(markers.dlpfc.t.1vAll[["Astro"]])
    ## For some reason all the results are in the second List entry (first is always empty)

    table(markers.dlpfc.t.1vAll[["Astro"]][[2]]$FDR<0.05) # 4892
    table(markers.dlpfc.t.1vAll[["Astro"]][[2]]$stats.0$log.FDR < log10(0.05))  # 6143
        ## Why are these different?  Though probably should take the former, because
         # there is only one test for each of these clusters and not sure
         # WHAT these are...

    head(markers.dlpfc.t.1vAll[["Oligo"]][[2]])
        ## Nice, MBP and PLP1 are in the top 6
    

sapply(markers.dlpfc.t.1vAll, function(x){
  table(x[[2]]$stats.0$log.FDR < log10(.001))
  })
    #       Oligo Astro Inhib.4 Excit.L4:5 Micro Inhib.6   OPC Excit.L2:3 Excit.ambig
    # FALSE 25366 23223   21034      19061 23999   20956 23802      22479       20869
    # TRUE   2745  4888    7077       9050  4112    7155  4309       5632        7242
    #       Excit.L3:4 Excit.L5:6 Excit.L6.broad Inhib.5 Excit.L5 Inhib.1 Inhib.2
    # FALSE      23916      22320          21501   23246    23878   26699   26540
    # TRUE        4195       5791           6610    4865     4233    1412    1571
    #       Inhib.3
    # FALSE   25766
    # TRUE     2345

# Replace that empty slot with the entry with the actul stats
markers.dlpfc.t.1vAll <- lapply(markers.dlpfc.t.1vAll, function(x){ x[[2]] })
# Same for that with std.lfc
temp.1vAll <- lapply(temp.1vAll, function(x){ x[[2]] })

# Now just pull from the 'stats.0' DataFrame column
markers.dlpfc.t.1vAll <- lapply(markers.dlpfc.t.1vAll, function(x){ x$stats.0 })
temp.1vAll <- lapply(temp.1vAll, function(x){ x$stats.0 })

# Re-name std.lfc column and add to the first result
for(i in names(temp.1vAll)){
  colnames(temp.1vAll[[i]])[1] <- "std.logFC"
  markers.dlpfc.t.1vAll[[i]] <- cbind(markers.dlpfc.t.1vAll[[i]], temp.1vAll[[i]]$std.logFC)
  # Oh the colname is kept weird
  colnames(markers.dlpfc.t.1vAll[[i]])[4] <- "std.logFC"
  # Then re-organize
  markers.dlpfc.t.1vAll[[i]] <- markers.dlpfc.t.1vAll[[i]][ ,c("logFC","std.logFC","log.p.value","log.FDR")]
}


  
  
## Let's save this along with the previous pairwise results
load("rdas/markers-stats_DLPFC_n2_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
save(markers.dlpfc.t.1vAll, markers.dlpfc.t.design,
     file="rdas/markers-stats_DLPFC_n2_findMarkers-SN-LEVEL_MNTApr2020.rda")




## Print these to pngs
markerList.t.1vAll <- lapply(markers.dlpfc.t.1vAll, function(x){
  rownames(x)[x$log.FDR < log10(0.000001)]
  }
)
genes.top40.t <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/DLPFC/DLPFC_t-sn-level_1vALL_top40markers-",gsub(":",".",i),"_logExprs_Apr2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.dlpfc.st, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:17], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers: single-nucleus-level p.w. t-tests, cluster-vs-all"))
  )
  dev.off()
}


## How do they intersect?
markerList.t.pw <- lapply(markers.dlpfc.t.design, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)

# From pairwise t-tests, FDR < 0.05
lengths(markerList.t.pw)
    #         Astro    Excit.ambig     Excit.L2:3     Excit.L3:4     Excit.L4:5       Excit.L5
    #           234             77              2             96             25            228
    #    Excit.L5:6 Excit.L6.broad        Inhib.1        Inhib.2        Inhib.3        Inhib.4
    #            48             38            245            200            138             27
    #       Inhib.5        Inhib.6          Micro          Oligo            OPC
    #            30             28            540            153            156


# From cluster-vs-all others, FDR < 1e6
lengths(markerList.t.1vAll)
    # thousands

# Intersection
sapply(names(markerList.t.pw), function(c){
  length(intersect(markerList.t.pw[[c]],
                   markerList.t.1vAll[[c]]))
})



## Write these top 40 lists to a csv
names(markerList.t.pw) <- paste0(names(markerList.t.pw),"_pw")
names(markerList.t.1vAll) <- paste0(names(markerList.t.1vAll),"_1vAll")

# Many of the PW results don't have 40 markers:
extend.idx <- names(which(lengths(markerList.t.pw) < 40))
for(i in extend.idx){
  markerList.t.pw[[i]] <- c(markerList.t.pw[[i]], rep("", 40-length(markerList.t.pw[[i]])))
}

top40genes <- cbind(sapply(markerList.t.pw, function(x) head(x, n=40)),
                    sapply(markerList.t.1vAll, function(y) head(y, n=40)))
top40genes <- top40genes[ ,sort(colnames(top40genes))]

write.csv(top40genes, file="tables/top40genesLists_DLPFC-n2_cellType.split_SN-LEVEL-tests_May2020.csv",
          row.names=FALSE)




### MNT add 02May2021 =================================
  # -> What if add param/requirement that for any given subcluster, median expression has to > 0?
load("rdas/markers-stats_DLPFC_n2_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.dlpfc.t.1vAll, markers.dlpfc.t.design

## Load SCE 
load("/rdas/regionSpecific_DLPFC-n2_SCE_cellTypesSplit-fromST_Apr2020.rda", verbose=T)
    # sce.dlpfc.st, clusterRefTab.dlpfc, chosen.hvgs.dlpfc, ref.sampleInfo

table(sce.dlpfc.st$cellType.split)

# First drop "Ambig.lowNtrxts" (168 nuclei)
sce.dlpfc.st <- sce.dlpfc.st[ ,sce.dlpfc.st$cellType.split != "Ambig.lowNtrxts"]
sce.dlpfc.st$cellType.split <- droplevels(sce.dlpfc.st$cellType.split)

# Remove 0 genes across all nuclei
sce.dlpfc.st <- sce.dlpfc.st[!rowSums(assay(sce.dlpfc.st, "counts"))==0, ]


## Make list of Boolean param / cell subtype ===
cellSubtype.idx <- splitit(sce.dlpfc.st$cellType.split)
medianNon0.idx <- lapply(cellSubtype.idx, function(x){
  apply(as.matrix(assay(sce.dlpfc.st, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

lengths(medianNon0.idx)
sapply(medianNon0.idx, head)

# Add these to the stats for each set of markers
for(i in names(markers.dlpfc.t.1vAll)){
  markers.dlpfc.t.1vAll[[i]] <- cbind(markers.dlpfc.t.1vAll[[i]],
                                    medianNon0.idx[[i]][match(rownames(markers.dlpfc.t.1vAll[[i]]),
                                                              names(medianNon0.idx[[i]]))])
  colnames(markers.dlpfc.t.1vAll[[i]])[5] <- "non0median"
}


## Use these restrictions to print (to png) a refined top 40, as before ===
markerList.t.1vAll <- lapply(markers.dlpfc.t.1vAll, function(x){
  rownames(x)[x$log.FDR < log10(0.000001) & x$non0median==TRUE]
  }
)
lengths(markerList.t.1vAll)     # ( **without $non0median==TRUE restriction )
    #         Oligo          Astro        Inhib.4     Excit.L4:5          Micro        Inhib.6 
    #          2337           3712           5051           6623           3300           4935 
    #           OPC     Excit.L2:3    Excit.ambig     Excit.L3:4     Excit.L5:6 Excit.L6.broad 
    #          3198           4015           5068           2774           4020           4577 
    #       Inhib.5       Excit.L5        Inhib.1        Inhib.2        Inhib.3 
    #          3451           2807            763            907           1469

lengths(markerList.t.1vAll)
    #         Oligo          Astro        Inhib.4     Excit.L4:5          Micro        Inhib.6 
    #           837            612           2832           3436            627           2311 
    #           OPC     Excit.L2:3    Excit.ambig     Excit.L3:4     Excit.L5:6 Excit.L6.broad 
    #          1059           2030           2775           1949           2289           2850 
    #       Inhib.5       Excit.L5        Inhib.1        Inhib.2        Inhib.3 
    #          1432           1756            433            507            986 

    ## (Don't bother printing top 40 markers with this restriction - will change with revision analysis)

## Do the same with the pairwise result ('markers.dlpfc.t.design') === === ===
# Add these to the stats for each set of markers
for(i in names(markers.dlpfc.t.design)){
  markers.dlpfc.t.design[[i]] <- cbind(markers.dlpfc.t.design[[i]],
                                     medianNon0.idx[[i]][match(rownames(markers.dlpfc.t.design[[i]]),
                                                               names(medianNon0.idx[[i]]))])
  colnames(markers.dlpfc.t.design[[i]])[19] <- "non0median"
}

markerList.t <- lapply(markers.dlpfc.t.design, function(x){
  rownames(x)[x$FDR < 0.05 & x$non0median==TRUE]
  }
)
# lengths(markerList.t)     # ( **without $non0median==TRUE restriction )
    #         Astro    Excit.ambig     Excit.L2:3     Excit.L3:4     Excit.L4:5       Excit.L5 
    #           234             77              2             96             25            228 
    #    Excit.L5:6 Excit.L6.broad        Inhib.1        Inhib.2        Inhib.3        Inhib.4 
    #            48             38            245            200            138             27 
    #       Inhib.5        Inhib.6          Micro          Oligo            OPC 
    #            30             28            540            153            156 

lengths(markerList.t)
    #         Astro    Excit.ambig     Excit.L2:3     Excit.L3:4     Excit.L4:5       Excit.L5 
    #           162             52              1             27             18            100 
    #    Excit.L5:6 Excit.L6.broad        Inhib.1        Inhib.2        Inhib.3        Inhib.4 
    #            26             26             57             23             62             21 
    #       Inhib.5        Inhib.6          Micro          Oligo            OPC 
    #            15              9            240            153            118 


    ## (Don't bother printing top 40 markers with this restriction - will change with revision analysis)


## Aside: add in 't.stat' as in 'step04' analyses to save for LoHu/LeCo ===
for(s in names(markers.dlpfc.t.1vAll)){
  markers.dlpfc.t.1vAll[[s]]$t.stat <- markers.dlpfc.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.dlpfc.st))
}

Readme <- "Making this iteration of marker/SCE objects so that can test MAGMA results with the 'non0median' filter, but for the preprint data/clusters (this was originally skipped for DLPFC) -MNT 02May2021"

save(markers.dlpfc.t.1vAll, markers.dlpfc.t.design, sce.dlpfc.st, Readme,
     file="rdas/markerStats-and-SCE_DLPFC-n2_sn-level_cleaned_MNTMay2021.rda")







