### MNT 10x snRNA-seq workflow: step 03 - marker detection
###   **Region-specific analyses**
### Initiated MNT 12Feb2020
### For revision 2021: (5x) amygdala samples
#####################################################################

library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(batchelor)
library(DropletUtils)
library(jaffelab)
library(limma)

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
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_Amyg-n5_cleaned-combined_SCE_MNT2021.rda",
     verbose=T)
    # sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo, annotationTab.amy, cell_colors.amy

table(sce.amy$cellType)
    # Astro_A       Astro_B drop.lowNTx_A drop.lowNTx_B          Endo       Excit_A
    #    1555            83          1067            71            31           344
    # Excit_B       Excit_C       Inhib_A       Inhib_B       Inhib_C       Inhib_D
    #      44            55           728           541           525           555
    # Inhib_E       Inhib_F       Inhib_G       Inhib_H         Micro         Mural
    #     414           216            86            52          1201            39
    #   Oligo           OPC         Tcell
    #    6080          1459            31


# First drop "drop.lowNTx_" (1138 nuclei)
sce.amy <- sce.amy[ ,-grep("drop.",sce.amy$cellType)]
sce.amy$cellType <- droplevels(sce.amy$cellType)

# Remove 0 genes across all nuclei
sce.amy <- sce.amy[!rowSums(assay(sce.amy, "counts"))==0, ]  # keeps same 29371 genes


## Re-create 'logcounts' (don't want to use 'multiBatchNorm's down-scaling across donor 'batches')
 # First 'hold' the MBN 'logcounts' for printing
sce.hold <- sce.amy

assay(sce.amy, "logcounts") <- NULL
sizeFactors(sce.amy) <- NULL
sce.amy <- logNormCounts(sce.amy)


### First make a list of Boolean param / cell subtype ===
# Will use this to assess more 'valid', non-noise-driving markers
cellSubtype.idx <- splitit(sce.amy$cellType)
medianNon0.amy <- lapply(cellSubtype.idx, function(x){
  apply(as.matrix(assay(sce.amy, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

sapply(medianNon0.amy, table)


## Traditional t-test implementation ===
mod <- with(colData(sce.amy), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

# Run pairwise t-tests
markers.amy.t.pw <- findMarkers(sce.amy, groups=sce.amy$cellType,
                                assay.type="logcounts", design=mod, test="t",
                                direction="up", pval.type="all", full.stats=T)

sapply(markers.amy.t.pw, function(x){table(x$FDR<0.05)})
    #       Astro_A Astro_B  Endo Excit_A Excit_B Excit_C Inhib_A Inhib_B Inhib_C Inhib_D
    # FALSE   28862   29311 28824   29228   28801   29061   29338   29202   29338   29187
    # TRUE      509      60   547     143     570     310      33     169      33     184
    #       Inhib_E Inhib_F Inhib_G Inhib_H Micro Mural Oligo   OPC Tcell
    # FALSE   29095   29265   29283   29223 28703 28807 28987 29109 28764
    # TRUE      276     106      88     148   668   564   384   262   607


## WMW: Blocking on donor (this test doesn't take 'design=' argument) ===
markers.amy.wilcox.block <- findMarkers(sce.amy, groups=sce.amy$cellType,
                                        assay.type="logcounts", block=sce.amy$donor, test="wilcox",
                                        direction="up", pval.type="all", full.stats=T)

# WMW FDR<0.05
sapply(markers.amy.wilcox.block, function(x){table(x$FDR<0.05)["TRUE"]})
    #Astro_A.TRUE   Astro_B.NA      Endo.NA   Excit_A.NA   Excit_B.NA   Excit_C.NA
    #         199           NA           NA           NA           NA           NA
    #  Inhib_A.NA Inhib_B.TRUE Inhib_C.TRUE Inhib_D.TRUE   Inhib_E.NA Inhib_F.TRUE
    #          NA           31            2           60           NA            4
    #  Inhib_G.NA   Inhib_H.NA   Micro.TRUE     Mural.NA   Oligo.TRUE     OPC.TRUE
    #          NA           NA           65           NA          245          125
    #  Tcell.TRUE
    #           4


## Binomial ===
markers.amy.binom.block <- findMarkers(sce.amy, groups=sce.amy$cellType,
                                       assay.type="logcounts", block=sce.amy$donor, test="binom",
                                       direction="up", pval.type="all", full.stats=T)

sapply(markers.amy.binom.block, function(x){table(x$FDR<0.05)["TRUE"]})
    # only a few hits for 'Astro_A'- disregard these
    #Astro_A.TRUE   Astro_B.NA      Endo.NA   Excit_A.NA   Excit_B.NA   Excit_C.NA
    #           6           NA           NA           NA           NA           NA
    #  Inhib_A.NA   Inhib_B.NA   Inhib_C.NA   Inhib_D.NA   Inhib_E.NA   Inhib_F.NA
    #          NA           NA           NA           NA           NA           NA
    #  Inhib_G.NA   Inhib_H.NA     Micro.NA     Mural.NA     Oligo.NA       OPC.NA
    #          NA           NA           NA           NA           NA           NA
    #    Tcell.NA
    #          NA


# Add respective 'non0median' column to the stats for each set of markers
for(i in names(markers.amy.t.pw)){
  markers.amy.t.pw[[i]] <- cbind(markers.amy.t.pw[[i]],
                                 medianNon0.amy[[i]][match(rownames(markers.amy.t.pw[[i]]),
                                                           names(medianNon0.amy[[i]]))])
  colnames(markers.amy.t.pw[[i]])[22] <- "non0median"
}

sapply(markers.amy.t.pw, function(x){table(x$FDR<0.05 & x$non0median == TRUE)["TRUE"]})
    # Astro_A.TRUE Astro_B.TRUE    Endo.TRUE Excit_A.TRUE Excit_B.TRUE Excit_C.TRUE
    #         326            3          180           59          173          142
    # Inhib_A.TRUE Inhib_B.TRUE Inhib_C.TRUE Inhib_D.TRUE Inhib_E.TRUE Inhib_F.TRUE
    #           26           88           14          131          147           43
    # Inhib_G.TRUE Inhib_H.TRUE   Micro.TRUE   Mural.TRUE   Oligo.TRUE     OPC.TRUE
    #            4           79          284          172          306          172
    #   Tcell.TRUE
    #          112


## Save all these for future reference ===
save(markers.amy.t.pw, markers.amy.wilcox.block, medianNon0.amy, #markers.amy.binom.block,
     file="rdas/revision/markers-stats_Amyg-n5_findMarkers-SN-LEVEL_MNT2021.rda")
        
# As needed:
#load("rdas/revision/markers-stats_Amyg-n5_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.amy.t.pw, markers.amy.wilcox.block, medianNon0.amy

# Print these to pngs
markerList.t.pw <- lapply(markers.amy.t.pw, function(x){
  rownames(x)[x$FDR < 0.05 & x$non0median == TRUE]
  }
)
genes.top40.t <- lapply(markerList.t.pw, function(x){head(x, n=40)})

#dir.create("pdfs/revision/Amyg/")
smaller.set <- names(genes.top40.t)[lengths(genes.top40.t) <= 20]
left.set <- setdiff(names(genes.top40.t), smaller.set)

# Smaller graphical window
for(i in smaller.set){
  png(paste0("pdfs/revision/Amyg/Amyg_t_pairwise_top40markers-", i, "_logExprs_MNT2021.png"), height=950, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         features = genes.top40.t[[i]], 
                         features_name = i,
                         anno_name = "cellType",
                         ncol=5, point_alpha=0.4) +
      scale_color_manual(values = cell_colors.amy) +  
      ggtitle(label=paste0("Amyg ", i, " top markers: single-nucleus-level p.w. t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}

# 20-40 markers
for(i in left.set){
  png(paste0("pdfs/revision/Amyg/Amyg_t_pairwise_top40markers-", i, "_logExprs_MNT2021.png"), height=1900, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         features = genes.top40.t[[i]], 
                         features_name = i,
                         anno_name = "cellType",
                         ncol=5, point_alpha=0.4) +
      scale_color_manual(values = cell_colors.amy) +  
      ggtitle(label=paste0("Amyg ", i, " top markers: single-nucleus-level p.w. t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}



### Cluster-vs-all single-nucleus-level iteration ========================
## Load SCE with new info
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_Amyg-n5_cleaned-combined_SCE_MNT2021.rda",
     verbose=T)
    # sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo, annotationTab.amy, cell_colors.amy


# First drop "drop.lowNTx_" (1138 nuclei)
sce.amy <- sce.amy[ ,-grep("drop.",sce.amy$cellType)]
sce.amy$cellType <- droplevels(sce.amy$cellType)

# Remove 0 genes across all nuclei
sce.amy <- sce.amy[!rowSums(assay(sce.amy, "counts"))==0, ]  # keeps same 29371 genes


## Re-create 'logcounts' (don't want to use 'multiBatchNorm's down-scaling across donor 'batches')
# First 'hold' the MBN 'logcounts' for printing
sce.hold <- sce.amy

assay(sce.amy, "logcounts") <- NULL
sizeFactors(sce.amy) <- NULL
sce.amy <- logNormCounts(sce.amy)

## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.amy), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

markers.amy.t.1vAll <- list()
for(i in levels(sce.amy$cellType)){
  # Make temporary contrast
  sce.amy$contrast <- ifelse(sce.amy$cellType==i, 1, 0)
  # Test cluster vs. all others
  markers.amy.t.1vAll[[i]] <- findMarkers(sce.amy, groups=sce.amy$contrast,
                                          assay.type="logcounts", design=mod, test="t",
                                          std.lfc=TRUE,
                                          direction="up", pval.type="all", full.stats=T)
}
    ## Since all other stats are the same, and don't really use the non-standardized
    #    logFC, just generate one object, unlike before

class(markers.amy.t.1vAll[["Oligo"]])
    # a SimpleList of length 2, named "0" and "1" (from the temporary 'contrast')
    # -> we want the second entry, named "1"
    #    (for other purposes, might be interesting to look into that "0" entry, which
    #     is basically what genes are depleted in the cell type of interest)


sapply(markers.amy.t.1vAll, function(x){
  table(x[["1"]]$stats.0$log.FDR < log(.001))
})
    #       Astro_A Astro_B  Endo Excit_A Excit_B Excit_C Inhib_A Inhib_B Inhib_C
    # FALSE   23522   28810 28028   21690   27214   26898   23342   23325   24620
    # TRUE     5849     561  1343    7681    2157    2473    6029    6046    4751
    
    #       Inhib_D Inhib_E Inhib_F Inhib_G Inhib_H Micro Mural Oligo   OPC Tcell
    # FALSE   23992   25883   25815   28348   28189 24995 27850 26614 25021 27921
    # TRUE     5379    3488    3556    1023    1182  4376  1521  2757  4350  1450

# Do some reorganizing
markers.amy.t.1vAll <- lapply(markers.amy.t.1vAll, function(x){
  # Basically take the 'stats.[1 or 0]' since is redundant with the 'summary'-level stats
  lapply(x, function(y){ y[ ,4] }) 
})

# Re-name std.lfc column and the entries; add non-0-median info
for(i in names(markers.amy.t.1vAll)){
  colnames(markers.amy.t.1vAll[[i]][["0"]])[1] <- "std.logFC"
  colnames(markers.amy.t.1vAll[[i]][["1"]])[1] <- "std.logFC"
  # Add non0median Boolean - might be informative for both sets of stats
  markers.amy.t.1vAll[[i]][["0"]] <- cbind(markers.amy.t.1vAll[[i]][["0"]],
                                           medianNon0.amy[[i]][match(rownames(markers.amy.t.1vAll[[i]][["0"]]),
                                                                     names(medianNon0.amy[[i]]))])
  colnames(markers.amy.t.1vAll[[i]][["0"]])[4] <- "non0median"
  
  # "1" aka 'enriched'
  markers.amy.t.1vAll[[i]][["1"]] <- cbind(markers.amy.t.1vAll[[i]][["1"]],
                                           medianNon0.amy[[i]][match(rownames(markers.amy.t.1vAll[[i]][["1"]]),
                                                                     names(medianNon0.amy[[i]]))])
  colnames(markers.amy.t.1vAll[[i]][["1"]])[4] <- "non0median"
  
  # Then re-name the entries to more interpretable, because we'll keeping both contrasts
  names(markers.amy.t.1vAll[[i]]) <- paste0(i,c("_depleted", "_enriched"))
}


## Let's save this along with the previous pairwise results
save(markers.amy.t.pw, markers.amy.wilcox.block, markers.amy.t.1vAll, medianNon0.amy,
     file="rdas/revision/markers-stats_Amyg-n5_findMarkers-SN-LEVEL_MNT2021.rda")


## Print these to pngs
markerList.t.1vAll <- lapply(markers.amy.t.1vAll, function(x){
  rownames(x[[2]])[ x[[2]]$log.FDR < log(0.05) & x[[2]]$non0median==TRUE ]
}
)
genes.top40.t <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0("pdfs/revision/Amyg/Amyg_t_1vALL_top40markers-",i,"_logExprs_MNT2021.png"), height=1900, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         features = genes.top40.t[[i]], 
                         features_name = i,
                         anno_name = "cellType",
                         ncol=5, point_alpha=0.4) +
      scale_color_manual(values = cell_colors.amy) +  
      ggtitle(label=paste0("Amyg ", i, " top markers: 'cluster-vs-all-others' t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}


## How do they intersect?
markerList.t.pw <- lapply(markers.amy.t.pw, function(x){
  rownames(x)[ x$FDR < 0.05 & x$non0median==TRUE ]
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
        # Astro_A Astro_B    Endo Excit_A Excit_B Excit_C Inhib_A Inhib_B Inhib_C Inhib_D 
        #      26       3      34       7      33      28       9      28      12      27 
        # Inhib_E Inhib_F Inhib_G Inhib_H   Micro   Mural   Oligo     OPC   Tcell 
        #      25      20       4      24      24      32      30      28      33


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

write.csv(top40genes, file="tables/revision/top40genesLists_Amyg-n5_cellType_SN-LEVEL-tests_MNT2021.csv",
          row.names=FALSE)





### For fig: Plot some top markers in vlnplot array (12Jun2020) === === === ===

load("rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
load("rdas/markers-stats_Amyg-n2_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.amy.t.1vAll, markers.amy.t.pw, markers.amy.wilcox.block
    # Focus on the pairwise result (".pw") bc more specific
    rm(markers.amy.t.1vAll, markers.amy.wilcox.block)

# First drop "ambig.lowNtrxts" (50 nuclei)
sce.amy <- sce.amy[ ,sce.amy$cellType != "Ambig.lowNtrxts"]
sce.amy$cellType <- droplevels(sce.amy$cellType)


# Take top two for broad glia
topToPrint <- as.data.frame(sapply(markers.amy.t.pw, function(x) {head(rownames(x),n=2)}))
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
                 x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=1,
                 add_legend=F) + stat_summary(fun = median, fun.min = median, fun.max = median,
                                              geom = "crossbar", width = 0.3,
                                              colour=rep(tableau20[1:12], length(topToPrint))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12), plot.title = element_text(size = 25)) +
    ggtitle(label="AMY top marker array (with glial)")
)
dev.off()

# Neuronal markers only (highlighted in paper)
topToPrint <- c("NRN1", "NPTX1", "SLC30A3", # Excit_A
  "SLC17A6", "VCAN", #Excit_B
  "MCHR2", "CDH22", # Excit_C
  "PENK", "ACTN2", # Inhib_A (& Inhib_H)
  "CCK", "CALB2", "KIT", # Inhib_B/D
  "CRH", # Inhib_B
  "NPFFR2", "TLL1", # Inhib_C
  "SYT2", "ONECUT2", # Inhib_E
  "LHX6", "ELAVL2", # Inhib_F
  "DIO2", # Inhib_G
  "NPY", "PRLR") # Inhib_H

#pdf("pdfs/revision/pubFigures/regionSpecific_Amyg-n5_topMarkers-ARRAY_MNT2021_1stHalf.pdf", height=11, width=3.8)
pdf("pdfs/revision/pubFigures/regionSpecific_Amyg-n5_topMarkers-ARRAY_MNT2021_2ndHalf.pdf", height=11, width=3.8)
print(
#  plotExpressionCustom(sce.amy, features=topToPrint[1:11], features_name="manually-selected",
  plotExpressionCustom(sce.amy, features=topToPrint[12:22], features_name="manually-selected",
                       point_alpha=0.5, point_size=.7, ncol=1,
                       scales="free_y") +
    scale_color_manual(values = cell_colors.amy) +
    xlab("") +
    ylab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(angle = 90, size = 14),
          plot.title = element_text(size = 12),
          panel.grid.major=element_line(colour="grey95", size=0.8),
          panel.grid.minor=element_line(colour="grey95", size=0.4)) +
    ggtitle(label="AMY top neuronal markers array")
)
dev.off()





# ## Aside: add in 't.stat' as in 'step04' analyses to save for LoHu/LeCo ===
# for(s in names(markers.amy.t.1vAll)){
#   markers.amy.t.1vAll[[s]]$t.stat <- markers.amy.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.amy))
# }
# 
# save(markers.amy.t.1vAll, markers.amy.t.pw, sce.amy,
#      file="rdas/markerStats-and-SCE_AMY-n2_sn-level_cleaned_MNTNov2020.rda")


### Session info for 02Jun2021 ============
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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices datasets  utils     methods  
# [9] base     
# 
# other attached packages:
#   [1] ggrepel_0.9.1               dynamicTreeCut_1.63-1       dendextend_1.14.0          
# [4] jaffelab_0.99.30            rafalib_1.0.0               DropletUtils_1.10.3        
# [7] batchelor_1.6.2             scran_1.18.5                scater_1.18.6              
# [10] ggplot2_3.3.3               EnsDb.Hsapiens.v86_2.99.0   ensembldb_2.14.1           
# [13] AnnotationFilter_1.14.0     GenomicFeatures_1.42.3      AnnotationDbi_1.52.0       
# [16] SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0 Biobase_2.50.0             
# [19] GenomicRanges_1.42.0        GenomeInfoDb_1.26.7         IRanges_2.24.1             
# [22] S4Vectors_0.28.1            BiocGenerics_0.36.1         MatrixGenerics_1.2.1       
# [25] matrixStats_0.58.0         
# 
# loaded via a namespace (and not attached):
#   [1] googledrive_1.0.1         ggbeeswarm_0.6.0          colorspace_2.0-0         
# [4] ellipsis_0.3.2            scuttle_1.0.4             bluster_1.0.0            
# [7] XVector_0.30.0            BiocNeighbors_1.8.2       rstudioapi_0.13          
# [10] farver_2.1.0              bit64_4.0.5               fansi_0.4.2              
# [13] xml2_1.3.2                splines_4.0.4             R.methodsS3_1.8.1        
# [16] sparseMatrixStats_1.2.1   cachem_1.0.4              Rsamtools_2.6.0          
# [19] ResidualMatrix_1.0.0      dbplyr_2.1.1              R.oo_1.24.0              
# [22] HDF5Array_1.18.1          compiler_4.0.4            httr_1.4.2               
# [25] dqrng_0.2.1               assertthat_0.2.1          Matrix_1.3-2             
# [28] fastmap_1.1.0             lazyeval_0.2.2            limma_3.46.0             
# [31] BiocSingular_1.6.0        prettyunits_1.1.1         tools_4.0.4              
# [34] rsvd_1.0.3                igraph_1.2.6              gtable_0.3.0             
# [37] glue_1.4.2                GenomeInfoDbData_1.2.4    dplyr_1.0.5              
# [40] rappdirs_0.3.3            Rcpp_1.0.6                vctrs_0.3.6              
# [43] Biostrings_2.58.0         rhdf5filters_1.2.0        rtracklayer_1.50.0       
# [46] DelayedMatrixStats_1.12.3 stringr_1.4.0             beachmat_2.6.4           
# [49] lifecycle_1.0.0           irlba_2.3.3               statmod_1.4.35           
# [52] XML_3.99-0.6              edgeR_3.32.1              zlibbioc_1.36.0          
# [55] scales_1.1.1              hms_1.0.0                 ProtGenerics_1.22.0      
# [58] rhdf5_2.34.0              RColorBrewer_1.1-2        curl_4.3                 
# [61] memoise_2.0.0             gridExtra_2.3             segmented_1.3-3          
# [64] biomaRt_2.46.3            stringi_1.5.3             RSQLite_2.2.7            
# [67] BiocParallel_1.24.1       rlang_0.4.10              pkgconfig_2.0.3          
# [70] bitops_1.0-7              lattice_0.20-41           purrr_0.3.4              
# [73] Rhdf5lib_1.12.1           labeling_0.4.2            GenomicAlignments_1.26.0 
# [76] cowplot_1.1.1             bit_4.0.4                 tidyselect_1.1.1         
# [79] magrittr_2.0.1            R6_2.5.0                  generics_0.1.0           
# [82] DelayedArray_0.16.3       DBI_1.1.1                 pillar_1.6.0             
# [85] withr_2.4.2               RCurl_1.98-1.3            tibble_3.1.1             
# [88] crayon_1.4.1              utf8_1.2.1                BiocFileCache_1.14.0     
# [91] viridis_0.6.0             progress_1.2.2            locfit_1.5-9.4           
# [94] grid_4.0.4                blob_1.2.1                digest_0.6.27            
# [97] R.utils_2.10.1            openssl_1.4.3             munsell_0.5.0            
# [100] beeswarm_0.3.1            viridisLite_0.4.0         vipor_0.4.5              
# [103] askpass_1.1

