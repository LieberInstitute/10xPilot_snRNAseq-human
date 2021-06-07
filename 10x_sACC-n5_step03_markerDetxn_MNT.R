### MNT 10x snRNA-seq workflow: step 03 - marker detection
###   **Region-specific analyses**
###     - (2x) sACC samples from: Br5161 & Br5212
### Initiated MNT 12Feb2020
### MNT 24May2021: add expansion samples (n=3, incl'g 2 female)
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
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_sACC-n5_cleaned-combined_SCE_MNT2021.rda",
     verbose=T)
    ## sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc, ref.sampleInfo

table(sce.sacc$cellType)
    # Astro_A        Astro_B   drop.doublet    drop.lowNTx        Excit_A        Excit_B
    #     747            160             28            298            856            575
    # Excit_C        Excit_D        Excit_E        Excit_F        Excit_G        Inhib_A
    #    1735            311            428            228             30            842
    # Inhib_B        Inhib_C        Inhib_D        Inhib_E        Inhib_F        Inhib_G
    #     912            465            384            330            521            206
    # Inhib_H        Inhib_I        Inhib_J        Inhib_K          Micro Neu_FAT2.CDH15
    #     208             39             42             25            784             20
    # Oligo_A        Oligo_B            OPC
    #    4389            195            911

# First drop "drop." lowNTx & doublet clusters (298; 28, respectively)
sce.sacc <- sce.sacc[ ,-grep("drop.", sce.sacc$cellType)]
sce.sacc$cellType <- droplevels(sce.sacc$cellType)

# Remove 0 genes across all nuclei
sce.sacc <- sce.sacc[!rowSums(assay(sce.sacc, "counts"))==0, ]  # keeps same 29583 genes

## Re-create 'logcounts' (don't want to use 'multiBatchNorm's down-scaling across donor 'batches')
 # First 'hold' the MBN 'logcounts' for printing
sce.hold <- sce.sacc

assay(sce.sacc, "logcounts") <- NULL
sizeFactors(sce.sacc) <- NULL
sce.sacc <- logNormCounts(sce.sacc)


### First make a list of Boolean param / cell subtype ===
# Will use this to assess more 'valid', non-noise-driving markers
cellSubtype.idx <- splitit(sce.sacc$cellType)
medianNon0.sacc <- lapply(cellSubtype.idx, function(x){
  apply(as.matrix(assay(sce.sacc, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

sapply(medianNon0.sacc, table)



## Traditional t-test implementation ===
mod <- with(colData(sce.sacc), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`


# Run pairwise t-tests
markers.sacc.t.pw <- findMarkers(sce.sacc, groups=sce.sacc$cellType,
                                    assay.type="logcounts", design=mod, test="t",
                                    direction="up", pval.type="all", full.stats=T)

sapply(markers.sacc.t.pw, function(x){table(x$FDR<0.05)})
    #       Astro Excit.1 Excit.2 Excit.3 Excit.4 Inhib.1 Inhib.2 Micro Oligo   OPC
    # FALSE 27821   28246   28493   28378   27967   28455   28282 27319 28059 28272
    # TRUE    953     528     281     396     807     319     492  1455   715   502


## WMW: Blocking on donor (this test doesn't take 'design=' argument) ===
markers.sacc.wilcox.block <- findMarkers(sce.sacc, groups=sce.sacc$cellType,
                                        assay.type="logcounts", block=sce.sacc$donor, test="wilcox",
                                        direction="up", pval.type="all", full.stats=T)

sapply(markers.sacc.wilcox.block, function(x){table(x$FDR<0.05)["TRUE"]})
    # Astro_A.TRUE        Astro_B.NA      Excit_A.TRUE      Excit_B.TRUE      Excit_C.TRUE
    #           76                NA                 7                18                27
    # Excit_D.TRUE      Excit_E.TRUE      Excit_F.TRUE      Excit_G.TRUE        Inhib_A.NA
    #           11                 2                22                 5                NA
    # Inhib_B.TRUE      Inhib_C.TRUE      Inhib_D.TRUE        Inhib_E.NA      Inhib_F.TRUE
    #           10                 3                 8                NA                 9
    # Inhib_G.TRUE      Inhib_H.TRUE        Inhib_I.NA        Inhib_J.NA      Inhib_K.TRUE
    #            6                33                NA                NA                 1
    #   Micro.TRUE Neu_FAT2.CDH15.NA      Oligo_A.TRUE        Oligo_B.NA          OPC.TRUE
    #           86                NA                46                NA                39


## Binomial ===
markers.sacc.binom.block <- findMarkers(sce.sacc, groups=sce.sacc$cellType,
                                       assay.type="logcounts", block=sce.sacc$donor, test="binom",
                                       direction="up", pval.type="all", full.stats=T)

sapply(markers.sacc.binom.block, function(x){table(x$FDR<0.05)["TRUE"]})
    # no results for any subpops



# Add respective 'non0median' column to the stats for each set of markers
for(i in names(markers.sacc.t.pw)){
  markers.sacc.t.pw[[i]] <- cbind(markers.sacc.t.pw[[i]],
                                 medianNon0.sacc[[i]][match(rownames(markers.sacc.t.pw[[i]]),
                                                           names(medianNon0.sacc[[i]]))])
  colnames(markers.sacc.t.pw[[i]])[28] <- "non0median"
}

sapply(markers.sacc.t.pw, function(x){table(x$FDR<0.05 & x$non0median == TRUE)["TRUE"]})
    #       Astro_A Astro_B Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Excit_G Inhib_A
    # FALSE   29348   29541   29567   29554   29537   29558   29551   29522   29520   29573
    # TRUE      235      42      16      29      46      25      32      61      63      10

    #       Inhib_B Inhib_C Inhib_D Inhib_E Inhib_F Inhib_G Inhib_H Inhib_I Inhib_J Inhib_K
    # FALSE   29547   29565   29562   29575   29569   29562   29515   29563   29555   29558
    # TRUE       36      18      21       8      14      21      68      20      28      25

    #       Micro Neu_FAT2.CDH15 Oligo_A Oligo_B   OPC
    # FALSE 29276          29471   29312   29579 29450
    # TRUE    307            112     271       4   133


## Save all these for future reference
save(markers.sacc.t.pw, markers.sacc.wilcox.block, medianNon0.sacc, #markers.sacc.binom.block,
     file="rdas/revision/markers-stats_sACC-n5_findMarkers-SN-LEVEL_MNT2021.rda")


# # As needed:
# load("rdas/revision/markers-stats_sACC-n5_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
#     # markers.sacc.t.pw, markers.sacc.wilcox.block, medianNon0.sacc

# Print these to pngs
markerList.t.pw <- lapply(markers.sacc.t.pw, function(x){
  rownames(x)[x$FDR < 0.05 & x$non0median == TRUE]
  }
)

genes.top40.t <- lapply(markerList.t.pw, function(x){head(x, n=40)})


#dir.create("pdfs/revision/sACC/")
smaller.set <- names(genes.top40.t)[lengths(genes.top40.t) <= 20]
left.set <- setdiff(names(genes.top40.t), smaller.set)

# Smaller graphical window
for(i in smaller.set){
  png(paste0("pdfs/revision/sACC/sACC_t_pairwise_top40markers-", i, "_logExprs_MNT2021.png"), height=950, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         features = genes.top40.t[[i]], 
                         features_name = i,
                         anno_name = "cellType",
                         ncol=5, point_alpha=0.4) +
      scale_color_manual(values = cell_colors.sacc) +  
      ggtitle(label=paste0("sACC ", i, " top markers: single-nucleus-level p.w. t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}

# 20-40 markers
for(i in left.set){
  png(paste0("pdfs/revision/sACC/sACC_t_pairwise_top40markers-", i, "_logExprs_MNT2021.png"), height=1900, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         features = genes.top40.t[[i]], 
                         features_name = i,
                         anno_name = "cellType",
                         ncol=5, point_alpha=0.4) +
      scale_color_manual(values = cell_colors.sacc) +  
      ggtitle(label=paste0("sACC ", i, " top markers: single-nucleus-level p.w. t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}





### Cluster-vs-all single-nucleus-level iteration ======

## Load SCE with new info
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_sACC-n5_cleaned-combined_SCE_MNT2021.rda",
     verbose=T)
    ## sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc, ref.sampleInfo

table(sce.sacc$cellType)

# First drop "drop." lowNTx & doublet clusters (298; 28, respectively)
sce.sacc <- sce.sacc[ ,-grep("drop.", sce.sacc$cellType)]
sce.sacc$cellType <- droplevels(sce.sacc$cellType)

# Remove 0 genes across all nuclei
sce.sacc <- sce.sacc[!rowSums(assay(sce.sacc, "counts"))==0, ]  # keeps same 29583 genes

## Re-create 'logcounts' (don't want to use 'multiBatchNorm's down-scaling across donor 'batches')
# First 'hold' the MBN 'logcounts' for printing
sce.hold <- sce.sacc

assay(sce.sacc, "logcounts") <- NULL
sizeFactors(sce.sacc) <- NULL
sce.sacc <- logNormCounts(sce.sacc)


# Traditional t-test ===
mod <- with(colData(sce.sacc), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

markers.sacc.t.1vAll <- list()
for(i in levels(sce.sacc$cellType)){
  # Make temporary contrast
  sce.sacc$contrast <- ifelse(sce.sacc$cellType==i, 1, 0)
  # Test cluster vs. all others
  markers.sacc.t.1vAll[[i]] <- findMarkers(sce.sacc, groups=sce.sacc$contrast,
                                          assay.type="logcounts", design=mod, test="t",
                                          std.lfc=TRUE,
                                          direction="up", pval.type="all", full.stats=T)
}
    ## Since all other stats are the same, and don't really use the non-standardized
    #    logFC, just generate one object, unlike before

class(markers.sacc.t.1vAll[["Micro"]])
    # a SimpleList of length 2, named "0" and "1" (from the temporary 'contrast')
    # -> we want the second entry, named "1"
    #    (for other purposes, might be interesting to look into that "0" entry, which
    #     is basically what genes are depleted in the cell type of interest)


sapply(markers.sacc.t.1vAll, function(x){
  table(x[["1"]]$stats.0$log.FDR < log(.001))
})


# Do some reorganizing
markers.sacc.t.1vAll <- lapply(markers.sacc.t.1vAll, function(x){
  # Basically take the 'stats.[1 or 0]' since is redundant with the 'summary'-level stats
  lapply(x, function(y){ y[ ,4] }) 
})

# Re-name std.lfc column and the entries; add non-0-median info
for(i in names(markers.sacc.t.1vAll)){
  colnames(markers.sacc.t.1vAll[[i]][["0"]])[1] <- "std.logFC"
  colnames(markers.sacc.t.1vAll[[i]][["1"]])[1] <- "std.logFC"
  # Add non0median Boolean - might be informative for both sets of stats
  markers.sacc.t.1vAll[[i]][["0"]] <- cbind(markers.sacc.t.1vAll[[i]][["0"]],
                                           medianNon0.sacc[[i]][match(rownames(markers.sacc.t.1vAll[[i]][["0"]]),
                                                                     names(medianNon0.sacc[[i]]))])
  colnames(markers.sacc.t.1vAll[[i]][["0"]])[4] <- "non0median"
  
  # "1" aka 'enriched'
  markers.sacc.t.1vAll[[i]][["1"]] <- cbind(markers.sacc.t.1vAll[[i]][["1"]],
                                           medianNon0.sacc[[i]][match(rownames(markers.sacc.t.1vAll[[i]][["1"]]),
                                                                     names(medianNon0.sacc[[i]]))])
  colnames(markers.sacc.t.1vAll[[i]][["1"]])[4] <- "non0median"
  
  # Then re-name the entries to more interpretable, because we'll keeping both contrasts
  names(markers.sacc.t.1vAll[[i]]) <- paste0(i,c("_depleted", "_enriched"))
}


## Let's save this along with the previous pairwise results
save(markers.sacc.t.pw, markers.sacc.wilcox.block, markers.sacc.t.1vAll, medianNon0.sacc,
     file="rdas/revision/markers-stats_sACC-n5_findMarkers-SN-LEVEL_MNT2021.rda")


## Print these to pngs
markerList.t.1vAll <- lapply(markers.sacc.t.1vAll, function(x){
  rownames(x[[2]])[ x[[2]]$log.FDR < log(0.05) & x[[2]]$non0median==TRUE ]
 }
)
genes.top40.t <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0("pdfs/revision/sACC/sACC_t_1vALL_top40markers-",i,"_logExprs_MNT2021.png"), height=1900, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         features = genes.top40.t[[i]], 
                         features_name = i,
                         anno_name = "cellType",
                         ncol=5, point_alpha=0.4) +
      scale_color_manual(values = cell_colors.sacc) +  
      ggtitle(label=paste0("sACC ", i, " top markers: 'cluster-vs-all-others' t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}


## How do they intersect?
markerList.t.pw <- lapply(markers.sacc.t.pw, function(x){
  rownames(x)[ x$FDR < 0.05 & x$non0median==TRUE ]
 }
)

# From pairwise t-tests, FDR < 0.05
lengths(markerList.t.pw)
    #   Astro_A        Astro_B        Excit_A        Excit_B        Excit_C 
    #       235             42             16             29             46 
    #   Excit_D        Excit_E        Excit_F        Excit_G        Inhib_A 
    #        25             32             61             63             10 
    #   Inhib_B        Inhib_C        Inhib_D        Inhib_E        Inhib_F 
    #        36             18             21              8             14 
    #   Inhib_G        Inhib_H        Inhib_I        Inhib_J        Inhib_K 
    #        21             68             20             28             25 
    #     Micro Neu_FAT2.CDH15        Oligo_A        Oligo_B            OPC 
    #       307            112            271              4            133

# From cluster-vs-all others, FDR < 0.05
lengths(markerList.t.1vAll)
    #   Astro_A        Astro_B        Excit_A        Excit_B        Excit_C 
    #      1050            373           4523           4226           4613 
    #   Excit_D        Excit_E        Excit_F        Excit_G        Inhib_A 
    #      4087           3969           3835            966           2507 
    #   Inhib_B        Inhib_C        Inhib_D        Inhib_E        Inhib_F 
    #      3707           2268           2590           1490           1784 
    #   Inhib_G        Inhib_H        Inhib_I        Inhib_J        Inhib_K 
    #      2035           2476            929            612            645 
    #     Micro Neu_FAT2.CDH15        Oligo_A        Oligo_B            OPC 
    #       719            472            797            152           1214 

# Intersection
sapply(names(markerList.t.pw), function(c){
  length(intersect(markerList.t.pw[[c]],
                   markerList.t.1vAll[[c]]))
})
    #   Astro_A        Astro_B        Excit_A        Excit_B        Excit_C 
    #       235             42             16             29             46 
    #   Excit_D        Excit_E        Excit_F        Excit_G        Inhib_A 
    #        25             32             61             63             10 
    #   Inhib_B        Inhib_C        Inhib_D        Inhib_E        Inhib_F 
    #        36             18             21              8             14 
    #   Inhib_G        Inhib_H        Inhib_I        Inhib_J        Inhib_K 
    #        21             68             20             28             25 
    #     Micro Neu_FAT2.CDH15        Oligo_A        Oligo_B            OPC 
    #       307            112            271              4            133

    # Of top 40's:
    sapply(names(markerList.t.pw), function(c){
      length(intersect(lapply(markerList.t.pw, function(l){head(l,n=40)})[[c]],
                       lapply(markerList.t.1vAll, function(l){head(l,n=40)})[[c]]
      ))
    })
        #   Astro_A        Astro_B        Excit_A        Excit_B        Excit_C 
        #        28             26             11             18             24 
        #   Excit_D        Excit_E        Excit_F        Excit_G        Inhib_A 
        #        10             16             29             25              9 
        #   Inhib_B        Inhib_C        Inhib_D        Inhib_E        Inhib_F 
        #        17             14             14              8             11 
        #   Inhib_G        Inhib_H        Inhib_I        Inhib_J        Inhib_K 
        #        13             27             16             23             15 
        #     Micro Neu_FAT2.CDH15        Oligo_A        Oligo_B            OPC 
        #        36             31             25              4             27


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

write.csv(top40genes, file="tables/revision/top40genesLists_sACC-n5_cellType_SN-LEVEL-tests_MNT2021.csv",
          row.names=FALSE)




# ## Aside: add in 't.stat' as in 'step04' analyses to save for LoHu/LeCo ===
# for(s in names(markers.sacc.t.1vAll)){
#   markers.sacc.t.1vAll[[s]]$t.stat <- markers.sacc.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.sacc))
# }
# 
# save(markers.sacc.t.1vAll, markers.sacc.t.pw, sce.sacc,
#      file="rdas/markerStats-and-SCE_sACC-n2_sn-level_cleaned_MNTNov2020.rda")


### Session info for 07Jun2021 ============
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
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
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

