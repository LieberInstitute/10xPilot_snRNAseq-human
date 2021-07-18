### MNT 10x snRNA-seq workflow: step 04 - downstream comparisons
###   **Pan-brain analyses**
###     - n=24 samples from 5 regions
###   * Cross-region analysis/correlation and comp. to other datasets
#####################################################################

library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
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

load("rdas/revision/markers-stats_DLPFC-n3_findMarkers-SN-LEVEL_MNT_v2_2021.rda", verbose=T)
    # markers.dlpfc.t.1vAll, medianNon0.dlpfc
    rm(medianNon0.dlpfc)

load("rdas/revision/markers-stats_HPC-n3_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.hpc.t.pw, markers.hpc.t.1vAll, medianNon0.hpc
    rm(markers.hpc.t.pw, medianNon0.hpc)

load("rdas/revision/markers-stats_NAc-n8_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.nac.t.pw, markers.nac.t.1vAll, medianNon0.nac
    rm(markers.nac.t.pw, medianNon0.nac)

load("rdas/revision/markers-stats_Amyg-n5_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.amy.t.pw, markers.amy.wilcox.block, markers.amy.t.1vAll, medianNon0.amy
    rm(markers.amy.t.pw, markers.amy.wilcox.block, medianNon0.amy)

load("rdas/revision/markers-stats_sACC-n5_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.sacc.t.pw, markers.sacc.wilcox.block, markers.sacc.t.1vAll, medianNon0.sacc
    rm(markers.sacc.t.pw, markers.sacc.wilcox.block, medianNon0.sacc)


# Re-order rows in each list entry for each set of stats
expressedGenes.list <- list(amy=rownames(markers.amy.t.1vAll[["Astro_A"]][[2]]),
                            dlpfc=rownames(markers.dlpfc.t.1vAll[["Astro"]][[2]]),
                            hpc=rownames(markers.hpc.t.1vAll[["Astro_A"]][[2]]),
                            nac=rownames(markers.nac.t.1vAll[["Astro_A"]][[2]]),
                            sacc=rownames(markers.sacc.t.1vAll[["Astro_A"]][[2]])
                            )

expressedGenes <- unique(unlist(expressedGenes.list))

expressedGenes <- expressedGenes[expressedGenes %in% expressedGenes.list[["amy"]] &
                                   expressedGenes %in% expressedGenes.list[["dlpfc"]] &
                                   expressedGenes %in% expressedGenes.list[["nac"]] &
                                   expressedGenes %in% expressedGenes.list[["hpc"]] &
                                   expressedGenes %in% expressedGenes.list[["sacc"]]
                                 ]
length(expressedGenes)  # 27875

# Store each set of stats into a list
FMstats.list <- list(amy=lapply(markers.amy.t.1vAll,function(x){x[[2]]}),
                     dlpfc=lapply(markers.dlpfc.t.1vAll,function(x){x[[2]]}),
                     hpc=lapply(markers.hpc.t.1vAll,function(x){x[[2]]}),
                     nac=lapply(markers.nac.t.1vAll,function(x){x[[2]]}),
                     sacc=lapply(markers.sacc.t.1vAll,function(x){x[[2]]}))

    # How many subclusters in each region?
    sapply(FMstats.list, length)
        #  amy dlpfc   hpc   nac  sacc
        #   19    19    20    24    25
    
    sapply(FMstats.list, function(x){nrow(x[[1]])})
        #  amy dlpfc   hpc   nac  sacc 
        #29371 29310 28764 29680 29583

# Subset and re-order for those intersecting genes across all regions
for(x in names(FMstats.list)){
  for(s in names(FMstats.list[[x]])){
    FMstats.list[[x]][[s]] <- FMstats.list[[x]][[s]][expressedGenes, ]
  }
}


sapply(FMstats.list, function(x){nrow(x[[1]])})
    # good
sapply(FMstats.list, function(x){head(x[[1]], n=3)})


### Get n Nuclei numbers for each region so can compute t-statistics ===
  # This can be done with Cohen's D (the 'std.lfc'), as d = t/sqrt(N)

load("rdas/revision/all-n24-samples_across-regions-analyses_forFigOnly_MNT2021.rda", verbose=T)
    # sce.allRegions, chosen.hvgs.union, ref.sampleInfo, Readme

sce.allRegions
    # class: SingleCellExperiment
    # dim: 33538 70497
table(sce.allRegions$region)
    #  amy dlpfc   hpc   nac  sacc 
    #14039 11202 10139 19892 15343

table(sce.allRegions$cellType)

sampleNumNuclei <- table(sce.allRegions$region)

## Calculate and add t-statistic (= std.logFC * sqrt(N))
for(x in names(FMstats.list)){
  for(s in names(FMstats.list[[x]])){
    FMstats.list[[x]][[s]]$t.stat <- FMstats.list[[x]][[s]]$std.logFC * sqrt(sampleNumNuclei[x])
  }
}


## Let's save these
readme.mnt <- "These stats are from region-specific specificity modeling (cluster-vs-all-others) at the single-nucleus level with 'scran::findMarkers()'. The t-statistic is computed by sqrt(N.nuclei) * std.logFC."
save(FMstats.list, sampleNumNuclei, readme.mnt, ref.sampleInfo,
     file="rdas/revision/markers-stats_all-regions-combined_SN-LEVEL-1vAll_MNT2021.rda")


# (If needed)
load("rdas/revision/markers-stats_all-regions-combined_SN-LEVEL-1vAll_MNT2021.rda", verbose=T)
    # FMstats.list, sampleNumNuclei, readme.mnt, ref.sampleInfo
readme.mnt

## Create matrix of t's with region:subcluster identifiers
ts.list <- lapply(FMstats.list, function(x){
  sapply(x, function(y){y$t.stat})
  }
)
# Add back in region suffix
for(i in names(ts.list)){
  colnames(ts.list[[i]]) <- paste0(colnames(ts.list[[i]]), "_", i)
}
# Cbind
ts.fullMat <- do.call(cbind, ts.list)


## Correlation; first shorten names
colnames(ts.fullMat) <- gsub("Excit", "Ex", colnames(ts.fullMat))
colnames(ts.fullMat) <- gsub("Inhib", "In", colnames(ts.fullMat))
colnames(ts.fullMat) <- gsub("Astro", "As", colnames(ts.fullMat))
colnames(ts.fullMat)[colnames(ts.fullMat)=="Neu_FAT2.CDH15_sacc"] <- "Neu_ambig_sacc"

# Perform in cluster-specific gene space, as with across-species comparisons
    clus_specific_indices = mapply(function(t) {
      oo = order(t, decreasing = TRUE)[1:100]
      },
    as.data.frame(ts.fullMat)
    )
    clus_ind = unique(as.numeric(clus_specific_indices))
    length(clus_ind)  # so of up to 10200 (100 x 102 cellType), 3715 unique
    
    ts.defined <- ts.fullMat[clus_ind, ]


cor_t_xRegions <- cor(ts.fullMat)
cor_t_defined <- cor(ts.defined)


### Heatmap - typically use levelplot (e.g. below), but will want pheatmap bc can cluster cols/rows
theSeq.all = seq(-1, 1, by = 0.025)
my.col.all <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq.all)-1)


pdf("pdfs/revision/acrossRegions_correlation_region-specific-subcluster-ts_MNT2021.pdf")
pheatmap(cor_t_xRegions,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=4.5, fontsize_col=4.5,
         main="Correlation of cluster-specific t's from all regions \n (all shared expressed genes)")
pheatmap(cor_t_defined,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=4.5, fontsize_col=4.5,
         main="Correlation of cluster-specific t's from all regions \n (top 100 cluster genes space)")
dev.off()


## Subset on neuronal subcluster t's and check
ts.fullMat.neu <- ts.fullMat
ts.defined.neu <- ts.defined
for(i in c("As", "Micro", "Endo", "Mural","Oligo", "OPC", "Tcell", "Macro")){
  ts.fullMat.neu <- ts.fullMat.neu[ ,-grep(i, colnames(ts.fullMat.neu))]
  ts.defined.neu <- ts.defined.neu[ ,-grep(i, colnames(ts.defined.neu))]
}

cor_t_xRegions.neu <- cor(ts.fullMat.neu)
cor_t_defined.neu <- cor(ts.defined.neu)


# Add some cluster info for add'l heatmap annotations
clusterInfo <- data.frame(region=ss(colnames(ts.fullMat.neu), "_",3))
rownames(clusterInfo) <- colnames(ts.fullMat.neu)

# Region cols to be consistent with the TSNE
clusterCols <- list(region=tableau10medium[1:5])
names(clusterCols[["region"]]) <- levels(as.factor(clusterInfo$region))

# Print
pdf("pdfs/revision/acrossRegions_correlation_region-specific-NeuronalSubcluster-ts_MNT2021.pdf",width=9, height=9)
# All genes
pheatmap(cor_t_xRegions.neu,
         annotation_col=clusterInfo,
         annotation_colors=clusterCols,
         #show_colnames=FALSE,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=6.2, fontsize_col=6.2,
         main="Correlation of neuronal cluster-specific t's from all regions \n (all shared expressed genes)")
# With numbers
pheatmap(cor_t_xRegions.neu,
         annotation_col=clusterInfo,
         annotation_colors=clusterCols,
         #show_colnames=FALSE,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=6.2, fontsize_col=6.2,
         display_numbers=TRUE, fontsize_number=2.6,
         main="Correlation of neuronal cluster-specific t's from all regions \n (all shared expressed genes)")

# Top 100 cluster genes space
pheatmap(cor_t_defined.neu,
         annotation_col=clusterInfo,
         annotation_colors=clusterCols,
         #show_colnames=FALSE,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=6.2, fontsize_col=6.2,
         main="Correlation of neuronal cluster-specific t's from all regions \n (top 100 cluster genes space, incl'g glial)")
# With numbers
pheatmap(cor_t_defined.neu,
         annotation_col=clusterInfo,
         annotation_colors=clusterCols,
         #show_colnames=FALSE,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=6.2, fontsize_col=6.2,
         display_numbers=TRUE, fontsize_number=2.6,
         main="Correlation of neuronal cluster-specific t's from all regions \n (top 100 cluster genes space, incl'g glial)")
dev.off()



## Non-neuronal set for supplement === === ===
ts.fullMat.non <- ts.fullMat
ts.defined.non <- ts.defined
glia.idx <- NA
for(i in c("As", "Micro", "Endo", "Mural","Oligo", "OPC", "Tcell", "Macro")){
  glia.idx <- c(glia.idx, grep(i, colnames(ts.fullMat.non)))
}
# Rm the empty NA
glia.idx <- glia.idx[-1]
ts.fullMat.non <- ts.fullMat.non[ ,glia.idx]
ts.defined.non <- ts.defined.non[ ,glia.idx]

cor_t_xRegions.non <- cor(ts.fullMat.non)
cor_t_defined.non <- cor(ts.defined.non)


# Add some cluster info for add'l heatmap annotations
clusterInfo.glia <- data.frame(region=ifelse(is.na(ss(colnames(ts.fullMat.non), "_",3)),
                                             ss(colnames(ts.fullMat.non), "_",2),
                                             ss(colnames(ts.fullMat.non), "_",3))
                               )

rownames(clusterInfo.glia) <- colnames(ts.fullMat.non)

# Print
pdf("pdfs/revision/acrossRegions_correlation_region-specific-NON-NeuronalSubcluster-ts_MNT2021.pdf",width=9)
pheatmap(cor_t_xRegions.non,
         annotation_col=clusterInfo.glia,
         annotation_colors=clusterCols,
         #show_colnames=FALSE,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=7, fontsize_col=7,
         display_numbers=TRUE, fontsize_number=4,
         main="Correlation of glia/other cluster-specific t's from all regions \n (all shared expressed genes)")
pheatmap(cor_t_defined.non,
         annotation_col=clusterInfo.glia,
         annotation_colors=clusterCols,
         #show_colnames=FALSE,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=7, fontsize_col=7,
         display_numbers=TRUE, fontsize_number=4,
         main="Correlation of glia/other cluster-specific t's from all regions \n (top 100 cluster genes space, incl'g neuronal)")
dev.off()


## For main section on AMY 'Astro_B' ===
table(droplevels(sce.allRegions$cellType[grep("Astro", sce.allRegions$cellType)])) 
    # amy_Astro_A  amy_Astro_B  dlpfc_Astro  hpc_Astro_A  hpc_Astro_B  nac_Astro_A 
    #        1555           83          782          936          234           99 
    # nac_Astro_B sacc_Astro_A sacc_Astro_B 
    #        1000          747          160

sce.astro <- sce.allRegions[ ,grep("Astro", sce.allRegions$cellType)]
sce.astro$cellType <- droplevels(sce.astro$cellType)
    #                br5161 br5207 br5212 br5276 br5287 br5400 br5701
    # amy_Astro_A     484      0    350    230      0    111    380
    # amy_Astro_B       7      0     10     49      0     12      5
    # dlpfc_Astro     371    274    137      0      0      0      0
    # hpc_Astro_A     424      0    375      0    137      0      0
    # hpc_Astro_B      83      0    125      0     26      0      0
    # nac_Astro_A      27      0      5      8      3     56      0
    # nac_Astro_B     115      0    377    173      8    294     33
    # sacc_Astro_A     87      0    390      8      0    224     38
    # sacc_Astro_B     85      0     19     23      0     28      5
    #       Note: br5182 not represented bc it was only used for an NAc-NeuN sample

# As in the step03's, re-create 'logcounts'
sce.astro.hold <- sce.astro
assay(sce.astro, "logcounts") <- NULL
sizeFactors(sce.astro) <- NULL
sce.astro <- logNormCounts(sce.astro)


## PW markers?
mod <- with(colData(sce.astro), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

# Run pairwise t-tests
markers.astro.t.pw <- findMarkers(sce.astro, groups=sce.astro$cellType,
                                assay.type="logcounts", design=mod, test="t",
                                direction="up", pval.type="all", full.stats=T)

sapply(markers.astro.t.pw, function(x){table(x$FDR<0.05)})
    #       amy_Astro_A amy_Astro_B dlpfc_Astro hpc_Astro_A hpc_Astro_B nac_Astro_A
    # FALSE       33520       33403       33534       33534       33368       33104
    # TRUE           18         135           4           4         170         434
    #       nac_Astro_B sacc_Astro_A sacc_Astro_B
    # FALSE       33260        33432        33532
    # TRUE          278          106            6

markerList.astro <- lapply(markers.astro.t.pw, function(x){
  rownames(x)[x$FDR < 0.05]
  })

# non-0-median
amy_astro_B <- which(sce.astro$cellType == "amy_Astro_B")
non0median.amy.As_B <- apply(as.matrix(assay(sce.astro, "logcounts")), 1, function(y){
  median(y[amy_astro_B]) > 0
})

table(non0median.amy.As_B)  # 146
markers.amy.As_B <- markers.astro.t.pw[["amy_Astro_B"]]
markers.amy.As_B <- cbind(markers.amy.As_B, non0median.amy.As_B[match(rownames(markers.amy.As_B),
                                                                      names(non0median.amy.As_B))])
colnames(markers.amy.As_B)[12] <- "non0median"
table(markers.amy.As_B$FDR < 0.05 & markers.amy.As_B$non0median==TRUE)  # just 4
rownames(markers.amy.As_B)[markers.amy.As_B$FDR < 0.05 & markers.amy.As_B$non0median==TRUE]
    # [1] "DST"     "COL19A1" "MACF1"   "RBFOX1"

# Check out these
plotExpressionCustom(sce.astro.hold, anno_name="cellType", features_name="Astro sub-class",
                     features=c("DST", "COL19A1", "MACF1", "RBFOX1"), ncol=2)

sce.astro$prelimCluster <- droplevels(sce.astro$prelimCluster)
table(sce.astro$cellType, sce.astro$prelimCluster)

sapply(splitit(sce.astro$region), function(x){table(droplevels(sce.astro$cellType[x]),
                                                    droplevels(sce.astro$prelimCluster[x]))})
    # $amy
    #               8  17  18  27  38  52  55
    # amy_Astro_A 836 131 107 347  90   0  44
    # amy_Astro_B   0   0   0   0   0  83   0
    # 
    # $dlpfc  * Note this:
    #              10  32  49  64  65  78  88  98
    # dlpfc_Astro  64  32 109  65 205 230  47  30
    # 
    # $hpc
    #               2   7  13  15  29  31  35  38
    # hpc_Astro_A 302 112   0 168  32 231  91   0
    # hpc_Astro_B   0   0 117   0   0   0   0 117
    # 
    # $nac
    #               16   17
    # nac_Astro_A   99    0
    # nac_Astro_B    0 1000
    # 
    # $sacc
    #               14  28  47
    # sacc_Astro_A 641   0 106
    # sacc_Astro_B   0 160   0



## AMY 'Inhib_B' vs DLPFC 'Inhib_A' - r=0.86 ===
sce.dlpfc <- sce.allRegions[ ,sce.allRegions$region=="dlpfc"]
sce.dlpfc$cellType <- droplevels(sce.dlpfc$cellType)

# Plot some AMY 'Inhib_B' markers (seen in Louise's DLPFC 'Inhib_A' lists)
plotExpressionCustom(sce.dlpfc, anno_name="cellType", features_name="some AMY 'Inhib_B'",
                     features=c("VIP", "CALB2", "CRH", "PTHLH"), ncol=2)



# Session info for 18Jul2021 ==============
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
#   [1] pheatmap_1.0.12             RColorBrewer_1.1-2          lattice_0.20-41            
# [4] limma_3.46.0                jaffelab_0.99.30            rafalib_1.0.0              
# [7] DropletUtils_1.10.3         batchelor_1.6.3             scran_1.18.7               
# [10] scater_1.18.6               ggplot2_3.3.3               org.Hs.eg.db_3.12.0        
# [13] EnsDb.Hsapiens.v86_2.99.0   ensembldb_2.14.1            AnnotationFilter_1.14.0    
# [16] GenomicFeatures_1.42.3      AnnotationDbi_1.52.0        SingleCellExperiment_1.12.0
# [19] SummarizedExperiment_1.20.0 Biobase_2.50.0              GenomicRanges_1.42.0       
# [22] GenomeInfoDb_1.26.7         IRanges_2.24.1              S4Vectors_0.28.1           
# [25] BiocGenerics_0.36.1         MatrixGenerics_1.2.1        matrixStats_0.58.0         
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
# [25] dqrng_0.3.0               assertthat_0.2.1          Matrix_1.3-4             
# [28] fastmap_1.1.0             lazyeval_0.2.2            BiocSingular_1.6.0       
# [31] prettyunits_1.1.1         tools_4.0.4               rsvd_1.0.5               
# [34] igraph_1.2.6              gtable_0.3.0              glue_1.4.2               
# [37] GenomeInfoDbData_1.2.4    dplyr_1.0.5               rappdirs_0.3.3           
# [40] Rcpp_1.0.6                vctrs_0.3.8               Biostrings_2.58.0        
# [43] rhdf5filters_1.2.0        rtracklayer_1.50.0        DelayedMatrixStats_1.12.3
# [46] stringr_1.4.0             beachmat_2.6.4            lifecycle_1.0.0          
# [49] irlba_2.3.3               statmod_1.4.35            XML_3.99-0.6             
# [52] edgeR_3.32.1              zlibbioc_1.36.0           scales_1.1.1             
# [55] hms_1.0.0                 ProtGenerics_1.22.0       rhdf5_2.34.0             
# [58] curl_4.3                  memoise_2.0.0             gridExtra_2.3            
# [61] segmented_1.3-4           biomaRt_2.46.3            stringi_1.5.3            
# [64] RSQLite_2.2.7             BiocParallel_1.24.1       rlang_0.4.11             
# [67] pkgconfig_2.0.3           bitops_1.0-7              purrr_0.3.4              
# [70] Rhdf5lib_1.12.1           GenomicAlignments_1.26.0  bit_4.0.4                
# [73] tidyselect_1.1.1          magrittr_2.0.1            R6_2.5.0                 
# [76] generics_0.1.0            DelayedArray_0.16.3       DBI_1.1.1                
# [79] pillar_1.6.0              withr_2.4.2               RCurl_1.98-1.3           
# [82] tibble_3.1.1              crayon_1.4.1              utf8_1.2.1               
# [85] BiocFileCache_1.14.0      viridis_0.6.0             progress_1.2.2           
# [88] locfit_1.5-9.4            grid_4.0.4                blob_1.2.1               
# [91] R.utils_2.10.1            openssl_1.4.3             munsell_0.5.0            
# [94] beeswarm_0.4.0            viridisLite_0.4.0         vipor_0.4.5              
# [97] askpass_1.1 

