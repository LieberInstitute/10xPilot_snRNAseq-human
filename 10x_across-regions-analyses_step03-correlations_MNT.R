### MNT 10x snRNA-seq workflow: step 04 - downstream comparisons
###   **Pan-brain analyses**
###     - n=12 samples from 5 regions, up to three donors
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

load("rdas/revision/markers-stats_DLPFC-n3_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.t.pw, markers.wilcox.block, markers.dlpfc.t.1vAll, medianNon0.dlpfc
    rm(markers.t.pw, markers.wilcox.block, medianNon0.dlpfc)

load("rdas/revision/markers-stats_HPC-n3_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.hpc.t.pw, markers.hpc.t.1vAll, medianNon0.hpc
    rm(markers.hpc.t.pw, medianNon0.hpc)

load("rdas/revision/markers-stats_NAc-n8_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.nac.t.pw, markers.nac.t.1vAll, medianNon0.nac
    rm(markers.nac.t.pw, medianNon0.nac)

load("rdas/revision/markers-stats_Amyg-n5_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.amy.t.pw, markers.amy.wilcox.block, markers.amy.t.1vAll, medianNon0.amy
    rm(markers.amy.t.pw, markers.amy.wilcox.block, medianNon0.amy)

load("rdas/markers-stats_sACC-n2_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
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
        #   19    18    20    24    25
    
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


## Clean up those decided a priori to remove from reported set
FMstats.list[["hpc"]][["OPC_COP"]] <- NULL
FMstats.list[["nac"]][["OPC_COP"]] <- NULL
FMstats.list[["nac"]][["Macro_infilt"]] <- NULL
FMstats.list[["nac"]][["Micro_resting"]] <- NULL
sapply(FMstats.list, length)  # good


### Get n Nuclei numbers for each region so can compute t-statistics ===
  # This can be done with Cohen's D (the 'std.lfc'), as d = t/sqrt(N)

load("rdas/revision/all-n24-samples_across-regions-analyses_forFigOnly_MNT2021.rda", verbose=T)
    # sce.allRegions, chosen.hvgs.union, ref.sampleInfo, Readme

sce.allRegions
    # class: SingleCellExperiment
    # dim: 33538 70497
table(sce.allRegions$region)
    #  amy dlpfc   hpc   nac  sacc 
    #14039 11202 10124 19789 15343

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
# Add back in row names and region suffix
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
cor_t_xRegions <- cor(ts.fullMat)



### Heatmap - typically use levelplot (e.g. below), but will want pheatmap bc can cluster cols/rows
theSeq.all = seq(-1, 1, by = 0.025)
my.col.all <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq.all)-1)


pdf("pdfs/revision/acrossRegions_correlation_region-specific-subcluster-ts_MNT2021.pdf")
pheatmap(cor_t_xRegions,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=6.5, fontsize_col=6.5,
         main="Correlation of cluster-specific t's from all regions \n (all shared expressed genes)")
dev.off()


## Subset on neuronal subcluster t's and check
ts.fullMat.neu <- ts.fullMat
for(i in c("As", "Micro", "Endo", "Mural","Oligo", "OPC", "Tcell")){
  ts.fullMat.neu <- ts.fullMat.neu[ ,-grep(i, colnames(ts.fullMat.neu))]
}

cor_t_xRegions.neu <- cor(ts.fullMat.neu)


### Heatmap - typically use levelplot (e.g. below), but will want pheatmap bc can cluster cols/rows
theSeq.all = seq(-1, 1, by = 0.025)
my.col.all <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq.all)-1)

## Updated 31Aug2020 - for paper === === ===
# Add some cluster info for add'l heatmap annotations
clusterInfo <- data.frame(region=ss(colnames(ts.fullMat.neu), "_",3))
rownames(clusterInfo) <- colnames(ts.fullMat.neu)


# Print
pdf("pdfs/revision/acrossRegions_correlation_region-specific-NeuronalSubcluster-ts_MNT2021.pdf")
pheatmap(cor_t_xRegions.neu,
         annotation_col=clusterInfo,
         #show_colnames=FALSE,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=5, fontsize_col=5,
         main="Correlation of neuronal cluster-specific t's from all regions \n (all shared expressed genes)")
dev.off()



## Non-neuronal set for supplement === === ===
ts.fullMat.non <- ts.fullMat
glia.idx <- NA
for(i in c("As", "Micro", "Endo", "Mural","Oligo", "OPC", "Tcell")){
  glia.idx <- c(glia.idx, grep(i, colnames(ts.fullMat.non)))
}
# Rm the empty NA
glia.idx <- glia.idx[-1]
ts.fullMat.non <- ts.fullMat.non[ ,glia.idx]

cor_t_xRegions.non <- cor(ts.fullMat.non)


### Heatmap ===
theSeq.all = seq(-1, 1, by = 0.025)
my.col.all <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq.all)-1)

# Add some cluster info for add'l heatmap annotations
clusterInfo.glia <- data.frame(region=ifelse(is.na(ss(colnames(ts.fullMat.non), "_",3)),
                                             ss(colnames(ts.fullMat.non), "_",2),
                                             ss(colnames(ts.fullMat.non), "_",3))
                               )

rownames(clusterInfo.glia) <- colnames(ts.fullMat.non)

# annotColors <- list(class = tableau10medium[1:5][factor(unique(clusterInfo$class))],
#                     region = tableau10medium[6:10][factor(unique(clusterInfo$region))])
# 
# names(annotColors[["class"]]) <- unique(clusterInfo$class)
# names(annotColors[["region"]]) <- unique(clusterInfo$region)


# Print
pdf("pdfs/revision/acrossRegions_correlation_region-specific-NON-NeuronalSubcluster-ts_MNT2021.pdf",width=9)
pheatmap(cor_t_xRegions.non,
         annotation_col=clusterInfo.glia,
         # annotation_colors=annotColors,
         #show_colnames=FALSE,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=6, fontsize_col=6,
         display_numbers=TRUE, fontsize_number=5,
         main="Correlation of glia/other cluster-specific t's from all regions \n (all shared expressed genes)")
dev.off()


# Session info for 13Jul2021 ==============
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
#   [1] lattice_0.20-41             pheatmap_1.0.12             RColorBrewer_1.1-2         
# [4] gridExtra_2.3               dynamicTreeCut_1.63-1       dendextend_1.14.0          
# [7] jaffelab_0.99.30            rafalib_1.0.0               DropletUtils_1.10.3        
# [10] batchelor_1.6.3             scran_1.18.7                scater_1.18.6              
# [13] ggplot2_3.3.3               EnsDb.Hsapiens.v86_2.99.0   ensembldb_2.14.1           
# [16] AnnotationFilter_1.14.0     GenomicFeatures_1.42.3      AnnotationDbi_1.52.0       
# [19] SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0 Biobase_2.50.0             
# [22] GenomicRanges_1.42.0        GenomeInfoDb_1.26.7         IRanges_2.24.1             
# [25] S4Vectors_0.28.1            BiocGenerics_0.36.1         MatrixGenerics_1.2.1       
# [28] matrixStats_0.58.0         
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
# [28] fastmap_1.1.0             lazyeval_0.2.2            limma_3.46.0             
# [31] BiocSingular_1.6.0        prettyunits_1.1.1         tools_4.0.4              
# [34] rsvd_1.0.5                igraph_1.2.6              gtable_0.3.0             
# [37] glue_1.4.2                GenomeInfoDbData_1.2.4    dplyr_1.0.5              
# [40] rappdirs_0.3.3            Rcpp_1.0.6                vctrs_0.3.8              
# [43] Biostrings_2.58.0         rhdf5filters_1.2.0        rtracklayer_1.50.0       
# [46] DelayedMatrixStats_1.12.3 stringr_1.4.0             beachmat_2.6.4           
# [49] lifecycle_1.0.0           irlba_2.3.3               statmod_1.4.35           
# [52] XML_3.99-0.6              edgeR_3.32.1              zlibbioc_1.36.0          
# [55] scales_1.1.1              hms_1.0.0                 ProtGenerics_1.22.0      
# [58] rhdf5_2.34.0              curl_4.3                  memoise_2.0.0            
# [61] segmented_1.3-4           biomaRt_2.46.3            stringi_1.5.3            
# [64] RSQLite_2.2.7             BiocParallel_1.24.1       rlang_0.4.11             
# [67] pkgconfig_2.0.3           bitops_1.0-7              purrr_0.3.4              
# [70] Rhdf5lib_1.12.1           labeling_0.4.2            GenomicAlignments_1.26.0 
# [73] cowplot_1.1.1             bit_4.0.4                 tidyselect_1.1.1         
# [76] magrittr_2.0.1            R6_2.5.0                  generics_0.1.0           
# [79] DelayedArray_0.16.3       DBI_1.1.1                 pillar_1.6.0             
# [82] withr_2.4.2               RCurl_1.98-1.3            tibble_3.1.1             
# [85] crayon_1.4.1              utf8_1.2.1                BiocFileCache_1.14.0     
# [88] viridis_0.6.0             progress_1.2.2            locfit_1.5-9.4           
# [91] grid_4.0.4                blob_1.2.1                digest_0.6.27            
# [94] R.utils_2.10.1            openssl_1.4.3             munsell_0.5.0            
# [97] beeswarm_0.4.0            viridisLite_0.4.0         vipor_0.4.5              
# [100] askpass_1.1

