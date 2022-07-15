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
load("/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_HPC-n3_cleaned-combined_SCE_MNT2021.rda")
    # sce.hpc, clusterRefTab.hpc, chosen.hvgs.hpc, ref.sampleInfo

table(sce.hpc$cellType)
    #  Astro_A       Astro_B  drop.doublet drop.lowNTx_A drop.lowNTx_B 
    #      936           234             5           105            19 
    #  Excit_A       Excit_B       Excit_C       Excit_D       Excit_E 
    #       87           421             6            35             6 
    #  Excit_F       Excit_G       Excit_H       Inhib_A       Inhib_B 
    #       29             6            33           300            30 
    #  Inhib_C       Inhib_D         Micro         Mural         Oligo 
    #        5            31          1161            43          5912 
    #      OPC       OPC_COP         Tcell 
    #      823            15            26 

# First drop decided "drop." clusters (129 nuclei)
sce.hpc <- sce.hpc[ ,-grep("drop.", sce.hpc$cellType)]
sce.hpc$cellType <- droplevels(sce.hpc$cellType)

# Remove 0 genes across all nuclei
sce.hpc <- sce.hpc[!rowSums(assay(sce.hpc, "counts"))==0, ]  # keeps same 28764 genes

## Re-create 'logcounts' (don't want to use 'multiBatchNorm's down-scaling across donor 'batches')

# First 'hold' the MBN 'logcounts' for printing
sce.hold <- sce.hpc

assay(sce.hpc, "logcounts") <- NULL
sizeFactors(sce.hpc) <- NULL
sce.hpc <- logNormCounts(sce.hpc)


### First make a list of Boolean param / cell subtype ===
  # Will use this to assess more 'valid', non-noise-driving markers
cellSubtype.idx <- splitit(sce.hpc$cellType)
medianNon0.hpc <- lapply(cellSubtype.idx, function(x){
  apply(as.matrix(assay(sce.hpc, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

sapply(medianNon0.hpc, table)


## Traditional t-test implementation ===
mod <- with(colData(sce.hpc), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`


# Run pairwise t-tests
markers.hpc.t.pw <- findMarkers(sce.hpc, groups=sce.hpc$cellType,
                                assay.type="logcounts", design=mod, test="t",
                                direction="up", pval.type="all", full.stats=T)

sapply(markers.hpc.t.pw, function(x){table(x$FDR<0.05)})
    #       Astro_A Astro_B Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Excit_G
    # FALSE   28622   28634   28701   28741   28534   28685   28653   28605   28576
    # TRUE      142     130      63      23     230      79     111     159     188
    #       Excit_H Inhib_A Inhib_B Inhib_C Inhib_D Micro Mural Oligo   OPC OPC_COP
    # FALSE   28607   28757   28664   28666   28690 28474 28383 28673 28701   28561
    # TRUE      157       7     100      98      74   290   381    91    63     203
    #         Tcell
    #   FALSE 28294
    #   TRUE    470


## WMW: Blocking on donor (this test doesn't take 'design=' argument) ===
markers.hpc.wilcox.block <- findMarkers(sce.hpc, groups=sce.hpc$cellType,
                                          assay.type="logcounts", block=sce.hpc$donor, test="wilcox",
                                          direction="up", pval.type="all", full.stats=T)

sapply(markers.hpc.wilcox.block, function(x){table(x$FDR<0.05)})
    # No results... disregard these


## Binomial ===
markers.hpc.binom.block <- findMarkers(sce.hpc, groups=sce.hpc$cellType,
                                         assay.type="logcounts", block=sce.hpc$donor, test="binom",
                                         direction="up", pval.type="all", full.stats=T)

sapply(markers.hpc.binom.block, function(x){table(x$FDR<0.05)})
    # Also no results... disregard these 


# Add respective 'non0median' column to the stats for each set of markers
for(i in names(markers.hpc.t.pw)){
  markers.hpc.t.pw[[i]] <- cbind(markers.hpc.t.pw[[i]],
                                 medianNon0.hpc[[i]][match(rownames(markers.hpc.t.pw[[i]]),
                                                           names(medianNon0.hpc[[i]]))])
  colnames(markers.hpc.t.pw[[i]])[23] <- "non0median"
}

sapply(markers.hpc.t.pw, function(x){table(x$FDR<0.05 & x$non0median == TRUE)["TRUE"]})
    # Astro_A.TRUE Astro_B.TRUE Excit_A.TRUE Excit_B.TRUE Excit_C.TRUE Excit_D.TRUE 
    #          124           83           46           13           57           40 
    # Excit_E.TRUE Excit_F.TRUE Excit_G.TRUE Excit_H.TRUE Inhib_A.TRUE Inhib_B.TRUE 
    #           48           61           55           52            1           44 
    # Inhib_C.TRUE Inhib_D.TRUE   Micro.TRUE   Mural.TRUE   Oligo.TRUE     OPC.TRUE 
    #           30           27          193           59           91           53 
    # OPC_COP.TRUE   Tcell.TRUE 
    #          101          114



## Save all these for future reference ===
save(markers.hpc.t.pw, #markers.hpc.wilcox.block, #markers.hpc.binom.block,
     medianNon0.hpc,
     file="rdas/revision/markers-stats_HPC-n3_findMarkers-SN-LEVEL_MNT2021.rda")



# Print these to pngs
markerList.t.pw <- lapply(markers.hpc.t.pw, function(x){
  rownames(x)[x$FDR < 0.05 & x$non0median == TRUE]
  }
)

genes.top40.t <- lapply(markerList.t.pw, function(x){head(x, n=40)})


#dir.create("pdfs/revision/HPC/")
smaller.set <- names(genes.top40.t)[lengths(genes.top40.t) <= 20]
left.set <- setdiff(names(genes.top40.t), smaller.set)

# Smaller graphical window
for(i in smaller.set){
  png(paste0("pdfs/revision/HPC/HPC_t_pairwise_top40markers-", i, "_logExprs_MNT2021.png"), height=950, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         features = genes.top40.t[[i]], 
                         features_name = i,
                         anno_name = "cellType",
                         ncol=5, point_alpha=0.4) +
      scale_color_manual(values = cell_colors.hpc) +  
      ggtitle(label=paste0(i, " top markers: single-nucleus-level p.w. t-tests (FDR<0.05)"))
  )
  dev.off()
}

# 20-40 markers
for(i in left.set){
  png(paste0("pdfs/revision/HPC/HPC_t_pairwise_top40markers-", i, "_logExprs_MNT2021.png"), height=1900, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         features = genes.top40.t[[i]], 
                         features_name = i,
                         anno_name = "cellType",
                         ncol=5, point_alpha=0.4) +
      scale_color_manual(values = cell_colors.hpc) +  
      ggtitle(label=paste0(i, " top markers: single-nucleus-level p.w. t-tests (FDR<0.05)"))
  )
  dev.off()
}


#source('plotExpressionCustom.R')


### Cluster-vs-all single-nucleus-level iteration ================================

## Load SCE with new info
load("/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_HPC-n3_cleaned-combined_SCE_MNT2021.rda")
# sce.hpc, clusterRefTab.hpc, chosen.hvgs.hpc, ref.sampleInfo

# First drop decided "drop." clusters (129 nuclei)
sce.hpc <- sce.hpc[ ,-grep("drop.", sce.hpc$cellType)]
sce.hpc$cellType <- droplevels(sce.hpc$cellType)

# Remove 0 genes across all nuclei
sce.hpc <- sce.hpc[!rowSums(assay(sce.hpc, "counts"))==0, ]  # keeps same 28764 genes

## Re-create 'logcounts' (don't want to use 'multiBatchNorm's down-scaling across donor 'batches')
# First 'hold' the MBN 'logcounts' for printing
sce.hold <- sce.hpc

assay(sce.hpc, "logcounts") <- NULL
sizeFactors(sce.hpc) <- NULL
sce.hpc <- logNormCounts(sce.hpc)


## Load pw marker stats .rda with the non0median Booleans/cellType
load("rdas/revision/markers-stats_HPC-n3_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.hpc.t.pw, medianNon0.hpc

## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.hpc), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

markers.hpc.t.1vAll <- list()
for(i in levels(sce.hpc$cellType)){
  # Make temporary contrast
  sce.hpc$contrast <- ifelse(sce.hpc$cellType==i, 1, 0)
  # Test cluster vs. all others
  markers.hpc.t.1vAll[[i]] <- findMarkers(sce.hpc, groups=sce.hpc$contrast,
                                          assay.type="logcounts", design=mod, test="t",
                                          std.lfc=TRUE,
                                          direction="up", pval.type="all", full.stats=T)
}
    ## Since all other stats are the same, and don't really use the non-standardized
     #    logFC, just generate one object, unlike before

class(markers.hpc.t.1vAll[["Oligo"]])
    # a SimpleList of length 2, named "0" and "1" (from the temporary 'contrast')
    # -> we want the second entry, named "1"
    #    (for other purposes, might be interesting to look into that "0" entry, which
    #     is basically what genes are depleted in the cell type of interest)


sapply(markers.hpc.t.1vAll, function(x){
  table(x[["1"]]$stats.0$log.FDR < log(.001))
})
    #      Astro_A Astro_B Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Excit_G Excit_H
    # FALSE   23674   25280   24676   21295   27759   27031   27943   25816   28162   26815
    # TRUE     5090    3484    4088    7469    1005    1733     821    2948     602    1949
    
    #       Inhib_A Inhib_B Inhib_C Inhib_D Micro Mural Oligo   OPC OPC_COP Tcell
    # FALSE   21767   26520   28193   26797 23995 26934 25713 24731   28086 27391
    # TRUE     6997    2244     571    1967  4769  1830  3051  4033     678  1373

# Do some reorganizing
markers.hpc.t.1vAll <- lapply(markers.hpc.t.1vAll, function(x){
  # Basically take the 'stats.[1 or 0]' since is redundant with the 'summary'-level stats
  lapply(x, function(y){ y[ ,4] }) 
  })

# Re-name std.lfc column and the entries; add non-0-median info
for(i in names(markers.hpc.t.1vAll)){
  colnames(markers.hpc.t.1vAll[[i]][["0"]])[1] <- "std.logFC"
  colnames(markers.hpc.t.1vAll[[i]][["1"]])[1] <- "std.logFC"
  # Add non0median Boolean - might be informative for both sets of stats
  markers.hpc.t.1vAll[[i]][["0"]] <- cbind(markers.hpc.t.1vAll[[i]][["0"]],
                                 medianNon0.hpc[[i]][match(rownames(markers.hpc.t.1vAll[[i]][["0"]]),
                                                           names(medianNon0.hpc[[i]]))])
  colnames(markers.hpc.t.1vAll[[i]][["0"]])[4] <- "non0median"
  
  # "1" aka 'enriched'
  markers.hpc.t.1vAll[[i]][["1"]] <- cbind(markers.hpc.t.1vAll[[i]][["1"]],
                                           medianNon0.hpc[[i]][match(rownames(markers.hpc.t.1vAll[[i]][["1"]]),
                                                                     names(medianNon0.hpc[[i]]))])
  colnames(markers.hpc.t.1vAll[[i]][["1"]])[4] <- "non0median"
  
  # Then re-name the entries to more interpretable, because we'll keeping both contrasts
  names(markers.hpc.t.1vAll[[i]]) <- paste0(i,c("_depleted", "_enriched"))
}


## Let's save this along with the previous pairwise results
save(markers.hpc.t.pw, markers.hpc.t.1vAll, medianNon0.hpc,
     file="rdas/revision/markers-stats_HPC-n3_findMarkers-SN-LEVEL_MNT2021.rda")


## Print these to pngs
markerList.t.1vAll <- lapply(markers.hpc.t.1vAll, function(x){
  rownames(x[[2]])[ x[[2]]$log.FDR < log(0.05) & x[[2]]$non0median==TRUE ]
 }
)
genes.top40.t <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0("pdfs/revision/HPC/HPC_t_1vALL_top40markers-",i,"_logExprs_MNT2021.png"), height=1900, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         features = genes.top40.t[[i]], 
                         features_name = i,
                         anno_name = "cellType",
                         ncol=5, point_alpha=0.4) +
      scale_color_manual(values = cell_colors.hpc) +  
      ggtitle(label=paste0(i, " top markers: 'cluster-vs-all-others' t-tests (FDR<0.05)"))
  )
  dev.off()
}


## How do they intersect?
markerList.t.pw <- lapply(markers.hpc.t.pw, function(x){
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
    #Astro_A Astro_B Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Excit_G Excit_H Inhib_A 
    #     24      30      16       7      20      21      31      23      23      28       1 
    #Inhib_B Inhib_C Inhib_D   Micro   Mural   Oligo     OPC OPC_COP   Tcell 
    #     22      21      16      31      30      26      17      32      37 


    
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

write.csv(top40genes, file="tables/revision/top40genesLists_HPC-n3_cellType_SN-LEVEL-tests_MNT2021.csv",
          row.names=FALSE)





## Aside: add in 't.stat' as in 'step04' analyses to save for LoHu/LeCo ===
# for(s in names(markers.hpc.t.1vAll)){
#   markers.hpc.t.1vAll[[s]]$t.stat <- markers.hpc.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.hpc))
# }
# 
# save(markers.hpc.t.1vAll, markers.hpc.t.pw, sce.hpc,
#      file="rdas/markerStats-and-SCE_HPC-n3_sn-level_cleaned_MNTNov2020.rda")



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
    #   [1] limma_3.46.0                jaffelab_0.99.30            rafalib_1.0.0              
    # [4] DropletUtils_1.10.3         batchelor_1.6.2             scran_1.18.5               
    # [7] scater_1.18.6               ggplot2_3.3.3               EnsDb.Hsapiens.v86_2.99.0  
    # [10] ensembldb_2.14.1            AnnotationFilter_1.14.0     GenomicFeatures_1.42.3     
    # [13] AnnotationDbi_1.52.0        SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0
    # [16] Biobase_2.50.0              GenomicRanges_1.42.0        GenomeInfoDb_1.26.7        
    # [19] IRanges_2.24.1              S4Vectors_0.28.1            BiocGenerics_0.36.1        
    # [22] MatrixGenerics_1.2.1        matrixStats_0.58.0         
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
    # [28] fastmap_1.1.0             lazyeval_0.2.2            BiocSingular_1.6.0       
    # [31] prettyunits_1.1.1         tools_4.0.4               rsvd_1.0.3               
    # [34] igraph_1.2.6              gtable_0.3.0              glue_1.4.2               
    # [37] GenomeInfoDbData_1.2.4    dplyr_1.0.5               rappdirs_0.3.3           
    # [40] Rcpp_1.0.6                vctrs_0.3.6               Biostrings_2.58.0        
    # [43] rhdf5filters_1.2.0        rtracklayer_1.50.0        DelayedMatrixStats_1.12.3
    # [46] stringr_1.4.0             beachmat_2.6.4            lifecycle_1.0.0          
    # [49] irlba_2.3.3               statmod_1.4.35            XML_3.99-0.6             
    # [52] edgeR_3.32.1              zlibbioc_1.36.0           scales_1.1.1             
    # [55] hms_1.0.0                 ProtGenerics_1.22.0       rhdf5_2.34.0             
    # [58] RColorBrewer_1.1-2        curl_4.3                  memoise_2.0.0            
    # [61] gridExtra_2.3             segmented_1.3-3           biomaRt_2.46.3           
    # [64] stringi_1.5.3             RSQLite_2.2.7             BiocParallel_1.24.1      
    # [67] rlang_0.4.10              pkgconfig_2.0.3           bitops_1.0-7             
    # [70] lattice_0.20-41           purrr_0.3.4               Rhdf5lib_1.12.1          
    # [73] labeling_0.4.2            GenomicAlignments_1.26.0  cowplot_1.1.1            
    # [76] bit_4.0.4                 tidyselect_1.1.1          magrittr_2.0.1           
    # [79] R6_2.5.0                  generics_0.1.0            DelayedArray_0.16.3      
    # [82] DBI_1.1.1                 pillar_1.6.0              withr_2.4.2              
    # [85] RCurl_1.98-1.3            tibble_3.1.1              crayon_1.4.1             
    # [88] utf8_1.2.1                BiocFileCache_1.14.0      viridis_0.6.0            
    # [91] progress_1.2.2            locfit_1.5-9.4            grid_4.0.4               
    # [94] blob_1.2.1                digest_0.6.27             R.utils_2.10.1           
    # [97] openssl_1.4.3             munsell_0.5.0             beeswarm_0.3.1           
    # [100] viridisLite_0.4.0         vipor_0.4.5               askpass_1.1 
