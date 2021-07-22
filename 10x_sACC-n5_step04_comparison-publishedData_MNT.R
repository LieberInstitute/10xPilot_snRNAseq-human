### MNT 10x snRNA-seq workflow: step 04
###   **Region-specific analyses**
###     - (5x) sACC samples (M & F donors)
###     - Comparison to Velmeshev, et al (Science 2019)
#####################################################################

library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(batchelor)
library(DropletUtils)
library(jaffelab)
library(limma)
library(pheatmap)
library(RColorBrewer)

### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

# ===

### Comparison to Velmeshev, et al (PFC & ACC) ========
## Load within-ACC statistics
load("/dcl01/ajaffe/data/lab/singleCell/velmeshev2019/analysis_MNT/markers-stats_velmeshev-et-al_ASD-cortex-withinRegion_findMarkers-SN-LEVEL_MNTAug2020.rda",
     verbose=T)
    # markers.asdVelm.t.pfc, markers.asdVelm.t.acc
    #rm(markers.asdVelm.t.pfc)

load("/dcl01/ajaffe/data/lab/singleCell/velmeshev2019/analysis_MNT/SCE_asd-velmeshev-etal_MNT.rda", verbose=T)
    # sce.asd, sce.asd.pfc, sce.asd.acc  
sce.asd
    # class: SingleCellExperiment 
    # dim: 36501 104559 
    # metadata(0):
    # assays(1): logcounts
    # rownames(36501): ENSG00000227232 ENSG00000243485 ... ENSG00000210195
    # ENSG00000210196
    # rowData names(0):
    # colnames(104559): AAACCTGGTACGCACC-1_1823_BA24
    #   AAACGGGCACCAGATT-1_1823_BA24 ... TTTGTCATCCCAAGTA-1_6033_BA9
    #   TTTGTCATCGTTACGA-1_6033_BA9
    # colData names(17): cell cluster ... BAregion contrast
    # reducedDimNames(1): TSNE.prov
    # altExpNames(0):

sce.asd.pfc <- sce.asd[ ,sce.asd$region=="PFC"]
sce.asd.acc <- sce.asd[ ,sce.asd$region=="ACC"]
rm(sce.asd)

# Need to convert Symbol in sce.dlpfc > EnsemblID, and also use n nuclei for t.stat
load("rdas/revision/regionSpecific_sACC-n5_cleaned-combined_SCE_MNT2021.rda", verbose=T)
    # sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc, ref.sampleInfo, annotationTab.sacc, cell_colors.sacc
sce.sacc
    # class: SingleCellExperiment 
    # dim: 33538 15669 
    # metadata(2): merge.info pca.info
    # assays(2): counts logcounts
    # rownames(33538): MIR1302-2HG FAM138A ... AC213203.1 FAM231C
    # rowData names(6): gene_id gene_version ... gene_biotype Symbol.uniq
    # colnames(15669): AAACCCAAGAGTCTTC-1 AAACCCAAGGCCCGTT-1 ...
    #   TTTGACTGTATCGTGT-1 TTTGGTTAGCAGCACA-1
    # colData names(20): Sample Barcode ... collapsedCluster cellType
    # reducedDimNames(4): PCA_corrected PCA_opt TSNE UMAP
    # altExpNames(0):

table(sce.sacc$cellType)
    #   Astro_A        Astro_B   drop.doublet    drop.lowNTx        Excit_A        Excit_B 
    #       747            160             28            298            856            575 
    #   Excit_C        Excit_D        Excit_E        Excit_F        Excit_G        Inhib_A 
    #      1735            311            428            228             30            842 
    #   Inhib_B        Inhib_C        Inhib_D        Inhib_E        Inhib_F        Inhib_G 
    #       912            465            384            330            521            206 
    #   Inhib_H        Inhib_I        Inhib_J        Inhib_K          Micro Neu_FAT2.CDH15 
    #       208             39             42             25            784             20 
    #   Oligo_A        Oligo_B            OPC 
    #      4389            195            911

# Drop the flagged "drop." clusters 
sce.sacc <- sce.sacc[ ,-grep("drop.", sce.sacc$cellType)]
sce.sacc$cellType <- droplevels(sce.sacc$cellType)


## Load LIBD sACC stats (just need the "1vAll" result)
load("rdas/revision/markers-stats_sACC-n5_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.sacc.t.pw, markers.sacc.wilcox.block, markers.sacc.t.1vAll, medianNon0.sacc
    rm(markers.sacc.t.pw, markers.sacc.wilcox.block, markers.sacc.binom.block)

markers.sacc.enriched <- lapply(markers.sacc.t.1vAll, function(x){x[[2]]})
    
for(i in names(markers.sacc.enriched)){
  rownames(markers.sacc.enriched[[i]]) <- rowData(sce.sacc)$gene_id[match(rownames(markers.sacc.enriched[[i]]),
                                                                    rownames(sce.sacc))]
}


## Calculate and add t-statistic (= std.logFC * sqrt(N)) from contrasts
#      and fix row order to the first entry "Astro"
fixTo <- rownames(markers.sacc.enriched[["Astro_A"]])
for(s in names(markers.sacc.enriched)){
  markers.sacc.enriched[[s]]$t.stat <- markers.sacc.enriched[[s]]$std.logFC * sqrt(ncol(sce.sacc))
  markers.sacc.enriched[[s]] <- markers.sacc.enriched[[s]][fixTo, ]
}

# Pull out the t's
ts.sacc <- sapply(markers.sacc.enriched, function(x){x$t.stat})
rownames(ts.sacc) <- fixTo



## Then for Velmeshev et al. - fix row order to the first entry "AST-FB"
fixTo <- rownames(markers.asdVelm.t.acc[["AST-FB"]])

for(s in names(markers.asdVelm.t.acc)){
  markers.asdVelm.t.acc[[s]]$t.stat <- markers.asdVelm.t.acc[[s]]$std.logFC * sqrt(ncol(sce.asd.acc))
  markers.asdVelm.t.acc[[s]] <- markers.asdVelm.t.acc[[s]][fixTo, ]
}

# Pull out the t's
ts.velmeshev.acc <- sapply(markers.asdVelm.t.acc, function(x){x$t.stat})
rownames(ts.velmeshev.acc) <- fixTo


## Take intersecting between two and subset/reorder
sharedGenes <- intersect(rownames(ts.velmeshev.acc), rownames(ts.sacc))
length(sharedGenes) # 27,442

ts.velmeshev.acc <- ts.velmeshev.acc[sharedGenes, ]
ts.sacc <- ts.sacc[sharedGenes, ]


cor_t_sacc <- cor(ts.sacc, ts.velmeshev.acc)
rownames(cor_t_sacc) = paste0(rownames(cor_t_sacc),"_","libd")
colnames(cor_t_sacc) = paste0(colnames(cor_t_sacc),"_","asd.acc")
range(cor_t_sacc)


## Heatmap
theSeq.all = seq(-.95, .95, by = 0.01)
my.col.all <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.all)-1)

pdf("pdfs/exploration/Velmeshev-ASD_pfc-acc/overlap-velmeshev-ASD-acc_with_LIBD-10x-sACC_Aug2020.pdf")
pheatmap(cor_t_sacc,
         color=my.col.all,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.all,
         fontsize=10, fontsize_row=11, fontsize_col=11,
         display_numbers=T, number_format="%.2f", fontsize_number=6.0,
         legend_breaks=c(seq(-0.95,0.95,by=0.475)),
         main="Correlation of cluster-specific t's between LIBD sACC to \n ACC from (Velmeshev et al. Science 2019)")
dev.off()




### What if compared between both the .acc set of stats vs the .pfc?? =============

## Set up PFC t's
fixTo <- rownames(markers.asdVelm.t.pfc[["AST-FB"]])

for(s in names(markers.asdVelm.t.pfc)){
  markers.asdVelm.t.pfc[[s]]$t.stat <- markers.asdVelm.t.pfc[[s]]$std.logFC * sqrt(ncol(sce.asd.pfc))
  markers.asdVelm.t.pfc[[s]] <- markers.asdVelm.t.pfc[[s]][fixTo, ]
}

# Pull out the t's
ts.velmeshev.pfc <- sapply(markers.asdVelm.t.pfc, function(x){x$t.stat})
rownames(ts.velmeshev.pfc) <- fixTo

sharedGenes.all <- intersect(rownames(ts.velmeshev.acc), rownames(ts.sacc))
sharedGenes.all <- intersect(sharedGenes.all, rownames(ts.velmeshev.pfc))
    # of length 27,890

# Subset/order
ts.sacc <- ts.sacc[sharedGenes.all, ]
ts.velmeshev.pfc <- ts.velmeshev.pfc[sharedGenes.all, ]
ts.velmeshev.acc <- ts.velmeshev.acc[sharedGenes.all, ]

colnames(ts.velmeshev.pfc) <- paste0(colnames(ts.velmeshev.pfc),"_pfc")
colnames(ts.velmeshev.acc) <- paste0(colnames(ts.velmeshev.acc),"_acc")

ts.velmeshev.full <- cbind(ts.velmeshev.pfc, ts.velmeshev.acc)

cor_t_sacc.asd <- cor(ts.sacc, ts.velmeshev.full)
range(cor_t_sacc.asd)


## Heatmap
# Add some cluster info for add'l heatmap annotations
regionInfo <- data.frame(region=ss(colnames(ts.velmeshev.full), "_",2))
rownames(regionInfo) <- colnames(ts.velmeshev.full)

# Re-name that "Neu_FAT2.CDH15" to "Neu_ambig", like for Fig 3
rownames(cor_t_sacc.asd)[rownames(cor_t_sacc.asd)=="Neu_FAT2.CDH15"] <- "Neu_ambig"

theSeq.all = seq(-.95, .95, by = 0.01)
my.col.all <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.all)-1)

pdf("pdfs/revision/exploration/overlap-velmeshev-ASD-bothRegions_with_LIBD-10x-sACC_MNT2021.pdf", width=10)
pheatmap(cor_t_sacc.asd,
         color=my.col.all,
         annotation_col=regionInfo,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.all,
         fontsize=10.5, fontsize_row=11, fontsize_col=10,
         display_numbers=T, number_format="%.2f", fontsize_number=4.5,
         legend_breaks=c(seq(-0.95,0.95,by=0.475)),
         main="Correlation of cluster-specific t's between LIBD sACC to \n ACC & PFC from (Velmeshev et al. Science 2019)")
dev.off()


### Session info for 22Jul2021 ==================================================
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
#   [1] RColorBrewer_1.1-2          pheatmap_1.0.12             limma_3.46.0               
# [4] jaffelab_0.99.30            rafalib_1.0.0               DropletUtils_1.10.3        
# [7] batchelor_1.6.3             scran_1.18.7                scater_1.18.6              
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
# [67] pkgconfig_2.0.3           bitops_1.0-7              lattice_0.20-41          
# [70] purrr_0.3.4               Rhdf5lib_1.12.1           GenomicAlignments_1.26.0 
# [73] bit_4.0.4                 tidyselect_1.1.1          magrittr_2.0.1           
# [76] R6_2.5.0                  generics_0.1.0            DelayedArray_0.16.3      
# [79] DBI_1.1.1                 pillar_1.6.0              withr_2.4.2              
# [82] RCurl_1.98-1.3            tibble_3.1.1              crayon_1.4.1             
# [85] utf8_1.2.1                BiocFileCache_1.14.0      viridis_0.6.0            
# [88] progress_1.2.2            locfit_1.5-9.4            grid_4.0.4               
# [91] blob_1.2.1                R.utils_2.10.1            openssl_1.4.3            
# [94] munsell_0.5.0             beeswarm_0.4.0            viridisLite_0.4.0        
# [97] vipor_0.4.5               askpass_1.1 
