### MNT 10x snRNA-seq workflow: step 02
###   ** Across-regions analyses **
###     - (n=24) all regions from up to 8 donors:
###     - Amyg, DLPFC, HPC, NAc, and sACC
### Initiated MNT 07Feb2020
### 
#####################################################################

library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(batchelor)
library(DropletUtils)
library(jaffelab)
library(dendextend)
library(dynamicTreeCut)
library(gridExtra)

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



### Read in region-specific SCEs ===

## Amyg
load("rdas/revision/regionSpecific_Amyg-n5_cleaned-combined_SCE_MNT2021.rda", verbose=T)
    # sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo, annotationTab.amy, cell_colors.amy
    rm(pc.choice.amy, clusterRefTab.amy, annotationTab.amy, cell_colors.amy)

    
## DLPFC
load("rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda", verbose=T)
    # sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc, ref.sampleInfo, annotationTab.dlpfc, cell_colors
    rm(pc.choice.dlpfc, clusterRefTab.dlpfc, annotationTab.dlpfc, cell_colors)  

    
## HPC
load("rdas/revision/regionSpecific_HPC-n3_cleaned-combined_SCE_MNT2021.rda", verbose=T)
    # sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo, annotationTab.hpc, cell_colors.hpc
    rm(pc.choice.hpc, clusterRefTab.hpc, annotationTab.hpc, cell_colors.hpc)

    
## sACC
load("rdas/revision/regionSpecific_sACC-n5_cleaned-combined_SCE_MNT2021.rda", verbose=T)
    # sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc, ref.sampleInfo, annotationTab.sacc, cell_colors.sacc
    rm(pc.choice.sacc, clusterRefTab.sacc, annotationTab.sacc, cell_colors.sacc)

    
## NAc
load("rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda", verbose=T)
    # sce.nac, chosen.hvgs.nac, pc.choice.nac, ref.sampleInfo, annotationTab.nac, cell_colors.nac
    rm(pc.choice.nac, clusterRefTab.nac, annotationTab.nac, cell_colors.nac)
    

    
### Create a new 'across-regions' n=24 SCE ===
table(rownames(sce.sacc) == rownames(sce.hpc))  # all good

# Remove 'collapsedCluster' and various reducedDims, metadata; add "region_" to $cellType
sce.amy$collapsedCluster <- NULL
sce.hpc$collapsedCluster <- NULL
sce.sacc$collapsedCluster <- NULL
sce.dlpfc$collapsedCluster <- NULL

reducedDims(sce.nac) <- NULL
sizeFactors(sce.nac) <- NULL
metadata(sce.nac) <- list(NULL)
sce.nac$cellType <- paste0("nac_", sce.nac$cellType)

reducedDims(sce.amy) <- NULL
sizeFactors(sce.amy) <- NULL
metadata(sce.amy) <- list(NULL)
sce.amy$cellType <- paste0("amy_", sce.amy$cellType)

reducedDims(sce.hpc) <- NULL
sizeFactors(sce.hpc) <- NULL
metadata(sce.hpc) <- list(NULL)
sce.hpc$cellType <- paste0("hpc_", sce.hpc$cellType)

reducedDims(sce.sacc) <- NULL
sizeFactors(sce.sacc) <- NULL
metadata(sce.sacc) <- list(NULL)
sce.sacc$cellType <- paste0("sacc_", sce.sacc$cellType)

reducedDims(sce.dlpfc) <- NULL
sizeFactors(sce.dlpfc) <- NULL
metadata(sce.dlpfc) <- list(NULL)
sce.dlpfc$cellType <- paste0("dlpfc_", sce.dlpfc$cellType)


## cbind() them
sce.allRegions <- cbind(sce.nac, sce.amy, sce.hpc, sce.sacc, sce.dlpfc)
sce.allRegions
    # class: SingleCellExperiment 
    # dim: 33538 72887 
    # metadata(5): '' '' '' '' ''
    # assays(2): counts logcounts
    # rownames(33538): MIR1302-2HG FAM138A ... AC213203.1 FAM231C
    # rowData names(6): gene_id gene_version ... gene_biotype Symbol.uniq
    # colnames(72887): AAACCCACATCGAACT-1 AAACCCATCCAACCAA-1 ...
    # TTTGTTGGTGACCGAA-1 TTTGTTGTCTCGAACA-1
    # colData names(18): Sample Barcode ... prelimCluster cellType
    # reducedDimNames(0):
    #   altExpNames(0):

table(sce.allRegions$sampleID)
    #  br5161.amy     br5161.dlpfc       br5161.hpc       br5161.nac 
    #        3294             4215             4421             2055 
    # br5161.sacc  br5182.nac.neun     br5207.dlpfc  br5207.nac.neun 
    #        3174             4256             5294             4425 
    #  br5212.amy     br5212.dlpfc       br5212.hpc       br5212.nac 
    #        3259             1693             3977             1773 
    # br5212.sacc  br5276.amy.neun       br5276.nac br5276.sacc.neun 
    #        3880             2465             2626              851 
    #  br5287.hpc       br5287.nac  br5400.amy.neun       br5400.nac 
    #        1870              681             2635             4108 
    # br5400.sacc       br5701.amy  br5701.nac.neun br5701.sacc.neun 
    #        3959             3524              647             3805

sce.allRegions$cellType <- factor(sce.allRegions$cellType)
# Take union of 'chosen.hvgs'
chosen.hvgs.union <- chosen.hvgs.nac | chosen.hvgs.amy | chosen.hvgs.hpc | chosen.hvgs.sacc | chosen.hvgs.dlpfc

## Save this
save(sce.allRegions, chosen.hvgs.union, ref.sampleInfo, 
     file="rdas/revision/all-n24-samples_across-regions-analyses_MNT2021.rda")


### (Optional:) Dimensionality reduction =========================================

# First remove those which are not being reported in the main Results section
#   i.e. 'drop.' and the 'OPC_COP', 'Macro_infilt', 'Micro_resting'

remove.idx <- c(grep("drop.", sce.allRegions$cellType),
                grep("OPC_COP", sce.allRegions$cellType),
                grep("Macro_infilt", sce.allRegions$cellType),
                grep("Micro_resting", sce.allRegions$cellType))

sce.allRegions <- sce.allRegions[ ,-remove.idx]
sce.allRegions$cellType <- droplevels(sce.allRegions$cellType)

## Run `fastMNN` (internally uses `multiBatchPCA`), taking default top 50
 #  (Note: won't cluster in these dimensions - just for TSNE/UMAP)
set.seed(109)
mnn.hold <-  fastMNN(sce.allRegions, batch=sce.allRegions$donor,
                     merge.order=c("br5161","br5212","br5400","br5207",
                                   "br5701","br5276","br5182","br5287"),
                     subset.row=chosen.hvgs.union, d=50,
                     correct.all=TRUE, get.variance=TRUE,
                     BSPARAM=BiocSingular::IrlbaParam())
    # This temp file just used for getting batch-corrected components (drops a variety of entries)
    date()

table(colnames(mnn.hold) == colnames(sce.allRegions))  # all TRUE
table(mnn.hold$batch == sce.allRegions$donor)

# Add them to the SCE, as well as the metadata (though the latter might not be so usefl)
reducedDim(sce.allRegions, "PCA_corrected_50") <- reducedDim(mnn.hold, "corrected")
metadata(sce.allRegions) <- metadata(mnn.hold)
names(metadata(sce.allRegions)) <- paste0(names(metadata(sce.allRegions)),"_50")


## For options, re-run but computing 100 corrected PCs ===
    set.seed(109)
    mnn.hold <-  fastMNN(sce.allRegions, batch=sce.allRegions$donor,
                         merge.order=c("br5161","br5212","br5400","br5207",
                                       "br5701","br5276","br5182","br5287"),
                         subset.row=chosen.hvgs.union, d=100,
                         correct.all=TRUE, get.variance=TRUE,
                         BSPARAM=BiocSingular::IrlbaParam())
    date()

    dim(reducedDim(mnn.hold))
    table(colnames(mnn.hold) == colnames(sce.allRegions))  # all TRUE
    table(mnn.hold$batch == sce.allRegions$donor)
    
    # Add them to the SCE, as well as the metadata (though the latter might not be so usefl)
    reducedDim(sce.allRegions, "PCA_corrected_100") <- reducedDim(mnn.hold, "corrected")
    metadata(sce.allRegions)[["merge.info_100"]] <- metadata(mnn.hold)[[1]]
    metadata(sce.allRegions)[["pca.info_100"]] <- metadata(mnn.hold)[[2]]
    
# Save for now
Readme <- "This SCE is just for MNT processing for main Fig only. Has 'drop.' clusters, in addition to those rarer/unexpected-to-be-normal-brain-niche cell types removed."
save(sce.allRegions, chosen.hvgs.union, ref.sampleInfo, Readme,
     file="rdas/revision/all-n24-samples_across-regions-analyses_forFigOnly_MNT2021.rda")


## t-SNE
set.seed(109)
sce.allRegions <- runTSNE(sce.allRegions, dimred="PCA_corrected_50")

## UMAP
set.seed(109)
sce.allRegions <- runUMAP(sce.allRegions, dimred="PCA_corrected_100")


# Save
save(sce.allRegions, chosen.hvgs.union, ref.sampleInfo, Readme,
     file="rdas/revision/all-n24-samples_across-regions-analyses_forFigOnly_MNT2021.rda")



## For main Fig3: Facet some different iterations of the best tSNE by region

# pdf("pdfs/revision/pubFigures/across-regions-n24_tSNEon50PCs_faceted_MNT2021.pdf", width=9)
# plotTSNE(sce.allRegions, colour_by="region", point_alpha=0.5, point_size=4.0, theme_size=22) +
#   facet_wrap(~ sce.allRegions$region)
#   ggtitle("t-SNE on top 50 corrected PCs") + theme(plot.title = element_text(size=18))
# dev.off()
    
    ## More manually to have shadow of those for each region ======
    sce.temp <- sce.allRegions
    
    ## DLPFC
    # Reorder to plot the region nuclei last
    sce.temp <- cbind(sce.temp[ ,sce.temp$region!="dlpfc"],
                      sce.temp[ ,sce.temp$region=="dlpfc"])
    # Color by region-specific annot
    sce.temp$annot.temp <- ifelse(sce.temp$region=="dlpfc", ss(as.character(sce.temp$cellType),"dlpfc_",2), NA)
    
    p.dlpfc <- plotTSNE(sce.temp, colour_by="annot.temp", point_alpha=0.6, point_size=3.5, theme_size=15,
                        add_legend=FALSE) + #, text_by="annot.temp", text_size=4) +
      ggtitle("DLPFC") + scale_color_manual(values=cell_colors.dlpfc) +
      theme(plot.title = element_text(size=30),
            axis.title = element_text(size=0),
            axis.text = element_text(size=20))
    
    ## HPC
    # Reorder to plot the region nuclei last
    sce.temp <- cbind(sce.temp[ ,sce.temp$region!="hpc"],
                      sce.temp[ ,sce.temp$region=="hpc"])
    # Color by region-specific annot
    sce.temp$annot.temp <- ifelse(sce.temp$region=="hpc", ss(as.character(sce.temp$cellType),"hpc_",2), NA)
    
    p.hpc <- plotTSNE(sce.temp, colour_by="annot.temp", point_alpha=0.6, point_size=3.5, theme_size=15,
                      add_legend=FALSE) + #, text_by="annot.temp", text_size=4) +
      ggtitle("HPC") + scale_color_manual(values=cell_colors.hpc) +
      theme(plot.title = element_text(size=30),
            axis.title = element_text(size=0),
            axis.text = element_text(size=20))
    
    ## sACC
    # Reorder to plot the region nuclei last
    sce.temp <- cbind(sce.temp[ ,sce.temp$region!="sacc"],
                      sce.temp[ ,sce.temp$region=="sacc"])
    # Color by region-specific annot
    sce.temp$annot.temp <- ifelse(sce.temp$region=="sacc", ss(as.character(sce.temp$cellType),"sacc_",2), NA)
    
    p.sacc <- plotTSNE(sce.temp, colour_by="annot.temp", point_alpha=0.6, point_size=3.5, theme_size=15,
                       add_legend=FALSE) + #, text_by="annot.temp", text_size=4) +
      ggtitle("sACC") + scale_color_manual(values=cell_colors.sacc) +
      theme(plot.title = element_text(size=30),
            axis.title = element_text(size=0),
            axis.text = element_text(size=20))
    
    ## AMY
    # Reorder to plot the region nuclei last
    sce.temp <- cbind(sce.temp[ ,sce.temp$region!="amy"],
                      sce.temp[ ,sce.temp$region=="amy"])
    # Color by region-specific annot
    sce.temp$annot.temp <- ifelse(sce.temp$region=="amy", ss(as.character(sce.temp$cellType),"amy_",2), NA)
    
    p.amy <- plotTSNE(sce.temp, colour_by="annot.temp", point_alpha=0.6, point_size=3.5, theme_size=15,
                      add_legend=FALSE) + #, text_by="annot.temp", text_size=4) +
      ggtitle("AMY") + scale_color_manual(values=cell_colors.amy) +
      theme(plot.title = element_text(size=30),
            axis.title = element_text(size=0),
            axis.text = element_text(size=20))
    
    ## NAc
    # Reorder to plot the region nuclei last
    sce.temp <- cbind(sce.temp[ ,sce.temp$region!="nac"],
                      sce.temp[ ,sce.temp$region=="nac"])
    # Color by region-specific annot
    sce.temp$annot.temp <- ifelse(sce.temp$region=="nac", ss(as.character(sce.temp$cellType),"nac_",2), NA)
    
    p.nac <- plotTSNE(sce.temp, colour_by="annot.temp", point_alpha=0.6, point_size=3.5, theme_size=15,
                      add_legend=FALSE) + #, text_by="annot.temp", text_size=4) +
      ggtitle("NAc") + scale_color_manual(values=cell_colors.nac) +
      theme(plot.title = element_text(size=30),
            axis.title = element_text(size=0),
            axis.text = element_text(size=20))
    
    ## end region-colored t-SNEs ========
    
    
    ## All nuclei (by region) ===
    p.full <- plotTSNE(sce.allRegions, colour_by="region", point_alpha=0.2, point_size=4.0,
                       text_size=8, theme_size=24) +
      ggtitle("t-SNE on top 50 corrected PCs") + theme(plot.title = element_text(size=28))
    
    lay <- rbind(c(1,1,2),
                 c(1,1,3),
                 c(6,5,4))
    
    pdf("pdfs/revision/pubFigures/across-regions-n24_tSNEon50PCs_rs-cellClasses_faceted_MNT2021.pdf", width=13.5, height=12.5)
    #pdf("pdfs/revision/pubFigures/across-regions-n24_tSNEon50PCs_rs-cellClasses_faceted_labeled_MNT2021.pdf", width=13.5, height=12.5)
    grid.arrange(grobs=list(p.full,
                            p.nac,
                            p.amy,
                            p.hpc,
                            p.dlpfc,
                            p.sacc),
                 layout_matrix=lay)
    dev.off()
    
    

        
## Print broad marker heatmap of pan-brain-defined clusters === === ===
load("rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda", verbose=T)
    #sce.all.n12, chosen.hvgs.all.n12, pc.choice.n12, ref.sampleInfo, clusterRefTab.all.n12

# As decided for marker detection, remove the clusters that won't focus on:
#     'Ambig.hiVCAN' & 'Excit.4' & those .RS-annot'd as 'Ambig.lowNtrxts'
sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType.RS != "Ambig.lowNtrxts"] # 445
sce.all.n12$cellType.RS <- droplevels(sce.all.n12$cellType.RS)

sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType != "Ambig.hiVCAN"] # 32 nuclei
sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType != "Excit.4"]  # 33 nuclei
sce.all.n12$cellType <- droplevels(sce.all.n12$cellType)

cell.idx <- splitit(sce.all.n12$cellType)
dat <- as.matrix(assay(sce.all.n12, "logcounts"))


genes <- c('SNAP25','SLC17A6','SLC17A7','GAD1','GAD2','AQP4','GFAP','C3','CD74','MBP','PDGFRA','VCAN','CLDN5','FLT1','SKAP1','TRAC')
current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMeans(dat[genes, ii])))

current_dat <- current_dat[ ,c(3:14, 1,2, 15:17)]

pdf("pdfs/pubFigures/heatmap-geneExprs_panBrain-annot_mean-broadMarkers_MNT.pdf")
pheatmap(current_dat, cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(0.02, 4, length.out = 101),
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "OrRd"))(100),
         fontsize_row = 17.5, fontsize_col = 17.5,
         main="Broad cell type marker expression (pan-brain annot.)", fontsize=11)
dev.off()



### Session info for 12-13Jul2021 ============
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
#   [1] dynamicTreeCut_1.63-1       dendextend_1.14.0           jaffelab_0.99.30           
# [4] rafalib_1.0.0               DropletUtils_1.10.3         batchelor_1.6.3            
# [7] scran_1.18.7                scater_1.18.6               ggplot2_3.3.3              
# [10] EnsDb.Hsapiens.v86_2.99.0   ensembldb_2.14.1            AnnotationFilter_1.14.0    
# [13] GenomicFeatures_1.42.3      AnnotationDbi_1.52.0        SingleCellExperiment_1.12.0
# [16] SummarizedExperiment_1.20.0 Biobase_2.50.0              GenomicRanges_1.42.0       
# [19] GenomeInfoDb_1.26.7         IRanges_2.24.1              S4Vectors_0.28.1           
# [22] BiocGenerics_0.36.1         MatrixGenerics_1.2.1        matrixStats_0.58.0         
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15                googledrive_1.0.1         ggbeeswarm_0.6.0         
# [4] colorspace_2.0-0          ellipsis_0.3.2            scuttle_1.0.4            
# [7] bluster_1.0.0             XVector_0.30.0            BiocNeighbors_1.8.2      
# [10] rstudioapi_0.13           bit64_4.0.5               RSpectra_0.16-0          
# [13] fansi_0.4.2               xml2_1.3.2                codetools_0.2-18         
# [16] splines_4.0.4             R.methodsS3_1.8.1         sparseMatrixStats_1.2.1  
# [19] cachem_1.0.4              Rsamtools_2.6.0           ResidualMatrix_1.0.0     
# [22] dbplyr_2.1.1              R.oo_1.24.0               uwot_0.1.10              
# [25] HDF5Array_1.18.1          compiler_4.0.4            httr_1.4.2               
# [28] dqrng_0.3.0               assertthat_0.2.1          Matrix_1.3-4             
# [31] fastmap_1.1.0             lazyeval_0.2.2            limma_3.46.0             
# [34] BiocSingular_1.6.0        prettyunits_1.1.1         tools_4.0.4              
# [37] rsvd_1.0.5                igraph_1.2.6              gtable_0.3.0             
# [40] glue_1.4.2                GenomeInfoDbData_1.2.4    dplyr_1.0.5              
# [43] rappdirs_0.3.3            Rcpp_1.0.6                vctrs_0.3.8              
# [46] Biostrings_2.58.0         rhdf5filters_1.2.0        rtracklayer_1.50.0       
# [49] DelayedMatrixStats_1.12.3 stringr_1.4.0             beachmat_2.6.4           
# [52] lifecycle_1.0.0           irlba_2.3.3               statmod_1.4.35           
# [55] XML_3.99-0.6              edgeR_3.32.1              zlibbioc_1.36.0          
# [58] scales_1.1.1              hms_1.0.0                 ProtGenerics_1.22.0      
# [61] rhdf5_2.34.0              RColorBrewer_1.1-2        curl_4.3                 
# [64] memoise_2.0.0             gridExtra_2.3             segmented_1.3-4          
# [67] biomaRt_2.46.3            stringi_1.5.3             RSQLite_2.2.7            
# [70] BiocParallel_1.24.1       rlang_0.4.11              pkgconfig_2.0.3          
# [73] bitops_1.0-7              lattice_0.20-41           purrr_0.3.4              
# [76] Rhdf5lib_1.12.1           GenomicAlignments_1.26.0  bit_4.0.4                
# [79] tidyselect_1.1.1          RcppAnnoy_0.0.18          magrittr_2.0.1           
# [82] R6_2.5.0                  generics_0.1.0            DelayedArray_0.16.3      
# [85] DBI_1.1.1                 pillar_1.6.0              withr_2.4.2              
# [88] RCurl_1.98-1.3            tibble_3.1.1              crayon_1.4.1             
# [91] utf8_1.2.1                BiocFileCache_1.14.0      viridis_0.6.0            
# [94] progress_1.2.2            locfit_1.5-9.4            grid_4.0.4               
# [97] blob_1.2.1                R.utils_2.10.1            openssl_1.4.3            
# [100] munsell_0.5.0             beeswarm_0.4.0            viridisLite_0.4.0        
# [103] vipor_0.4.5               askpass_1.1 

