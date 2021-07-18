### LAH 10x snRNA-seq workflow: step 03 - marker detection
###   **Region-specific analyses**
###     - (3x) DLPFC samples
### 25May2021
#####################################################################

library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(batchelor)
library(DropletUtils)
library(jaffelab)
library(limma)
library(DeconvoBuddies)
library(here)

source("plotExpressionCustom.R")

load(here("rdas","revision","tableau_colors.rda"), verbose = TRUE)

## Load SCE with new info
load(here("rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda"), verbose=T)
# sce.dlpfc
# chosen.hvgs.dlpfc
# pc.choice.dlpfc
# clusterRefTab.dlpfc
# ref.sampleInfo
# annotationTab.dlpfc
# cell_colors

table(sce.dlpfc$cellType)
# Astro Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Inhib_A Inhib_B Inhib_C Inhib_D Inhib_E Inhib_F   Micro   Mural 
# 782     529     773     524     132     187     243     333     454     365     413       7       8     388      18 
# Oligo     OPC   Tcell 
# 5455     572      19 

dim(sce.dlpfc)
# [1] 33538 11202

# Remove 0 genes across all nuclei - there are None
sce.dlpfc <- sce.dlpfc[!rowSums(assay(sce.dlpfc, "counts"))==0, ]  # 29310

# First 'hold' the MBN 'logcounts' for printing
sce.hold <- sce.dlpfc

assay(sce.dlpfc, "logcounts") <- NULL
sizeFactors(sce.dlpfc) <- NULL
sce.dlpfc <- logNormCounts(sce.dlpfc)

### First make a list of Boolean param / cell subtype ===
# Will use this to assess more 'valid', non-noise-driving markers
cellSubtype.idx <- splitit(sce.dlpfc$cellType)
medianNon0 <- lapply(cellSubtype.idx, function(x){
  apply(as.matrix(assay(sce.dlpfc, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

sapply(medianNon0, table)
#       Astro Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Inhib_A Inhib_B Inhib_C Inhib_D Inhib_E Inhib_F Micro Mural
# FALSE 28052   23023   24201   22110   23293   22397   23313   25080   25774   24528   23874   25636   24170 28439 28287
# TRUE   1258    6287    5109    7200    6017    6913    5997    4230    3536    4782    5436    3674    5140   871  1023
#       Oligo   OPC Tcell
# FALSE 27985 27262 28735
# TRUE   1325  2048   575

## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.dlpfc), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`


# Run pairwise t-tests
markers.t.pw <- findMarkers(sce.dlpfc, groups=sce.dlpfc$cellType,
                                      assay.type="logcounts", design=mod, test="t",
                                      direction="up", pval.type="all", full.stats=T)

sapply(markers.t.pw, function(x){table(x$FDR<0.05)})
#       Astro Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Inhib_A Inhib_B Inhib_C Inhib_D Inhib_E Inhib_F Micro Mural
# FALSE 28997   29273   29233   29272   29204   29278   29254   29252   29304   29300   29237   29104   29152 28860 28819
# TRUE    313      37      77      38     106      32      56      58       6      10      73     206     158   450   491
#       Oligo   OPC Tcell
# FALSE 29085 29139 28882
# TRUE    225   171   428


## WMW: Blocking on donor (this test doesn't take 'design=' argument) ===
markers.wilcox.block <- findMarkers(sce.dlpfc, groups=sce.dlpfc$cellType,
                                          assay.type="logcounts", block=sce.dlpfc$donor, test="wilcox",
                                          direction="up", pval.type="all", full.stats=T)


sapply(markers.wilcox.block, function(x){table(x$FDR<0.05)})
      # Actually some decent results but many subclusters with 0 hits
# $Micro
# FALSE  TRUE 
# 29268    42 
# $Oligo
# FALSE  TRUE 
# 29224    86 
# $OPC
# FALSE  TRUE 
# 29293    17

## Binomial ===
markers.binom.block <- findMarkers(sce.dlpfc, groups=sce.dlpfc$cellType,
                                   assay.type="logcounts", block=sce.dlpfc$donor, test="binom",
                                   direction="up", pval.type="all", full.stats=T)

sapply(markers.binom.block, function(x){table(x$FDR<0.05)})
    # All FALSE

# ## Save all these for future reference
# save(markers.t.pw, markers.wilcox.block, #markers.binom.block,
#      file="rdas/revision/markers-stats_DLPFC-n3_findMarkers-SN-LEVEL_LAHMay2021.rda")


## Get top 40 pw genes
markerList.t <- lapply(markers.t.pw, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)

genes.top40.t <- lapply(markerList.t, function(x){head(x, n=40)})

# Print these to pngs
#dir.create("pdfs/revision/DLPFC/")
for(i in names(genes.top40.t)){
  png(paste0("pdfs/revision/DLPFC/DLPFC_t-sn-level_pairwise_top40markers-", i, "_logExprs_LAH2021.png"), height=1900, width=1200)
  print(
    plotExpressionCustom(sce = sce.dlpfc,
                         features = genes.top40.t[[i]], 
                         features_name = i, 
                         anno_name = "cellType") +
      scale_color_manual(values = cell_colors)
  )
  dev.off()
}


#### 1vALL test ####
markers.t.1vAll <- list()
for(i in levels(sce.dlpfc$cellType)){
  # Make temporary contrast
  sce.dlpfc$contrast <- ifelse(sce.dlpfc$cellType==i, 1, 0)
  # Test cluster vs. all
  markers.t.1vAll[[i]] <- findMarkers(sce.dlpfc, groups=sce.dlpfc$contrast,
                                            assay.type="logcounts", design=mod, test="t",
                                            direction="up", pval.type="all", full.stats=T)
}


    ## Then, temp set of stats to get the standardized logFC
    temp.1vAll <- list()
    for(i in levels(sce.dlpfc$cellType)){
      # Make temporary contrast
      sce.dlpfc$contrast <- ifelse(sce.dlpfc$cellType==i, 1, 0)
      # Test cluster vs. all
      temp.1vAll[[i]] <- findMarkers(sce.dlpfc, groups=sce.dlpfc$contrast,
                                     assay.type="logcounts", design=mod, test="t",
                                     std.lfc=TRUE,
                                     direction="up", pval.type="all", full.stats=T)
    }



    ## As with DLPFC, for some reason all the results are in the
     #    second List entry (first is always empty)

head(markers.t.1vAll[["Oligo"]][[2]])
    ## Nice, MBP and PLP1 are again in the top 6

markers.t.1vAll.db <- findMarkers_1vAll(sce.dlpfc, assay_name = "logcounts")

sapply(markers.t.1vAll, function(x){
  table(x[[2]]$stats.0$log.FDR < log10(.001))
})
#       Astro Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Inhib_A Inhib_B Inhib_C Inhib_D Inhib_E Inhib_F Micro Mural
# FALSE 24313   18625   18603   19054   22991   22044   21319   21276   22932   21041   20449   28245   28333 25473 27693
# TRUE   4997   10685   10707   10256    6319    7266    7991    8034    6378    8269    8861    1065     977  3837  1617
#       Oligo   OPC Tcell
# FALSE 26128 25312 27579
# TRUE   3182  3998  1731

# Replace that empty slot with the entry with the actul stats
markers.t.1vAll <- lapply(markers.t.1vAll, function(x){ x[[2]] })
# Same for that with std.lfc
temp.1vAll <- lapply(temp.1vAll, function(x){ x[[2]] })

# Now just pull from the 'stats.0' DataFrame column
markers.t.1vAll <- lapply(markers.t.1vAll, function(x){ x$stats.0 })
temp.1vAll <- lapply(temp.1vAll, function(x){ x$stats.0 })

# Re-name std.lfc column and add to the first result
for(i in names(temp.1vAll)){
  colnames(temp.1vAll[[i]])[1] <- "std.logFC"
  markers.t.1vAll[[i]] <- cbind(markers.t.1vAll[[i]], temp.1vAll[[i]]$std.logFC)
  # Oh the colname is kept weird
  colnames(markers.t.1vAll[[i]])[4] <- "std.logFC"
  # Then re-organize
  markers.t.1vAll[[i]] <- markers.t.1vAll[[i]][ ,c("logFC","std.logFC","log.p.value","log.FDR")]
}





## Let's save this along with the previous pairwise results
save(markers.t.1vAll, markers.t.1vAll.db, markers.t.pw, markers.wilcox.block,
     file="rdas/revision/markers-stats_DLPFC-n3_findMarkers-SN-LEVEL_LAHMay2021.rda")


## Print these to pngs
markerList.t.1vAll <- lapply(markers.t.1vAll, function(x){
  rownames(x)[x$log.FDR < log10(0.000001)]
 }
)
genes.top40.t <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(here("pdfs/revision/DLPFC", paste0("DLPFC_t-sn-level_1vALL_top40markers-",i,"_logExprs_Apr2020.png")), height=1900, width=1200)
  print(
      plotExpressionCustom(sce = sce.dlpfc,
                           features = genes.top40.t[[i]], 
                           features_name = paste0(i, " top 40 markers: single-nucleus-level p.w. t-tests, cluster-vs-all"), 
                           anno_name = "cellType") +
        scale_color_manual(values = cell_colors)
  )
  dev.off()
}


## How do they intersect?
markerList.t.pw <- lapply(markers.t.pw, function(x){
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
# Astro Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Inhib_A Inhib_B Inhib_C Inhib_D Inhib_E Inhib_F   Micro   Mural 
# 313      37      77      38     106      32      56      58       6      10      73     206     158     450     491 
# Oligo     OPC   Tcell 
# 225     171     428 
    # Of top 40's:
    sapply(names(markerList.t.pw), function(c){
      length(intersect(lapply(markerList.t.pw, function(l){head(l,n=40)})[[c]],
                       lapply(markerList.t.1vAll, function(l){head(l,n=40)})[[c]]
                       ))
    })
    # Astro Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Inhib_A Inhib_B Inhib_C Inhib_D Inhib_E Inhib_F   Micro   Mural 
    # 32      18      26      18      28      18      22      20       3       5      23      33      30      30      38 
    # Oligo     OPC   Tcell 
    # 31      34      38 


    
## Write these top 40 lists to a csv
names(markerList.t.pw) <- paste0(names(markerList.t.pw),"_pw")
names(markerList.t.1vAll) <- paste0(names(markerList.t.1vAll),"_1vAll")


## Add empty string
pad <- rep("",40 - min(sapply(markerList.t.pw, length)))
markerList.t.pw <- sapply(markerList.t.pw, function(x) c(x, pad))

top40genes <- cbind(sapply(markerList.t.pw, function(x) head(x, n=40)),
                    sapply(markerList.t.1vAll, function(y) head(y, n=40)))
top40genes <- top40genes[ ,sort(colnames(top40genes))]

## fix this
write.csv(top40genes, file=here("tables/revision/top40genesLists_DLPFC-n3_cellType_SN-LEVEL-tests_May2020.csv"),
          row.names=FALSE)


### MNT add 09Jul2021 =========
  # Another way ('cluster-vs-all-others' method used in other regions):

## Post-hoc: 'prelimCluster' 101 are T cells; 90 look like macrophages;
 #           (from interactively exploration) -> Edit/make a copy of this SCE for MNT work
load("rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda", verbose=T)

table(droplevels(sce.dlpfc$prelimCluster[sce.dlpfc$cellType == "Tcell"]))
    #  90 101 
    #  10   9

# First convert to 'character' class
sce.dlpfc$cellType <- as.character(sce.dlpfc$cellType)
sce.dlpfc$cellType[sce.dlpfc$prelimCluster == "90"] <- "Macrophage"
# Re-factor
sce.dlpfc$cellType <- factor(sce.dlpfc$cellType)

# Add new color
cell_colors["Macrophage"] <- setdiff(tableau20, cell_colors)[1]

# For reference
    annotationTab.dlpfc$cellType[annotationTab.dlpfc$cellType=="Tcell"] <- "Tcell_Macrophage"
    clusterRefTab.dlpfc$annot.MNT <- annotationTab.dlpfc$cellType[match(clusterRefTab.dlpfc$merged,
                                                                        annotationTab.dlpfc$collapsedCluster)]
    clusterRefTab.dlpfc$annot.MNT[clusterRefTab.dlpfc$origClust=="90"] <- "Macrophage"
    clusterRefTab.dlpfc$annot.MNT[clusterRefTab.dlpfc$origClust=="101"] <- "Tcell"
    
# Check
plotTSNE(sce.dlpfc, colour_by="cellType", point_alpha=0.5, text_by="cellType") +
  scale_color_manual(values = cell_colors) + labs(colour="Cell type")

save(sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc, ref.sampleInfo, annotationTab.dlpfc, cell_colors, 
     file="rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_MNT2021.rda")

## (Filter all-0 genes; set up `logNormCounts()`, as above)

## Re-create list of Boolean param / cell subtype (will append/save this info):
cellSubtype.idx <- splitit(sce.dlpfc$cellType)
medianNon0.dlpfc <- lapply(cellSubtype.idx, function(x){
  apply(as.matrix(assay(sce.dlpfc, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

sapply(medianNon0.dlpfc, table) # see above
    #      Astro Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Inhib_A Inhib_B Inhib_C
    # FALSE 28052   23023   24201   22110   23293   22397   23313   25080   25774   24528
    # TRUE   1258    6287    5109    7200    6017    6913    5997    4230    3536    4782
    #       Inhib_D Inhib_E Inhib_F   Macrophage Micro Mural Oligo   OPC Tcell
    # FALSE   23874   25636   24170        28283 28439 28287 27985 27262 28655
    # TRUE     5436    3674    5140         1027   871  1023  1325  2048   655
    #     - now we can see that 'Tcell's have more non-0-median-expressing genes

    # Confirm with some consistent T cell markers seen in other regionos:
    c("SKAP1","ITK","CD247") %in% names(medianNon0.dlpfc[["Tcell"]][medianNon0.dlpfc[["Tcell"]]==T])
    # and similarly (to NAc's 'Macrophage')
    c("CD163","MRC1","SIGLEC1") %in% names(medianNon0.dlpfc[["Macrophage"]][medianNon0.dlpfc[["Macrophage"]]==T])

    
mod <- with(colData(sce.dlpfc), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

markers.dlpfc.t.1vAll <- list()
for(i in levels(sce.dlpfc$cellType)){
  # Make temporary contrast
  sce.dlpfc$contrast <- ifelse(sce.dlpfc$cellType==i, 1, 0)
  # Test cluster vs. all others
  markers.dlpfc.t.1vAll[[i]] <- findMarkers(sce.dlpfc, groups=sce.dlpfc$contrast,
                                          assay.type="logcounts", design=mod, test="t",
                                          std.lfc=TRUE,
                                          direction="up", pval.type="all", full.stats=T)
}
    ## Since all other stats are the same, and don't really use the non-standardized
    #    logFC, just generate one object, unlike before

class(markers.dlpfc.t.1vAll[["Oligo"]])
    # a SimpleList of length 2, named "0" and "1" (from the temporary 'contrast')
    # -> we want the second entry, named "1"
        #    (for other purposes, might be interesting to look into that "0" entry, which
        #     is basically what genes are depleted in the cell type of interest)


# Do some reorganizing
markers.dlpfc.t.1vAll <- lapply(markers.dlpfc.t.1vAll, function(x){
  # Basically take the 'stats.[1 or 0]' since is redundant with the 'summary'-level stats
  lapply(x, function(y){ y[ ,4] }) 
})

# Re-name std.lfc column and the entries; add non-0-median info
for(i in names(markers.dlpfc.t.1vAll)){
  colnames(markers.dlpfc.t.1vAll[[i]][["0"]])[1] <- "std.logFC"
  colnames(markers.dlpfc.t.1vAll[[i]][["1"]])[1] <- "std.logFC"
  # Add non0median Boolean - might be informative for both sets of stats
  markers.dlpfc.t.1vAll[[i]][["0"]] <- cbind(markers.dlpfc.t.1vAll[[i]][["0"]],
                                           medianNon0.dlpfc[[i]][match(rownames(markers.dlpfc.t.1vAll[[i]][["0"]]),
                                                                     names(medianNon0.dlpfc[[i]]))])
  colnames(markers.dlpfc.t.1vAll[[i]][["0"]])[4] <- "non0median"
  
  # "1" aka 'enriched'
  markers.dlpfc.t.1vAll[[i]][["1"]] <- cbind(markers.dlpfc.t.1vAll[[i]][["1"]],
                                           medianNon0.dlpfc[[i]][match(rownames(markers.dlpfc.t.1vAll[[i]][["1"]]),
                                                                     names(medianNon0.dlpfc[[i]]))])
  colnames(markers.dlpfc.t.1vAll[[i]][["1"]])[4] <- "non0median"
  
  # Then re-name the entries to more interpretable, because we'll keeping both contrasts
  names(markers.dlpfc.t.1vAll[[i]]) <- paste0(i,c("_depleted", "_enriched"))
}

## Some interactive exploration of Inhib_E / Inhib_F ===
    # More believable markers numbers
    markerList.t.1vAll <- lapply(markers.dlpfc.t.1vAll, function(x){
      rownames(x[[2]])[ x[[2]]$log.FDR < log(0.05) & x[[2]]$non0median==TRUE ]
      }
    )
    lengths(markerList.t.1vAll)
        #   Astro Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Inhib_A Inhib_B Inhib_C 
        #     769    4500    3534    4868    3241    4122    3720    2631    2104    3085 
        # Inhib_D Inhib_E Inhib_F   Micro   Mural   Oligo     OPC   Tcell 
        #    3692     552     491     649     305     903    1129     250
            # With splitting the 'Tcell' into 'Tcell' & 'Macrophage', only difference:
            # Tcell   Macrophage
            #   260          429

    # Macrophage
    plotExpressionCustom(sce.dlpfc, anno_name="cellType",features_name="Check: Macrophage",
                         features=head(markerList.t.1vAll[["Macrophage"]])) +
      scale_color_manual(values = cell_colors)
    
    # Tcell
    plotExpressionCustom(sce.dlpfc, anno_name="cellType",features_name="Check: Tcell",
                         features=head(markerList.t.1vAll[["Tcell"]])) +
      scale_color_manual(values = cell_colors)

        # Save this into a separate iteration of .rda
        save(markers.dlpfc.t.1vAll, medianNon0.dlpfc,
             file="rdas/revision/markers-stats_DLPFC-n3_findMarkers-SN-LEVEL_MNT_v2_2021.rda")
        

## Load previous results for reference
load("rdas/revision/markers-stats_DLPFC-n3_findMarkers-SN-LEVEL_LAHMay2021.rda", verbose=T)
    # markers.t.1vAll, markers.t.1vAll.db, markers.t.pw, markers.wilcox.block

    # ** Another observation: These are interesting
    table(rownames(markers.t.1vAll[["Inhib_C"]]) ==
            rownames(markers.dlpfc.t.1vAll[["Inhib_C"]][["Inhib_C_enriched"]]))
        # 97 FALSE (and this varies on the cell class tested)


# Save back into a 'duplicate'/MNT copy, with the new objects
save(markers.t.pw, markers.wilcox.block,
     markers.dlpfc.t.1vAll, medianNon0.dlpfc,
     file="rdas/revision/markers-stats_DLPFC-n3_findMarkers-SN-LEVEL_MNT2021.rda")



### Session info for 03Jun2021 ============
sessionInfo()

# R version 4.1.0 Patched (2021-05-18 r80330)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS:   /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/lib/libRblas.so
# LAPACK: /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/lib/libRlapack.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices datasets  utils     methods   base     
# 
# other attached packages:
#   [1] here_1.0.1                  DeconvoBuddies_0.99.0       limma_3.48.0                jaffelab_0.99.31           
# [5] rafalib_1.0.0               DropletUtils_1.12.0         batchelor_1.8.0             scran_1.20.1               
# [9] scater_1.20.0               ggplot2_3.3.3               scuttle_1.2.0               EnsDb.Hsapiens.v86_2.99.0  
# [13] ensembldb_2.16.0            AnnotationFilter_1.16.0     GenomicFeatures_1.44.0      AnnotationDbi_1.54.0       
# [17] SingleCellExperiment_1.14.1 SummarizedExperiment_1.22.0 Biobase_2.52.0              GenomicRanges_1.44.0       
# [21] GenomeInfoDb_1.28.0         IRanges_2.26.0              S4Vectors_0.30.0            BiocGenerics_0.38.0        
# [25] MatrixGenerics_1.4.0        matrixStats_0.59.0          colorout_1.2-2             
# 
# loaded via a namespace (and not attached):
#   [1] googledrive_1.0.1         ggbeeswarm_0.6.0          colorspace_2.0-1          rjson_0.2.20             
# [5] ellipsis_0.3.2            rprojroot_2.0.2           bluster_1.2.1             XVector_0.32.0           
# [9] BiocNeighbors_1.10.0      rstudioapi_0.13           bit64_4.0.5               fansi_0.5.0              
# [13] splines_4.1.0             R.methodsS3_1.8.1         sparseMatrixStats_1.4.0   cachem_1.0.5             
# [17] Rsamtools_2.8.0           ResidualMatrix_1.2.0      cluster_2.1.2             dbplyr_2.1.1             
# [21] R.oo_1.24.0               png_0.1-7                 HDF5Array_1.20.0          compiler_4.1.0           
# [25] httr_1.4.2                dqrng_0.3.0               assertthat_0.2.1          Matrix_1.3-4             
# [29] fastmap_1.1.0             lazyeval_0.2.2            BiocSingular_1.8.0        prettyunits_1.1.1        
# [33] tools_4.1.0               rsvd_1.0.5                igraph_1.2.6              gtable_0.3.0             
# [37] glue_1.4.2                GenomeInfoDbData_1.2.6    dplyr_1.0.6               rappdirs_0.3.3           
# [41] Rcpp_1.0.6                vctrs_0.3.8               Biostrings_2.60.0         rhdf5filters_1.4.0       
# [45] rtracklayer_1.52.0        DelayedMatrixStats_1.14.0 stringr_1.4.0             beachmat_2.8.0           
# [49] lifecycle_1.0.0           irlba_2.3.3               restfulr_0.0.13           statmod_1.4.36           
# [53] XML_3.99-0.6              edgeR_3.34.0              zlibbioc_1.38.0           scales_1.1.1             
# [57] hms_1.1.0                 ProtGenerics_1.24.0       rhdf5_2.36.0              RColorBrewer_1.1-2       
# [61] yaml_2.2.1                curl_4.3.1                memoise_2.0.0             gridExtra_2.3            
# [65] segmented_1.3-4           biomaRt_2.48.0            stringi_1.6.2             RSQLite_2.2.7            
# [69] BiocIO_1.2.0              ScaledMatrix_1.0.0        filelock_1.0.2            BiocParallel_1.26.0      
# [73] rlang_0.4.11              pkgconfig_2.0.3           bitops_1.0-7              lattice_0.20-44          
# [77] Rhdf5lib_1.14.0           purrr_0.3.4               GenomicAlignments_1.28.0  bit_4.0.4                
# [81] tidyselect_1.1.1          magrittr_2.0.1            R6_2.5.0                  generics_0.1.0           
# [85] metapod_1.0.0             DelayedArray_0.18.0       DBI_1.1.1                 pillar_1.6.1             
# [89] withr_2.4.2               KEGGREST_1.32.0           RCurl_1.98-1.3            tibble_3.1.2             
# [93] crayon_1.4.1              utf8_1.2.1                BiocFileCache_2.0.0       viridis_0.6.1            
# [97] progress_1.2.2            locfit_1.5-9.4            grid_4.1.0                blob_1.2.1               
# [101] digest_0.6.27             R.utils_2.10.1            munsell_0.5.0             beeswarm_0.3.1           
# [105] viridisLite_0.4.0         vipor_0.4.5    
