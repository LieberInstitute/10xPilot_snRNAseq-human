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

# non-0-median - run them all
#amy_astro_B <- which(sce.astro$cellType == "amy_Astro_B")
astro.idx <- splitit(sce.astro$cellType)
medianNon0.astro <- lapply(astro.idx, function(x){
  apply(as.matrix(assay(sce.astro, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})
# non0median.amy.As_B <- apply(as.matrix(assay(sce.astro, "logcounts")), 1, function(y){
#   median(y[amy_astro_B]) > 0
# })
#table(non0median.amy.As_B)  # 146

sapply(medianNon0.astro, table)
    #       amy_Astro_A amy_Astro_B dlpfc_Astro hpc_Astro_A hpc_Astro_B nac_Astro_A
    # FALSE       31172       33392       32280       31669       32205       32155
    # TRUE         2366        *146        1258        1869        1333        1383
    #       nac_Astro_B sacc_Astro_A sacc_Astro_B
    # FALSE       31122        31761        32827
    # TRUE         2416         1777          711

    sapply(astro.idx, function(x){quantile(sce.astro$sum[x])})
        #      amy_Astro_A amy_Astro_B dlpfc_Astro hpc_Astro_A hpc_Astro_B nac_Astro_A
        # 0%        1182.0       188.0      884.00      1940.0        1127       798.0
        # 25%       7471.5       850.5     4008.75      5631.5        3609      4049.5
        # 50%      11275.0     *1156.0     5737.00      7996.5        5767      6524.0
        # 75%      16469.5      1842.5     7953.00     11802.5        8420      9175.5
        # 100%     35728.0      7739.0    26618.00     30085.0       20088     17601.0
        #      nac_Astro_B sacc_Astro_A sacc_Astro_B
        # 0%        635.00        174.0       102.00
        # 25%      7457.75       5781.0      1977.25
        # 50%     10811.50       7575.0      3872.50
        # 75%     15545.00      10077.5      6176.25
        # 100%    37265.00      23974.0     14206.00


# Add respective 'non0median' column to the stats for each set of markers
for(i in names(markers.astro.t.pw)){
  markers.astro.t.pw[[i]] <- cbind(markers.astro.t.pw[[i]],
                                 medianNon0.astro[[i]][match(rownames(markers.astro.t.pw[[i]]),
                                                           names(medianNon0.astro[[i]]))])
  colnames(markers.astro.t.pw[[i]])[12] <- "non0median"
}

sapply(markers.astro.t.pw, function(x){table(x$FDR<0.05 & x$non0median == TRUE)["TRUE"]})
    # amy_Astro_A.TRUE  amy_Astro_B.TRUE  dlpfc_Astro.TRUE  hpc_Astro_A.TRUE 
    #               14                 4                 3                 3 
    # hpc_Astro_B.TRUE  nac_Astro_A.TRUE  nac_Astro_B.TRUE sacc_Astro_A.TRUE 
    #               72                11               202                75 
    #  sacc_Astro_B.NA 
    #               NA 

markerList.astro <- lapply(markers.astro.t.pw, function(x){
  rownames(x)[x$FDR < 0.05 & x$non0median==TRUE]
})

markerList.astro[["amy_Astro_B"]]
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


## Make Astro class marker array from pw tests like with the NAc interneurons ===
markers.astro.t.pw.full <- markers.astro.t.pw
markers.astro.t.pw[["sacc_Astro_B"]] <- NULL
topToPrint <- as.data.frame(sapply(markers.astro.t.pw, function(x) {
  head(rownames(x)[x$FDR < 0.05 & x$non0median==TRUE], n=3)}))

table(unlist(topToPrint) %in% rownames(sce.astro)) # good

topToPrint

# Print
pdf("pdfs/revision/pubFigures/suppFig_across-regions_astros-marker-array_MNT2021.pdf", height=4, width=10)
print(
  plotExpressionCustom(sce.astro.hold, features=c(t(topToPrint)), features_name="",
                       anno_name="cellType", point_alpha=0.3, point_size=0.7, ncol=8, scales="free_y") +
    ggtitle(label="amy_Astro_A   amy_Astro_B   dlpfc_Astro      hpc_Astro_A    hpc_Astro_B    nac_Astro_A    nac_Astro_B    sacc_Astro_A") + xlab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
          axis.title.y = element_text(angle = 90, size = 12),
          axis.text.y = element_text(size = 9),
          plot.title = element_text(size = 12),
          panel.grid.major=element_line(colour="grey95", size=0.8),
          panel.grid.minor=element_line(colour="grey95", size=0.4))
)
dev.off()


### MNT revision add ========
  #   -> What are depleted in AMY's Astro_B, compared to the rest?

    # Same mod
    mod <- with(colData(sce.astro), model.matrix(~ donor))
    mod <- mod[ , -1, drop=F]
    
    markers.astro.t.1vAll <- list()
    for(i in levels(sce.astro$cellType)){
      # Make temporary contrast
      sce.astro$contrast <- ifelse(sce.astro$cellType==i, 1, 0)
      # Test cluster vs. all others
      markers.astro.t.1vAll[[i]] <- findMarkers(sce.astro, groups=sce.astro$contrast,
                                              assay.type="logcounts", design=mod, test="t",
                                              std.lfc=TRUE,
                                              direction="up", pval.type="all", full.stats=T)
    }
    
    
    # Re-organize like with other region-specific markers stats:
    markers.astro.t.1vAll <- lapply(markers.astro.t.1vAll, function(x){
      # Basically take the 'stats.[1 or 0]' since is redundant with the 'summary'-level stats
      lapply(x, function(y){ y[ ,4] }) 
    })
    
    #Re-name std.lfc column and the entries; add non-0-median info
    for(i in names(markers.astro.t.1vAll)){
      colnames(markers.astro.t.1vAll[[i]][["0"]])[1] <- "std.logFC"
      colnames(markers.astro.t.1vAll[[i]][["1"]])[1] <- "std.logFC"
      # Add non0median Boolean - might be informative for both sets of stats
      markers.astro.t.1vAll[[i]][["0"]] <- cbind(markers.astro.t.1vAll[[i]][["0"]],
                                               medianNon0.astro[[i]][match(rownames(markers.astro.t.1vAll[[i]][["0"]]),
                                                                         names(medianNon0.astro[[i]]))])
      colnames(markers.astro.t.1vAll[[i]][["0"]])[4] <- "non0median"
      
      # "1" aka 'enriched'
      markers.astro.t.1vAll[[i]][["1"]] <- cbind(markers.astro.t.1vAll[[i]][["1"]],
                                               medianNon0.astro[[i]][match(rownames(markers.astro.t.1vAll[[i]][["1"]]),
                                                                         names(medianNon0.astro[[i]]))])
      colnames(markers.astro.t.1vAll[[i]][["1"]])[4] <- "non0median"
      
      # Then re-name the entries to more interpretable, because we'll keeping both contrasts
      names(markers.astro.t.1vAll[[i]]) <- paste0(i,c("_depleted", "_enriched"))
    }
    
    
    ## Let's save this
    save(markers.astro.t.1vAll, medianNon0.astro,
         file="rdas/revision/markers-stats_acrossRegions-astros_findMarkers_MNT2021.rda")
    
    # Explore some expression
    plotExpressionCustom(sce.astro.hold, anno_name="cellType", features=head(rownames(markers.astro.t.1vAll[["amy_Astro_B"]][["amy_Astro_B_depleted"]]),n=12),
                         features_name="Downregulated in AMY Astro_B", ncol=4, scales="free_y") +
      xlab("") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
            axis.title.y = element_text(angle = 90, size = 12),
            axis.text.y = element_text(size = 10),
            plot.title = element_text(size = 12),
            panel.grid.major=element_line(colour="grey95", size=0.8),
            panel.grid.minor=element_line(colour="grey95", size=0.4))

    sapply(markers.astro.t.1vAll, function(x){table(head(x[[1]]$non0median,n=40))})
        #       amy_Astro_A amy_Astro_B dlpfc_Astro hpc_Astro_A hpc_Astro_B nac_Astro_A
        # FALSE          26          22          21          34          15           8
        # TRUE           14          18          19           6          25          32
        #       nac_Astro_B sacc_Astro_A sacc_Astro_B
        # FALSE          24           33            6
        # TRUE           16            7           34
            # So each population has its own set of genes specifically depleted

    # For AMY Astro_B, those are (median expression == 0):
    printThese <- head(rownames(markers.astro.t.1vAll[["amy_Astro_B"]][["amy_Astro_B_depleted"]])[
                      markers.astro.t.1vAll[["amy_Astro_B"]][["amy_Astro_B_depleted"]]$non0median==FALSE
                    ],n=22)

    head(which(markers.astro.t.1vAll[["amy_Astro_B"]][["amy_Astro_B_depleted"]]$non0median==FALSE),n=22)
        # PTPRZ1     GABRB1  LINC00461      NHSL1       TCF4      PDE7B       FUT9 
        #      7         12         16         19         22         23         24 
        #  SASH1       AKT3      MAST4       NFIB    ATP13A4        APC AC091826.2 
        #     25         27         28         29         30         31         32 
        #   SOX6    CARMIL1      PREX2    DENND1A       MBD5       GNAQ     ZFAND3 
        #     33         34         35         36         37         38         39 
        # STXBP5 
        #     40 
    
    plotExpressionCustom(sce.astro.hold, anno_name="cellType", features=printThese,
                         features_name="Downregulated in AMY Astro_B", ncol=4, scales="free_y") +
      xlab("") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
            axis.title.y = element_text(angle = 90, size = 12),
            axis.text.y = element_text(size = 10),
            plot.title = element_text(size = 12),
            panel.grid.major=element_line(colour="grey95", size=0.8),
            panel.grid.minor=element_line(colour="grey95", size=0.4))
    
    
    # end revision exploration ==========


# Plot the $sum densities
coldat <- as.data.frame(colData(sce.astro))
coldat$log10.sum <- log10(coldat$sum)

pdf("pdfs/revision/pubFigures/suppFig_across-regions_astros-totalUMIs-density_MNT2021.pdf", height=3.5, width=7.5)
ggplot(coldat, aes(x=log10.sum, color=cellType, fill=cellType)) +
  geom_density(alpha=0.15,size=1.2) +
  scale_color_manual(values=tableau10medium[1:9],
                     labels=paste0(levels(sce.astro$cellType)," (",table(sce.astro$cellType),")")) +
  labs(colour="Cell type") +
  scale_fill_manual(values=tableau10medium[1:9]) + guides(fill=FALSE) +
  xlab("log10(total.n.UMIs)") + ylab("Density") +
  ggtitle("Distribution of total UMIs captured per astrocyte cell class") +
  theme(axis.title.x = element_text(size = 14),
        axis.text.x = element_text(hjust = 1, size = 11),
        axis.title.y = element_text(angle = 90, size = 14),
        axis.text.y = element_text(size = 11),
        plot.title = element_text(size = 14))
dev.off()


## Do the same with 'Micro' to contrast to Astro:
sce.micro <- sce.allRegions[ ,grep("Micro", sce.allRegions$cellType)]
sce.micro$cellType <- droplevels(sce.micro$cellType)
table(sce.micro$cellType, sce.micro$donor)
    #                   br5161 br5207 br5212 br5276 br5287 br5400 br5701
    # amy_Micro            411      0    304     14      0    117    355
    # dlpfc_Micro          152    144     92      0      0      0      0
    # hpc_Micro            487      0    481      0    193      0      0
    # nac_Micro             66      0     59     33     34    222     15
    # nac_Micro_resting      3      0     33     22      0      5      0
    # sacc_Micro           232      0    243      3      0    292     14

# Broad distribution:
micro.idx <- splitit(sce.micro$cellType)
sapply(micro.idx, function(x){quantile(sce.micro$sum[x])})
    #      amy_Micro dlpfc_Micro hpc_Micro nac_Micro nac_Micro_resting sacc_Micro
    # 0%         435      879.00       454      1504             105.0     113.00
    # 25%       3952     3019.25      3580      4268             469.5    3198.50
    # 50%       5381     3883.50      4661      5545             866.0    4378.00
    # 75%       7055     4911.75      5966      6789            1327.5    5458.25
    # 100%     20500    11137.00     12463     12537            5239.0   11896.00

# Plot the $sum densities
coldat <- as.data.frame(colData(sce.micro))
coldat$log10.sum <- log10(coldat$sum)

pdf("pdfs/revision/pubFigures/suppFig_across-regions_micro-totalUMIs-density_MNT2021.pdf", height=3.5, width=7.5)
ggplot(coldat, aes(x=log10.sum, color=cellType, fill=cellType)) +
  geom_density(alpha=0.15,size=1.2) +
  scale_color_manual(values=tableau20[15:20],
                     labels=paste0(levels(sce.micro$cellType)," (",table(sce.micro$cellType),")")) +
  labs(colour="Cell type") +
  scale_fill_manual(values=tableau20[15:20]) + guides(fill=FALSE) +
  xlab("log10(total.n.UMIs)") + ylab("Density") +
  ggtitle("Distribution of total UMIs captured per microglia cell class") +
  theme(axis.title.x = element_text(size = 14),
        axis.text.x = element_text(hjust = 1, size = 11),
        axis.title.y = element_text(angle = 90, size = 14),
        axis.text.y = element_text(size = 11),
        plot.title = element_text(size = 14))
dev.off()






## AMY 'Inhib_B' vs DLPFC 'Inhib_A' - r=0.86 ===
sce.dlpfc <- sce.allRegions[ ,sce.allRegions$region=="dlpfc"]
sce.dlpfc$cellType <- droplevels(sce.dlpfc$cellType)

# Plot some AMY 'Inhib_B' markers (seen in Louise's DLPFC 'Inhib_A' lists)
plotExpressionCustom(sce.dlpfc, anno_name="cellType", features_name="some AMY 'Inhib_B'",
                     features=c("VIP", "CALB2", "CRH", "PTHLH"), ncol=2) +
  


### NAc 'excit' MSNs vs others?? =============================
  # 'MSN.D1_A', 'MSN.D1_D', 'MSN.D2_A', 'MSN.D2_B' vs the rest
load("rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda", verbose=T)

sce.test <- sce.nac[ ,-grep("drop.", sce.nac$cellType)]
sce.test$cellType <- droplevels(sce.test$cellType)

sce.test$contrast <- 0
# Make the 'excit' ones 1
msn.excit.idx <- c(grep("MSN.D1_A", sce.test$cellType),
                   grep("MSN.D1_D",sce.test$cellType),
                   grep("MSN.D2_A",sce.test$cellType),
                   grep("MSN.D2_B",sce.test$cellType))
sce.test$contrast[msn.excit.idx] <- 1
# And the others -1
msn.rest.idx <- setdiff(grep("MSN", sce.test$cellType), msn.excit.idx)
sce.test$contrast[msn.rest.idx] <- -1
  
table(sce.test$cellType, sce.test$contrast) # good

# Do the PC components correlate with this contrast?
apply(reducedDim(sce.test, "PCA_opt"), 2, function(x){cor(x, sce.test$contrast)})
    # [1]  0.5451128939  0.3584196778 -0.1526565503  0.4616036003  0.0412478196 -0.5014262238
    # [7]  0.2344336760  0.4066390607 -0.5440848773 -0.0247296291  0.1251926241 -0.2562537482
    # [13]  0.2353769426 -0.2995765135  0.1100234036 -0.0645767070  0.3322698995  0.0861300085
    # [19]  0.1999916772 -0.1836870541  0.0628782010 -0.0735522902 -0.1502806994  0.1031656553
    # [ etc. ]
which(abs(apply(reducedDim(sce.test, "PCA_opt"), 2, function(x){cor(x, sce.test$contrast)})) >= 0.4)
    #[1] 1 4 6 8 9

sce.test$cellType.con <- as.character(sce.test$cellType)
sce.test$cellType.con[msn.excit.idx] <- "MSN.excit"
sce.test$cellType.con[msn.rest.idx] <- "MSN.inhib"
sce.test$cellType.con <- factor(sce.test$cellType.con)

coldat <- cbind(reducedDim(sce.test, "PCA_opt"), colData(sce.test))
colnames(coldat) <- gsub("V", "PC", colnames(coldat))
coldat <- as.data.frame(coldat)
coldat$cellType.MSNs <- as.character(coldat$cellType)
coldat$cellType.MSNs[grep("MSN", coldat$cellType.MSNs)] <- "MSN.broad"
coldat$cellType.MSNs <- factor(coldat$cellType.MSNs)

# Trick some cell class colors - use blue & red
cell_colors.nac.hold <- cell_colors.nac
cell_colors.nac <- cell_colors.nac[names(cell_colors.nac) %in% levels(sce.test$cellType.con)]
cell_colors.nac["MSN.excit"] <- cell_colors.nac.hold["MSN.D1_A"]
cell_colors.nac["MSN.inhib"] <- cell_colors.nac.hold["drop.doublet_B"]
cell_colors.nac["MSN.broad"] <- cell_colors.nac.hold["drop.doublet_C"]

## Plot top 10 ===
lay <- rbind(c(1,1),
             c(2,2))

coldat$MSNsize <- ifelse(coldat$cellType.MSNs=="MSN.broad",0.5,0.15)

pdf("pdfs/revision/regionSpecific_top10PCs_cellClass_MSNsGrouped-or-split_MNT2021.pdf")
for(i in 1:10){
  # All MSNs combined:
  grouped <- ggplot(coldat, aes_string(x=colnames(coldat)[i], color="cellType.MSNs", fill="cellType.MSNs")) +
                geom_density(alpha=0.15,size=0.8) +
                scale_color_manual(values=cell_colors.nac) +
                labs(colour="Cell type") +
                scale_fill_manual(values=cell_colors.nac) + guides(fill=FALSE) +
                xlab(paste0("PC", i)) + ylab("Density") +
                ggtitle(paste0("Principal component ",i," by cell class; all MSNs grouped")) +
                theme(axis.title.x = element_text(size = 13),
                      axis.text.x = element_text(hjust = 1, size = 11),
                      axis.title.y = element_text(angle = 90, size = 14),
                      axis.text.y = element_text(size = 11),
                      plot.title = element_text(size = 12),
                      legend.text = element_text(size = 9),
                      legend.key.size = unit(0.4, "cm"))
  # Separated into 'Excit' & 'Inhib'
 separated <- ggplot(coldat, aes_string(x=colnames(coldat)[i], color="cellType.con", fill="cellType.con")) +
                geom_density(alpha=0.15,size=0.8) +
                scale_color_manual(values=cell_colors.nac) +
                labs(colour="Cell type") +
                scale_fill_manual(values=cell_colors.nac) + guides(fill=FALSE, size=FALSE) +
                xlab(paste0("PC", i)) + ylab("Density") +
                ggtitle(paste0("Principal component ",i," by cell class; MSNs separated by signature")) +
                theme(axis.title.x = element_text(size = 13),
                      axis.text.x = element_text(hjust = 1, size = 11),
                      axis.title.y = element_text(angle = 90, size = 14),
                      axis.text.y = element_text(size = 11),
                      plot.title = element_text(size = 12),
                      legend.text = element_text(size = 9),
                      legend.key.size = unit(0.4, "cm"))
 # Plot both per page:
 grid.arrange(grobs=list(grouped,
                         separated),
              layout_matrix=lay)
}
dev.off()


## Violin plot iteration =======
pdf("pdfs/revision/regionSpecific_top10PCs_cellClass_MSNsGrouped-or-split_violins_MNT2021.pdf")
for(i in 1:10){
grouped <- ggplot(coldat, aes_string(x="cellType.MSNs",y=colnames(coldat)[i],
                                     fill="cellType.MSNs",color="cellType.MSNs")) +
              geom_violin(alpha=0.4,size=1.2) +
              scale_color_manual(values=cell_colors.nac) +
              labs(colour="Cell class") +
              scale_fill_manual(values=cell_colors.nac) + guides(fill=FALSE) +
              ylab(paste0("PC", i)) + xlab("") +
              ggtitle(paste0("Principal component ",i," by cell class; all MSNs grouped")) +
              theme(axis.title.x = element_text(size = 13),
                    axis.text.x = element_text(angle=90, hjust = 1, size = 11),
                    axis.title.y = element_text(angle = 90, size = 14),
                    axis.text.y = element_text(size = 11),
                    plot.title = element_text(size = 12),
                    legend.text = element_text(size = 9),
                    legend.key.size = unit(0.4, "cm"))
# Separated into 'Excit' & 'Inhib'
separated <- ggplot(coldat, aes_string(x="cellType.MSNs",y=colnames(coldat)[i],
                                       fill="cellType.con",color="cellType.con")) +
                geom_violin(alpha=0.4,size=1.2) +
                scale_color_manual(values=cell_colors.nac) +
                labs(colour="Cell class") +
                scale_fill_manual(values=cell_colors.nac) + guides(fill=FALSE) +
                ylab(paste0("PC", i)) + xlab("") +
                ggtitle(paste0("Principal component ",i," by cell class; MSNs separated by signature")) +
                theme(axis.title.x = element_text(size = 13),
                      axis.text.x = element_text(angle=90, hjust = 1, size = 11),
                      axis.title.y = element_text(angle = 90, size = 14),
                      axis.text.y = element_text(size = 11),
                      plot.title = element_text(size = 12),
                      legend.text = element_text(size = 9),
                      legend.key.size = unit(0.4, "cm"))
# Plot both per page:
grid.arrange(grobs=list(grouped,
                        separated),
             layout_matrix=lay)
}
dev.off()


### 'Excitatory' MSNs [+ rest of excitatory 'branch'] vs the rest (or vs other MSNs at least) =========
## From above:
# Top 100 cluster genes space
cor_t_neu_cluster <- pheatmap(cor_t_defined.neu,
                              color=my.col.all,
                              breaks=theSeq.all,
                              fontsize_row=6.2, fontsize_col=6.2, cex=1,
                              clustering_distance_rows="euclidean",
                              clustering_distance_cols="euclidean",
                              clustering_method="complete")

rownames(cor_t_defined.neu[cor_t_neu_cluster$tree_row[["order"]],]) # This is the order we want
    # [1] "MSN.D1_D_nac"  "MSN.D2_B_nac"  "Inhib_E_amy"   "MSN.D1_A_nac"  "MSN.D2_A_nac" 
    # [6] "Excit_A_sacc"  "Excit_E_sacc"  "Excit_F_hpc"   "Excit_A_dlpfc" "Excit_D_dlpfc"
    # [11] "Excit_F_sacc"  "Excit_C_sacc"  "Excit_E_dlpfc" "Excit_D_sacc"  "Excit_C_amy"  
    # [16] "Excit_B_dlpfc" "Excit_C_dlpfc" "Excit_A_hpc"   "Excit_D_hpc"   "Excit_A_amy"  
    # [21] "Excit_B_hpc"   "Excit_F_dlpfc" "Excit_B_sacc"  "Inhib_C_amy"   "Inhib_B_hpc"  
    # [26] "Excit_C_hpc"   "Inhib_G_amy"   "Inhib_A_amy"   "Excit_E_hpc"  
excit.classes <- rownames(cor_t_defined.neu[cor_t_neu_cluster$tree_row[["order"]],])[1:29]
# Clean this up (needs to look like 'region_Class')
excit.classes <- paste0(ss(excit.classes,"_",3),
                        "_",
                        ss(excit.classes,"_",1),
                        "_",
                        ss(excit.classes,"_",2))

# Wanna test specifically this 'excitatory branch' vs the other MSNs
sub.classes <- c(excit.classes, "nac_MSN.D1_B","nac_MSN.D1_C","nac_MSN.D1_E","nac_MSN.D1_F",
                 "nac_MSN.D2_C","nac_MSN.D2_D")

# Load and subset those for some contrasts
load("rdas/revision/all-n24-samples_across-regions-analyses_forFigOnly_MNT2021.rda", verbose=T)

sce.sub <- sce.allRegions[ ,sce.allRegions$cellType %in% sub.classes]
sce.sub$cellType <- droplevels(sce.sub$cellType)
table(sce.sub$cellType)
# Re-create 'logcounts'
sizeFactors(sce.sub) <- NULL
assay(sce.sub, "logcounts") <- NULL
sce.sub <- logNormCounts(sce.sub)

# Make contrast
sce.sub$excit.branch <- ifelse(sce.sub$cellType %in% excit.classes, 1, 0)

# t-test: test 'excit' vs. rest of MSNs
#         since binarizing, can use the '1vAll' approach to keep both contrasts
mod <- with(colData(sce.sub), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

markers.excit.branch <- findMarkers(sce.sub, groups=sce.sub$excit.branch,
                                    assay.type="logcounts", design=mod, test="t",
                                    std.lfc=TRUE,
                                    direction="up", pval.type="all", full.stats=T)

# Only one contrast, so take that column
markers.excit.branch <- lapply(markers.excit.branch, function(x){ x[ ,4] })
names(markers.excit.branch) <- c("Excit_depleted", "Excit_enriched")

sapply(markers.excit.branch, function(x){table(x$log.FDR <= log(0.05))})
    #       Excit_depleted Excit_enriched
    # FALSE          28108          27182
    # TRUE            5430           6356


Readme <- "This test is comparing all 'excitatory branch' classes from the across-regions (including four MSN D1/D2 classes) analysis vs those remaining 'inhibitory' MSNs"
save(markers.excit.branch, excit.classes, Readme,
     file="rdas/revision/markers-stats_x-regions_excitBranch-vs-restOfMSNs_MNT2021.rda")

head(rownames(markers.excit.branch[["Excit_enriched"]]), n=40)
    # [1] "NRG3"       "RGS7"       "GRM5"       "ARAP2"      "RGS6"       "ZMAT4"     
    # [7] "AK5"        "CSMD1"      "ARNT2"      "UTRN"       "KCNMA1"     "LDLRAD4"   
    # [13] "GPC5"       "EPB41L2"    "LIMCH1"     "PRKCB"      "MAP3K5"     "COCH"      
    # [19] "AC073050.1" "SLC35F3"    "SLC8A1"     "TJP1"       "ADAMTS19"   "FRMD5"     
    # [25] "FMNL2"      "MAN1A1"     "EFNA5"      "HTR2C"      "XKR4"       "AC008574.1"
    # [31] "LDB2"       "KIRREL3"    "ADK"        "SGCZ"       "ERC2"       "PDZD2"     
    # [37] "HS6ST3"     "ABLIM1"     "PDE3B"      "DOCK4"    

head(rownames(markers.excit.branch[["Excit_depleted"]]), n=40)
    # [1] "ADARB2"     "ALK"        "LUZP2"      "CASZ1"      "KIRREL1"    "RXRG"      
    # [7] "KCNT2"      "SORCS2"     "CACNG5"     "BACH2"      "PMEPA1"     "EPS8"      
    # [13] "TACR1"      "SEMA5B"     "SEMA6D"     "DLX6-AS1"   "TSHZ1"      "GDNF-AS1"  
    # [19] "CRHR2"      "GABRG3"     "SLC35F1"    "KIT"        "FOXP2"      "AC068722.1"
    # [25] "TRHDE"      "PLCB1"      "ARHGAP18"   "PLCXD3"     "AL590867.1" "FHOD3"     
    # [31] "FAM155A"    "MTSS1"      "AC022126.1" "PDLIM5"     "CPNE4"      "TLL1"      
    # [37] "PCSK2"      "AL589740.1" "SPON1"      "LOXHD1"   




## What about just a within-NAc MSNs contrast? =====================
load("rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda", verbose=T)
    # sce.nac, chosen.hvgs.nac, pc.choice.nac, ref.sampleInfo, annotationTab.nac, cell_colors.nac

sce.test <- sce.nac[ ,grep("MSN.", sce.nac$cellType)]
sce.test$cellType <- droplevels(sce.test$cellType)
# For plotting:
sce.hold <- sce.test

# Split into 'MSN.excit' or 'MSN.inhib'
msn.excit.idx <- c(grep("MSN.D1_A", sce.test$cellType),
                   grep("MSN.D1_D",sce.test$cellType),
                   grep("MSN.D2_A",sce.test$cellType),
                   grep("MSN.D2_B",sce.test$cellType))
sce.test$cellType.sig <- "MSN.inhib"
sce.test$cellType.sig[msn.excit.idx] <- "MSN.excit"
table(sce.test$cellType, sce.test$cellType.sig)
    #           MSN.excit MSN.inhib
    # MSN.D1_A      3927         0
    # MSN.D1_B         0       239
    # MSN.D1_C         0       283
    # MSN.D1_D       718         0
    # MSN.D1_E         0       638
    # MSN.D1_F         0        86
    # MSN.D2_A      4262         0
    # MSN.D2_B       285         0
    # MSN.D2_C         0       314
    # MSN.D2_D         0        58

# Re-create 'logcounts'
sizeFactors(sce.test) <- NULL
assay(sce.test, "logcounts") <- NULL
sce.test <- logNormCounts(sce.test)

mod <- with(colData(sce.test), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

markers.excit.msns <- findMarkers(sce.test, groups=sce.test$cellType.sig,
                                  assay.type="logcounts", design=mod, test="t",
                                  std.lfc=TRUE,
                                  direction="up", pval.type="all", full.stats=T)

# Only one contrast, so take that column
markers.excit.msns <- lapply(markers.excit.msns, function(x){ x[ ,4] })
names(markers.excit.msns)
    #[1] "MSN.excit" "MSN.inhib"

Readme <- "This test is comparing the excitatory-signature (4) MSN classes vs those (6) more 'inhibitory' MSNs"
save(markers.excit.msns, Readme,
     file="rdas/revision/markers-stats_excitMSNs-vs-inhibMSNs_MNT2021.rda")
lapply(markers.excit.msns, function(x){table(x$log.FDR<log(0.05))})
    #       MSN.excit MSN.inhib
    # FALSE     26937     28396
    # TRUE       6601      5142

head(rownames(markers.excit.msns[["MSN.excit"]]), n=40)
    # [1] "HTR2C"      "SLC35F3"    "UTRN"       "RGS6"       "ARAP2"      "COCH"      
    # [7] "RGS7"       "SGCZ"       "PDZD2"      "ZMAT4"      "GPC5"       "LDLRAD4"   
    # [13] "MAP3K5"     "AC008574.1" "NRG3"       "EPB41L2"    "FOXO1"      "MAN1A1"    
    # [19] "PLPPR1"     "ADK"        "HTR4"       "ARNT2"      "ADAMTS19"   "AC073050.1"
    # [25] "SGK3"       "DCC"        "GRM5"       "AK5"        "TENM2"      "LRRC4C"    
    # [31] "PRKCB"      "PSD3"       "CSMD1"      "KCTD8"      "CDH12"      "DIAPH2"    
    # [37] "ANKS1B"     "MCTP1"      "SLC8A1"     "LINC01322"  

head(rownames(markers.excit.msns[["MSN.inhib"]]), n=40)
    # [1] "LUZP2"   "ADARB2"  "KIRREL1" "TRHDE"   "CPNE4"   "FAM155A" "ALK"     "CASZ1"  
    # [9] "SEMA5B"  "KCNT2"   "CNTN5"   "RXFP1"   "FHOD3"   "ILDR2"   "EPS8"    "ZNF804B"
    # [17] "DPP10"   "RXRG"    "CACNG5"  "CRHR2"   "SLIT1"   "GABRG3"  "SORCS2"  "BACH2"  
    # [25] "SPECC1"  "FMN1"    "PLCXD3"  "SLC35F1" "PLCB1"   "DIRAS2"  "XYLT1"   "CAMK2D" 
    # [33] "ADAM19"  "OPCML"   "TLL1"    "PPM1E"   "CACNA1C" "NRXN3"   "GRIN2A"  "KCNJ3"  


# Re-order and plot top 6 per grouping
sce.hold$cellType <- factor(sce.hold$cellType,
                            levels=c("MSN.D1_A","MSN.D1_D","MSN.D2_A","MSN.D2_B",
                                     "MSN.D1_B", "MSN.D1_C","MSN.D1_E",
                                     "MSN.D1_F","MSN.D2_C", "MSN.D2_D"))

pdf("pdfs/revision/pubFigures/suppFig_NAc_excit-vs-inhib-MSN-markers_MNT2021.pdf", height=4, width=6)
plotExpressionCustom(sce.hold, anno_name="cellType", features=c(head(rownames(markers.excit.msns[["MSN.excit"]]), n=6)),
                     features_name="Upregulated in 'excitatory'-signature MSN", ncol=3, scales="free_y") +
  scale_color_manual(values=cell_colors.nac) + xlab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.title.y = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 12),
        panel.grid.major=element_line(colour="grey95", size=0.8),
        panel.grid.minor=element_line(colour="grey95", size=0.4))

plotExpressionCustom(sce.hold, anno_name="cellType", features=c(head(rownames(markers.excit.msns[["MSN.inhib"]]), n=6)),
                     features_name="Upregulated in 'inhibitory'-signature MSN", ncol=3, scales="free_y") +
  scale_color_manual(values=cell_colors.nac) + xlab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.title.y = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 12),
        panel.grid.major=element_line(colour="grey95", size=0.8),
        panel.grid.minor=element_line(colour="grey95", size=0.4))
dev.off()


## Gene ontology of these broad 'excit' / 'inhib' MSN markers? ===========
#BiocManager::install("clusterProfiler")
library(clusterProfiler)

# (load markers, above)

# Make marker lists - be pretty strict here (take FDR < 0.0001)
markerList.msns <- lapply(markers.excit.msns, function(x){rownames(x)[x$log.FDR < log(0.0001)]})
lengths(markerList.msns)
    #MSN.excit MSN.inhib 
    #     4187      3174
length(intersect(markerList.msns[["MSN.excit"]], markerList.msns[["MSN.inhib"]])) # none

# Input list
inputList.msns <- list(MSN.excit=rowData(sce.nac)$gene_id[match(markerList.msns[["MSN.excit"]],
                                                           rownames(sce.nac))],
                  MSN.inhib=rowData(sce.nac)$gene_id[match(markerList.msns[["MSN.inhib"]],
                                                           rownames(sce.nac))]
                  )

# Convert and remove probes without conversion (Entrez) ID's
inputList.msns <- lapply(inputList.msns, function(x) bitr(x, fromType="ENSEMBL",
                                                toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = TRUE)[,2])
    # Warning messages:
    #   1: In bitr(x, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db",  :
    #                13.47% of input gene IDs are fail to map...
    #              2: In bitr(x, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db",  :
    #                           8.79% of input gene IDs are fail to map...
    # Resulting lengths
    # MSN.excit MSN.inhib 
    #      3654      2913 

# Set gene universe
geneMap <- rowRanges(sce.nac)
geneUniverse = bitr(geneMap$gene_id, fromType="ENSEMBL",
                    toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = TRUE)[,2]
    #Warning message:
    # In bitr(geneMap$gene_id, fromType = "ENSEMBL", toType = "ENTREZID",  :
    #           30.36% of input gene IDs are fail to map...


# GO & KEGG enrichment test:
goBP_msns <- compareCluster(inputList.msns, fun = "enrichGO",
                                 universe = geneUniverse, OrgDb = org.Hs.eg.db,
                                 ont = "BP", pvalueCutoff = 0.05, readable= TRUE)
goMF_msns <- compareCluster(inputList.msns, fun = "enrichGO",
                                 universe = geneUniverse, OrgDb = org.Hs.eg.db,
                                 ont = "MF", pvalueCutoff = 0.05, readable= TRUE)
goCC_msns <- compareCluster(inputList.msns, fun = "enrichGO",
                                 universe = geneUniverse, OrgDb = org.Hs.eg.db,
                                 ont = "CC", pvalueCutoff = 0.05, readable= TRUE)
KEGG_msns = compareCluster(inputList.msns, fun ='enrichKEGG', universe = geneUniverse,
                                organism = "hsa", pvalueCutoff=0.05)

save(goBP_msns, goMF_msns, goCC_msns, KEGG_msns,
     file="rdas/revision/clusterGOs_excit-vs-inhib-MSNs_fdr-0.0001_MNT2021.rda")

# Print some enriched pathways
pdf("pdfs/revision/pubFigures/suppFig_NAc_excit-vs-inhib-MSN_GO-KEGG_MNT2021.pdf", width=6)
dotplot(goBP_msns, font.size=9, title = "MSN 'excit' vs. 'inhib':\nBiological Process GO", showCategory=12)
dotplot(goMF_msns, font.size=9, title = "MSN 'excit' vs. 'inhib':\nMolecular Function GO", showCategory=12)
dotplot(goCC_msns, font.size=9, title = "MSN 'excit' vs. 'inhib':\nCellular Component GO", showCategory=12)
dotplot(KEGG_msns, font.size=9, title = "MSN 'excit' vs. 'inhib':\nKEGG Process", showCategory=12)
dev.off()




# Session info for 28Jul2021 ==============
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
#   [1] clusterProfiler_3.18.1      gridExtra_2.3               pheatmap_1.0.12            
# [4] RColorBrewer_1.1-2          lattice_0.20-41             limma_3.46.0               
# [7] jaffelab_0.99.30            rafalib_1.0.0               DropletUtils_1.10.3        
# [10] batchelor_1.6.3             scran_1.18.7                org.Hs.eg.db_3.12.0        
# [13] EnsDb.Hsapiens.v86_2.99.0   ensembldb_2.14.1            AnnotationFilter_1.14.0    
# [16] GenomicFeatures_1.42.3      AnnotationDbi_1.52.0        scater_1.18.6              
# [19] ggplot2_3.3.3               SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0
# [22] Biobase_2.50.0              GenomicRanges_1.42.0        GenomeInfoDb_1.26.7        
# [25] IRanges_2.24.1              S4Vectors_0.28.1            BiocGenerics_0.36.1        
# [28] MatrixGenerics_1.2.1        matrixStats_0.58.0         
# 
# loaded via a namespace (and not attached):
#   [1] shadowtext_0.0.8          fastmatch_1.1-3           BiocFileCache_1.14.0     
# [4] plyr_1.8.6                igraph_1.2.6              lazyeval_0.2.2           
# [7] splines_4.0.4             BiocParallel_1.24.1       digest_0.6.27            
# [10] GOSemSim_2.16.1           viridis_0.6.0             GO.db_3.12.1             
# [13] fansi_0.4.2               magrittr_2.0.1            memoise_2.0.0            
# [16] graphlayouts_0.7.1        Biostrings_2.58.0         R.utils_2.10.1           
# [19] askpass_1.1               enrichplot_1.10.2         prettyunits_1.1.1        
# [22] colorspace_2.0-0          ggrepel_0.9.1             blob_1.2.1               
# [25] rappdirs_0.3.3            dplyr_1.0.5               crayon_1.4.1             
# [28] RCurl_1.98-1.3            scatterpie_0.1.6          glue_1.4.2               
# [31] polyclip_1.10-0           gtable_0.3.0              zlibbioc_1.36.0          
# [34] XVector_0.30.0            DelayedArray_0.16.3       BiocSingular_1.6.0       
# [37] Rhdf5lib_1.12.1           HDF5Array_1.18.1          scales_1.1.1             
# [40] DOSE_3.16.0               DBI_1.1.1                 edgeR_3.32.1             
# [43] Rcpp_1.0.6                viridisLite_0.4.0         progress_1.2.2           
# [46] dqrng_0.3.0               bit_4.0.4                 rsvd_1.0.5               
# [49] ResidualMatrix_1.0.0      httr_1.4.2                fgsea_1.16.0             
# [52] ellipsis_0.3.2            farver_2.1.0              pkgconfig_2.0.3          
# [55] XML_3.99-0.6              R.methodsS3_1.8.1         scuttle_1.0.4            
# [58] dbplyr_2.1.1              locfit_1.5-9.4            utf8_1.2.1               
# [61] labeling_0.4.2            reshape2_1.4.4            tidyselect_1.1.1         
# [64] rlang_0.4.11              munsell_0.5.0             tools_4.0.4              
# [67] cachem_1.0.4              downloader_0.4            generics_0.1.0           
# [70] RSQLite_2.2.7             stringr_1.4.0             fastmap_1.1.0            
# [73] bit64_4.0.5               tidygraph_1.2.0           purrr_0.3.4              
# [76] ggraph_2.0.5              sparseMatrixStats_1.2.1   R.oo_1.24.0              
# [79] DO.db_2.9                 xml2_1.3.2                biomaRt_2.46.3           
# [82] compiler_4.0.4            rstudioapi_0.13           beeswarm_0.4.0           
# [85] curl_4.3                  tweenr_1.0.2              tibble_3.1.1             
# [88] statmod_1.4.35            stringi_1.5.3             bluster_1.0.0            
# [91] ProtGenerics_1.22.0       Matrix_1.3-4              vctrs_0.3.8              
# [94] pillar_1.6.0              lifecycle_1.0.0           rhdf5filters_1.2.0       
# [97] BiocManager_1.30.12       BiocNeighbors_1.8.2       cowplot_1.1.1            
# [100] data.table_1.14.0         bitops_1.0-7              irlba_2.3.3              
# [103] qvalue_2.22.0             rtracklayer_1.50.0        R6_2.5.0                 
# [106] vipor_0.4.5               MASS_7.3-53.1             assertthat_0.2.1         
# [109] rhdf5_2.34.0              openssl_1.4.3             withr_2.4.2              
# [112] GenomicAlignments_1.26.0  Rsamtools_2.6.0           GenomeInfoDbData_1.2.4   
# [115] hms_1.0.0                 grid_4.0.4                beachmat_2.6.4           
# [118] tidyr_1.1.3               rvcheck_0.1.8             DelayedMatrixStats_1.12.3
# [121] segmented_1.3-4           googledrive_2.0.0         ggforce_0.3.3            
# [124] ggbeeswarm_0.6.0 

