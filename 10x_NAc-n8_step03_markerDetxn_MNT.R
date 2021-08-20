### MNT 10x snRNA-seq workflow: step 03 - marker detection
###   **Region-specific analyses: nucleus accumbens (NAc)**
###     - Preprint: (3x) un-selected samples + (2x) NeuN-sorted samples
###     - Revision: (3x) samples (2 female, 1 NeuN-sorted)
### MNT 25Jun2021
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
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda",
     verbose=T)
    # sce.nac, chosen.hvgs.nac, pc.choice.nac, ref.sampleInfo, annotationTab.nac, cell_colors.nac

table(sce.nac$cellType)
    #  ambig.glial_A  ambig.glial_B        Astro_A        Astro_B drop.doublet_A drop.doublet_B 
    #             63             22             99           1000             36             52 
    # drop.doublet_C drop.doublet_D    drop.lowNTx        Inhib_A        Inhib_B        Inhib_C 
    #             41             21            529            251             40             98 
    #        Inhib_D        Inhib_E          Micro       MSN.D1_A       MSN.D1_B       MSN.D1_C 
    #            240             37            429           3927            239            283 
    #       MSN.D1_D       MSN.D1_E       MSN.D1_F       MSN.D2_A       MSN.D2_B       MSN.D2_C 
    #            718            638             86           4262            285            314 
    #       MSN.D2_D        Oligo_A        Oligo_B            OPC        OPC_COP 
    #             58            988           5146            651             18 

# First drop the "drop." clusters (669 nuclei, total)
sce.nac <- sce.nac[ ,-grep("drop.", sce.nac$cellType)]
sce.nac$cellType <- droplevels(sce.nac$cellType)
    
# Drop genes with all 0's
sce.nac <- sce.nac[!rowSums(assay(sce.nac, "counts"))==0, ]
    ## keeps 29680 genes


## Re-create 'logcounts' (don't want to use 'multiBatchNorm's down-scaling across donor 'batches')
# First 'hold' the MBN 'logcounts' for printing
sce.hold <- sce.nac

assay(sce.nac, "logcounts") <- NULL
sizeFactors(sce.nac) <- NULL
sce.nac <- logNormCounts(sce.nac)


### First make a list of Boolean param / cell subtype ===
# Will use this to assess more 'valid', non-noise-driving markers
cellSubtype.idx <- splitit(sce.nac$cellType)
medianNon0.nac <- lapply(cellSubtype.idx, function(x){
  apply(as.matrix(assay(sce.nac, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

sapply(medianNon0.nac, table)


## Traditional t-test implementation ===
mod <- with(colData(sce.nac), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

# Run pairwise t-tests
markers.nac.t.pw <- findMarkers(sce.nac, groups=sce.nac$cellType,
                                assay.type="logcounts", design=mod, test="t",
                                direction="up", pval.type="all", full.stats=T)

sapply(markers.nac.t.pw, function(x){table(x$FDR<0.05)})


## WMW: Blocking on sample (this test doesn't take 'design=' argument) ===
markers.nac.wilcox.block <- findMarkers(sce.nac, groups=sce.nac$cellType,
                                        assay.type="logcounts", block=sce.nac$donor, test="wilcox",
                                        direction="up", pval.type="all", full.stats=T)


sapply(markers.nac.wilcox.block, function(x){table(x$FDR<0.05)})
    ## none for about ~2/3 of these; ignore


## Binomial ===
markers.nac.binom.block <- findMarkers(sce.nac, groups=sce.nac$cellType,
                                       assay.type="logcounts", block=sce.nac$donor, test="binom",
                                       direction="up", pval.type="all", full.stats=T)


sapply(markers.nac.binom.block, function(x){table(x$FDR<0.05)})
    ## none; ignore

# Add respective 'non0median' column to the stats for each set of markers
for(i in names(markers.nac.t.pw)){
  markers.nac.t.pw[[i]] <- cbind(markers.nac.t.pw[[i]],
                                 medianNon0.nac[[i]][match(rownames(markers.nac.t.pw[[i]]),
                                                           names(medianNon0.nac[[i]]))])
  colnames(markers.nac.t.pw[[i]])[27] <- "non0median"
}

sapply(markers.nac.t.pw, function(x){table(x$FDR<0.05 & x$non0median == TRUE)["TRUE"]})
    #  Astro_A.TRUE      Astro_B.TRUE      Inhib_A.TRUE      Inhib_B.TRUE      Inhib_C.TRUE
    #            97               250                78                40                69
    #  Inhib_D.TRUE      Inhib_E.TRUE Macro_infilt.TRUE        Micro.TRUE  Micro_resting.NA
    #            40                70               154               227                NA
    # MSN.D1_A.TRUE     MSN.D1_B.TRUE     MSN.D1_C.TRUE     MSN.D1_D.TRUE     MSN.D1_E.TRUE
    #            12                16                15                21                38
    # MSN.D1_F.TRUE     MSN.D2_A.TRUE     MSN.D2_B.TRUE     MSN.D2_C.TRUE     MSN.D2_D.TRUE
    #            33                16                 6                99                41
    #  Oligo_A.TRUE      Oligo_B.TRUE          OPC.TRUE      OPC_COP.TRUE
    #           213                59                67                23

## Save all these for future reference
save(markers.nac.t.pw, medianNon0.nac,#markers.nac.wilcox.block, markers.nac.binom.block,
     file="rdas/revision/markers-stats_NAc-n8_findMarkers-SN-LEVEL_MNT2021.rda")


# # As needed:
# load("rdas/revision/markers-stats_NAc-n8_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
#     # markers.nac.t.pw, medianNon0.nac

# Print these to pngs
markerList.t.pw <- lapply(markers.nac.t.pw, function(x){
  rownames(x)[x$FDR < 0.05 & x$non0median == TRUE]
  }
)
genes.top40.t <- lapply(markerList.t.pw, function(x){head(x, n=40)})

#dir.create("pdfs/revision/NAc/")
smaller.set <- names(genes.top40.t)[lengths(genes.top40.t) <= 20 & lengths(genes.top40.t) != 0]
left.set <- names(genes.top40.t)[lengths(genes.top40.t) > 20]

# Smaller graphical window
for(i in smaller.set){
  png(paste0("pdfs/revision/NAc/NAc_t_pairwise_top40markers-", i, "_logExprs_MNT2021.png"), height=950, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         features = genes.top40.t[[i]],
                         features_name = i,
                         anno_name = "cellType",
                         ncol=5, point_alpha=0.4) +
      scale_color_manual(values = cell_colors.nac) +
      ggtitle(label=paste0("NAc ", i, " top markers: single-nucleus-level p.w. t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}

# 20-40 markers
for(i in left.set){
  png(paste0("pdfs/revision/NAc/NAc_t_pairwise_top40markers-", i, "_logExprs_MNT2021.png"), height=1900, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         features = genes.top40.t[[i]],
                         features_name = i,
                         anno_name = "cellType",
                         ncol=5, point_alpha=0.4) +
      scale_color_manual(values = cell_colors.nac) +
      ggtitle(label=paste0("NAc ", i, " top markers: single-nucleus-level p.w. t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}



### Cluster-vs-all single-nucleus-level iteration ========================
## Load SCE with new info
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_NAc-n8_cleaned-combined_SCE_MNT2021.rda",
     verbose=T)
# sce.nac, chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac, ref.sampleInfo, annotationTab.nac, cell_colors.nac


# First drop "drop.lowNTx_" (669 nuclei)
sce.nac <- sce.nac[ ,-grep("drop.",sce.nac$cellType)]
sce.nac$cellType <- droplevels(sce.nac$cellType)

# Remove 0 genes across all nuclei
sce.nac <- sce.nac[!rowSums(assay(sce.nac, "counts"))==0, ]  # keeps same 29680 genes


## Re-create 'logcounts' (don't want to use 'multiBatchNorm's down-scaling across donor 'batches')
# First 'hold' the MBN 'logcounts' for printing
sce.hold <- sce.nac

assay(sce.nac, "logcounts") <- NULL
sizeFactors(sce.nac) <- NULL
sce.nac <- logNormCounts(sce.nac)

## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.nac), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

markers.nac.t.1vAll <- list()
for(i in levels(sce.nac$cellType)){
  # Make temporary contrast
  sce.nac$contrast <- ifelse(sce.nac$cellType==i, 1, 0)
  # Test cluster vs. all others
  markers.nac.t.1vAll[[i]] <- findMarkers(sce.nac, groups=sce.nac$contrast,
                                          assay.type="logcounts", design=mod, test="t",
                                          std.lfc=TRUE,
                                          direction="up", pval.type="all", full.stats=T)
}
    ## Since all other stats are the same, and don't really use the non-standardized
    #    logFC, just generate one object, unlike before

class(markers.nac.t.1vAll[["Oligo_A"]])
    # a SimpleList of length 2, named "0" and "1" (from the temporary 'contrast')
    # -> we want the second entry, named "1"
    #    (for other purposes, might be interesting to look into that "0" entry, which
    #     is basically what genes are depleted in the cell type of interest)


sapply(markers.nac.t.1vAll, function(x){
  table(x[["1"]]$stats.0$log.FDR < log(.001))
})
    #

# Do some reorganizing
markers.nac.t.1vAll <- lapply(markers.nac.t.1vAll, function(x){
  # Basically take the 'stats.[1 or 0]' since is redundant with the 'summary'-level stats
  lapply(x, function(y){ y[ ,4] }) 
})

# Re-name std.lfc column and the entries; add non-0-median info
for(i in names(markers.nac.t.1vAll)){
  colnames(markers.nac.t.1vAll[[i]][["0"]])[1] <- "std.logFC"
  colnames(markers.nac.t.1vAll[[i]][["1"]])[1] <- "std.logFC"
  # Add non0median Boolean - might be informative for both sets of stats
  markers.nac.t.1vAll[[i]][["0"]] <- cbind(markers.nac.t.1vAll[[i]][["0"]],
                                           medianNon0.nac[[i]][match(rownames(markers.nac.t.1vAll[[i]][["0"]]),
                                                                     names(medianNon0.nac[[i]]))])
  colnames(markers.nac.t.1vAll[[i]][["0"]])[4] <- "non0median"
  
  # "1" aka 'enriched'
  markers.nac.t.1vAll[[i]][["1"]] <- cbind(markers.nac.t.1vAll[[i]][["1"]],
                                           medianNon0.nac[[i]][match(rownames(markers.nac.t.1vAll[[i]][["1"]]),
                                                                     names(medianNon0.nac[[i]]))])
  colnames(markers.nac.t.1vAll[[i]][["1"]])[4] <- "non0median"
  
  # Then re-name the entries to more interpretable, because we'll keeping both contrasts
  names(markers.nac.t.1vAll[[i]]) <- paste0(i,c("_depleted", "_enriched"))
}


## Let's save this along with the previous pairwise results
save(markers.nac.t.pw, markers.nac.t.1vAll, medianNon0.nac,
     file="rdas/revision/markers-stats_NAc-n8_findMarkers-SN-LEVEL_MNT2021.rda")

## MNT add 17Jul: change 'Macro_infilt' annotation to 'Macrophage':
load("rdas/revision/markers-stats_NAc-n8_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)

names(medianNon0.nac)[names(medianNon0.nac)=="Macro_infilt"] <- "Macrophage"
names(markers.nac.t.pw)[names(markers.nac.t.pw)=="Macro_infilt"] <- "Macrophage"
names(markers.nac.t.1vAll)[names(markers.nac.t.1vAll)=="Macro_infilt"] <- "Macrophage"
    # and for the _enriched vs _depleted for this one:
    names(markers.nac.t.1vAll[["Macrophage"]]) <- gsub("Macro_infilt", "Macrophage",
                                                       names(markers.nac.t.1vAll[["Macrophage"]]))

save(markers.nac.t.pw, markers.nac.t.1vAll, medianNon0.nac,
     file="rdas/revision/markers-stats_NAc-n8_findMarkers-SN-LEVEL_MNT2021.rda")
    # -> since this is pretty negligible, don't worry about re-printing the top 40 marker plots
    #    (just re-print the top 40 marker lists, below)

## Print these to pngs
markerList.t.1vAll <- lapply(markers.nac.t.1vAll, function(x){
  rownames(x[[2]])[ x[[2]]$log.FDR < log(0.05) & x[[2]]$non0median==TRUE ]
  }
)
genes.top40.t <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0("pdfs/revision/NAc/NAc_t_1vALL_top40markers-",i,"_logExprs_MNT2021.png"), height=1900, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         features = genes.top40.t[[i]], 
                         features_name = i,
                         anno_name = "cellType",
                         ncol=5, point_alpha=0.4) +
      scale_color_manual(values = cell_colors.nac) +  
      ggtitle(label=paste0("NAc ", i, " top markers: 'cluster-vs-all-others' t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}


## How do they intersect?
markerList.t.pw <- lapply(markers.nac.t.pw, function(x){
  rownames(x)[ x$FDR < 0.05 & x$non0median==TRUE ]
  }
)

# From pairwise t-tests, FDR < 0.05 (& median restriction)
lengths(markerList.t.pw)

# From cluster-vs-all others, FDR < 0.05 (& median restriction)
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
    #


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

write.csv(top40genes, file="tables/revision/top40genesLists_NAc-n8_cellType_SN-LEVEL-tests_MNT2021.csv",
          row.names=FALSE)



## Make marker array for Supp figure (MNT suggested panel A) ===========
load("rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda", verbose=T)
    # sce.nac, chosen.hvgs.nac, pc.choice.nac, ref.sampleInfo, annotationTab.nac, cell_colors.nac

load("rdas/markers-stats_NAc-n5_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.nac.t.pw, markers.nac.t.1vAll, medianNon0.nac

# Clean up
sce.nac <- sce.nac[ ,-grep("drop.",sce.nac$cellType)]
sce.nac$cellType <- droplevels(sce.nac$cellType)

# First make interneuron subset
sce.nac.int <- sce.nac[ ,grep("Inhib_", sce.nac$cellType)]
sce.nac.int$cellType <- droplevels(sce.nac.int$cellType)

# Take top four for 4 inhib. interneuron pops
# Rm the 'Micro_resting' for now, bc there are no pairwise markers for that cell type
markers.nac.t.pw[["Micro_resting"]] <- NULL
topToPrint <- as.data.frame(sapply(markers.nac.t.pw, function(x) {
  head(rownames(x)[x$FDR < 0.05 & x$non0median==TRUE], n=4)}))
topToPrint <- topToPrint[grep("Inhib_", names(topToPrint))]

table(unlist(topToPrint) %in% rownames(sce.nac.int)) # good

topToPrint

# Print
pdf("pdfs/revision/pubFigures/suppFig_NAc_interneuron-marker-array_MNT2021.pdf", height=8, width=6.5)
print(
  plotExpressionCustom(sce.nac.int, features=c(t(topToPrint)), features_name="",
                 anno_name="cellType", point_alpha=0.6, point_size=1.2, ncol=5, scales="free_y") +
    scale_color_manual(values = cell_colors.nac) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12), plot.title = element_text(size = 25)) +
    ggtitle(label="Inhib_A       Inhib_B         Inhib_C        Inhib_D        Inhib_E") + xlab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 11),
          axis.title.y = element_text(angle = 90, size = 16),
          plot.title = element_text(size = 15),
          panel.grid.major=element_line(colour="grey95", size=0.8),
          panel.grid.minor=element_line(colour="grey95", size=0.4))
)
dev.off()


## Supp Fig 5: Other markers pointed out in text ===
genes2print <- c("DRD1", "DRD2", "CASZ1", "GPR6", "EBF1", "GRM8")

pdf("pdfs/revision/pubFigures/suppFig_NAc_other-MSN-markers_MNT2021.pdf", height=4, width=7.5)
print(
  plotExpressionCustom(sce.nac, features=genes2print, features_name="",
                       anno_name="cellType", point_alpha=0.3, point_size=1, ncol=2, scales="free_y") +
    scale_color_manual(values = cell_colors.nac) +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
          axis.title.y = element_text(angle = 90, size = 10),
          panel.grid.major=element_line(colour="grey95", size=0.8),
          panel.grid.minor=element_line(colour="grey95", size=0.4))
)
dev.off()



## Top markers for D1_A / D2_A often co-expressed ===
load("rdas/markers-stats_NAc-n5_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.nac.t.pw, markers.nac.t.1vAll

cell.idx <- splitit(sce.nac$cellType)
dat <- as.matrix(assay(sce.nac, "logcounts"))
markerList.t.1vAll <- lapply(markers.nac.t.1vAll, function(x){
  rownames(x[[2]])[ x[[2]]$log.FDR < log(0.05) & x[[2]]$non0median==TRUE ]
  }
)
genes <- unique(c(head(markerList.t.1vAll[["MSN.D1_A"]], n=20),
                head(markerList.t.1vAll[["MSN.D2_A"]], n=20)))

#pdf('pdfs/revision/pubFigures/suppFigure_heatmap-Exprs_NAc-n8_D1or2_A-markers-1vAlltest_MNT2021.pdf',
pdf('pdfs/revision/pubFigures/suppFigure_heatmap-medians-Exprs_NAc-n8_D1or2_A-markers-1vAlltest_MNT2021.pdf',
    useDingbats=TRUE, height=6, width=10)
#current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMeans(dat[genes, ii])))
    # or medians:
    current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMedians(dat[genes, ii])))
    # For some reason rownames aren't kept:
    rownames(current_dat) <- genes
pheatmap(t(current_dat), cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(0.02, 4, length.out = 101),
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "OrRd"))(100),
         fontsize_row = 15, fontsize_col=14)
grid::grid.text(label="log2-\nExprs", x=0.975, y=0.60, gp=grid::gpar(fontsize=10))
dev.off()



## For revision - can show more evidence for overlapping another way? =====
load("rdas/revision/markers-stats_NAc-n8_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.nac.t.pw, markers.nac.t.1vAll, medianNon0.nac

colnames(markers.nac.t.pw[["MSN.D1_A"]])
    # [1] "p.value"             "FDR"                 "summary.stats"       "stats.Astro_A"      
    # [5] "stats.Astro_B"       "stats.Inhib_A"       "stats.Inhib_B"       "stats.Inhib_C"      
    # [9] "stats.Inhib_D"       "stats.Inhib_E"       "stats.Macro_infilt"  "stats.Micro"        
    # [13] "stats.Micro_resting" "stats.MSN.D1_B"      "stats.MSN.D1_C"      "stats.MSN.D1_D"     
    # [17] "stats.MSN.D1_E"      "stats.MSN.D1_F"      "stats.MSN.D2_A"      "stats.MSN.D2_B"     
    # [21] "stats.MSN.D2_C"      "stats.MSN.D2_D"      "stats.Oligo_A"       "stats.Oligo_B"      
    # [25] "stats.OPC"           "stats.OPC_COP"       "non0median"

# E.g. pairwise result:
head(markers.nac.t.pw[["MSN.D1_A"]][ ,4])
    # DataFrame with 6 rows and 3 columns
    #                logFC log.p.value   log.FDR
    #            <numeric>   <numeric> <numeric>
    # AGBL1       2.738819   -218.2578 -214.0532
    # ADAMTS19    2.967084   -393.2135 -388.4605
    # AP001993.1  0.997718   -140.8300 -136.9631
    # PRKCH       2.513853   -577.9395 -572.8060
    # AC068875.1  1.263544   -141.9177 -138.0492
    # ADAMTSL1    0.992064    -73.2992  -70.0266

sapply(colnames(markers.nac.t.pw[["MSN.D1_A"]])[c(4:26)], function(x){
  table(markers.nac.t.pw[["MSN.D1_A"]][ ,x]$log.FDR < log(0.05))
  })
    # D1_A closest to _D (which is probably expected) and the next closest (aka least # DEGs)
    #      is to D2_A

# With non0median restriction:
sapply(colnames(markers.nac.t.pw[["MSN.D1_A"]])[c(4:26)], function(x){
  table(markers.nac.t.pw[["MSN.D1_A"]][ ,x]$log.FDR < log(0.05) & 
          markers.nac.t.pw[["MSN.D1_A"]]$non0median == TRUE)
})
    # same pattern


# And pairwise result from D2_A:
sapply(colnames(markers.nac.t.pw[["MSN.D2_A"]])[c(4:26)], function(x){
  table(markers.nac.t.pw[["MSN.D2_A"]][ ,x]$log.FDR < log(0.05) & 
          markers.nac.t.pw[["MSN.D2_A"]]$non0median == TRUE)["TRUE"]
})
    # stats.Astro_A.TRUE       stats.Astro_B.TRUE       stats.Inhib_A.TRUE       stats.Inhib_B.TRUE 
    #               3944                     3726                     2030                     1521 
    # stats.Inhib_C.TRUE       stats.Inhib_D.TRUE       stats.Inhib_E.TRUE  stats.Macro_infilt.TRUE 
    #               1614                     2165                     1221                     2928 
    #   stats.Micro.TRUE stats.Micro_resting.TRUE      stats.MSN.D1_A.TRUE      stats.MSN.D1_B.TRUE 
    #               3909                     4725                     *776                     1577 
    #stats.MSN.D1_C.TRUE      stats.MSN.D1_D.TRUE      stats.MSN.D1_E.TRUE      stats.MSN.D1_F.TRUE 
    #               1935                     *791                     2296                     1565 
    #stats.MSN.D2_B.TRUE      stats.MSN.D2_C.TRUE      stats.MSN.D2_D.TRUE       stats.Oligo_A.TRUE 
    #               1146                     2310                     1196                     3743 
    # stats.Oligo_B.TRUE           stats.OPC.TRUE       stats.OPC_COP.TRUE 
    #               4297                     3859                     1857




### (From preprint exploration)
# ## Markers pulled from Gokce, et al (doi: 10.1016/j.celrep.2016.06.059) =========
# markers.gokce <- list(
#   "D1.MSN" = c("Tac1","Drd1","Asic4","Slc35d3","Pdyn","Sfxn1","Nrxn1"),
#                       # ^ edited from 'Drd1a', Accn4
#   "D2.MSN" = c("Penk","Adora2a","Drd2","Gpr6","Grik3","Gpr52","Gnas"),
#                       # ^ edited from 'A2a'
#   "D1.Pcdh8" = c("Pcdh8","Adarb2","Tacr1","Tac1" ,"Nrxn2","Sema3e","Sema4a","Sema5a","Sema5b",
#                  "Sema6d","Pcdh7","Ptprg","Ptprm","Ptpro","Ptpru","TAC3","Elavl4",
#                                                                   # ^ edited from 'Tac2'
#                  "Khdrbs3","Rbm20","Aff2","Lrpprc","Celf4",
#                  # Depleted set:
#                  "Nlgn1", "Calb1"),
#   "D1.Foxp1" = c("Foxp1","Camk4"),
#   "D2.Htr7" = c("Htr7","AGTR1","Penk","Tac1","Ptprt","Ngfr","Grik3","Cacng5",
#                         # ^ edited from 'Agtr1a'
#                 "Tmeff2","Sox9","Sp8","Runx1","Mafb","Litaf",
#                 # Depleted set:
#                 "Cacna2d3","Synpr"),
#   "D2.Synpr" = c("Synpr"),
#   "gradient" = c("Dner","Cxcl14","Tnnt1","Meis2","Cartpt","Kcnip1","Calb1",
#                  "Crym","Cnr1","Nnat","Gfra1","Wfs1","Th")
# )
# 
# markers.gokce <- lapply(markers.gokce, toupper)
# 
# # Load SCE
# load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda",
#      verbose=T)
#     # sce.nac, chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac, ref.sampleInfo
# 
# 
# # First drop "ambig.lowNtrxts" (93 nuclei)
# sce.nac <- sce.nac[ ,sce.nac$cellType != "ambig.lowNtrxts"]
# sce.nac$cellType <- droplevels(sce.nac$cellType)
# 
# # Which are there?
# sapply(markers.gokce, function(x){x %in% rownames(sce.nac)})  # most of them
# 
# lapply(sapply(markers.gokce, function(x){x %in% rownames(sce.nac)}),  # most of them
#        function(n){which(n==FALSE)})
#     # So 'Drd1a', 'Accn4', 'A2a',       'Tac2',       'Agtr1a'
# 
# 
#     # Exploring/identifying homologous gene names ====
#     load("/dcl01/ajaffe/data/lab/singleCell/day_rat_snRNAseq/SCE_rat-NAc-PBd_w_matchingHsap-NAc-PBd_HomoloGene.IDs_MNT.rda", verbose=T)
#     # sce.rat.PBsub, sce.hsap.PBsub, Readme
#     
#         # 'sce.hsap.PBsub' has 'HomoloGene.ID'
# 
#     hom = read.delim("http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt",
#                      as.is=TRUE)
#     
#     hom_mm <- hom[hom$Common.Organism.Name == "mouse, laboratory", ]
#     
#     c('Drd1a', 'Accn4', 'A2a','Tac2','Agtr1a') %in% hom_mm$Symbol
#     # FALSE     FALSE   FALSE   TRUE    TRUE
#     # Drd1      Asic4  Adora2a
#     
#     # Then for 'Tac2'
#     hom_mm$HomoloGene.ID[which(hom_mm$Symbol=="Tac2")]  # 7560
#     rowData(sce.hsap.PBsub)$Symbol[rowData(sce.hsap.PBsub)$HomoloGene.ID==7560] # none...
#     
#     hom_hs <- hom[hom$Common.Organism.Name == "human", ]
#     hom_hs[hom_hs$HomoloGene.ID==7560, ]  # symbol is TAC3
#     'TAC3' %in% rowData(sce.nac)$Symbol # TRUE
#         # ahhh so this one just didn't have a shared homolog with rat, I guess
#     
#     
#     hom_mm$HomoloGene.ID[which(hom_mm$Symbol=="Agtr1a")]
#     rowData(sce.hsap.PBsub)$Symbol[rowData(sce.hsap.PBsub)$HomoloGene.ID==3556]
#         # AGTR1
#     # end find synonyms =======
# 
# 
# 
# 
# ## Let's make a new dir and files for these graphics, since these are of various length
# #dir.create("pdfs/exploration/gokce-etal_markers/")
# 
# for(i in names(markers.gokce)){
#   pdf(paste0("./pdfs/exploration/gokce-etal_markers/",i,"-mouseStriatum-markers_human-NAcExpression_Apr2020.pdf"), height=2.6, width=3)
#   # Print each gene's expression in its own page of the pdf
#   for(g in 1:length(markers.gokce[[i]])){
#     print(
#       plotExpression(sce.nac, exprs_values = "logcounts", features=c(markers.gokce[[i]][g]),
#                      x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7,
#                      add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
#                                                   geom = "crossbar", width = 0.3,
#                                                   colour=tableau20[1:14]) +
#         theme(axis.text.x = element_text(angle=90, hjust=1, size=5.5), axis.text.y = element_text(size=7.5),
#               plot.title = element_text(size=7)) +  
#         ggtitle(label=paste0(i, " markers in human NAc subclusters: ", markers.gokce[[i]][g]))
#     )
#   }
#   dev.off()
# }

rm(list=ls())
### Session info for 30Jun2021 ==================================================
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
#   [1] parallel  stats4    stats     graphics  grDevices datasets  utils    
# [8] methods   base     
# 
# other attached packages:
#   [1] pheatmap_1.0.12             RColorBrewer_1.1-2         
# [3] lattice_0.20-41             limma_3.46.0               
# [5] jaffelab_0.99.30            rafalib_1.0.0              
# [7] DropletUtils_1.10.3         batchelor_1.6.2            
# [9] scran_1.18.5                scater_1.18.6              
# [11] ggplot2_3.3.3               EnsDb.Hsapiens.v86_2.99.0  
# [13] ensembldb_2.14.1            AnnotationFilter_1.14.0    
# [15] GenomicFeatures_1.42.3      AnnotationDbi_1.52.0       
# [17] SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0
# [19] Biobase_2.50.0              GenomicRanges_1.42.0       
# [21] GenomeInfoDb_1.26.7         IRanges_2.24.1             
# [23] S4Vectors_0.28.1            BiocGenerics_0.36.1        
# [25] MatrixGenerics_1.2.1        matrixStats_0.58.0         
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
# [57] rhdf5_2.34.0              curl_4.3                 
# [59] memoise_2.0.0             gridExtra_2.3            
# [61] segmented_1.3-3           biomaRt_2.46.3           
# [63] stringi_1.5.3             RSQLite_2.2.7            
# [65] BiocParallel_1.24.1       rlang_0.4.10             
# [67] pkgconfig_2.0.3           bitops_1.0-7             
# [69] purrr_0.3.4               Rhdf5lib_1.12.1          
# [71] labeling_0.4.2            GenomicAlignments_1.26.0 
# [73] cowplot_1.1.1             bit_4.0.4                
# [75] tidyselect_1.1.1          magrittr_2.0.1           
# [77] R6_2.5.0                  generics_0.1.0           
# [79] DelayedArray_0.16.3       DBI_1.1.1                
# [81] pillar_1.6.0              withr_2.4.2              
# [83] RCurl_1.98-1.3            tibble_3.1.1             
# [85] crayon_1.4.1              utf8_1.2.1               
# [87] BiocFileCache_1.14.0      viridis_0.6.0            
# [89] progress_1.2.2            locfit_1.5-9.4           
# [91] grid_4.0.4                blob_1.2.1               
# [93] digest_0.6.27             R.utils_2.10.1           
# [95] openssl_1.4.3             munsell_0.5.0            
# [97] beeswarm_0.3.1            viridisLite_0.4.0        
# [99] vipor_0.4.5               askpass_1.1  
