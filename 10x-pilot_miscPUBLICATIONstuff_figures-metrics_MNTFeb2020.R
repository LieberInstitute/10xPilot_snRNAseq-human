### MNT 10x snRNA-seq workflow: (step 04?:)
### Miscellaneous lookings-into / finalizing graphics for manuscript
###   - Brief neuron-specific clustering for DLPFC (Maynard-Collado-Torres et al.)
###   - 10x pilot snRNA-seq paper (Tran-Maynard et al.)
##################################################################################

library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(batchelor)
library(DropletUtils)
library(jaffelab)


### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

#VIRIDFUN(name = colour_by_name, discrete = TRUE)

### DLPFC ========================================================================
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_DLPFC-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    #   sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc, ref.sampleInfo

#sce.dlpfc$cellType <- factor(sce.dlpfc$cellType,
#                             levels=setdiff(unique(sce.dlpfc$cellType), "Ambig.lowNtrxts"))

sce.dlpfc.noAmbig <- sce.dlpfc[ ,!sce.dlpfc$cellType=="Ambig.lowNtrxts"]
sce.dlpfc.noAmbig$prelimCluster.labels <- factor(paste0(sce.dlpfc.noAmbig$prelimCluster," (",sce.dlpfc.noAmbig$cellType,")"))
levels(sce.dlpfc.noAmbig$prelimCluster.labels) <- levels(sce.dlpfc.noAmbig$prelimCluster.labels)[order(as.numeric(ss(levels(sce.dlpfc.noAmbig$prelimCluster.labels)," ",1)))]


#dir.create("pdfs/pubFigures/")
pdf("pdfs/pubFigures/tSNE-DLPFC-n2_collapsedClusters-cellType_MNTFeb2020.pdf", height=5, width=5)
plotTSNE(sce.dlpfc.noAmbig, colour_by="cellType", point_alpha=0.6, text_by="cellType",
         text_size=7, point_size=5.0, add_legend=F, theme_size=20)
dev.off()

# From stack overflow haha
# https://stackoverflow.com/questions/52297978/decrease-overal-legend-size-elements-and-text
addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 10, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
  }

#pdf("pdfs/pubFigures/tSNE-DLPFC-n2_prelimClusters-withCellType_MNTFeb2020.pdf", height=5, width=7)
plotb <- plotTSNE(sce.dlpfc.noAmbig, point_alpha=0.6, text_by="cellType", colour_by="prelimCluster.labels",
         text_size=7, point_size=5.0, add_legend=T, theme_size=20) + guides(fill=guide_legend(title=NULL))
plotb

# then
addSmallLegend(plotb)
#dev.off()




## Visual sub-structure of DLPFC neurons? (for ST) ============
sce.dlpfc.neu <- sce.dlpfc[ ,sce.dlpfc$cellType %in% c("Excit", "Inhib")]
sizeFactors(sce.dlpfc.neu) <- NULL
assay(sce.dlpfc.neu, "logcounts") <- NULL
sce.dlpfc.neu <- logNormCounts(sce.dlpfc.neu)
geneVar.neuron <- modelGeneVar(sce.dlpfc.neu)
chosen.hvgs.neuron <- geneVar.neuron$bio > 0
    ## of sum 11701

# remove reducedDims
reducedDim(sce.dlpfc.neu, "PCA") <- NULL
reducedDim(sce.dlpfc.neu, "PCA_opt") <- NULL
reducedDim(sce.dlpfc.neu, "TSNE") <- NULL
reducedDim(sce.dlpfc.neu, "UMAP") <- NULL

# run PCA (just take 50 PCs)
set.seed(109)
sce.dlpfc.neu <- runPCA(sce.dlpfc.neu, subset_row=chosen.hvgs.neuron, ncomponents=50,
                        BSPARAM=BiocSingular::RandomParam())

sum(attr(reducedDim(sce.dlpfc.neu, "PCA"), "percentVar"))
#[1] 39.63105

## t-SNE
set.seed(109)
sce.dlpfc.neu <- runTSNE(sce.dlpfc.neu, dimred="PCA")



sce.dlpfc.neu$prelimCluster.labels <- factor(paste0(sce.dlpfc.neu$prelimCluster," (",sce.dlpfc.neu$cellType,")"))
levels(sce.dlpfc.neu$prelimCluster.labels) <- levels(sce.dlpfc.neu$prelimCluster.labels)[order(as.numeric(ss(levels(sce.dlpfc.neu$prelimCluster.labels)," ",1)))]


# Plot
pdf("pdfs/pubFigures/tSNE-DLPFC-n2_neuronalSub-exploration_MNTFeb2020.pdf", height=5, width=5)
plotc <- plotTSNE(sce.dlpfc.neu, colour_by="prelimCluster.labels", point_size=5.0,
                  text_by="prelimCluster.labels", text_size=3.0, theme_size=15) +
  guides(fill=guide_legend(title=NULL))
addSmallLegend(plotc, textSize=6)
plotTSNE(sce.dlpfc.neu, colour_by="sample", point_size=5.0, theme_size=15)
dev.off()



# Save work done to these
save(sce.dlpfc, sce.dlpfc.neu, file="rdas/regionSpecific_DLPFC-n2_neuronalSCE-and-moreExploration.rda")







### METRICS ====================================================

## NAc, all n=5
load("rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda", verbose=T)
    #sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all, ref.sampleInfo

## Distribution of n transcripts
# 'Depth' per sample
nac.samples <- splitit(sce.nac.all$sample)
sapply(nac.samples, function(x){quantile(sce.nac.all[ ,x]$sum)})
# 'Capture' per cell type
cellType_idx <- splitit(sce.nac.all$cellType)
sapply(cellType_idx, function(x){quantile(sce.nac.all[ ,x]$sum)})

# n Nuclei from each sample
table(sce.nac.all$sample)
    # nac.5161      nac.5212      nac.5287 nac.neun.5182 nac.neun.5207
    #     2067          1774           707          4267          4426


# Make diff object for 'ambig.lowNtrxts' removed
sce.nac.noLNT <- sce.nac.all[ ,!sce.nac.all$cellType=="ambig.lowNtrxts"]
sce.nac.noLNT$cellType.final <- droplevels(sce.nac.noLNT$cellType.final)

## Proportion of cell types
table(sce.nac.noLNT$cellType.final)
prop.table(table(sce.nac.noLNT$cellType.final))
    #       Astro     Inhib.1     Inhib.2     Inhib.3     Inhib.4       Micro
    # 0.041451171 0.001901430 0.004259203 0.020991786 0.013918467 0.013766352
    #    MSN.D1.1    MSN.D1.2    MSN.D1.3    MSN.D1.4    MSN.D2.1    MSN.D2.2
    # 0.010039550 0.022893216 0.054685123 0.285442653 0.022817159 0.276163675
    #       Oligo         OPC
    # 0.213492546 0.018177670

# How many MSNs in our dataset?
sum(prop.table(table(sce.nac.noLNT$cellType.final))[grep("MSN", levels(sce.nac.noLNT$cellType.final))])
    # However this is including two samples where we enriched on neurons


table(sce.nac.noLNT$cellType.final, sce.nac.noLNT$sample)
    # 
    #          nac.5161 nac.5212 nac.5287 nac.neun.5182 nac.neun.5207
    # Astro         149      384       12             0             0
    # Inhib.1         1        3        0            16             5
    # Inhib.2         1        1        1            42            11
    # Inhib.3         7        7        9            86           167
    # Inhib.4         9        8        4           104            58
    # Micro          72       72       37             0             0
    # MSN.D1.1        2        0        0           117            13
    # MSN.D1.2       10        3        0           285             3
    # MSN.D1.3       17        8        6           369           319
    # MSN.D1.4      178      169       72          1505          1829
    # MSN.D2.1        9        6        3           134           148
    # MSN.D2.2       41      113        5          1602          1870
    # Oligo        1454      854      499             0             0
    # OPC            98      104       37             0             0

## Proportions of cell type by sample
apply(table(sce.nac.noLNT$cellType.final, sce.nac.noLNT$sample), 2, function(x){round(prop.table(x),3)})
    #          nac.5161 nac.5212 nac.5287 nac.neun.5182 nac.neun.5207
    # Astro       0.073    0.222    0.018         0.000         0.000
    # Inhib.1     0.000    0.002    0.000         0.004         0.001
    # Inhib.2     0.000    0.001    0.001         0.010         0.002
    # Inhib.3     0.003    0.004    0.013         0.020         0.038
    # Inhib.4     0.004    0.005    0.006         0.024         0.013
    # Micro       0.035    0.042    0.054         0.000         0.000
    # MSN.D1.1    0.001    0.000    0.000         0.027         0.003
    # MSN.D1.2    0.005    0.002    0.000         0.067         0.001
    # MSN.D1.3    0.008    0.005    0.009         0.087         0.072
    # MSN.D1.4    0.087    0.098    0.105         0.353         0.414
    # MSN.D2.1    0.004    0.003    0.004         0.031         0.033
    # MSN.D2.2    0.020    0.065    0.007         0.376         0.423
    # Oligo       0.710    0.493    0.728         0.000         0.000
    # OPC         0.048    0.060    0.054         0.000         0.000

    # Of neurons, how many are MSNs?  Let's just look at the NeuN-sorted samples
    
    table(ss(as.character(sce.nac.noLNT$cellType.final[grep("neun", sce.nac.noLNT$sample)]),
             "\\.",1),
          sce.nac.noLNT$sample[grep("neun", sce.nac.noLNT$sample)])
        #        nac.neun.5182 nac.neun.5207
        #  Inhib           248           241
        #  MSN            4012          4182
    apply(table(ss(as.character(sce.nac.noLNT$cellType.final[grep("neun", sce.nac.noLNT$sample)]),
                   "\\.",1),
                sce.nac.noLNT$sample[grep("neun", sce.nac.noLNT$sample)]),2,prop.table)
        #        nac.neun.5182 nac.neun.5207
        #  Inhib    0.05821596     0.0544879
        #  MSN      0.94178404     0.9455121

    
## Numbers for markers ===
load("rdas/markers-stats_NAc-n5_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.nac.t.design, markers.nac.t.1vAll

sapply(markers.nac.t.design, function(x){sum(x$FDR < 0.05)})
    #    Astro  Inhib.1  Inhib.2  Inhib.3  Inhib.4    Micro MSN.D1.1 MSN.D1.2
    #      831      340      179      256      283     1417       82      295
    # MSN.D1.3 MSN.D1.4 MSN.D2.1 MSN.D2.2    Oligo      OPC
    #       87       61      230       66      493      372


## Are DRD1/DRD2 ever expressed together? ===
plotExpression(sce.nac.noLNT, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7, ncol=5,
               add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                            geom = "crossbar", width = 0.3,
                                            colour=rep(tableau20[1:14], 2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25))

geneExprs <- assay(sce.nac.noLNT, "logcounts")

sce.nac.noLNT$MSN.D1 <- NA
sce.nac.noLNT$MSN.D1[grep("MSN.D1", sce.nac.noLNT$cellType.final)] <- TRUE

sce.nac.noLNT$MSN.D1.0 <- NA
sce.nac.noLNT$MSN.D1.0 <- sce.nac.noLNT$MSN.D1 & geneExprs["DRD1", ]==0

    # Now plot & color by this - is DRD2 expressed?
    plotExpression(sce.nac.noLNT, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
                   x="cellType.final", colour_by="MSN.D1.0", point_alpha=0.5, point_size=1.2, ncol=5,
                   add_legend=T) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25))

    
sce.nac.noLNT$MSN.D2 <- NA
sce.nac.noLNT$MSN.D2[grep("MSN.D2", sce.nac.noLNT$cellType.final)] <- TRUE
    
sce.nac.noLNT$MSN.D2.0 <- NA
sce.nac.noLNT$MSN.D2.0 <- sce.nac.noLNT$MSN.D2 & geneExprs["DRD2", ]==0

# Now plot & color by this - is DRD1 expressed?
plotExpression(sce.nac.noLNT, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType.final", colour_by="MSN.D2.0", point_alpha=0.5, point_size=1.2, ncol=5,
               add_legend=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25))




# non-0 expression of BOTH? ===
table(geneExprs["DRD1", union(which(sce.nac.noLNT$MSN.D1), which(sce.nac.noLNT$MSN.D2))] > 0 &
        geneExprs["DRD2", union(which(sce.nac.noLNT$MSN.D1), which(sce.nac.noLNT$MSN.D2))] > 0)
    #FALSE  TRUE
    # 8017   819

sce.nac.noLNT$MSN.broad <- NA
sce.nac.noLNT$MSN.broad[grep("MSN", sce.nac.noLNT$cellType.final)] <- TRUE
    # the above is the same as
table(geneExprs["DRD1", which(sce.nac.noLNT$MSN.broad)] > 0 &
        geneExprs["DRD2", which(sce.nac.noLNT$MSN.broad)] > 0)


geneExprs.mat <- as.matrix(geneExprs)
sce.nac.noLNT$bothD1D2 <- ifelse(geneExprs.mat["DRD1", sce.nac.noLNT$MSN.broad] > 0 &
                                   geneExprs.mat["DRD2", sce.nac.noLNT$MSN.broad] > 0,
                                 TRUE, FALSE)


# Now plot & color by this - are these high-expressing in some??
plotExpression(sce.nac.noLNT, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType.final", colour_by="bothD1D2", point_alpha=0.5, point_size=1.2, ncol=5,
               add_legend=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25))


## Plot these three subsets: ===
pdf("pdfs/exploration/DRD1-vs-DRD2-exclusivity_all-NAc-n5_May202.pdf", height=6)
# D1-MSN but 0 DRD1
plotExpression(sce.nac.noLNT, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType.final", colour_by="MSN.D1.0", point_alpha=0.5, point_size=1.3, ncol=5,
               add_legend=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 12)) +
  ggtitle(label="Classified with D1 subclusters but 0 DRD1 expression")

# D2-MSN but 0 DRD2
plotExpression(sce.nac.noLNT, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType.final", colour_by="MSN.D2.0", point_alpha=0.5, point_size=1.3, ncol=5,
               add_legend=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 12)) +
  ggtitle(label="Classified with D2 subclusters but 0 DRD2 expression")

# MSN but non-0 expression for BOTH
plotExpression(sce.nac.noLNT, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType.final", colour_by="bothD1D2", point_alpha=0.5, point_size=1.3, ncol=5,
               add_legend=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 12)) +
  ggtitle(label="Classified as MSN (D1 or D2) but (+) expression of BOTH DRD1/DRD2")
dev.off()


# More exploration of these dual-expressing..??
table(sce.nac.noLNT$prelimCluster, sce.nac.noLNT$bothD1D2)
    # don't seem to make up or be enriched in their own cluster

table(sce.nac.noLNT$cellType.final, sce.nac.noLNT$bothD1D2)
    #          FALSE TRUE
    # Astro        0    0
    # Inhib.1      0    0
    # Inhib.2      0    0
    # Inhib.3      0    0
    # Inhib.4      0    0
    # Micro        0    0
    # MSN.D1.1   127    5
    # MSN.D1.2   292    9
    # MSN.D1.3   480  239   # however these ARE a little bit more enriched in this group
    # MSN.D1.4  3351  402
    # MSN.D2.1   294    6
    # MSN.D2.2  3473  158
    # Oligo        0    0
    # OPC          0    0     

table(sce.nac.noLNT$MSN.broad, sce.nac.noLNT$bothD1D2)
    #      FALSE TRUE
    # TRUE  8017  819





