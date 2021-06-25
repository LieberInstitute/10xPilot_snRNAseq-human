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

    ### Exploration: What are the 'ambig.glial's? ====
    ## For identification/annotation of those 'ambig.glial' clusters, check out some top markers:
    prelimMarkers <- lapply(markers.nac.t.pw, function(x){rownames(x)[x$FDR<0.05]})
    lengths(prelimMarkers)
    # ambig.glial_A ambig.glial_B       Astro_A       Astro_B       Inhib_A       Inhib_B 
    #            43           509           444           419           183           164 
    #       Inhib_C       Inhib_D       Inhib_E         Micro      MSN.D1_A      MSN.D1_B 
    #           137            80           206           547            12            22 
    #      MSN.D1_C      MSN.D1_D      MSN.D1_E      MSN.D1_F      MSN.D2_A      MSN.D2_B 
    #            26            32            55            76            21            15 
    #      MSN.D2_C      MSN.D2_D       Oligo_A       Oligo_B           OPC       OPC_COP 
    #           156           127           256            62           102           169

    # (ran interactively just for these 'ambig.glial' clusters)
    sapply(medianNon0.nac, table)
    #       ambig.glial_A ambig.glial_B
    # FALSE         29581         28549
    # TRUE             99          1131
    
    head(prelimMarkers[["ambig.glial_A"]],n=40)
        # [1] "AC063949.2" "AC002480.3" "LINC00958"  "ASGR2"      "AC008687.4" "AC048383.1"
        # [7] "AL139351.1" "AC007952.7" "FOSL1"      "AL359636.2" "OR3A2"      "SNX20"     
        # [13] "SMG8"       "TGIF2"      "AC080078.1" "SH2D2A"     "HIST1H1B"   "FABP12"    
        # [19] "LINCMD1"    "AC093627.7" "ZNF503"     "SENCR"      "TSPEAR-AS2" "ADGRF1"    
        # [25] "IL6"        "SAGE1"      "OSM"        "AL359979.1" "CHAT"       "AP000977.1"
        # [31] "AC006974.2" "LINC02422"  "C2orf91"    "ANKRD22"    "AL590764.1" "URAD"      
        # [37] "Z84488.2"   "GRK1"       "FOXN4"      "KLK15"     
            # None of these 99 non-0-median-expressing genes are pairwise-FDR<0.05 markers for this cluster...
            i <- "ambig.glial_A"
            table(markers.nac.t.pw[[i]]$FDR < 0.5 & markers.nac.t.pw[[i]]$non0median==TRUE)
                # Even none here...
            plotExpressionCustom(sce = sce.hold,
                                 features = head(prelimMarkers[[i]][-c(grep("^AC", prelimMarkers[[i]]),
                                                                       grep("^AL", prelimMarkers[[i]]))],n=12), 
                                 features_name = i,
                                 anno_name = "cellType",
                                 ncol=4, point_alpha=0.4) +
              scale_color_manual(values = cell_colors.nac)
            
      head(rownames(markers.nac.t.pw[[i]])[markers.nac.t.pw[[i]]$non0median==TRUE],n=12)
          # [1] "AC244021.1" "FP236383.1" "FP671120.1" "FRMD4A"     "SYNDIG1"    "BACH1"     
          # [7] "SSH2"       "SLC8A1"     "LDLRAD4"    "MAML2"      "PLXDC2"     "MAML3"
              # Plot these:
            plotExpressionCustom(sce = sce.hold,
                                 features = head(rownames(markers.nac.t.pw[[i]])[markers.nac.t.pw[[i]]$non0median==TRUE],n=12), 
                                 features_name = i,
                                 anno_name = "cellType",
                                 ncol=4, point_alpha=0.4) +
              scale_color_manual(values = cell_colors.nac)
            
              # -> Since it's pretty unclear from the literature as to what these might be, and since
              #    they're quite low in transcriptional activity, let's call 'Micro_resting'

            
    i <- "ambig.glial_B"
    head(prelimMarkers[[i]],n=40)
        # [1] "CD163"      "F13A1"      "MS4A4E"     "MS4A4A"     "LILRB5"     "TGFBI"     
        # [7] "MS4A6A"     "CD209"      "SIGLEC1"    "MRC1"       "MARCO"      "STAB1"     
        # [13] "HOXB-AS1"   "CD28"       "IQGAP2"     "FCGR2B"     "MS4A14"     "IFI44L"    
        # [19] "ADGRG6"     "CCL8"       "AP000812.1" "CASP4"      "MAFB"       "LINC01839" 
        # [25] "VSIG4"      "AP001636.3" "LINC00278"  "KIR2DL4"    "CD163L1"    "LYVE1"     
        # [31] "LY96"       "JAML"       "MPEG1"      "MMP2-AS1"   "GYPC"       "FPR3"      
        # [37] "MYO7A"      "CD200R1"    "TNFRSF14"   "FCMR"
    plotExpressionCustom(sce = sce.hold,
                         features = c("AIF1","SIGLEC1", "CD44","EGR3", "MPEG1","CD28",
                                      # Broad microglial markers (used for annotation)
                                      "CD74", "CSF1R"), 
                         features_name = i,
                         anno_name = "cellType",
                         ncol=4, point_alpha=0.4) +
      scale_color_manual(values = cell_colors.nac)
        # --> These are definitely a macrophage or microglia subtype; former: SIGLEC1 (aka CD169)
        #     - is an infiltrating-Mphage-specific marker: https://doi.org/10.1038/s41583-018-0057-5
        #     -> Call 'Macro_infilt'
    # End prelim marker exploration ======
    

## WMW: Blocking on sample (this test doesn't take 'design=' argument) ===
markers.nac.wilcox.block <- findMarkers(sce.nac, groups=sce.nac$cellType,
                                        assay.type="logcounts", block=sce.nac$donor, test="wilcox",
                                        direction="up", pval.type="all", full.stats=T)


sapply(markers.nac.wilcox.block, function(x){table(x$FDR<0.05)})
    ## none for about ~1/3 of these... none for Micros or Oligos...


## Binomial ===
markers.nac.binom.block <- findMarkers(sce.nac, groups=sce.nac$cellType,
                                       assay.type="logcounts", block=sce.nac$donor, test="binom",
                                       direction="up", pval.type="all", full.stats=T)


sapply(markers.nac.binom.block, function(x){table(x$FDR<0.05)})
    ## even worse than WMW... basically none

# Add respective 'non0median' column to the stats for each set of markers
for(i in names(markers.nac.t.pw)){
  markers.nac.t.pw[[i]] <- cbind(markers.nac.t.pw[[i]],
                                 medianNon0.nac[[i]][match(rownames(markers.nac.t.pw[[i]]),
                                                           names(medianNon0.nac[[i]]))])
  colnames(markers.nac.t.pw[[i]])[27] <- "non0median"
}

sapply(markers.nac.t.pw, function(x){table(x$FDR<0.05 & x$non0median == TRUE)["TRUE"]})


## Save all these for future reference
save(markers.nac.t.pw, #markers.nac.wilcox.block, markers.nac.binom.block,
     file="rdas/markers-stats_NAc-n5_findMarkers-SN-LEVEL_MNTApr2020.rda")


## Print these to PNGs
markerList.t <- lapply(markers.nac.t.pw, function(x){
    rownames(x)[x$FDR < 0.05]
  }
)

# Take top 40
genes.top40.t <- lapply(markerList.t, function(x){head(x, n=40)})

#dir.create("pdfs/exploration/NAc-n5-markers/")
for(i in names(genes.top40.t)){
  png(paste0("pdfs/exploration/NAc-n5-markers/NAc-all-n5_t-sn-level-top40markers-",i,"_logExprs_Apr2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.nac, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:14], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers: single-nucleus-level p.w. t-tests"))
  )
  dev.off()
}



### Cluster-vs-all single-nucleus-level iteration ======

## Load SCE with new info
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda",
     verbose=T)
    # sce.nac, chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac, ref.sampleInfo

# First drop "ambig.lowNtrxts" (93 nuclei)
sce.nac <- sce.nac[ ,sce.nac$cellType != "ambig.lowNtrxts"]
sce.nac$cellType <- droplevels(sce.nac$cellType)

# Drop genes with all 0's
sce.nac <- sce.nac[!rowSums(assay(sce.nac, "counts"))==0, ]  ## keeps 29236 genes


## Traditional t-test with design as in PB'd/limma approach ===
# mod <- with(colData(sce.nac), model.matrix(~ donor))
# mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`


#mod <- with(colData(sce.nac), model.matrix(~ processDate + donor))
    # Error in .ranksafe_qr(full.design) : design matrix is not of full rank
    # -> just try processDate bc it describes more var (at least at PB-pan-brain level)
mod <- with(colData(sce.nac), model.matrix(~ processDate))
mod <- mod[ ,-1]


markers.nac.t.1vAll <- list()
for(i in levels(sce.nac$cellType)){
  # Make temporary contrast
  sce.nac$contrast <- ifelse(sce.nac$cellType==i, 1, 0)
  # Test cluster vs. all
  markers.nac.t.1vAll[[i]] <- findMarkers(sce.nac, groups=sce.nac$contrast,
                                          assay.type="logcounts", design=mod, test="t",
                                          direction="up", pval.type="all", full.stats=T)
}

    ## Then, temp set of stats to get the standardized logFC
    temp.1vAll <- list()
    for(i in levels(sce.nac$cellType)){
      # Make temporary contrast
      sce.nac$contrast <- ifelse(sce.nac$cellType==i, 1, 0)
      # Test cluster vs. all
      temp.1vAll[[i]] <- findMarkers(sce.nac, groups=sce.nac$contrast,
                                     assay.type="logcounts", design=mod, test="t",
                                     std.lfc=TRUE,
                                     direction="up", pval.type="all", full.stats=T)
    }


## For some reason all the results are in the second List entry (first is always empty)

# Replace that empty slot with the entry with the actul stats
markers.nac.t.1vAll <- lapply(markers.nac.t.1vAll, function(x){ x[[2]] })
# Same for that with std.lfc
temp.1vAll <- lapply(temp.1vAll, function(x){ x[[2]] })

# Now just pull from the 'stats.0' DataFrame column
markers.nac.t.1vAll <- lapply(markers.nac.t.1vAll, function(x){ x$stats.0 })
temp.1vAll <- lapply(temp.1vAll, function(x){ x$stats.0 })

# Re-name std.lfc column and add to the first result
for(i in names(temp.1vAll)){
  colnames(temp.1vAll[[i]])[1] <- "std.logFC"
  markers.nac.t.1vAll[[i]] <- cbind(markers.nac.t.1vAll[[i]], temp.1vAll[[i]]$std.logFC)
  # Oh the colname is kept weird
  colnames(markers.nac.t.1vAll[[i]])[4] <- "std.logFC"
  # Then re-organize
  markers.nac.t.1vAll[[i]] <- markers.nac.t.1vAll[[i]][ ,c("logFC","std.logFC","log.p.value","log.FDR")]
}


## Let's save this along with the previous pairwise results
save(markers.nac.t.pw, markers.nac.t.1vAll,
     file="rdas/markers-stats_NAc-n5_findMarkers-SN-LEVEL_MNTApr2020.rda")



## Print these to pngs
markerList.t.1vAll <- lapply(markers.nac.t.1vAll, function(x){
    rownames(x)[x[ ,"log.FDR"] < log10(0.001)]
  }
)
genes.top40.t.1vAll <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t.1vAll)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/NAc-n5-markers/NAc-all-n5_t-sn-level_1vALL_top40markers-",i,"_logExprs_May2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.nac, exprs_values = "logcounts", features=genes.top40.t.1vAll[[i]],
                   x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:14], length(genes.top40.t.1vAll[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers: single-nucleus-level t-tests, cluster-vs-all-others"))
  )
  dev.off()
}



## How do these top 40 intersect? ===
sapply(names(genes.top40.t), function(c){
  length(intersect(genes.top40.t[[c]],
                   genes.top40.t.1vAll[[c]]))
})
    #     Astro  Inhib.1  Inhib.2  Inhib.3  Inhib.4    Micro MSN.D1.1 MSN.D1.2
    #        35       30       16       26       22       40       17       24
    #  MSN.D1.3 MSN.D1.4 MSN.D2.1 MSN.D2.2    Oligo      OPC
    #         7        9       23       16       32       29


## Write these top 40 lists to a csv
names(markerList.t) <- paste0(names(markerList.t),"_pw")
names(markerList.t.1vAll) <- paste0(names(markerList.t.1vAll),"_1vAll")

top40genes <- cbind(sapply(markerList.t, function(x) head(x, n=40)),
                    sapply(markerList.t.1vAll, function(y) head(y, n=40)))
top40genes <- top40genes[ ,sort(colnames(top40genes))]

write.csv(top40genes, file="tables/top40genesLists_NAc-n5_cellType.final_SN-LEVEL-tests_May2020.csv",
          row.names=FALSE)



## Make marker array for Supp figure (MNT suggested panel A) ===========
load("rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda", verbose=T)
    #sce.nac, chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac, ref.sampleInfo

load("rdas/markers-stats_NAc-n5_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.nac.t.pw, markers.nac.t.1vAll

# First make interneuron subset
sce.nac.int <- sce.nac[ ,grep("Inhib.", sce.nac$cellType)]
sce.nac.int$cellType <- droplevels(sce.nac.int$cellType)


# Take top four for 4 inhib. interneuron pops
topToPrint <- as.data.frame(sapply(markers.nac.t.1vAll, function(x) {head(rownames(x),n=4)}))
topToPrint <- topToPrint[grep("Inhib.", names(topToPrint))]

# Manual assignment, bc 'Inhib.2' is mostly driven by noise (but 0 median)
topToPrint["Inhib.2"] <- c("KCNJ6", "SDK1", "LRFN2","NR2F1-AS1")

table(unlist(topToPrint) %in% rownames(sce.nac.int)) # good

# Print
pdf("pdfs/pubFigures/suppFig_NAc_interneuron-marker-array_MNTSep2020.pdf", height=8, width=5.5)
print(
  plotExpression(sce.nac.int, exprs_values = "logcounts", features=c(t(topToPrint)),
                 x="cellType.final", colour_by="cellType.final", point_alpha=0.6, point_size=1.5, ncol=4,
                 add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                              geom = "crossbar", width = 0.3,
                                              colour=rep(tableau10medium[1:4], length(unlist(topToPrint)))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12), plot.title = element_text(size = 25)) +  
    ggtitle(label="Inhib.1          Inhib.2           Inhib.3          Inhib.4") + xlab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 13),
          axis.title.y = element_text(angle = 90, size = 16),
          plot.title = element_text(size = 15),
          panel.grid.major=element_line(colour="grey95", size=0.8),
          panel.grid.minor=element_line(colour="grey95", size=0.4))
)
dev.off()


## For MNT version panel B
pdf("pdfs/pubFigures/suppFig_NAc_interneuron-experiment-panelB_MNTSep2020.pdf", height=3.2, width=4.5)
print(
  plotExpression(sce.nac.int, exprs_values = "logcounts", features=c("GAD1", "KIT", "PTHLH", "PVALB"),
                 x="cellType.final", colour_by="cellType.final", point_alpha=0.7, point_size=1.2, ncol=2,
                 add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                              geom = "crossbar", width = 0.3,
                                              colour=rep(tableau10medium[1:4], 4)) +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9.5),
          axis.title.y = element_text(angle = 90, size = 10),
          panel.grid.major=element_line(colour="grey95", size=0.8),
          panel.grid.minor=element_line(colour="grey95", size=0.4))
)
dev.off()


## Supp Fig 5: Other markers pointed out in text ===
genes2print <- c("DRD1", "DRD2", "CASZ1", "GPR6", "EBF1", "GRM8")

# First drop "ambig.lowNtrxts" (93 nuclei)
sce.nac <- sce.nac[ ,sce.nac$cellType != "ambig.lowNtrxts"]
sce.nac$cellType <- droplevels(sce.nac$cellType)

pdf("pdfs/pubFigures/suppFig_NAc_other-MSN-markers_MNTSep2020.pdf", height=4.5, width=6.5)
print(
  plotExpression(sce.nac, exprs_values = "logcounts", features=genes2print,
                 x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=1.0, ncol=2,
                 add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                              geom = "crossbar", width = 0.3,
                                              colour=rep(tableau20[1:14], length(genes2print))) +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9.5),
          axis.title.y = element_text(angle = 90, size = 10),
          panel.grid.major=element_line(colour="grey95", size=0.8),
          panel.grid.minor=element_line(colour="grey95", size=0.4))
)
dev.off()



## Top markers for D1.4 / D2.2 often co-expressed ===
load("rdas/markers-stats_NAc-n5_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.nac.t.pw, markers.nac.t.1vAll

cell.idx <- splitit(sce.nac$cellType)
dat <- as.matrix(assay(sce.nac, "logcounts"))
genes <- head(rownames(markers.nac.t.1vAll[["MSN.D1.4"]]), n=40)

pdf('pdfs/pubFigures/suppFigure_heatmap-Exprs_NAc-n5_top40-D1.4markers-1vAlltest_MNTSep2020.pdf',
    useDingbats=TRUE, height=6, width=10)
current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMeans(dat[genes, ii])))
pheatmap(t(current_dat), cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(0.02, 5, length.out = 101),
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "OrRd"))(100),
         fontsize_row = 20, fontsize_col=15)
dev.off()






### MNT add 18Nov2020 =================================
# -> What if add param/requirement that for any given subcluster, median expression has to > 0?
load("rdas/markers-stats_NAc-n5_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.nac.t.pw, markers.nac.t.1vAll

## Load SCE 
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda",
     verbose=T)
    # sce.nac, chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac, ref.sampleInfo

table(sce.nac$cellType)

# First drop "ambig.lowNtrxts" (93 nuclei)
sce.nac <- sce.nac[ ,sce.nac$cellType != "ambig.lowNtrxts"]
sce.nac$cellType <- droplevels(sce.nac$cellType)

# Remove 0 genes across all nuclei
sce.nac <- sce.nac[!rowSums(assay(sce.nac, "counts"))==0, ]


## Make list of Boolean param / cell subtype ===
cellSubtype.idx <- splitit(sce.nac$cellType)
medianNon0.idx <- lapply(cellSubtype.idx, function(x){
  apply(as.matrix(assay(sce.nac, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

lengths(medianNon0.idx)
sapply(medianNon0.idx, head)

# Add these to the stats for each set of markers
for(i in names(markers.nac.t.1vAll)){
  markers.nac.t.1vAll[[i]] <- cbind(markers.nac.t.1vAll[[i]],
                                    medianNon0.idx[[i]][match(rownames(markers.nac.t.1vAll[[i]]),
                                                              names(medianNon0.idx[[i]]))])
  colnames(markers.nac.t.1vAll[[i]])[5] <- "non0median"
}


## Use these restrictions to print (to png) a refined top 40, as before ===
markerList.t.1vAll <- lapply(markers.nac.t.1vAll, function(x){
  rownames(x)[x$log.FDR < log10(0.000001) & x$non0median==TRUE]
  }
)
    # lengths(markerList.t.1vAll)     # ( **without $non0median==TRUE restriction )
        #    Astro  Inhib.1  Inhib.2  Inhib.3  Inhib.4    Micro MSN.D1.1 MSN.D1.2
        #     4956     1295     1589     4917     4719     3656     1474     2468
        # MSN.D1.3 MSN.D1.4 MSN.D2.1 MSN.D2.2    Oligo      OPC
        #     2434     5893     2027     3106     2248     3259

lengths(markerList.t.1vAll)
    #    Astro  Inhib.1  Inhib.2  Inhib.3  Inhib.4    Micro MSN.D1.1 MSN.D1.2
    #     1214      690      561     1953     1716      763      700      989
    # MSN.D1.3 MSN.D1.4 MSN.D2.1 MSN.D2.2    Oligo      OPC
    #     1173     3054     1069     1912      700     1086

genes.top40.t <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/NAc-n5-markers/NAc_t-sn-level_1vALL_top40markers-REFINED-",gsub(":",".",i),"_logExprs_Nov2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.nac, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:14], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers, refined: single-nucleus-level p.w. t-tests, cluster-vs-all"))
  )
  dev.off()
}



## Do the same with the pairwise result ('markers.nac.t.pw') === === ===
# Add these to the stats for each set of markers
for(i in names(markers.nac.t.pw)){
  markers.nac.t.pw[[i]] <- cbind(markers.nac.t.pw[[i]],
                                     medianNon0.idx[[i]][match(rownames(markers.nac.t.pw[[i]]),
                                                               names(medianNon0.idx[[i]]))])
  colnames(markers.nac.t.pw[[i]])[16] <- "non0median"
}

markerList.t <- lapply(markers.nac.t.pw, function(x){
  rownames(x)[x$FDR < 0.05 & x$non0median==TRUE]
  }
)
    # lengths(markerList.t)     # ( **without $non0median==TRUE restriction )
        #   Astro  Inhib.1  Inhib.2  Inhib.3  Inhib.4    Micro MSN.D1.1 MSN.D1.2
        #     831      340      179      256      283     1417       82      295
        #MSN.D1.3 MSN.D1.4 MSN.D2.1 MSN.D2.2    Oligo      OPC
        #      87       61      230       66      493      372

lengths(markerList.t)
    #   Astro  Inhib.1  Inhib.2  Inhib.3  Inhib.4    Micro MSN.D1.1 MSN.D1.2
    #     378      110        9      145      105      402       45      201
    #MSN.D1.3 MSN.D1.4 MSN.D2.1 MSN.D2.2    Oligo      OPC
    #      64       54      143       57      333      187


genes.top40.t <- lapply(markerList.t, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/NAc-n5-markers/NAc_t-sn-level_pairwise_top40markers-REFINED-", i, "_logExprs_Nov2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.nac, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:14], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers, refined: single-nucleus-level p.w. t-tests"))
  )
  dev.off()
}

## Then write a new CSV of these refined top 40 genes ===
names(markerList.t) <- paste0(names(markerList.t),"_pw")
names(markerList.t.1vAll) <- paste0(names(markerList.t.1vAll),"_1vAll")

# Many of the PW results don't have at least 40 markers:
extend.idx <- names(which(lengths(markerList.t) < 40))
for(i in extend.idx){
  markerList.t[[i]] <- c(markerList.t[[i]], rep("", 40-length(markerList.t[[i]])))
}

top40genes <- cbind(sapply(markerList.t, function(x) head(x, n=40)),
                    sapply(markerList.t.1vAll, function(y) head(y, n=40)))
top40genes <- top40genes[ ,sort(colnames(top40genes))]

write.csv(top40genes, file="tables/top40genesLists-REFINED_NAc-n5_cellType.final_Nov2020.csv",
          row.names=FALSE)




## Aside: add in 't.stat' as in 'step04' analyses to save for LoHu/LeCo ===
for(s in names(markers.nac.t.1vAll)){
  markers.nac.t.1vAll[[s]]$t.stat <- markers.nac.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.nac))
}

save(markers.nac.t.1vAll, markers.nac.t.pw, sce.nac,
     file="rdas/markerStats-and-SCE_NAc-n5_sn-level_cleaned_MNTNov2020.rda")


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




