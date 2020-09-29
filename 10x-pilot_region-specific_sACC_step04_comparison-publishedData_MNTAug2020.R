### MNT 10x snRNA-seq workflow: step 04
###   **Region-specific analyses**
###     - (2x) sACC samples from: Br5161 & Br5212
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

sce.asd.pfc <- sce.asd[ ,sce.asd$region=="PFC"]
sce.asd.acc <- sce.asd[ ,sce.asd$region=="ACC"]


# Need to convert Symbol in sce.dlpfc > EnsemblID, and also use n nuclei for t.stat
load("rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    # sce.sacc, clusterRefTab.sacc, pc.choice.sacc, chosen.hvgs.sacc, ref.sampleInfo

# Drop "Ambig.lowNtrxts" (43 nuclei)
sce.sacc <- sce.sacc[ ,sce.sacc$cellType != "Ambig.lowNtrxts"]
sce.sacc$cellType <- droplevels(sce.sacc$cellType)

## Load LIBD sACC stats (just need the "1vAll" result)
load("rdas/markers-stats_sACC-n2_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.sacc.t.design, markers.sacc.wilcox.block, markers.sacc.binom.block, markers.sacc.t.1vAll
    rm(markers.sacc.t.design, markers.sacc.wilcox.block, markers.sacc.binom.block)

for(i in names(markers.sacc.t.1vAll)){
  rownames(markers.sacc.t.1vAll[[i]]) <- rowData(sce.sacc)$ID[match(rownames(markers.sacc.t.1vAll[[i]]),
                                                                    rownames(sce.sacc))]
}


## Calculate and add t-statistic (= std.logFC * sqrt(N)) from contrasts
#      and fix row order to the first entry "Astro"
fixTo <- rownames(markers.sacc.t.1vAll[["Astro"]])
for(s in names(markers.sacc.t.1vAll)){
  markers.sacc.t.1vAll[[s]]$t.stat <- markers.sacc.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.sacc))
  markers.sacc.t.1vAll[[s]] <- markers.sacc.t.1vAll[[s]][fixTo, ]
}

# Pull out the t's
ts.sacc <- sapply(markers.sacc.t.1vAll, function(x){x$t.stat})
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
    # of length 27,422

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

theSeq.all = seq(-.95, .95, by = 0.01)
my.col.all <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.all)-1)

pdf("pdfs/exploration/Velmeshev-ASD_pfc-acc/overlap-velmeshev-ASD-bothRegions_with_LIBD-10x-sACC_Aug2020.pdf", width=10)
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


