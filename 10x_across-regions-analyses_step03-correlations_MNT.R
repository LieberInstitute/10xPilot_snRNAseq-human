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

load("rdas/markers-stats_DLPFC_n2_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.dlpfc.t.1vAll, markers.dlpfc.t.design

load("rdas/markers-stats_HPC-n3_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.hpc.t.1vAll, markers.hpc.t.design, markers.hpc.wilcox.block

load("rdas/markers-stats_NAc-n5_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.nac.t.design, markers.nac.t.1vAll

load("rdas/markers-stats_Amyg-n2_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.amy.t.1vAll, markers.amy.t.design, markers.amy.wilcox.block

load("rdas/markers-stats_sACC-n2_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.sacc.t.design, markers.sacc.wilcox.block, markers.sacc.binom.block, markers.sacc.t.1vAll

    # Remove all pairwise stats
    rm(list=setdiff(ls(), ls(pattern=".1vAll")))


# Re-order rows in each list entry for each set of stats
expressedGenes.list <- list(amy=rownames(markers.amy.t.1vAll[["Astro"]]),
                            dlpfc=rownames(markers.dlpfc.t.1vAll[["Astro"]]),
                            hpc=rownames(markers.hpc.t.1vAll[["Astro"]]),
                            nac=rownames(markers.nac.t.1vAll[["Astro"]]),
                            sacc=rownames(markers.sacc.t.1vAll[["Astro"]])
                            )

expressedGenes <- unique(unlist(expressedGenes.list))

expressedGenes <- expressedGenes[expressedGenes %in% expressedGenes.list[["amy"]] &
                                   expressedGenes %in% expressedGenes.list[["dlpfc"]] &
                                   expressedGenes %in% expressedGenes.list[["nac"]] &
                                   expressedGenes %in% expressedGenes.list[["hpc"]] &
                                   expressedGenes %in% expressedGenes.list[["sacc"]]
                                 ]
length(expressedGenes)  # 26888

# Store each set of stats into a list
FMstats.list <- list(amy=markers.amy.t.1vAll,
                     dlpfc=markers.dlpfc.t.1vAll,
                     hpc=markers.hpc.t.1vAll,
                     nac=markers.nac.t.1vAll,
                     sacc=markers.sacc.t.1vAll)

    # How many subclusters in each region?
    sapply(FMstats.list, length)
        # amy dlpfc   hpc   nac  sacc
        #  12    17    15    14    10
    
    sapply(FMstats.list, function(x){nrow(x[["Astro"]])})
        #  amy dlpfc   hpc   nac  sacc
        #28464 28111 28757 29236 28774

# Subset and re-order for those intersecting genes across all regions
for(x in names(FMstats.list)){
  for(s in names(FMstats.list[[x]])){
    FMstats.list[[x]][[s]] <- FMstats.list[[x]][[s]][expressedGenes, ]
  }
}


sapply(FMstats.list, function(x){nrow(x[["Astro"]])})
    # good
sapply(FMstats.list, function(x){head(x[["Astro"]], n=3)})







### Get n Nuclei numbers for each region so can compute t-statistics ===
  # This can be done with Cohen's D (the 'std.lfc'), as d = t/sqrt(N)

load("rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda", verbose=T)
    # sce.all.n12, chosen.hvgs.all.n12, pc.choice.n12, ref.sampleInfo, clusterRefTab.all.n12
    rm(chosen.hvgs.all.n12, pc.choice.n12, clusterRefTab.all.n12, ref.sampleInfo)

sce.all.n12
    # class: SingleCellExperiment
    # dim: 33538 34070    
table(sce.all.n12$cellType.RS)

# Get rid of those which were 'Ambig.lowNtrxts' bc removed before cluster-specificity tests run
sce.all.n12 <- sce.all.n12[ ,!sce.all.n12$cellType.RS=="Ambig.lowNtrxts"]
sampleNumNuclei <- table(sce.all.n12$region)
# Remove NAc because these clusters were defined from all-n5-sample region analysis
sampleNumNuclei <- sampleNumNuclei[c("amy","dlpfc","hpc","sacc")]

# Get N nuclei from all-n=5-NAc dataset
load("rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda", verbose=T)
    # sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all, ref.sampleInfo
    rm(chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all)

sampleNumNuclei.nac <- table(sce.nac.all$region) - table(sce.nac.all$cellType.final=="ambig.lowNtrxts")["TRUE"]
    #   nac
    # 13148

sampleNumNuclei <- c(sampleNumNuclei, sampleNumNuclei.nac)
    #  amy dlpfc   hpc  sacc   nac
    # 6582  5231 10343  7004 13148


## Calculate and add t-statistic (= std.logFC * sqrt(N))
for(x in names(FMstats.list)){
  for(s in names(FMstats.list[[x]])){
    FMstats.list[[x]][[s]]$t.stat <- FMstats.list[[x]][[s]]$std.logFC * sqrt(sampleNumNuclei[x])
  }
}


## Let's save these
readme.mnt <- "These stats are from region-specific specificity modeling (cluster-vs-all-others) at the single-nucleus level with 'scran::findMarkers()'. The t-statistic is computed by sqrt(N.nuclei) * std.logFC."
save(FMstats.list, sampleNumNuclei, readme.mnt, ref.sampleInfo,
     file="rdas/markers-stats_SN-LEVEL_1vAll_all-regions-combined_May2020.rda")


# (If needed)
load("rdas/markers-stats_SN-LEVEL_1vAll_all-regions-combined_May2020.rda", verbose=T)
    # FMstats.list, sampleNumNuclei, readme.mnt, ref.sampleInfo
readme.mnt

## Create matrix of t's with region:subcluster identifiers
ts.list <- lapply(FMstats.list, function(x){
  sapply(x, function(y){y$t.stat})
  }
)
# Add back in row names and region suffix
for(i in names(ts.list)){
  rownames(ts.list[[i]]) <- rownames(FMstats.list[[1]][[1]])
  colnames(ts.list[[i]]) <- paste0(colnames(ts.list[[i]]), "_", i)
}
# Cbind
ts.fullMat <- do.call(cbind, ts.list)


## Correlation; first shorten names
colnames(ts.fullMat) <- gsub("Excit", "Ex", colnames(ts.fullMat))
colnames(ts.fullMat) <- gsub("Inhib", "In", colnames(ts.fullMat))
cor_t_panBrain <- cor(ts.fullMat)



### Heatmap - typically use levelplot (e.g. below), but will want pheatmap bc can cluster cols/rows
theSeq.all = seq(-1, 1, by = 0.025)
my.col.all <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq.all)-1)


pdf("pdfs/exploration/pan-Brain_correlation_region-specific-subcluster-ts_May2020.pdf")
pheatmap(cor_t_panBrain,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=6.5, fontsize_col=6.5,
         main="Correlation of cluster-specific t's from all regions \n (all shared expressed genes)")
dev.off()


## Subset on neuronal subcluster t's and check
ts.fullMat.neu <- ts.fullMat
for(i in c("Astro", "Micro", "Oligo", "OPC", "Tcell")){
  ts.fullMat.neu <- ts.fullMat.neu[ ,-grep(i, colnames(ts.fullMat.neu))]
}


colnames(ts.fullMat.neu) <- gsub("Excit", "Ex", colnames(ts.fullMat.neu))
colnames(ts.fullMat.neu) <- gsub("Inhib", "In", colnames(ts.fullMat.neu))
cor_t_panBrain.neu <- cor(ts.fullMat.neu)


### Heatmap - typically use levelplot (e.g. below), but will want pheatmap bc can cluster cols/rows
theSeq.all = seq(-1, 1, by = 0.025)
my.col.all <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq.all)-1)

# pdf("pdfs/exploration/pan-Brain_correlation_region-specific-NeuronalSubcluster-ts_May2020.pdf")
# pheatmap(cor_t_panBrain.neu,
#          color=my.col.all,
#          breaks=theSeq.all,
#          fontsize_row=6.5, fontsize_col=6.5,
#          main="Correlation of neuronal cluster-specific t's from all regions \n (all shared expressed genes)")
# dev.off()

## Updated 31Aug2020 - for paper === === ===
# Add some cluster info for add'l heatmap annotations
clusterInfo <- data.frame(class=ss(colnames(ts.fullMat.neu), "\\.", 1),
                          region=ss(colnames(ts.fullMat.neu), "_",2))
rownames(clusterInfo) <- colnames(ts.fullMat.neu)


# Print
pdf("pdfs/pubFigures/pan-Brain_correlation_region-specific-NeuronalSubcluster-ts_Aug2020.pdf")
pheatmap(cor_t_panBrain.neu,
         annotation_col=clusterInfo,
         show_colnames=FALSE,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=6.5, #fontsize_col=6.5,
         main="Correlation of neuronal cluster-specific t's from all regions \n (all shared expressed genes)")
dev.off()



## Non-neuronal set for supplement === === ===
load("rdas/markers-stats_SN-LEVEL_1vAll_all-regions-combined_May2020.rda", verbose=T)
    # FMstats.list, sampleNumNuclei, readme.mnt, ref.sampleInfo

ts.list <- lapply(FMstats.list, function(x){
  sapply(x, function(y){y$t.stat})
  }
)
# Add back in row names and region suffix
for(i in names(ts.list)){
  rownames(ts.list[[i]]) <- rownames(FMstats.list[[1]][[1]])
  colnames(ts.list[[i]]) <- paste0(colnames(ts.list[[i]]), "_", i)
}
# Cbind
ts.fullMat <- do.call(cbind, ts.list)


## Subset OUT neuronal subcluster t's and check
ts.fullMat.non <- ts.fullMat
#keepThese <- c("Astro", "Micro", "Oligo", "OPC", "Tcell")
ts.fullMat.non <- ts.fullMat.non[ ,-grep("Excit.", colnames(ts.fullMat.non))]
ts.fullMat.non <- ts.fullMat.non[ ,-grep("Inhib.", colnames(ts.fullMat.non))]
ts.fullMat.non <- ts.fullMat.non[ ,-grep("MSN.", colnames(ts.fullMat.non))]

cor_t_panBrain.non <- cor(ts.fullMat.non)


### Heatmap ===
theSeq.all = seq(-1, 1, by = 0.025)
my.col.all <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq.all)-1)

# Add some cluster info for add'l heatmap annotations
clusterInfo <- data.frame(class=ss(colnames(ts.fullMat.non), "\\_", 1),
                          region=ss(colnames(ts.fullMat.non), "_",2))
rownames(clusterInfo) <- colnames(ts.fullMat.non)

annotColors <- list(class = tableau10medium[1:5][factor(unique(clusterInfo$class))],
                    region = tableau10medium[6:10][factor(unique(clusterInfo$region))])

names(annotColors[["class"]]) <- unique(clusterInfo$class)
names(annotColors[["region"]]) <- unique(clusterInfo$region)



# Print
pdf("pdfs/pubFigures/pan-Brain_correlation_region-specific-NON-NeuronalSubcluster-ts_Sep2020.pdf",width=9)
pheatmap(cor_t_panBrain.non,
         annotation_col=clusterInfo,
         annotation_colors=annotColors,
         show_colnames=FALSE,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=10, #fontsize_col=6.5,
         display_numbers=TRUE, fontsize_number=6,
         main="Correlation of glia/other cluster-specific t's from all regions \n (all shared expressed genes)")
dev.off()


