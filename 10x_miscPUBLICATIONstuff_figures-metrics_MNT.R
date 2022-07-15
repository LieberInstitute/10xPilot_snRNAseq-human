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

source("plotExpressionCustom.R")

### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

#VIRIDFUN(name = colour_by_name, discrete = TRUE)

### DLPFC (for ST) ========================================================================
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


## Re-make neuron-specific t-SNE - 17Jun2020 (for R01 Visium atlas grant) ==============
# ## Aside - for BoG poster
# pdf("pdfs/pubFigures/DLPFC-n2_tSNE-cellType.split_MNTApr2020.pdf", width=9)
# plotTSNE(sce.dlpfc.st, colour_by="cellType.split", point_size=5.0, point_alpha=0.5,
#          text_by="cellType", text_size=7, theme_size=18)
# dev.off()

## Load SCE with new '.ST' annotations
load("rdas/regionSpecific_DLPFC-n2_SCE_cellTypesSplit-fromST_Apr2020.rda", verbose=T)
# sce.dlpfc.st, clusterRefTab.dlpfc, chosen.hvgs.dlpfc, ref.sampleInfo

# Add those annotations to the corresponding barcodes in the .neu sub
sce.dlpfc.neu$cellType.split <- sce.dlpfc.st$cellType.split[match(colnames(sce.dlpfc.neu),
                                                                  colnames(sce.dlpfc.st))]

sce.dlpfc.neu$cellType.split <- droplevels(sce.dlpfc.neu$cellType.split)

# Plot tSNE that had generated above
plotTSNE(sce.dlpfc.neu, colour_by="cellType.split", point_size=4.0, point_alpha=0.5,
         text_by="cellType.split", text_size=6, theme_size=16)

# * Will need to adjust text labels due to some slight overplotting
# Adapted from scater::plotReducedDim():
# Hidden function needed for this chunk ====
.coerce_to_factor <- function(x, level.limit, msg) {
  if (!is.null(x)) {
    x <- as.factor(x)
    if (nlevels(x) > level.limit) {
      stop(sprintf("more than %i levels for '%s'", level.limit, msg))
    }
  }
  x
}
# ====

text_by <- "cellType.split"
text_out <- retrieveCellInfo(sce.dlpfc.neu, text_by, search="colData")
text_out$val <- .coerce_to_factor(text_out$val, level.limit=Inf)
## actually not necessary if the colData chosen (usually cellType[.etc] is factorized)
df_to_plot <- data.frame(reducedDim(sce.dlpfc.neu, "TSNE"))
by_text_x <- vapply(split(df_to_plot$X1, text_out$val), median, FUN.VALUE=0)
by_text_y <- vapply(split(df_to_plot$X2, text_out$val), median, FUN.VALUE=0)

# This should recreate the automatic:
plotTSNE(sce.dlpfc.neu, colour_by="cellType.split", point_size=4.5, point_alpha=0.5,
         text_size=8, theme_size=18) +
  annotate("text", x=by_text_x, y=by_text_y, 
           label=names(by_text_x), size=6) +
  ggtitle("t-SNE on DLPFC neuronal nuclei (n=961)")

# OR
sce.dlpfc.neu$labels <- ifelse(!duplicated(sce.dlpfc.neu$cellType.split), as.character(sce.dlpfc.neu$cellType.split), NA)
Labs.df <- data.frame(by_text_x, by_text_y, labs=names(by_text_x))

colDF <- data.frame(colData(sce.dlpfc.neu))
DFforLabs <- cbind(reducedDim(sce.dlpfc.neu,"TSNE"), data.frame(colDF$labels))
colnames(DFforLabs) <- c("X","Y","labels")

# -> can replace those X,Y with the median positions for those labels?

DFforLabs.edit <- DFforLabs
DFforLabs.edit$X[!is.na(DFforLabs$labels)] <- by_text_x[match(as.character(DFforLabs$labels[!is.na(DFforLabs$labels)]),
                                                              names(by_text_x))]
DFforLabs.edit$Y[!is.na(DFforLabs$labels)] <- by_text_y[match(as.character(DFforLabs$labels[!is.na(DFforLabs$labels)]),
                                                              names(by_text_y))]

## Finally print
library(ggrepel)

pdf("pdfs/pubFigures/DLPFC-n2_tSNE_neuronalBarcodes_STregisteredLabs_MNTJun2020.pdf", width=8)
set.seed(109)
plotTSNE(sce.dlpfc.neu, colour_by="cellType.split", point_size=6, point_alpha=0.5,
         theme_size=18) +
  geom_text_repel(data=DFforLabs.edit, size=6.0,
                  aes(label=labels)) +
  ggtitle("t-SNE on DLPFC neuronal nuclei (n=961)")
dev.off()



## =====







## For BrBa-AnJa, et al (PTSD) ===========
# print: c("NPY", "CORT", "CRHBP", "DLL3", "NXPH2", "SST")
#        c("FERMT3", "CRHBP", "FOLR2", "PTGS1", "SLCO1C1", "P2RY13", "GLT8D2")
load("rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
# sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo

printThese <- list('first' = c("NPY", "CORT", "CRHBP", "DLL3", "NXPH2", "SST"),
                   'second' = c("FERMT3", "CRHBP", "FOLR2", "PTGS1", "SLCO1C1", "P2RY13", "GLT8D2"))

# First remove "Ambig.lowNtrxts" (50 nuclei):
sce.amy <- sce.amy[ ,sce.amy$cellType.split != "Ambig.lowNtrxts"]
sce.amy$cellType.split <- droplevels(sce.amy$cellType.split)

pdf("pdfs/pubFigures/zforBrBa-AnJa_AMY-cellType.split_requestedExprsPlots_May2020.pdf", height=6, width=8)
for(i in 1:length(printThese)){
  print(
    plotExpression(sce.amy, exprs_values = "logcounts", features=c(printThese[[i]]),
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7,
                   ncol=3, add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                        geom = "crossbar", width = 0.3,
                                                        colour=rep(tableau20[1:12], length(printThese[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+  
    #ggtitle(label=paste0(names(markers.mathys.custom)[i], " markers"))
  )
}
dev.off()

pdf("pdfs/pubFigures/zforBrBa-AnJa_AMY-cellType.split_CORT-Exprs_May2020.pdf", height=4, width=6)
print(
  plotExpression(sce.amy, exprs_values = "logcounts", features="CORT",
                 x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7,
                 add_legend=F, theme_size=18) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                             geom = "crossbar", width = 0.3,
                                                             colour=rep(tableau20[1:12], 1)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) #+  
  #ggtitle(label=paste0(names(markers.mathys.custom)[i], " markers"))
)
dev.off()

sum(assay(sce.amy, "counts")["CORT", ]) # 0
as.data.frame(prop.table(table(sce.amy$cellType.split)))
#       Var1        Freq
# 1    Astro 0.129443938
# 2  Excit.1 0.050744455
# 3  Excit.2 0.006077180
# 4  Excit.3 0.008356123
# 5  Inhib.1 0.025979945
# 6  Inhib.2 0.016560316
# 7  Inhib.3 0.005317533
# 8  Inhib.4 0.003646308
# 9  Inhib.5 0.014889091
# 10   Micro 0.116074142
# 11   Oligo 0.527651170
# 12     OPC 0.095259799


## For BrBa-AnJa, et al (PTSD) REVISION ===========
# print: c("NPY", "CORT", "CRHBP", "DLL3", "NXPH2", "SST")
#        c("FERMT3", "CRHBP", "FOLR2", "PTGS1", "SLCO1C1", "P2RY13", "GLT8D2")
load("rdas/revision/regionSpecific_Amyg-n5_cleaned-combined_SCE_MNT2021.rda", verbose=T)
# sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo

sce.amy = sce.amy[,!grepl("drop", sce.amy$cellType)]
sce.amy$cellType = droplevels(sce.amy$cellType)

printThese <- list('first' = c("NPY", "CORT", "CRHBP", "DLL3", "NXPH2", "SST"),
                   'second' = c("FERMT3", "CRHBP", "FOLR2", "PTGS1", "SLCO1C1", "P2RY13", "GLT8D2"))

pdf("pdfs/pubFigures/zforBrBa-AnJa_AMY-cellType_requestedExprsPlots_Nov2021.pdf", height=6, width=8)
for(i in 1:length(printThese)){
  print(
    plotExpression(sce.amy, exprs_values = "logcounts", features=c(printThese[[i]]),
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,
                   ncol=3, add_legend=F) + stat_summary(fun = median, fun.min = median, fun.max = median,
                                                        geom = "crossbar", width = 0.3,
                                                        colour="black") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))  
  )
}
dev.off()

pdf("pdfs/pubFigures/zforBrBa-AnJa_AMY-cellType_CORT-Exprs_Nov2021.pdf", height=4, width=6)
print(
  plotExpression(sce.amy, exprs_values = "logcounts", features="CORT",
                 x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,
                 add_legend=F, theme_size=18) + stat_summary(fun = median, fun.min = median, fun.max = median,
                                                        geom = "crossbar", width = 0.3,
                                                        colour="black") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20))  
)
dev.off()

sum(assay(sce.amy, "counts")["CORT", ]) # 0


### METRICS ====================================================

## Metrics and stuff for paper, all n=14 samples that went into analyses ===
load("rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda", verbose=T)
# sce.all.n12, ...
load("rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda", verbose=T)
# sce.nac.all, ...

# Get only NeuNs to add to 'sce.all.n12'
sce.neuns <- sce.nac.all[ ,sce.nac.all$protocol == "Frank.NeuN"]
# Trick 'sce.neuns' to keep the column which has regionally-defined annotations
sce.neuns$cellType.RS.sub <- sce.neuns$cellType
# Clean up
sce.neuns$cellType.RS.sub <- gsub("Excit", "Ex", sce.neuns$cellType.RS.sub)
sce.neuns$cellType.RS.sub <- gsub("Inhib", "In", sce.neuns$cellType.RS.sub)
sce.neuns$cellType.RS.sub <- gsub("ambig", "Ambig", sce.neuns$cellType.RS.sub)
sce.neuns$cellType.RS.sub <- paste0(sce.neuns$cellType.RS.sub,"_nac")

intersectingMetrics <- intersect(colnames(colData(sce.neuns)), colnames(colData(sce.all.n12)))
colData(sce.neuns) <- colData(sce.neuns)[ ,intersectingMetrics]
colData(sce.all.n12) <- colData(sce.all.n12)[ ,intersectingMetrics]

# Also have to remove reducedDims to do this
reducedDims(sce.all.n12) <- NULL
reducedDims(sce.neuns) <- NULL

sce.n14 <- cbind(sce.neuns, sce.all.n12)
# 42,763 nuclei

# Ok what's the median UMI and gene coverage across all 
quantile(sce.n14$sum)
#    0%      25%      50%      75%     100%
# 101.0   5279.5   8747.0  19895.0 196431.0
quantile(sce.n14$detected)
# 0%   25%   50%   75%  100%
# 95  2224  3047  5359 11838

# What does this coverage look like across cell types?
cellType.idx <- splitit(sce.n14$cellType)
sapply(cellType.idx, function(x){quantile(sce.n14[ ,x]$detected)})
# it's similiar to what we see in total transcript capture ('$sum')
sapply(cellType.idx, function(x){quantile(sce.n14[ ,x]$sum)})



## Create compiled file for Cell Ranger metrics - MNT 29Jul2021 ====
n14.samples <- c("Br5161_Amy", "Br5182_NAc_NeuN_reseq", "Br5212_DLPFC", "Br5212_sACC", "Br5287_HPC",
                 "Br5161_DLPFC", "Br5161_NAc", "Br5207_NAc_NeuN", "Br5212_HPC", "Br5287_NAc",
                 "Br5161_HPC", "Br5161_sACC", "Br5212_Amy", "Br5212_NAc")

n14.metrics <- list()
for(i in n14.samples){
  tempFile <- paste0("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/",i,"/outs/metrics_summary.csv")
  n14.metrics[[i]] <- read.csv(tempFile, header=TRUE)
}

n14.metrics.collapsed <- n14.metrics[[1]]
sharedMetrics <- colnames(n14.metrics.collapsed)
for(i in c(2:14)){
  n14.metrics.collapsed <- rbind(n14.metrics.collapsed, n14.metrics[[i]][ ,sharedMetrics])
}
rownames(n14.metrics.collapsed) <- n14.samples


## And revision n=10
samples.revision <- read.table("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/FASTQ/Feb2021/samples.manifest",
                               sep="\t", header=F)$V1
    #  [1] "Br5276_sACC_neun" "Br5400_NAc"       "Br5276_NAc"       "Br5701_NAc_neun" 
    #  [5] "Br5701_sACC_neun" "Br5207_DLPFC"     "Br5276_Amy_neun"  "Br5400_Amy_neun" 
    #  [9] "Br5400_sACC"      "Br5701_Amy"
rev.metrics <- list()
for(i in samples.revision){
  tempFile <- paste0("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/Feb2021/",i,"/outs/metrics_summary.csv")
  rev.metrics[[i]] <- read.csv(tempFile, header=TRUE)
}

rev.metrics.collapsed <- rev.metrics[[1]]
sharedMetrics <- colnames(rev.metrics.collapsed)
for(i in c(2:10)){
  rev.metrics.collapsed <- rbind(rev.metrics.collapsed, rev.metrics[[i]][ ,sharedMetrics])
}
rownames(rev.metrics.collapsed) <- samples.revision

colnames(rev.metrics.collapsed) == colnames(n14.metrics.collapsed)  # all TRUE

cellRanger.metrics <- rbind(n14.metrics.collapsed, rev.metrics.collapsed)

# Write this out for future reference
write.table(cellRanger.metrics, row.names=T, col.names=T, sep="\t", quote=F,
            file="tables/revision/METRICS-n24-analyzed_CellRanger-premRNA-output_MNT2021.tsv")
# As semi-colon-SV
write.table(cellRanger.metrics, row.names=T, col.names=T, sep=";", quote=F, 
            file="tables/revision/METRICS-n24-analyzed_CellRanger-premRNA-output_MNT2021.csv")


# Sequencing depth
quantile(as.numeric(gsub(",","",cellRanger.metrics$Number.of.Reads)))
            #       0%       25%       50%       75%      100%
            #118752669 148651409 253012622 274925035 296198922  (formerly, n=14 preprint)
    # Revision n=24:
    #        0%       25%       50%       75%      100% 
    # 118752669 253695421 284257444 418976634 479112744


## Version with additional batch info
load("rdas/revision/all-n24-samples_across-regions-analyses_MNT2021.rda", verbose=T)
ref.sampleInfo <- ref.sampleInfo[-which(rownames(ref.sampleInfo)=="br5182.nac.neun"), ]

## Re-order, so can cbind
# Add those new 'de-identified' donor IDs
donorRef <- data.frame(donor=unique(sce.allRegions$donor),
                       pub.donorID=paste0("donor", 1:8))

ref.sampleInfo$DonorID.pub <- donorRef$pub.donorID[match(ref.sampleInfo$donor,
                                                         donorRef$donor)]
ref.sampleInfo <- ref.sampleInfo[match(gsub("_",".",tolower(rownames(cellRanger.metrics))),
                                               rownames(ref.sampleInfo)), ]
rownames(ref.sampleInfo) <- rownames(cellRanger.metrics)
colnames(ref.sampleInfo) <- c("SampleID","Region","DonorID.BrNum","Sex","Process.Batch",
                              "Protocol","Sequencer","DonorID.pub")
# Re-organize (don't need 'sampleID' bc redundant)
ref.sampleInfo <- ref.sampleInfo[ ,c("DonorID.pub","DonorID.BrNum","Sex","Region",
                                     "Process.Batch","Protocol","Sequencer")]

# Add FACS machine
ref.sampleInfo$Sorter <- "BD_FACSAriaII"
ref.sampleInfo$Sorter[rownames(ref.sampleInfo) %in% samples.revision] <- "BioRad_S3e"
ref.sampleInfo <- ref.sampleInfo[ ,c(1:6,8,7)]


sampleMetrics.full <- cbind(ref.sampleInfo, cellRanger.metrics)


# Write this out for future reference
write.table(sampleMetrics.full, row.names=T, col.names=T, sep="\t", quote=F,
            file="tables/revision/METRICS-n24_CellRanger-premRNA-output_wBatchInfo_MNT2021.tsv")
# As semi-colon-SV
write.table(sampleMetrics.full, row.names=T, col.names=T, sep=";", quote=F, 
            file="tables/revision/METRICS-n24_CellRanger-premRNA-output_wBatchInfo_MNT2021.csv")




## NAc, all n=8 ===
load("rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda", verbose=T)
    # sce.nac, chosen.hvgs.nac, pc.choice.nac, ref.sampleInfo, annotationTab.nac, cell_colors.nac

## Distribution of n transcripts
# 'Depth' per sample
nac.samples <- splitit(sce.nac$sample)
sapply(nac.samples, function(x){quantile(sce.nac[ ,x]$sum)})
# 'Capture' per cell type
cellType_idx <- splitit(sce.nac$cellType)
sapply(cellType_idx, function(x){quantile(sce.nac[ ,x]$sum)})
# Median n UMIs per nucleus
quantile(sce.nac$sum, probs=seq(0.1,1,by=0.1))
#   10%   20%   30%   40%   50%   60%   70%   80%   90%  100%
#  4472  7308 11360 16017 19249 21819 24411 27653 32519 85460

# n Nuclei from each sample
table(sce.nac$sample)
# nac.5161      nac.5212      nac.5287 nac.neun.5182 nac.neun.5207
#     2067          1774           707          4267          4426

signif(sum(assay(sce.nac, "counts")),3)
#[1] 2.51e+08



# Clean up '.drop' clusters (669 nuclei, total)
sce.nac <- sce.nac[ ,-grep("drop.", sce.nac$cellType)]
sce.nac$cellType <- droplevels(sce.nac$cellType)

## Proportion of cell types
table(sce.nac$cellType)
round(prop.table(table(sce.nac$cellType)),3)
    #  Astro_A       Astro_B       Inhib_A       Inhib_B       Inhib_C 
    #    0.005         0.050         0.013         0.002         0.005 
    #  Inhib_D       Inhib_E  Macro_infilt         Micro Micro_resting 
    #    0.012         0.002         0.001         0.022         0.003 
    # MSN.D1_A      MSN.D1_B      MSN.D1_C      MSN.D1_D      MSN.D1_E 
    #    0.197         0.012         0.014         0.036         0.032 
    # MSN.D1_F      MSN.D2_A      MSN.D2_B      MSN.D2_C      MSN.D2_D 
    #    0.004         0.214         0.014         0.016         0.003 
    #  Oligo_A       Oligo_B           OPC       OPC_COP 
    #    0.050         0.259         0.033         0.001 

# How many MSNs in our dataset?
sum(prop.table(table(sce.nac$cellType))[grep("MSN", levels(sce.nac$cellType))])
    #[1] 0.5434345

table(sce.nac$cellType, sce.nac$donor)
    #               br5161 br5182 br5207 br5212 br5276 br5287 br5400 br5701
    # Astro_A           27      0      0      5      8      3     56      0
    # Astro_B          115      0      0    377    173      8    294     33
    # Inhib_A           10    101     58      8     29      4     26     15
    # Inhib_B            0     23      6      0      3      0      4      4
    # Inhib_C            4     33     40      2      4      3     11      1
    # Inhib_D            4     56    128      5      4      6     24     13
    # Inhib_E            1     16      5      3      6      0      6      0
    # Macro_infilt       4      0      0      7      2      1      4      4
    # Micro             66      0      0     59     33     34    222     15
    # Micro_resting      3      0      0     33     22      0      5      0
    # MSN.D1_A          96   1269   1680    129    419     35    258     41
    # MSN.D1_B           4    161     64      2      1      0      0      7
    # MSN.D1_C           8    224      4     22     21      0      2      2
    # MSN.D1_D          27    267    173     26    108      5    111      1
    # MSN.D1_E          15    260    255      6     15      6     30     51
    # MSN.D1_F           2     61     10      0      2      0      2      9
    # MSN.D2_A          95   1421   1819    101    488     29    252     57
    # MSN.D2_B           9    134     29      7     58      5     36      7
    # MSN.D2_C           9    130    131      3     15      3     18      5
    # MSN.D2_D           3     53      0      0      2      0      0      0
    # Oligo_A          237      0      0     50    385     69    247      0
    # Oligo_B         1202      0      0    804    523    421   2186     10
    # OPC               98      0      0    104    200     34    209      6
    # OPC_COP            0      0      0      0      1      3     14      0
        # ** br5182, 5207 & 5701 are the neun-enriched samples


## Proportions of cell type by sample
apply(table(sce.nac$cellType, sce.nac$donor), 2, function(x){round(prop.table(x),3) * 100})
    #               br5161 br5182 br5207 br5212 br5276 br5287 br5400 br5701
    # Astro_A          1.3    0.0    0.0    0.3    0.3    0.4    1.4    0.0
    # Astro_B          5.6    0.0    0.0   21.5    6.9    1.2    7.3   11.7
    # Inhib_A          0.5    2.4    1.3    0.5    1.1    0.6    0.6    5.3
    # Inhib_B          0.0    0.5    0.1    0.0    0.1    0.0    0.1    1.4
    # Inhib_C          0.2    0.8    0.9    0.1    0.2    0.4    0.3    0.4
    # Inhib_D          0.2    1.3    2.9    0.3    0.2    0.9    0.6    4.6
    # Inhib_E          0.0    0.4    0.1    0.2    0.2    0.0    0.1    0.0
    # Macro_infilt     0.2    0.0    0.0    0.4    0.1    0.1    0.1    1.4
    # Micro            3.2    0.0    0.0    3.4    1.3    5.1    5.5    5.3
    # Micro_resting    0.1    0.0    0.0    1.9    0.9    0.0    0.1    0.0
    # MSN.D1_A         4.7   30.1   38.2    7.4   16.6    5.2    6.4   14.6
    # MSN.D1_B         0.2    3.8    1.5    0.1    0.0    0.0    0.0    2.5
    # MSN.D1_C         0.4    5.3    0.1    1.3    0.8    0.0    0.0    0.7
    # MSN.D1_D         1.3    6.3    3.9    1.5    4.3    0.7    2.8    0.4
    # MSN.D1_E         0.7    6.2    5.8    0.3    0.6    0.9    0.7   18.1
    # MSN.D1_F         0.1    1.4    0.2    0.0    0.1    0.0    0.0    3.2
    # MSN.D2_A         4.7   33.8   41.3    5.8   19.3    4.3    6.3   20.3
    # MSN.D2_B         0.4    3.2    0.7    0.4    2.3    0.7    0.9    2.5
    # MSN.D2_C         0.4    3.1    3.0    0.2    0.6    0.4    0.4    1.8
    # MSN.D2_D         0.1    1.3    0.0    0.0    0.1    0.0    0.0    0.0
    # Oligo_A         11.6    0.0    0.0    2.9   15.3   10.3    6.1    0.0
    # Oligo_B         59.0    0.0    0.0   45.9   20.7   62.9   54.4    3.6
    # OPC              4.8    0.0    0.0    5.9    7.9    5.1    5.2    2.1
    # OPC_COP          0.0    0.0    0.0    0.0    0.0    0.4    0.3    0.0

# Of neurons, how many are MSNs?  Let's just look at the NeuN-sorted samples

table(ss(as.character(sce.nac$cellType[grep("neun", sce.nac$sampleID)]),
         "_",1),
      sce.nac$sampleID[grep("neun", sce.nac$sampleID)])[c("Inhib", "MSN.D1", "MSN.D2"), ]
    #        br5182.nac.neun br5207.nac.neun br5701.nac.neun
    # Inhib              229             237              33
    # MSN.D1            2242            2186             111
    # MSN.D2            1738            1979              69
apply(table(ss(as.character(sce.nac$cellType[grep("neun", sce.nac$sampleID)]),
               "_",1),
            sce.nac$sampleID[grep("neun", sce.nac$sampleID)])[c("Inhib", "MSN.D1", "MSN.D2"), ],2,prop.table)
    #        br5182.nac.neun br5207.nac.neun br5701.nac.neun
    # Inhib       0.05440722      0.05383916       0.1549296
    # MSN.D1      0.53266809      0.49659246       0.5211268
    # MSN.D2      0.41292469      0.44956838       0.3239437
    
    # (MSN prop:) 0.9455928       0.9461608        0.8450705

# Proportions of MSN subpopulations across donors
round(apply(prop.table(table(sce.nac$cellType[grep("MSN", sce.nac$cellType)],
                             sce.nac$donor[grep("MSN", sce.nac$cellType)])),2,prop.table),3)
    #                br5161 br5182 br5207 br5212 br5276 br5287 br5400 br5701
    # [(0's)]
    # MSN.D1_A       0.358  0.319  0.403  0.436  0.371  0.422  0.364  0.228
    # MSN.D1_B       0.015  0.040  0.015  0.007  0.001  0.000  0.000  0.039
    # MSN.D1_C       0.030  0.056  0.001  0.074  0.019  0.000  0.003  0.011
    # MSN.D1_D       0.101  0.067  0.042  0.088  0.096  0.060  0.157  0.006
    # MSN.D1_E       0.056  0.065  0.061  0.020  0.013  0.072  0.042  0.283
    # MSN.D1_F       0.007  0.015  0.002  0.000  0.002  0.000  0.003  0.050
    # MSN.D2_A       0.354  0.357  0.437  0.341  0.432  0.349  0.355  0.317
    # MSN.D2_B       0.034  0.034  0.007  0.024  0.051  0.060  0.051  0.039
    # MSN.D2_C       0.034  0.033  0.031  0.010  0.013  0.036  0.025  0.028
    # MSN.D2_D       0.011  0.013  0.000  0.000  0.002  0.000  0.000  0.000
    # [(0's)]



## Numbers for markers ===
load("rdas/markers-stats_NAc-n5_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
# markers.nac.t.design, markers.nac.t.1vAll

sapply(markers.nac.t.design, function(x){sum(x$FDR < 0.05)})
#    Astro  Inhib.1  Inhib.2  Inhib.3  Inhib.4    Micro MSN.D1.1 MSN.D1.2
#      831      340      179      256      283     1417       82      295
# MSN.D1.3 MSN.D1.4 MSN.D2.1 MSN.D2.2    Oligo      OPC
#       87       61      230       66      493      372


## Are DRD1/DRD2 ever expressed together? ===
plotExpressionCustom(sce.nac, anno_name="cellType", features_name="DRDs",
                     features=c("DRD1", "DRD2", "DRD3")) +
  scale_color_manual(values = cell_colors.nac)

geneExprs <- assay(sce.nac, "logcounts")

sce.nac$MSN.D1 <- NA
sce.nac$MSN.D1[grep("MSN.D1", sce.nac$cellType)] <- TRUE

sce.nac$MSN.D1.0 <- NA
sce.nac$MSN.D1.0 <- sce.nac$MSN.D1 & geneExprs["DRD1", ]==0

# Now plot & color by this - is DRD2 expressed?
plotExpression(sce.nac, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType", colour_by="MSN.D1.0", point_alpha=0.5, point_size=1.2, ncol=5,
               add_legend=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25))


sce.nac$MSN.D2 <- NA
sce.nac$MSN.D2[grep("MSN.D2", sce.nac$cellType)] <- TRUE

sce.nac$MSN.D2.0 <- NA
sce.nac$MSN.D2.0 <- sce.nac$MSN.D2 & geneExprs["DRD2", ]==0

# Now plot & color by this - is DRD1 expressed?
plotExpression(sce.nac, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType", colour_by="MSN.D2.0", point_alpha=0.5, point_size=1.2, ncol=5,
               add_legend=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25))




# non-0 expression of BOTH? ===
table(geneExprs["DRD1", union(which(sce.nac$MSN.D1), which(sce.nac$MSN.D2))] > 0 &
        geneExprs["DRD2", union(which(sce.nac$MSN.D1), which(sce.nac$MSN.D2))] > 0)
    #FALSE  TRUE 
    # 9596  1214
    1214 / (1214+9596)  # 11.2%

sce.nac$MSN.broad <- NA
sce.nac$MSN.broad[grep("MSN", sce.nac$cellType)] <- TRUE
# the above is the same as
table(geneExprs["DRD1", which(sce.nac$MSN.broad)] > 0 &
        geneExprs["DRD2", which(sce.nac$MSN.broad)] > 0)


geneExprs.mat <- as.matrix(geneExprs)
sce.nac$bothD1D2 <- ifelse(geneExprs.mat["DRD1", sce.nac$MSN.broad] > 0 &
                                   geneExprs.mat["DRD2", sce.nac$MSN.broad] > 0,
                                 TRUE, FALSE)


# Now plot & color by this - are these high-expressing in some??
plotExpression(sce.nac, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType", colour_by="bothD1D2", point_alpha=0.5, point_size=1.2, ncol=5,
               add_legend=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25))


## Plot these three subsets: ===
pdf("pdfs/revision/regionSpecific_NAc-n8_DRD1-vs-DRD2-exclusivity_MNT2021.pdf", height=4)
# D1-MSN but 0 DRD1
plotExpression(sce.nac, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType", colour_by="MSN.D1.0", point_alpha=0.3, point_size=1.3, ncol=5,
               add_legend=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 12)) +
  ggtitle(label="Classified with D1 subclusters but 0 DRD1 expression")

# D2-MSN but 0 DRD2
plotExpression(sce.nac, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType", colour_by="MSN.D2.0", point_alpha=0.3, point_size=1.3, ncol=5,
               add_legend=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 12)) +
  ggtitle(label="Classified with D2 subclusters but 0 DRD2 expression")

# MSN but non-0 expression for BOTH
plotExpression(sce.nac, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType", colour_by="bothD1D2", point_alpha=0.3, point_size=1.3, ncol=5,
               add_legend=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 12)) +
  ggtitle(label="Classified as MSN (D1 or D2) but (+) expression of BOTH DRD1/DRD2")
dev.off()


# More exploration of these dual-expressing..??
table(sce.nac$cellType, sce.nac$bothD1D2)
    #               FALSE TRUE
    # MSN.D1_A       3360  567
    # MSN.D1_B        210   29
    # MSN.D1_C        271   12
    # MSN.D1_D        581  137
    # MSN.D1_E        361  277
    # MSN.D1_F         84    2
    # MSN.D2_A       4085  177
    # MSN.D2_B        278    7
    # MSN.D2_C        308    6
    # MSN.D2_D         58    0
t(round(apply(table(sce.nac$cellType, sce.nac$bothD1D2)[11:20, ],1,prop.table),3))
    #          FALSE  TRUE
    # MSN.D1_A 0.856 0.144
    # MSN.D1_B 0.879 0.121
    # MSN.D1_C 0.958 0.042
    # MSN.D1_D 0.809 0.191
    # MSN.D1_E 0.566 0.434
    # MSN.D1_F 0.977 0.023
    # MSN.D2_A 0.958 0.042
    # MSN.D2_B 0.975 0.025
    # MSN.D2_C 0.981 0.019
    # MSN.D2_D 1.000 0.000



## MNT 27Jul2021: Tabulating all regionally-defined clusters per donor ===
 #                * keep the 'drop.' clusters for transparency - kept in supp fig.
load("rdas/revision/all-n24-samples_across-regions-analyses_MNT2021.rda", verbose=T)
    # sce.allRegions, chosen.hvgs.union, ref.sampleInfo, cell_colors.amy, cell_colors.dlpfc, cell_colors.hpc, cell_colors.sacc, cell_colors.nac

# Add those new 'de-identified' donor IDs
donorRef <- data.frame(donor=unique(sce.allRegions$donor),
                       pub.donorID=paste0("donor", 1:8))
sce.allRegions$pub.donorID <- donorRef$pub.donorID[match(sce.allRegions$donor, donorRef$donor)]

RS.donorTab <- as.data.frame.matrix(table(sce.allRegions$cellType, sce.allRegions$pub.donorID))
# Trick some reordering by region
RS.donorTab$Region <- ss(rownames(RS.donorTab),"_",1)
RS.donorTab$CellType <- paste0(ss(rownames(RS.donorTab),"_",2), "_", ss(rownames(RS.donorTab),"_",3))
RS.donorTab$CellType <- gsub("_NA","",RS.donorTab$CellType)


# Reorder from most sampled donor to least
RS.donorTab <- RS.donorTab[ ,c("donor1","donor2","donor4","donor8",
                               "donor5","donor6","donor3","donor7",
                               "CellType","Region")]

# Write into table
#write.table(RS.donorTab, file="tables/revision/suppTable_107-RS-classes_donorStratified_MNT2021.tsv",
write.table(RS.donorTab, file="tables/revision/suppTable_107-RS-classes_12drops_donorStratified_MNT2021.tsv",
            sep="\t", quote=F, row.names=F,col.names=T)


### Prelim cluster-level tabulation ===
clusterRefTab <- list()
region.idx <- splitit(sce.allRegions$region)

for(i in names(region.idx)){
  clusterRefTab[[i]] <- as.matrix(table(droplevels(sce.allRegions$cellType[region.idx[[i]]]),
                                            droplevels(sce.allRegions$prelimCluster[region.idx[[i]]]))
                                  )
}
sapply(clusterRefTab, dim)
    #      amy dlpfc hpc nac sacc
    # [1,]  19    19  20  24   25
    # [2,]  57   107  45  24   65
    57 + 107 + 45 + 24 + 65 # 298 total [sub]clusters

for(i in names(clusterRefTab)){
  write.csv(clusterRefTab[[i]],
            file=paste0("tables/revision/suppTable_",i,"_cellClasses-by-prelimCluster_MNT2021.tsv"))
}





### RNAscope misc checks for NAc ===============================

load("rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda", verbose=T)
    # sce.nac, chosen.hvgs.nac, pc.choice.nac, ref.sampleInfo, annotationTab.nac, cell_colors.nac


# Clean up '.drop' clusters (669 nuclei, total)
sce.nac <- sce.nac[ ,-grep("drop.", sce.nac$cellType)]
sce.nac$cellType <- droplevels(sce.nac$cellType)

## Save logcounts of experimental genes for future use (with FULL SCE)
RNAscope.probes <- c("DRD1","TAC1","RXFP2","GABRQ","CRHR2","RXFP1",
                     "DRD2","PENK","HTR7", "PVALB", "GAD1", "PTHLH","KIT")

# Print these for revision
pdf("pdfs/revision/NAc/forRef_RNAscope-probes_NAc-n8_revision-cellClasses_MNT2021.pdf", width=10)
plotExpressionCustom(sce.nac, anno_name="cellType", features_name="select RNAscope",
                     features=RNAscope.probes, ncol=3) +
  scale_color_manual(values = cell_colors.nac)
dev.off()


## For RNAscope images ===
 #  -> For each experiment, print four probes (in the relevant subset of cell classes) in a single column
 #     ( will show side-by-side with the quantified RNAscope images )

## 'Circle' (main fig one) - D1 classes
sce.d1 <- sce.nac[ ,grep("MSN.D1", sce.nac$cellType)]
sce.d1$cellType <- droplevels(sce.d1$cellType)

pdf("./rnascope/circle_pdfs/MNT-revision_snRNA-seq-D1-classes_circleProbes.pdf", width=2)
plotExpressionCustom(sce.d1, anno_name="cellType", features=c("DRD1","TAC1","RXFP1","CRHR2"),
                     features_name="D1-classes differentiation", ncol=1, scales="free_y", xlab="") +
  scale_color_manual(values = cell_colors.nac) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(angle = 90, size = 13),
        plot.title = element_text(size = 10),
        panel.grid.major=element_line(colour="grey95", size=0.8),
        panel.grid.minor=element_line(colour="grey95", size=0.4))
dev.off()


## 'Swirl' - Inhib. interneuron classes
sce.inter <- sce.nac[ ,grep("Inhib", sce.nac$cellType)]
sce.inter$cellType <- droplevels(sce.inter$cellType)

pdf("./rnascope/swirl_pdfs/MNT-revision_snRNA-seq-interneuron-classes_swirlProbes.pdf", width=2)
plotExpressionCustom(sce.inter, anno_name="cellType", features=c("GAD1","KIT","PTHLH","PVALB"),
                     features_name="Interneuron-classes differentiation", ncol=1, scales="free_y",
                     xlab="", point_size=1.1, point_alpha=0.4) +
  scale_color_manual(values = cell_colors.nac) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(angle = 90, size = 13),
        plot.title = element_text(size = 10),
        panel.grid.major=element_line(colour="grey95", size=0.8),
        panel.grid.minor=element_line(colour="grey95", size=0.4))
dev.off()


## 'Triangle' - MSN classes - neuropeptides
sce.msn <- sce.nac[ ,grep("MSN", sce.nac$cellType)]
sce.msn$cellType <- droplevels(sce.msn$cellType)

pdf("./rnascope/triangle_pdfs/MNT-revision_snRNA-seq-MSN-classes-neuropeptides_triangleProbes.pdf", width=2.7)
plotExpressionCustom(sce.msn, anno_name="cellType", features=c("DRD1","DRD2","TAC1","PENK"),
                     features_name="MSN classes - neuropeptide", ncol=1, scales="free_y", xlab="") +
  scale_color_manual(values = cell_colors.nac) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(angle = 90, size = 13),
        plot.title = element_text(size = 10),
        panel.grid.major=element_line(colour="grey95", size=0.8),
        panel.grid.minor=element_line(colour="grey95", size=0.4))
dev.off()


## 'Star' - add'l D1 or D2 class misc markers 
pdf("./rnascope/star_pdfs/MNT-revision_snRNA-seq-MSN-classes-CRHR2-HTR7_starProbes.pdf", width=2.7)
plotExpressionCustom(sce.msn, anno_name="cellType", features=c("DRD1","DRD2","CRHR2","HTR7"),
                     features_name="MSN classes - misc", ncol=1, scales="free_y", xlab="") +
  scale_color_manual(values = cell_colors.nac) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(angle = 90, size = 13),
        plot.title = element_text(size = 10),
        panel.grid.major=element_line(colour="grey95", size=0.8),
        panel.grid.minor=element_line(colour="grey95", size=0.4))
dev.off()


## 'Square' - D1 classes add'l markers
sce.d1 <- sce.nac[ ,grep("MSN.D1", sce.nac$cellType)]
sce.d1$cellType <- droplevels(sce.d1$cellType)

pdf("./rnascope/square_pdfs/MNT-revision_snRNA-seq-D1-classes-addl_squareProbes.pdf", width=2)
plotExpressionCustom(sce.d1, anno_name="cellType", features=c("DRD1","TAC1","RXFP2","GABRQ"),
                     features_name="D1-classes add'l differentiation", ncol=1, scales="free_y", xlab="") +
  scale_color_manual(values = cell_colors.nac) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(angle = 90, size = 13),
        plot.title = element_text(size = 10),
        panel.grid.major=element_line(colour="grey95", size=0.8),
        panel.grid.minor=element_line(colour="grey95", size=0.4))
dev.off()




### Regionally-defined subcluster broad marker heatmaps =============
  # (for supplement: AMY, HIPPO, DLPFC, sACC)
library(pheatmap)

load("rdas/revision/all-n24-samples_across-regions-analyses_forFigOnly_MNT2021.rda", verbose=T)
    # sce.allRegions, chosen.hvgs.union, ref.sampleInfo, Readme
    #cell_colors.amy, cell_colors.dlpfc, cell_colors.hpc, cell_colors.sacc, cell_colors.nac


# First drop all "Ambig.lowNtrxts" & '_nac' because already have
sce.allRegions <- sce.allRegions[ ,-grep("nac_", sce.allRegions$cellType)]
sce.allRegions$cellType <- droplevels(sce.allRegions$cellType)

table(sce.allRegions$cellType)

        RSsubs.donorTab <- as.data.frame.matrix(table(sce.n14$cellType, sce.n14$donor))
        # Trick some reordering by region
        RSsubs.donorTab$Region <- ss(rownames(RSsubs.donorTab),"_",2)
        RSsubs.donorTab$CellType <- paste0(RSsubs.donorTab$Region, "_", ss(rownames(RSsubs.donorTab),"_",1))
        
        # Now basically re-order alphabetically
        RSsubs.donorTab <- RSsubs.donorTab[order(RSsubs.donorTab$CellType), ]


cell.idx <- splitit(sce.allRegions$cellType)
regionCell.idx <- list()

for(i in c("amy_","hpc_","dlpfc_","sacc_")){
  regionCell.idx[[i]] <- lapply(grep(i, names(cell.idx)), function(x){
      y <- cell.idx[[x]]
      y
  })
  
  names(regionCell.idx[[i]]) <- names(cell.idx)[grep(i, names(cell.idx))]
}


dat <- as.matrix(assay(sce.allRegions, "logcounts"))


# Print AMY first
pdf('pdfs/revision/pubFigures/heatmap_AMY-cellTypes_mean-broadCellMarkers_MNT.pdf', useDingbats=TRUE, height=6, width=6)
genes <- c('SNAP25','SLC17A6','SLC17A7','SLC17A8','GAD1','GAD2','AQP4','GFAP','CLDN5','FLT1',
           'CD163','SIGLEC1','C3','CD74','COL1A2','PDGFRB','MBP','PDGFRA','VCAN','SKAP1','CD247')
current_dat <- do.call(cbind, lapply(regionCell.idx[["amy_"]], function(ii) rowMeans(dat[genes, ii])))

neuron.pos <- c(grep("Excit", names(regionCell.idx[["amy_"]])),
                grep("Inhib", names(regionCell.idx[["amy_"]])),
                grep("Neu", names(regionCell.idx[["amy_"]])))
current_dat <- current_dat[ ,c(neuron.pos, setdiff(1:length(names(regionCell.idx[["amy_"]])),
                                                   neuron.pos))
                            ]
pheatmap(current_dat, cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(0.02, 4, length.out = 101),
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "OrRd"))(100),
         fontsize_row = 15.5, fontsize_col = 15.5)
grid::grid.text(label="log2-\nExprs", x=0.96, y=0.6, gp=grid::gpar(fontsize=10))
dev.off()


## Print rest
# First call sACC's 'Neu_FAT2.CDH15' 'Neu_ambig'
names(regionCell.idx[["sacc_"]])[names(regionCell.idx[["sacc_"]])=="sacc_Neu_FAT2.CDH15"] <- "sacc_Neu_ambig"

pdf('pdfs/revision/pubFigures/heatmap_HPC-DLPFC-sACC-cellTypes_mean-broadCellMarkers_MNT.pdf', useDingbats=TRUE, height=6, width=6)
genes <- c('SNAP25','SLC17A6','SLC17A7','SLC17A8','GAD1','GAD2','AQP4','GFAP','CLDN5','FLT1',
           'CD163','SIGLEC1','C3','CD74','COL1A2','PDGFRB','MBP','PDGFRA','VCAN','SKAP1','CD247')

for(i in c("hpc_","dlpfc_","sacc_")){
current_dat <- do.call(cbind, lapply(regionCell.idx[[i]], function(ii) rowMeans(dat[genes, ii])))
neuron.pos <- c(grep("Excit", names(regionCell.idx[[i]])),
                grep("Inhib", names(regionCell.idx[[i]])),
                grep("Neu", names(regionCell.idx[[i]])))
current_dat <- current_dat[ ,c(neuron.pos, setdiff(1:length(names(regionCell.idx[[i]])),
                                                   neuron.pos))
                            ]
pheatmap(current_dat, cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(0.02, 4, length.out = 101),
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "OrRd"))(100),
         fontsize_row = 14.5, fontsize_col = 14.5)
grid::grid.text(label="log2-\nExprs", x=0.96, y=0.6, gp=grid::gpar(fontsize=10))
}
dev.off()




## Supplementary Table 1 === ===
pheno <- read.csv("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/10x-snRNA-seq_Pilot_Subject_Info_rev2021.csv",
                   header=T, sep=",")
# Add more 'de-identifying' donor ID used in figs
donorRef <- data.frame(donor=unique(sce.allRegions$donor),
                       pub.donorID=paste0("donor", 1:8))
donorRef$donor <- gsub("br","Br",donorRef$donor)

pheno$pub.donorID <- donorRef$pub.donorID[match(pheno$brnum, donorRef$donor)]

# Re-order
pheno <- pheno[order(pheno$pub.donorID), c("pub.donorID", "brnum", "sex", "agedeath", "race", "primarydx",
                                           "pmi","pmi_confidence_level","smoking","codeine","morphine","bmi_calculated")]

# Cleaner names
colnames(pheno) <- c("DonorID.pub","BrNum", "Sex", "AgeDeath_yrs", "Race", "PrimaryDx",
                     "PMI_hrs", "PMI_confidence_level",
                     "Smoking","Codeine","Morphine","BMI")

write.table(pheno, row.names=F, col.names=T, sep="\t", quote=F,
            file="tables/revision/suppTab1_SampleInformation_MNT2021.tsv")




### Donor proportion across clusters ===========
library(reshape2)
load("rdas/revision/all-n24-samples_across-regions-analyses_MNT2021.rda", verbose=T)
    # sce.allRegions, chosen.hvgs.union, ref.sampleInfo, cell_colors.amy, cell_colors.dlpfc, cell_colors.hpc, cell_colors.sacc, cell_colors.nac

donor.cols <- tableau10medium[1:8]
names(donor.cols) <- paste0("donor", c("A","B","C","D","E","F","G","H"))
# or
names(donor.cols) <- paste0("donor", 1:8)

# Add those new 'de-identified' donor IDs
donorRef <- data.frame(donor=unique(sce.allRegions$donor),
                       pub.donorID=paste0("donor", 1:8))
sce.allRegions$pub.donorID <- donorRef$pub.donorID[match(sce.allRegions$donor, donorRef$donor)]


## Do by region
region.idx <- splitit(sce.allRegions$region)
region.names <- c("AMY", "DLPFC", "HPC", "NAc", "sACC")
names(region.names) <- names(region.idx)

pdf("pdfs/revision/pubFigures/suppFig_cellClass-donorDistrib_MNT2021.pdf", height=3)
for(i in names(region.idx)){
  sce.i <- sce.allRegions[ ,region.idx[[i]]]
  sce.i$cellType <- droplevels(sce.i$cellType)
  coldat <- as.data.frame(colData(sce.i))
  coldat$Donor <- coldat$pub.donorID
  coldat$Donor[grep("NeuN", coldat$protocol)] <- paste0(coldat$Donor[grep("NeuN", coldat$protocol)],
                                                        "_neun")
  coldat$Donor <- as.factor(coldat$Donor)
  
  # Temp colors bc the same donors aren't always NeuN-enriched
  cols.i <- donor.cols
  names(cols.i) <- ifelse(paste0(names(cols.i),"_neun") %in% levels(coldat$Donor),
                          paste0(names(cols.i),"_neun"),
                          names(cols.i))
  
  customLabs <- paste0(levels(coldat$Donor)," (",table(coldat$Donor),")")
  
  print(
    ggplot(coldat,aes(x=cellType, color=Donor, fill=Donor)) + geom_bar(position="fill") +
      scale_color_manual(values = cols.i,
                         labels=customLabs) +
      scale_fill_manual(values = cols.i,
                        labels=customLabs) +
      theme(axis.text.x = element_text(angle=90, hjust = 1, size=8),
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.4, "cm")) +
      ylab("Proportion") +
      xlab("") +
      scale_x_discrete(labels=paste0(levels(sce.i$cellType)," (",table(sce.i$cellType),")")) +
      ggtitle(paste0("Distribution of ",region.names[i]," cell classes by donor"))
  )
}
dev.off()


# How many total neurons?
length(c(grep("Inhib", sce.allRegions$cellType),
         grep("Excit", sce.allRegions$cellType),
         grep("Neu", sce.allRegions$cellType),
         grep("MSN",sce.allRegions$cellType)))
  # [1] 28150








### For LeCo/iSEE ==========================
load("rdas/ztemp_Amyg-n2_SCE-with-tSNEonOptPCs-minus-PC5_MNT.rda", verbose=T)
    #sce.amy.tsne.optb, Readme

## Percentage of ALL nuclei expressing non-0 amount of each gene
rowData(sce.amy.tsne.optb)$propNucleiExprs <- apply(assay(sce.amy.tsne.optb, "counts"), 1,
                                                    function(x){round(mean(x != 0), 4)})

# The above, by cell type ===
cellType.idx <- splitit(sce.amy.tsne.optb$cellType.split)

rowdat.sce <- rowData(sce.amy.tsne.optb)
for(i in names(cellType.idx)){
  rowdat.sce[ ,paste0("propExprsIn.",i)] <- apply(assay(sce.amy.tsne.optb, "counts")[ ,cellType.idx[[i]] ], 1,
                        function(x){round(mean(x != 0), 4)})
}
rowData(sce.amy.tsne.optb) <- rowdat.sce


### Clean region-specific SCEs Amazon hosting ====================
  # Updated MNT 2021 to clean up SCEs
  # Notes: Keeping 'drop.' clusters, since other investigators might want to explore
  #        such technical artifact-driven clusters / QC differently

# Assign 'publication' donor IDs
donor_map <- paste0("donor", seq_len(8))
names(donor_map) <- c("br5161", "br5212", "br5287", "br5400", "br5276", "br5207", "br5182", "br5701")

keepCols <-  c("Barcode","sum","detected","doubletScore","region","donor","sex",
               "processBatch","protocol","sequencer","sizeFactor",
               "prelimCluster","collapsedCluster","cellType")

## Adapted from Leo's `create_small_sce_2021`:
  create_unif_sce <- function(sce_original) {

      # Assign 'publication' donor IDs
      stopifnot(all(unique(sce_original$donor) %in% names(donor_map)))
      sce_original$donor <- unname(donor_map[sce_original$donor])
      
      message(Sys.time(), " Uniformizing the sce object")
      sce_small <- sce_original
      #assays(sce_small) <- assays(sce_small)["logcounts"] - keep this
      sce_small$region <- factor(sce_small$region)
      sce_small$donor <- factor(sce_small$donor)
      sce_small$sex <- factor(sce_small$sex)
      sce_small$processBatch <- factor(sce_small$processBatch)
      sce_small$protocol <- factor(sce_small$protocol)
      sce_small$sequencer <- factor(sce_small$sequencer)
      sce_small$cellType <- factor(sce_small$cellType)
      
      colData(sce_small) <- colData(sce_small)[ ,keepCols]
      
      
      ## Make the rows more browsable
      sce_small$Barcode <- make.names(sce_small$Barcode, unique = TRUE)
      colnames(sce_small) <- paste0(sce_small$donor, '_', sce_small$Barcode)
      #metadata(sce_small) <- list()  - keep this

          ## What proportion of [all] nuclei is gene X captured in? ===
          message(Sys.time(), " computing propNucleiExprs")
          
          rowData(sce_small)$propNucleiExprs <- apply(
            assay(sce_original, "counts"),
            1,
            function(x) {
              mean(x != 0)
            }
          )
          # The above, by cell type:
          cellType.idx <- splitit(sce_small$cellType)
          rowdat.sce <- rowData(sce_small)
          for(i in names(cellType.idx)){
            message(Sys.time(), " computing propNucleiExprs for ", i)
            rowdat.sce[, paste0("propExprsIn.", i)] <- apply(
              assay(sce_original, "counts")[, cellType.idx[[i]]],
              1,
              function(x){
                mean(x != 0)
              }
            )
          }
          rowData(sce_small) <- rowdat.sce
          
          print(pryr::object_size(sce_original))
          print(pryr::object_size(sce_small))
          
          return(sce_small)
        }


## NAc === ===
load("rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda", verbose=T)

# This one (only) didn't have a 'collapsedCluster' bc clustering step 2 wasn't needed
sce.nac$collapsedCluster <- sce.nac$prelimCluster
    
# Reduce
sce.nac.red <- create_unif_sce(sce.nac)
# Check
head(colData(sce.nac.red))
head(rowData(sce.nac.red))
plotTSNE(sce.nac.red, colour_by="cellType")

sce.nac.tran <- sce.nac.red
#dir.create("rdas/revision/forAmazonS3/")
save(sce.nac.tran, file="rdas/revision/forAmazonS3/SCE_NAc-n8_tran-etal.rda")

rm(list=ls(pattern="nac"))


## AMY ===
load("rdas/revision/regionSpecific_Amyg-n5_cleaned-combined_SCE_MNT2021.rda", verbose=T)

# Reduce
sce.amy.red <- create_unif_sce(sce.amy)
# Check
head(colData(sce.amy.red))
head(rowData(sce.amy.red))
plotTSNE(sce.amy.red, colour_by="cellType")

sce.amy.tran <- sce.amy.red
save(sce.amy.tran, file="rdas/revision/forAmazonS3/SCE_AMY-n5_tran-etal.rda")

rm(list=ls(pattern="amy"))


## sACC ===
load("rdas/revision/regionSpecific_sACC-n5_cleaned-combined_SCE_MNT2021.rda", verbose=T)

# Reduce
sce.sacc.red <- create_unif_sce(sce.sacc)
# Check
head(colData(sce.sacc.red))
head(rowData(sce.sacc.red))
plotTSNE(sce.sacc.red, colour_by="cellType")

sce.sacc.tran <- sce.sacc.red
save(sce.sacc.tran, file="rdas/revision/forAmazonS3/SCE_sACC-n5_tran-etal.rda")

rm(list=ls(pattern="sacc"))


## DLPFC ===
load("rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_MNT2021.rda", verbose=T)

# Reduce
sce.dlpfc.red <- create_unif_sce(sce.dlpfc)
# Check
head(colData(sce.dlpfc.red))
head(rowData(sce.dlpfc.red))
plotTSNE(sce.dlpfc.red, colour_by="cellType")

sce.dlpfc.tran <- sce.dlpfc.red
save(sce.dlpfc.tran, file="rdas/revision/forAmazonS3/SCE_DLPFC-n3_tran-etal.rda")

rm(list=ls(pattern="dlpfc"))


## HPC ===
load("rdas/revision/regionSpecific_HPC-n3_cleaned-combined_SCE_MNT2021.rda", verbose=T)

# Reduce
sce.hpc.red <- create_unif_sce(sce.hpc)
# Check
head(colData(sce.hpc.red))
head(rowData(sce.hpc.red))
plotTSNE(sce.hpc.red, colour_by="cellType")

sce.hpc.tran <- sce.hpc.red

save(sce.hpc.tran, file="rdas/revision/forAmazonS3/SCE_HPC-n3_tran-etal.rda")



### Session info for 27Jul2021 ====================================================
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
#   [1] reshape2_1.4.4              dplyr_1.0.5                 plotly_4.9.3               
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
# [16] sparseMatrixStats_1.2.1   cachem_1.0.4              jsonlite_1.7.2           
# [19] Rsamtools_2.6.0           ResidualMatrix_1.0.0      dbplyr_2.1.1             
# [22] R.oo_1.24.0               HDF5Array_1.18.1          compiler_4.0.4           
# [25] httr_1.4.2                dqrng_0.3.0               assertthat_0.2.1         
# [28] Matrix_1.3-4              fastmap_1.1.0             lazyeval_0.2.2           
# [31] limma_3.46.0              BiocSingular_1.6.0        htmltools_0.5.1.1        
# [34] prettyunits_1.1.1         tools_4.0.4               rsvd_1.0.5               
# [37] igraph_1.2.6              gtable_0.3.0              glue_1.4.2               
# [40] GenomeInfoDbData_1.2.4    rappdirs_0.3.3            Rcpp_1.0.6               
# [43] vctrs_0.3.8               Biostrings_2.58.0         rhdf5filters_1.2.0       
# [46] rtracklayer_1.50.0        DelayedMatrixStats_1.12.3 stringr_1.4.0            
# [49] beachmat_2.6.4            lifecycle_1.0.0           irlba_2.3.3              
# [52] statmod_1.4.35            XML_3.99-0.6              edgeR_3.32.1             
# [55] zlibbioc_1.36.0           scales_1.1.1              hms_1.0.0                
# [58] ProtGenerics_1.22.0       rhdf5_2.34.0              RColorBrewer_1.1-2       
# [61] curl_4.3                  memoise_2.0.0             gridExtra_2.3            
# [64] segmented_1.3-4           biomaRt_2.46.3            stringi_1.5.3            
# [67] RSQLite_2.2.7             BiocParallel_1.24.1       rlang_0.4.11             
# [70] pkgconfig_2.0.3           bitops_1.0-7              lattice_0.20-41          
# [73] purrr_0.3.4               Rhdf5lib_1.12.1           htmlwidgets_1.5.3        
# [76] labeling_0.4.2            GenomicAlignments_1.26.0  bit_4.0.4                
# [79] tidyselect_1.1.1          plyr_1.8.6                magrittr_2.0.1           
# [82] R6_2.5.0                  generics_0.1.0            DelayedArray_0.16.3      
# [85] DBI_1.1.1                 pillar_1.6.0              withr_2.4.2              
# [88] RCurl_1.98-1.3            tibble_3.1.1              crayon_1.4.1             
# [91] utf8_1.2.1                BiocFileCache_1.14.0      viridis_0.6.0            
# [94] progress_1.2.2            locfit_1.5-9.4            grid_4.0.4               
# [97] data.table_1.14.0         blob_1.2.1                digest_0.6.27            
# [100] tidyr_1.1.3               R.utils_2.10.1            openssl_1.4.3            
# [103] munsell_0.5.0             beeswarm_0.4.0            viridisLite_0.4.0        
# [106] vipor_0.4.5               askpass_1.1  


