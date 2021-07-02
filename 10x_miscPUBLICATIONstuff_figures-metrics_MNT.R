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

### METRICS ====================================================

## Metrics and stuff for paper, all n=14 samples that went into analyses ===
load("rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda", verbose=T)
# sce.all.n12, ...
load("rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda", verbose=T)
# sce.nac.all, ...

# Get only NeuNs to add to 'sce.all.n12'
sce.neuns <- sce.nac.all[ ,sce.nac.all$protocol == "Frank.NeuN"]
# Trick 'sce.neuns' to keep the column which has regionally-defined annotations
sce.neuns$cellType.RS.sub <- sce.neuns$cellType.final
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



## Create compiled file for Cell Ranger metrics - MNT 19Jun2020 ====
n14.samples <- c("Br5161_Amy", "Br5182_NAc_NeuN", "Br5212_DLPFC", "Br5212_sACC", "Br5287_HPC",
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

# Write this out for future reference
write.table(n14.metrics.collapsed, row.names=T, col.names=T, sep="\t", quote=F,
            file="tables/METRICS-n14-analyzed_CellRanger-premRNA-output_MNT.tsv")
# As semi-colon-SV
write.table(n14.metrics.collapsed, row.names=T, col.names=T, sep=";", quote=F, 
            file="tables/METRICS-n14-analyzed_CellRanger-premRNA-output_MNT.csv")


# Sequencing depth
quantile(as.numeric(gsub(",","",n14.metrics.collapsed$Number.of.Reads)))
#       0%       25%       50%       75%      100%
#118752669 148651409 253012622 274925035 296198922


## NAc, all n=5 ===
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



# Make diff object for 'ambig.lowNtrxts' removed
sce.nac <- sce.nac[ ,!sce.nac$cellType=="ambig.lowNtrxts"]
sce.nac$cellType.final <- droplevels(sce.nac$cellType.final)

## Proportion of cell types
table(sce.nac$cellType.final)
prop.table(table(sce.nac$cellType.final))
#       Astro     Inhib.1     Inhib.2     Inhib.3     Inhib.4       Micro
# 0.041451171 0.001901430 0.004259203 0.020991786 0.013918467 0.013766352
#    MSN.D1.1    MSN.D1.2    MSN.D1.3    MSN.D1.4    MSN.D2.1    MSN.D2.2
# 0.010039550 0.022893216 0.054685123 0.285442653 0.022817159 0.276163675
#       Oligo         OPC
# 0.213492546 0.018177670

# How many MSNs in our dataset?
sum(prop.table(table(sce.nac$cellType.final))[grep("MSN", levels(sce.nac$cellType.final))])
# However this is including two samples where we enriched on neurons


table(sce.nac$cellType.final, sce.nac$sample)
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
apply(table(sce.nac$cellType.final, sce.nac$sample), 2, function(x){round(prop.table(x),3)})
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

table(ss(as.character(sce.nac$cellType.final[grep("neun", sce.nac$sample)]),
         "\\.",1),
      sce.nac$sample[grep("neun", sce.nac$sample)])
#        nac.neun.5182 nac.neun.5207
#  Inhib           248           241
#  MSN            4012          4182
apply(table(ss(as.character(sce.nac$cellType.final[grep("neun", sce.nac$sample)]),
               "\\.",1),
            sce.nac$sample[grep("neun", sce.nac$sample)]),2,prop.table)
#        nac.neun.5182 nac.neun.5207
#  Inhib    0.05821596     0.0544879
#  MSN      0.94178404     0.9455121

# Proportions of MSN subpopulations across donors
round(apply(prop.table(table(sce.nac$cellType.final[grep("MSN", sce.nac$cellType.final)],
                             sce.nac$sample[grep("MSN", sce.nac$cellType.final)])),2,prop.table),3)
#          nac.5161 nac.5212 nac.5287 nac.neun.5182 nac.neun.5207
# Astro       0.000    0.000    0.000         0.000         0.000
# Inhib.1     0.000    0.000    0.000         0.000         0.000
# Inhib.2     0.000    0.000    0.000         0.000         0.000
# Inhib.3     0.000    0.000    0.000         0.000         0.000
# Inhib.4     0.000    0.000    0.000         0.000         0.000
# Micro       0.000    0.000    0.000         0.000         0.000
# MSN.D1.1    0.008    0.000    0.000         0.029         0.003
# MSN.D1.2    0.039    0.010    0.000         0.071         0.001
# MSN.D1.3    0.066    0.027    0.070         0.092         0.076
# MSN.D1.4    0.693    0.565    0.837         0.375         0.437
# MSN.D2.1    0.035    0.020    0.035         0.033         0.035
# MSN.D2.2    0.160    0.378    0.058         0.399         0.447
# Oligo       0.000    0.000    0.000         0.000         0.000
# OPC         0.000    0.000    0.000         0.000         0.000



## Numbers for markers ===
load("rdas/markers-stats_NAc-n5_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
# markers.nac.t.design, markers.nac.t.1vAll

sapply(markers.nac.t.design, function(x){sum(x$FDR < 0.05)})
#    Astro  Inhib.1  Inhib.2  Inhib.3  Inhib.4    Micro MSN.D1.1 MSN.D1.2
#      831      340      179      256      283     1417       82      295
# MSN.D1.3 MSN.D1.4 MSN.D2.1 MSN.D2.2    Oligo      OPC
#       87       61      230       66      493      372


## Are DRD1/DRD2 ever expressed together? ===
plotExpression(sce.nac, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType.final", colour_by="cellType.final", point_alpha=0.5, point_size=.7, ncol=5,
               add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                            geom = "crossbar", width = 0.3,
                                            colour=rep(tableau20[1:14], 2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25))

geneExprs <- assay(sce.nac, "logcounts")

sce.nac$MSN.D1 <- NA
sce.nac$MSN.D1[grep("MSN.D1", sce.nac$cellType.final)] <- TRUE

sce.nac$MSN.D1.0 <- NA
sce.nac$MSN.D1.0 <- sce.nac$MSN.D1 & geneExprs["DRD1", ]==0

# Now plot & color by this - is DRD2 expressed?
plotExpression(sce.nac, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType.final", colour_by="MSN.D1.0", point_alpha=0.5, point_size=1.2, ncol=5,
               add_legend=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25))


sce.nac$MSN.D2 <- NA
sce.nac$MSN.D2[grep("MSN.D2", sce.nac$cellType.final)] <- TRUE

sce.nac$MSN.D2.0 <- NA
sce.nac$MSN.D2.0 <- sce.nac$MSN.D2 & geneExprs["DRD2", ]==0

# Now plot & color by this - is DRD1 expressed?
plotExpression(sce.nac, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType.final", colour_by="MSN.D2.0", point_alpha=0.5, point_size=1.2, ncol=5,
               add_legend=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25))




# non-0 expression of BOTH? ===
table(geneExprs["DRD1", union(which(sce.nac$MSN.D1), which(sce.nac$MSN.D2))] > 0 &
        geneExprs["DRD2", union(which(sce.nac$MSN.D1), which(sce.nac$MSN.D2))] > 0)
#FALSE  TRUE
# 8017   819

sce.nac$MSN.broad <- NA
sce.nac$MSN.broad[grep("MSN", sce.nac$cellType.final)] <- TRUE
# the above is the same as
table(geneExprs["DRD1", which(sce.nac$MSN.broad)] > 0 &
        geneExprs["DRD2", which(sce.nac$MSN.broad)] > 0)


geneExprs.mat <- as.matrix(geneExprs)
sce.nac$bothD1D2 <- ifelse(geneExprs.mat["DRD1", sce.nac$MSN.broad] > 0 &
                                   geneExprs.mat["DRD2", sce.nac$MSN.broad] > 0,
                                 TRUE, FALSE)


# Now plot & color by this - are these high-expressing in some??
plotExpression(sce.nac, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType.final", colour_by="bothD1D2", point_alpha=0.5, point_size=1.2, ncol=5,
               add_legend=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25))


## Plot these three subsets: ===
pdf("pdfs/exploration/DRD1-vs-DRD2-exclusivity_all-NAc-n5_May202.pdf", height=6)
# D1-MSN but 0 DRD1
plotExpression(sce.nac, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType.final", colour_by="MSN.D1.0", point_alpha=0.5, point_size=1.3, ncol=5,
               add_legend=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 12)) +
  ggtitle(label="Classified with D1 subclusters but 0 DRD1 expression")

# D2-MSN but 0 DRD2
plotExpression(sce.nac, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType.final", colour_by="MSN.D2.0", point_alpha=0.5, point_size=1.3, ncol=5,
               add_legend=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 12)) +
  ggtitle(label="Classified with D2 subclusters but 0 DRD2 expression")

# MSN but non-0 expression for BOTH
plotExpression(sce.nac, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType.final", colour_by="bothD1D2", point_alpha=0.5, point_size=1.3, ncol=5,
               add_legend=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 12)) +
  ggtitle(label="Classified as MSN (D1 or D2) but (+) expression of BOTH DRD1/DRD2")
dev.off()


# More exploration of these dual-expressing..??
table(sce.nac$prelimCluster, sce.nac$bothD1D2)
# don't seem to make up or be enriched in their own cluster

table(sce.nac$cellType.final, sce.nac$bothD1D2)
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

table(sce.nac$MSN.broad, sce.nac$bothD1D2)
#      FALSE TRUE
# TRUE  8017  819



## MNT 23Sep2020: Tabulating all regionally-defined clusters per donor ===
# -> first make 'sce.n14', above
sce.n14 <- sce.n14[ ,-grep("Ambig.lowNtrxts", sce.n14$cellType.RS.sub)]
sce.n14$cellType.RS.sub <- droplevels(sce.n14$cellType.RS.sub)

RSsubs.donorTab <- as.data.frame.matrix(table(sce.n14$cellType.RS.sub, sce.n14$donor))
# Trick some reordering by region
RSsubs.donorTab$Region <- ss(rownames(RSsubs.donorTab),"_",2)
RSsubs.donorTab$CellType <- paste0(RSsubs.donorTab$Region, "_", ss(rownames(RSsubs.donorTab),"_",1))

# Now basically re-order alphabetically
RSsubs.donorTab <- RSsubs.donorTab[order(RSsubs.donorTab$CellType), ]
# And most sampled donor to least
RSsubs.donorTab <- RSsubs.donorTab[ ,c(1,4,5,2,3, 7,6)]

# Write into table
write.table(RSsubs.donorTab, file="tables/suppTable_68-RSsubclusters_donorStratified_MNT.tsv",
            sep="\t", quote=F, row.names=F,col.names=T)


    ## For reference - how the N's have changed after dropping 'ambig.lowNtrxts'
        table(sce.n14$region)
            #  amy dlpfc   hpc   nac  sacc
            # 6632  5399 10444 13241  7047
        
        table(sce.n14$region)
            #  amy dlpfc   hpc   nac  sacc
            # 6582  5231 10343 13148  7004
        
        # % Kept?
        round(table(sce.n14$region) / table(sce.n14$region),3)
            #   amy dlpfc   hpc   nac  sacc
            # 0.992 0.969 0.990 0.993 0.994

## Similarly for pan-brain clustering/annotation ===

# As for downstream analyses, remove the clusters that don't pursue:
#     'Ambig.hiVCAN' & 'Excit.4' & those .RS-annot'd as 'Ambig.lowNtrxts'
sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType.RS != "Ambig.lowNtrxts"] # 445
sce.all.n12$cellType.RS <- droplevels(sce.all.n12$cellType.RS)

sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType != "Ambig.hiVCAN"] # 32 nuclei
sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType != "Excit.4"]  # 33 nuclei
sce.all.n12$cellType <- droplevels(sce.all.n12$cellType)

cellType.rgnTab <- as.data.frame.matrix(table(sce.all.n12$cellType, sce.all.n12$region))

# Write into table
cellType.rgnTab <- cbind(rownames(cellType.rgnTab), cellType.rgnTab)
colnames(cellType.rgnTab)[1] <- "cellType.panBrain"
write.table(cellType.rgnTab, file="tables/suppTable_17panBrain-subclusters_regionStratified_MNT.tsv",
            sep="\t", quote=F, row.names=F,col.names=T)
  
### RNAscope misc checks for NAc ===============================
# MNT/AnJa 21Sep2020
load("rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda", verbose=T)
#sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all, ref.sampleInfo

# For AnJa: within DRD2's:
sce.nac.all$DRD2 <- FALSE
sce.nac.all$DRD2[grep("MSN.D2.", sce.nac.all$cellType.final)] <- TRUE
table(sce.nac.all$DRD2) # 3931


# Amongst DRD2's, could one define which are of which class based on non-0 expression of these markers?:
table(sce.nac.all$DRD2 & assay(sce.nac.all, "logcounts")["PENK", ] > 0,
      sce.nac.all$DRD2 & assay(sce.nac.all, "logcounts")["TAC1", ] > 0)
#       FALSE TRUE
# FALSE  9467  277
# TRUE   3116  381

(3116+381)/(3116+381+277)
#[1] 0.9266031

(277+381)/(3116+381+277)
#[1] 0.1743508


table(assay(sce.nac.all, "logcounts")["DRD2", ] > 0 & assay(sce.nac.all, "logcounts")["PENK", ] > 0,
      assay(sce.nac.all, "logcounts")["DRD2", ] > 0 & assay(sce.nac.all, "logcounts")["TAC1", ] > 0)
#        FALSE TRUE
# FALSE  8554  601
# TRUE   3385  701

## Without the 'data-driven'defined cluster
# ## cross-tab - first transpose and make data.frame 
# d_match <- t(as.data.frame(as.matrix(assay(sce.nac.all, "logcounts"))))
# 
# with(d_match[d_match[ ,"DRD2"] > 0, ],
#      table(PENK = PENK > 0 ,  TAC1 = TAC1 > 0))
# 
# 
# with(d_match[ ,d_match["DRD2",] > 0],
#      table(PENK = d_match["PENK", ] > 0 ,  TAC1 = d_match["TAC1", ] > 0))
# 
# with(d_match[ ,d_match["DRD2",] > 0 & d_match["DRD1",]==0],
#      table(PENK = "PENK" >0 ,  TAC1 = "TAC1"> 0))
#     ** not happy with trying to re-formate the dgCMatrix into a data.frame


## Try subsetting ===
sce.nac.d2 <- sce.nac.all[ ,assay(sce.nac.d2, "logcounts")["DRD2", ] > 0]

# Now tabulate
logExprs.DRD2 <-  assay(sce.nac.d2, "logcounts")
table(PENK = logExprs.DRD2["PENK", ] > 0,
      TAC1 = logExprs.DRD2["TAC1", ] > 0)
#       TAC1
# PENK    FALSE TRUE
#  FALSE  4987 2307
#  TRUE   3707 1958



d = as.data.frame(t(as.matrix(assay(sce.nac.all, "logcounts")[c("TAC1","PENK","DRD1","DRD2"), ])))
with(d[d$DRD2 > 0,],
     table(PENK = PENK >0 ,  TAC1 = TAC1> 0))
#       TAC1
# PENK    FALSE TRUE
#   FALSE   430  601
#   TRUE   3385  701
with(d[d$DRD2 > 0 & d$DRD1==0,],
     table(PENK = PENK >0 ,  TAC1 = TAC1> 0))
#       TAC1
# PENK    FALSE TRUE
#   FALSE   339  323
#   TRUE   3252  372


## Save logcounts of experimental genes for future use (with FULL SCE)
RNAscope.probes <- c("DRD1","TAC1","RXFP2","GABRQ","CRHR2","RXFP1",
                     "DRD2","PENK","HTR7", "PVALB", "GAD1", "PTHLH","KIT")

logExprs.rnascope <- as.data.frame(t(as.matrix(assay(sce.nac.all, "logcounts")[c(RNAscope.probes), ])))

save(logExprs.rnascope, file="rnascope/snRNA-seq-logExprs_13probes_MNT.rda")



# On-the-fly printing ===
sce.nac.all <- sce.nac.all[ ,!sce.nac.all$cellType=="ambig.lowNtrxts"]
sce.nac.all$cellType.final <- droplevels(sce.nac.all$cellType.final)

plotExpression(sce.nac.all, exprs_values = "logcounts", features=c("RXFP1"),
               x="cellType.final", colour_by="cellType.final", point_alpha=0.6, point_size=1.3, ncol=5,
               add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                            geom = "crossbar", width = 0.3,
                                            colour=rep(tableau20[1:14], 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, size=20), plot.title = element_text(size = 25))




### Regionally-defined subcluster broad marker heatmaps =============
  # (for supplement: AMY, HIPPO, DLPFC, sACC)
library(pheatmap)

load("rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda", verbose=T)
    # sce.all.n12, chosen.hvgs.all.n12, pc.choice.n12, ref.sampleInfo, clusterRefTab.all.n12


# First drop all "Ambig.lowNtrxts" & '_nac' because already have
sce.rest <- sce.all.n12[ ,-grep("Ambig.lowNtrxts", sce.all.n12$cellType.RS.sub)]
sce.rest <- sce.rest[ ,-grep("_nac", sce.rest$cellType.RS.sub)]
sce.rest$cellType.RS.sub <- droplevels(sce.rest$cellType.RS.sub)

table(sce.rest$cellType.RS.sub)

        RSsubs.donorTab <- as.data.frame.matrix(table(sce.n14$cellType.RS.sub, sce.n14$donor))
        # Trick some reordering by region
        RSsubs.donorTab$Region <- ss(rownames(RSsubs.donorTab),"_",2)
        RSsubs.donorTab$CellType <- paste0(RSsubs.donorTab$Region, "_", ss(rownames(RSsubs.donorTab),"_",1))
        
        # Now basically re-order alphabetically
        RSsubs.donorTab <- RSsubs.donorTab[order(RSsubs.donorTab$CellType), ]


cell.idx <- splitit(sce.rest$cellType.RS.sub)
regionCell.idx <- list()

for(i in c("_amy","_hpc","_dlpfc","_sacc")){
  regionCell.idx[[i]] <- lapply(grep(i, names(cell.idx)), function(x){
      y <- cell.idx[[x]]
      y
  })
  
  names(regionCell.idx[[i]]) <- names(cell.idx)[grep(i, names(cell.idx))]
}


dat <- as.matrix(assay(sce.rest, "logcounts"))


# Print AMY first
pdf('pdfs/pubFigures/heatmap-geneExprs_AMY_cellType.split_mean-broadMarkers_MNT.pdf', useDingbats=TRUE, height=6, width=6)
genes <- c('SNAP25','SLC17A6','SLC17A7','GAD1','GAD2','AQP4','GFAP','C3','CD74','MBP','PDGFRA','VCAN','CLDN5','FLT1','SKAP1','TRAC')
current_dat <- do.call(cbind, lapply(regionCell.idx[["_amy"]], function(ii) rowMeans(dat[genes, ii])))

neuron.pos <- grep("\\.", names(regionCell.idx[["_amy"]]))
current_dat <- current_dat[ ,c(neuron.pos, setdiff(1:length(names(regionCell.idx[["_amy"]])),
                                                   neuron.pos))
                            ]
pheatmap(current_dat, cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(0.02, 4, length.out = 101),
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "OrRd"))(100),
         fontsize_row = 17.5, fontsize_col = 17.5)
dev.off()


# Print rest
pdf('pdfs/pubFigures/heatmap-geneExprs_HIPPO-DLPFC-sACC_cellType.split_mean-broadMarkers_MNT.pdf', useDingbats=TRUE, height=6, width=6)
genes <- c('SNAP25','SLC17A6','SLC17A7','GAD1','GAD2','AQP4','GFAP','C3','CD74','MBP','PDGFRA','VCAN','CLDN5','FLT1','SKAP1','TRAC')

for(i in c("_hpc","_dlpfc","_sacc")){
current_dat <- do.call(cbind, lapply(regionCell.idx[[i]], function(ii) rowMeans(dat[genes, ii])))
neuron.pos <- grep("\\.", names(regionCell.idx[[i]]))
current_dat <- current_dat[ ,c(neuron.pos, setdiff(1:length(names(regionCell.idx[[i]])),
                                                   neuron.pos))
                            ]
pheatmap(current_dat, cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(0.02, 4, length.out = 101),
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "OrRd"))(100),
         fontsize_row = 17.5, fontsize_col = 17.5)
}
dev.off()




## Supplementary Table 1 === ===
sheet1 <- as.data.frame(readxl::read_excel("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/10xpilot_Subject_Info.xlsx", sheet=1))
sheet2 <- as.data.frame(readxl::read_excel("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/10xpilot_Subject_Info.xlsx", sheet=2))
sheet2 <- sheet2[ ,c("brnum", "sex", "race", "agedeath", "primarydx", "pmi")]
colnames(sheet2) <- c("BrNum", "Sex", "Race", "AgeDeath", "PrimaryDx", "PMI")

sheet1 <- sheet1[!duplicated(sheet1$BrNum), 1:3]
sheet2$BrNum == sheet1$BrNum # order good

suppTab1 <- cbind(sheet2, sheet1[ ,"BestRIN PFC"])
colnames(suppTab1)[7] <- "PFC_RIN"

write.table(suppTab1, row.names=F, col.names=T, sep="\t", quote=F,
            file="tables/suppTab1_SampleInformation.tsv")





## Brief exploration of AnJa's LDSC (ignore for first submission) ==========
dat.ldsc <- read.table("LDSC/suppTable_ldsc.tsv", sep="\t")
colnames(dat.ldsc) <- dat.ldsc[1, ]
dat.ldsc <- dat.ldsc[-1, ]

for(i in colnames(dat.ldsc)[3:11]){
  dat.ldsc[ ,i] <- signif(as.numeric(dat.ldsc[ ,i],2))
}

unique(dat.ldsc$Trait)
    # [1] "ADHD"                            "Alzheimers_disease"
    # [3] "Anorexia_nervosa"                "Anxiety_disorder"
    # [5] "Autism_spectrum_disorder"        "Bipolar_disorder"
    # [7] "BMI"                             "Childhood_cognitive_performance"
    # [9] "Cigarettes_per_day"              "College_attainment"
    # [11] "Conscientiousness"               "Coronary_artery_disease"
    # [13] "Crohns_disease"                  "Depressive_symptoms"
    # [15] "Epilepsy"                        "Ever_smoked"
    # [17] "Extraversion"                    "Focal_epilepsy"
    # [19] "Generalized_epilepsy"            "Height"
    # [21] "Intracarebral_hemorrhage"        "IQ"
    # [23] "Ischemic_stroke"                 "Major_depressive_disorder"
    # [25] "Neuroticism"                     "Openness"
    # [27] "PTSD"                            "Schizophrenia"
    # [29] "Subjective_well-being"           "Years_of_education"

# Bonferroni-significant
dat.ldsc[dat.ldsc$Trait=="Schizophrenia" & dat.ldsc$Enrichment_p < .05/nrow(dat.ldsc), ]$Category
dat.ldsc[dat.ldsc$Trait=="Bipolar_disorder" & dat.ldsc$Enrichment_p < .05/nrow(dat.ldsc), ]$Category
    # (a lot, but for AMY/NAc:)
        # [1] "amy:Excit.1"          "amy:Excit.2"          "amy:Excit.3"          "amy:Inhib.1"
        # [5] "amy:Inhib.2"          "amy:Inhib.3"          "amy:Inhib.5"          "amy:OPC"
        # [30] "nac:MSN.D1.3"
dat.ldsc[dat.ldsc$Trait=="Major_depressive_disorder" & dat.ldsc$Enrichment_p < .05/nrow(dat.ldsc), ]$Category
    # [1] "amy:Inhib.1"      "amy:Inhib.2"      "dlpfc:Excit.L3:4" "dlpfc:Inhib.5"    "dlpfc:Inhib.6"
    # [6] "hpc:Inhib.2"      "hpc:Inhib.3"      "hpc:Inhib.4"      "hpc:Inhib.5"      "sacc:Inhib.1"

# None for PTSD, Depressive_symptoms, Autism_spectrum_disorder at Bonferroni
dat.ldsc[dat.ldsc$Trait=="PTSD" & dat.ldsc$Enrichment_holm < .05, ]$Category
    # still none
dat.ldsc[dat.ldsc$Trait=="Autism_spectrum_disorder" & dat.ldsc$Enrichment_holm < .05, ]$Category
    # [1] "hpc:OPC"      "sacc:Inhib.1"






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
keepCols <-  c("Barcode","sum","detected","sample","region","donor","processDate","protocol","prelimCluster")

## Adapted from Leo's `create_small_sce`:
    create_small_sce <- function(sce_original, cell_var = "cellType.final") {

          sce_original <- sce_original[, tolower(colData(sce_original)[[cell_var]]) != tolower("ambig.lowNtrxts")]
          colData(sce_original)[[cell_var]] <- factor(colData(sce_original)[[cell_var]])
          colData(sce_original)[[cell_var]] <- factor(colData(sce_original)[[cell_var]])
          
          message(Sys.time(), " reducing the sce object")
          sce_small <- sce_original
          # Just provide raw counts, as opposed to logcounts used for shiny app
          assays(sce_small) <- assays(sce_small)["counts"]
          
          sce_small$sample <- factor(sce_small$sample)
          sce_small$region <- factor(sce_small$region)
          sce_small$donor <- factor(sce_small$donor)
          sce_small$processDate <- factor(paste0(sce_small$processDate,".2019"))
          sce_small$protocol <- factor(sce_small$protocol)
          colData(sce_small) <- colData(sce_small)[ ,colnames(colData(sce_small)) %in% keepCols]
          # Add back in standardized 'cell_type' naming:
          sce_small$cell_type <- factor(colData(sce_original)[[cell_var]])
          
          ## Make the rows more browsable
          colnames(sce_small) <- paste0(sce_small$sample, '_', sce_small$Barcode)
          sce_small$Barcode <- make.names(sce_small$Barcode, unique = TRUE)
          ## Can remove bc just points to location of raw count matrices:
          metadata(sce_small) <- list()
          
          ## It's all "Gene Expression", so we can remove it
          rowData(sce_small)$Type <- NULL
          
          message(Sys.time(), " computing propNucleiExprs")
          
          rowData(sce_small)$propNucleiExprs <- apply(
            assay(sce_original, "counts"),
            1,
            function(x) {
              mean(x != 0)
            }
          )
          # The above, by cell type ===
          cellType.idx <- splitit(colData(sce_original)[[cell_var]])
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
load("rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda", verbose=T)
# Publication tSNE
load("rdas/ztemp_NAc-ALL-n5_SCE-with-tSNEon15-10PCs_MNT.rda", verbose=T)
    rm(sce.all.tsne.10pcs)

# Already has 'ambig.lowNtrxts' removed, so remove these and just assign to $cellType
sce.nac <- sce.nac.all[ ,!sce.nac.all$cellType.final=="ambig.lowNtrxts"]

reducedDim(sce.nac, "TSNE") <- reducedDim(sce.all.tsne.15pcs, "TSNE")

# Reduce
sce.nac.red <- create_small_sce(sce.nac)
# Check
head(colData(sce.nac.red))
head(rowData(sce.nac.red))
plotTSNE(sce.nac.red, colour_by="cell_type")

sce.nac <- sce.nac.red

#dir.create("rdas/forAmazonS3/")
save(sce.nac, file="rdas/forAmazonS3/SCE_NAc_tran-etal.rda")


## AMY === ===
load("rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
# Publication tSNE
load("rdas/ztemp_Amyg-n2_SCE-with-tSNEonOptPCs-minus-PC5_MNT.rda", verbose=T)

# Already has 'Ambig.lowNtrxts' removed, so remove these and just assign to $cellType
sce.amy <- sce.amy[ ,!sce.amy$cellType.split=="Ambig.lowNtrxts"]

reducedDim(sce.amy, "TSNE") <- reducedDim(sce.amy.tsne.optb, "TSNE")

# Reduce
sce.amy.red <- create_small_sce(sce.amy, cell_var="cellType.split")
# Check
head(colData(sce.amy.red))
head(rowData(sce.amy.red))
plotTSNE(sce.amy.red, colour_by="cell_type")









