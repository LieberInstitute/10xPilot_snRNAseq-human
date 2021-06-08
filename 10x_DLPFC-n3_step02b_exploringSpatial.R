### MNT 10x snRNA-seq workflow: step 02b
###   **Region-specific analyses**
###     - (3x) DLPFC samples from: Br5161, Br5212, Br207
###   **Comparing to Spatial Transcriptomics data
### Taken/adapted from /dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses
### Initiated MNT 17Feb2020
### Updated LAH 08Jun2021
#####################################################################

library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(batchelor)
library(DropletUtils)
library(jaffelab)
library(here)
## Pasted from the reference script (see header), as of 10:45 17Feb2020

# library(jaffelab)
# library(Seurat)
# #library(scater)
# #library(DropletUtils)
# library(limma)
# library(lattice)
# library(RColorBrewer)
# library(pheatmap)



load(here("rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda"), verbose=T)
# sce.dlpfc
# chosen.hvgs.dlpfc
# pc.choice.dlpfc
# clusterRefTab.dlpfc
# ref.sampleInfo
# annotationTab.dlpfc
# cell_colors
    
########################
### broad clusters #####      - well actually 'prelim' clusters
########################

## get pseudobulk
sce.dlpfc$PseudoSample = paste0(sce.dlpfc$donor, ":", sce.dlpfc$prelimCluster)

cIndexes = splitit(sce.dlpfc$PseudoSample)
umiComb <- sapply(cIndexes, function(ii)
  rowSums(assays(sce.dlpfc)$counts[, ii, drop = FALSE]))

#phenoComb = colData(sce.dlpfc)[!duplicated(sce.dlpfc$PseudoSample), 13:15]
    phenoComb = colData(sce.dlpfc)[!duplicated(sce.dlpfc$PseudoSample), c(13:21)]

    # MTN comment: I think this is missing
    phenoComb$PseudoSample <- sce.dlpfc$PseudoSample[!duplicated(sce.dlpfc$PseudoSample)]

rownames(phenoComb) = phenoComb$PseudoSample
phenoComb = phenoComb[colnames(umiComb), ]
phenoComb = DataFrame(phenoComb)

sce_pseudobulk <-
  logNormCounts(SingleCellExperiment(
    list(counts = umiComb),
    colData = phenoComb,
    rowData = rowData(sce.dlpfc)
  ))
    
    # Check not treating LSFs weird:
    sce_raw <-
      SingleCellExperiment(
        list(counts = umiComb),
        colData = phenoComb,
        rowData = rowData(sce.dlpfc)
      )   # of mean(librarySizeFactors(sce_raw)) == 1, so this is ok

#save(sce_pseudobulk, file = "rda/dlpfc_snRNAseq_pseudobulked.Rdata")

###############################
## extract expression
mat <- assays(sce_pseudobulk)$logcounts

## Build a group model
mod <- with(colData(sce_pseudobulk),
            model.matrix(~ 0 + prelimCluster))
colnames(mod) <- gsub('prelimCluster', '', colnames(mod))

## get duplicate correlation
corfit <- duplicateCorrelation(mat, mod,
                               block = sce_pseudobulk$sample)
    corfit$consensus.correlation
        # 0.02392859

#save(corfit, file = "rda/dlpfc_snRNAseq_pseudobulked_dupCor.Rdata")

## Next for each layer test that layer vs the rest
cell_idx <- splitit(sce_pseudobulk$prelimCluster)

eb0_list_cell <- lapply(cell_idx, function(x) {
  res <- rep(0, ncol(sce_pseudobulk))
  res[x] <- 1
  m <- with(colData(sce_pseudobulk),
            model.matrix(~ res))
  eBayes(
    lmFit(
      mat,design = m,
      block = sce_pseudobulk$sample,
      correlation = corfit$consensus.correlation
    )
  )
})

    ## MNT addition: [\n] Warning messages:
     #1: Zero sample variances detected, have been offset away from zero
     #2: Zero sample variances detected, have been offset away from zero

#save(eb0_list_cell, file = "rda/dlpfc_snRNAseq_pseudobulked_specific_Ts.Rdata")

##########
## Extract the p-values
pvals0_contrasts_cell <- sapply(eb0_list_cell, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts_cell) = rownames(mat)

t0_contrasts_cell <- sapply(eb0_list_cell, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts_cell) = rownames(mat)
fdrs0_contrasts_cell = apply(pvals0_contrasts_cell, 2, p.adjust, 'fdr')

data.frame(
  'FDRsig' = colSums(fdrs0_contrasts_cell < 0.05 &
                       t0_contrasts_cell > 0),
  'Pval10-6sig' = colSums(pvals0_contrasts_cell < 1e-6 &
                            t0_contrasts_cell > 0),
  'Pval10-8sig' = colSums(pvals0_contrasts_cell < 1e-8 &
                            t0_contrasts_cell > 0)
)

# FDRsig Pval10.6sig Pval10.8sig
# 1    1646         770         557
# 2     346         147         103
# 3     315         153         127
# 4     306         118          66
# 5     846         288         165
# 6     247          92          48
# 7     965         281         170
# 8     308         183         142
# 9     302         136          94
# 10    381         159          96
# 11    218          91          51
# 12    356         174         123
# 13    307         177         129
# 14    292         133         104
# 15    207          96          48
# 16    312         114          65
# 17    251         179          31
# 18    198          94          65
# 19    426         149          98
# 20    192          57          17
# 21    173          56          27
# 22    229          86          39
# 23    192          69          39
# 24    609         177         105
# 25    408         153          84
# 26    236          71          28
# 27    274         111          63
# 28    261         119          81
# 29    178          64          33
# 30    318         106          60
# 31    396         157          96

############################
### correlate to layer?? ###
############################

###################
## load modeling outputs
#load("rda/eb_contrasts.Rdata")
load("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/rda/eb_contrasts.Rdata", verbose=T)
    # eb_contrasts

#load("rda/eb0_list.Rdata")
load("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/rda/eb0_list.Rdata", verbose=T)
    # eb0_list

## Extract the p-values
pvals0_contrasts <- sapply(eb0_list, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts) = rownames(eb_contrasts)
fdrs0_contrasts = apply(pvals0_contrasts, 2, p.adjust, "fdr")

## Extract the t-stats
t0_contrasts <- sapply(eb0_list, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts) = rownames(eb_contrasts)

############
# line up ##

mm = match(rownames(pvals0_contrasts), rowData(sce_pseudobulk)$ID)

    # Layer
pvals0_contrasts = pvals0_contrasts[!is.na(mm), ]
t0_contrasts = t0_contrasts[!is.na(mm), ]
fdrs0_contrasts = fdrs0_contrasts[!is.na(mm), ]

    # SN-PB'd
pvals0_contrasts_cell = pvals0_contrasts_cell[mm[!is.na(mm)], ]
t0_contrasts_cell = t0_contrasts_cell[mm[!is.na(mm)], ]
fdrs0_contrasts_cell = fdrs0_contrasts_cell[mm[!is.na(mm)], ]

cor_t = cor(t0_contrasts_cell, t0_contrasts)
signif(cor_t, 2)

### just layer specific genes from ones left
layer_specific_indices = mapply(function(t, p) {
  oo = order(t, decreasing = TRUE)[1:100]
},
as.data.frame(t0_contrasts),
as.data.frame(pvals0_contrasts))
layer_ind = unique(as.numeric(layer_specific_indices))

cor_t_layer = cor(t0_contrasts_cell[layer_ind, ],
                  t0_contrasts[layer_ind, ])
signif(cor_t_layer, 2)

### heatmap
theSeq = seq(-.81, .81, by = 0.01)
my.col <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq))

ct = colData(sce_pseudobulk)
ct = ct[!duplicated(sce_pseudobulk$prelimCluster),]
ct = ct[order(ct$cellType, ct$prelimCluster),]
ct$lab = paste0(ct$prelimCluster, " (", ct$cellType,")")

#cor_t_layer_toPlot = cor_t_layer[as.character(ct$prelimCluster), c(1, 7:2)]
    # MNT add, instead
    cor_t_layer_toPlot = cor_t_layer[setdiff(as.character(ct$prelimCluster),"24"), c(1, 7:2)]

#rownames(cor_t_layer_toPlot) = ct$lab
    # MNT add, instead
    rownames(cor_t_layer_toPlot) = ct$lab[-grep("Ambig.lowNtrxts", ct$lab)]
    

#pdf("pdf/dlpfc_snRNAseq_overlap_heatmap.pdf", width = 10)
#pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/overlap-ST-dlpfc_with_10x-snRNA-seq_n2_MNTannotations_17Feb2020.pdf",
pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/overlap-ST-dlpfc_with_10x-snRNA-seq_n2_MNTannotations_24Apr2020.pdf",
    height=5)
print(
  levelplot(
    cor_t_layer_toPlot,
    aspect = "fill",
    at = theSeq,
    col.regions = my.col,
    ylab = "",
    xlab = "",
    scales = list(x = list(rot = 90, cex = 0.9), y = list(cex = 1.5))
  )
)
dev.off()

#pdf("pdf/dlpfc_snRNAseq_overlap_pheatmap.pdf", width = 10)
print(
  pheatmap(
    cor_t_layer_toPlot,
    ylab = "",
    xlab = "",
  )
)
#dev.off()

#### gene expression
g = c("SNAP25", "CAMK2A", "GAD2", "SOX11",
      "FOXP2", "PDGFRA", "MBP", "PLP1",
      "AQP4", "GFAP", "CD74")
t0_contrasts_cell_markers = t0_contrasts_cell[g,]

cc_cell_layer = cor(t(t0_contrasts_cell_markers), cor_t_layer)
signif(cc_cell_layer,3)

### heatmap
theSeq2 = seq(-5, 5, by = 0.01)
my.col2 <- colorRampPalette(brewer.pal(7, "RdBu"))(length(theSeq2))
t0_contrasts_cell_markers_plot = t(t0_contrasts_cell_markers)
t0_contrasts_cell_markers_plot[t0_contrasts_cell_markers_plot > 5] = 5
t0_contrasts_cell_markers_plot = t0_contrasts_cell_markers_plot[,ncol(t0_contrasts_cell_markers_plot):1]

#pdf("pdf/dlpfc_snRNAseq_marker_heatmap.pdf", width = 10)
print(
  levelplot(
    t0_contrasts_cell_markers_plot,
    aspect = "fill",
    at = theSeq2,
    col.regions = my.col2,
    ylab = "",
    xlab = "",
    scales = list(x = list(rot = 90, cex = 1.5), y = list(cex = 1.5))
  )
)
#dev.off()

########################################
#### specific /collapsed clusters ######
########################################

      ## (Not sure if this was completed? - MNT 25Mar2020)



### MNT addition 05May2020 ===================================
  # What if compare markers of the annotations from the above, but 
  #     tested at the single-nucleus level (i.e. `findMarkers()`)?

load("rdas/markers-stats_DLPFC_n2_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.dlpfc.t.1vAll, markers.dlpfc.t.design
    # ^ t-statistics can be computed with the formula [Cohen's] d = t/sqrt(N)


## Following the above, but with:
    ############
    # line up ##
    
    mm = match(rownames(pvals0_contrasts), rowData(sce.dlpfc.st)$ID)
    
    # Layer
    pvals0_contrasts = pvals0_contrasts[!is.na(mm), ]
    t0_contrasts = t0_contrasts[!is.na(mm), ]
    fdrs0_contrasts = fdrs0_contrasts[!is.na(mm), ]
    
    
## Single-nucleus-level stats
    
# ** FIRST WILL HAVE TO FIX ROW ORDER IN 'markers.dlpfc.t.1vAll'
markers.dlpfc.t.1vAll <- lapply(markers.dlpfc.t.1vAll, function(x){
  x[rownames(sce.dlpfc.st), ]
})

# Then get t's, re-computing by multiplying Cohen's D (std.logFC) * N
t0_contrasts_cell <- sapply(markers.dlpfc.t.1vAll, function(x){
  x$std.logFC * ncol(sce.dlpfc.st)
})
rownames(t0_contrasts_cell) <- rowData(sce.dlpfc.st)$ID

#pvals0_contrasts_cell = pvals0_contrasts_cell[mm[!is.na(mm)], ]    - don't need
t0_contrasts_cell = t0_contrasts_cell[mm[!is.na(mm)], ]
#fdrs0_contrasts_cell = fdrs0_contrasts_cell[mm[!is.na(mm)], ]      - don't need

# Re-order col order - for some reason isn't maintained as in pairwise stats..
      ## ahhh this is the order of `unique(sce.dlpfc.st$cellType.split)` (iterative contrasts)
t0_contrasts_cell <- t0_contrasts_cell[ ,names(markers.dlpfc.t.design)]


cor_t = cor(t0_contrasts_cell, t0_contrasts)
signif(cor_t, 2)
    #                    WM Layer1 Layer2 Layer3 Layer4  Layer5 Layer6
    # Astro          -0.093  0.390  0.040  0.018 -0.075 -0.1100 -0.120
    # Excit.ambig    -0.520 -0.140  0.440  0.420  0.200  0.1300  0.090
    # Excit.L2:3     -0.510 -0.140  0.410  0.440  0.220  0.1300  0.075
    # Excit.L3:4     -0.510 -0.180  0.250  0.380  0.340  0.2500  0.057
    # Excit.L4:5     -0.530 -0.230  0.220  0.290  0.350  0.3900  0.140
    # Excit.L5       -0.380 -0.150  0.150  0.180  0.200  0.3100  0.130
    # Excit.L5:6     -0.440 -0.190  0.240  0.230  0.170  0.2600  0.260
    # Excit.L6.broad -0.500 -0.190  0.300  0.310  0.210  0.2300  0.240
    # Inhib.1        -0.310 -0.019  0.190  0.160  0.130  0.1500  0.035
    # Inhib.2        -0.330 -0.096  0.160  0.180  0.190  0.2000  0.069
    # Inhib.3        -0.360 -0.070  0.180  0.210  0.230  0.2000  0.025
    # Inhib.4        -0.460 -0.120  0.200  0.280  0.300  0.2700  0.048
    # Inhib.5        -0.440 -0.096  0.250  0.250  0.220  0.2500  0.072
    # Inhib.6        -0.460 -0.049  0.290  0.260  0.220  0.2200  0.036
    # Micro           0.039  0.150 -0.052 -0.040 -0.068 -0.0800 -0.059
    # Oligo           0.520 -0.088 -0.290 -0.300 -0.190 -0.1700 -0.031
    # OPC            -0.130  0.091  0.099  0.072  0.016  0.0059 -0.019

    ## of range -0.533  0.520;;  with layer-specific markers (690 genes): -0.593  0.653

### heatmap
theSeq = seq(-.71, .71, by = 0.01)
my.col <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq))

print(
  levelplot(
    cor_t,
    aspect = "fill",
    at = theSeq,
    col.regions = my.col,
    ylab = "",
    xlab = "",
    scales = list(x = list(rot = 90, cex = 0.9), y = list(cex = 1.5)),
    main="Single-nucleus-level marker (cluster-vs-all) t's vs layer t's \n (with spatial annotations)"
  )
)




### Test sACC stats vs DLPFC layers ===========================================
  # ** just try the 'collapsed' clusters (use $cellType), since these aren't 'broad'
  #    - also bc there are 46 prelim clusters
load("rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    # sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc, ref.sampleInfo

# First remove 'Ambig.lowNtrxts' (43 nuclei)
sce.sacc <- sce.sacc[ ,sce.sacc$cellType != "Ambig.lowNtrxts"]
sce.sacc$cellType <- droplevels(sce.sacc$cellType)

table(sce.sacc$cellType, sce.sacc$sample)
    #         sacc.5161 sacc.5212
    # Astro         188       435
    # Excit.1       187       399
    # Excit.2       151       259
    # Excit.3        68       117
    # Excit.4        38        47
    # Inhib.1       142       373
    # Inhib.2        82       233
    # Micro         235       267
    # Oligo        1852      1400
    # OPC           227       304

## get pseudobulk
sce.sacc$PseudoSample = paste0(sce.sacc$sample, ":", sce.sacc$cellType)

cIndexes = splitit(sce.sacc$PseudoSample)
umiComb <- sapply(cIndexes, function(ii)
  rowSums(assays(sce.sacc)$counts[, ii, drop = FALSE]))

phenoComb = colData(sce.sacc)[!duplicated(sce.sacc$PseudoSample), c(13:21)]

rownames(phenoComb) = phenoComb$PseudoSample
phenoComb = phenoComb[colnames(umiComb), ]
phenoComb = DataFrame(phenoComb)

sce_pseudobulk <-
  logNormCounts(SingleCellExperiment(
    list(counts = umiComb),
    colData = phenoComb,
    rowData = rowData(sce.sacc)
  ))


## extract expression
mat <- assays(sce_pseudobulk)$logcounts

## Build a group model
mod <- with(colData(sce_pseudobulk),
            model.matrix(~ 0 + cellType))
colnames(mod) <- gsub('cellType', '', colnames(mod))

## get duplicate correlation
corfit <- duplicateCorrelation(mat, mod,
                               block = sce_pseudobulk$sample)
corfit$consensus.correlation
    # 0.0113936

#save(corfit, file = "rda/dlpfc_snRNAseq_pseudobulked_dupCor.Rdata")

## Next for each layer test that layer vs the rest
cell_idx <- splitit(sce_pseudobulk$cellType)

eb0_list_cell <- lapply(cell_idx, function(x) {
  res <- rep(0, ncol(sce_pseudobulk))
  res[x] <- 1
  m <- with(colData(sce_pseudobulk),
            model.matrix(~ res))
  eBayes(
    lmFit(
      mat,design = m,
      block = sce_pseudobulk$sample,
      correlation = corfit$consensus.correlation
    )
  )
})

    ## Warning messages:
    #1: Zero sample variances detected, have been offset away from zero
    #2: Zero sample variances detected, have been offset away from zero


## Extract the p-values
pvals0_contrasts_cell <- sapply(eb0_list_cell, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts_cell) = rownames(mat)

t0_contrasts_cell <- sapply(eb0_list_cell, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts_cell) = rownames(mat)
fdrs0_contrasts_cell = apply(pvals0_contrasts_cell, 2, p.adjust, 'fdr')

data.frame(
  'FDRsig' = colSums(fdrs0_contrasts_cell < 0.05 &
                       t0_contrasts_cell > 0),
  'Pval10-6sig' = colSums(pvals0_contrasts_cell < 1e-6 &
                            t0_contrasts_cell > 0),
  'Pval10-8sig' = colSums(pvals0_contrasts_cell < 1e-8 &
                            t0_contrasts_cell > 0)
)

    #        FDRsig Pval10.6sig Pval10.8sig
    # Astro     1134         195          78
    # Excit.1    156          80          62
    # Excit.2    112          45          24
    # Excit.3     67          29          16
    # Excit.4    176          46          26
    # Inhib.1     91          37          23
    # Inhib.2    135          49          23
    # Micro     3900         786         396
    # Oligo     1271         218          86
    # OPC        329          62          39




### correlate to layer?? === === === ===

## load modeling outputs
#load("rda/eb_contrasts.Rdata")
load("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/rda/eb_contrasts.Rdata", verbose=T)
# eb_contrasts

#load("rda/eb0_list.Rdata")
load("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/rda/eb0_list.Rdata", verbose=T)
# eb0_list

## Extract the p-values
pvals0_contrasts <- sapply(eb0_list, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts) = rownames(eb_contrasts)
fdrs0_contrasts = apply(pvals0_contrasts, 2, p.adjust, "fdr")

## Extract the t-stats
t0_contrasts <- sapply(eb0_list, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts) = rownames(eb_contrasts)

############
# line up ##

mm = match(rownames(pvals0_contrasts), rowData(sce_pseudobulk)$ID)

# Layer
pvals0_contrasts = pvals0_contrasts[!is.na(mm), ]
t0_contrasts = t0_contrasts[!is.na(mm), ]
fdrs0_contrasts = fdrs0_contrasts[!is.na(mm), ]

# SN-PB'd
pvals0_contrasts_cell = pvals0_contrasts_cell[mm[!is.na(mm)], ]
t0_contrasts_cell = t0_contrasts_cell[mm[!is.na(mm)], ]
fdrs0_contrasts_cell = fdrs0_contrasts_cell[mm[!is.na(mm)], ]

cor_t = cor(t0_contrasts_cell, t0_contrasts)
signif(cor_t, 2)

### just layer specific genes from ones left
layer_specific_indices = mapply(function(t, p) {
  oo = order(t, decreasing = TRUE)[1:100]
},
as.data.frame(t0_contrasts),
as.data.frame(pvals0_contrasts))
layer_ind = unique(as.numeric(layer_specific_indices))

cor_t_layer = cor(t0_contrasts_cell[layer_ind, ],
                  t0_contrasts[layer_ind, ])
signif(cor_t_layer, 3)

### heatmap
theSeq = seq(-.81, .81, by = 0.01)
my.col <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq))

ct = colData(sce_pseudobulk)
ct = ct[!duplicated(sce_pseudobulk$cellType),]
ct = ct[order(ct$cellType, ct$cellType),]
ct$lab = paste0(ct$cellType, " (", ct$cellType,")")

cor_t_layer_toPlot = cor_t_layer[as.character(ct$cellType), c(1, 7:2)]
rownames(cor_t_layer_toPlot) = ct$lab


pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/overlap-ST-dlpfc_with_10x-snRNAseq_sACC_MNTannotations_25Mar2020.pdf")
print(
  levelplot(
    cor_t_layer_toPlot,
    aspect = "fill",
    at = theSeq,
    col.regions = my.col,
    ylab = "",
    xlab = "",
    scales = list(x = list(rot = 90, cex = 0.9), y = list(cex = 1.5))
  )
)
dev.off()











### MNT exploration for RNA-scope ==============================================
## Looking at requested genes from Stephanie Page: CHL1 & CNTN4 are enriched in any clusters
 #    * Using this framework to look at prelimClusters, since it seems some unique 'sub-'clusters
 #      in fact do pertain to particular layers

"UCHL1" %in% rownames(fdrs0_contrasts_cell)

fdrs0_contrasts_cell["UCHL1", ]
fdrs0_contrasts_cell["CNTN4", ]

plotExpression(sce.dlpfc, x="cellType", features="UCHL1")
plotExpression(sce.dlpfc, x="cellType", features="CNTN4")


print(
  plotExpression(sce.dlpfc, exprs_values = "logcounts", features=c("UCHL1", "CNTN4"),
                 x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,
                 add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                              geom = "crossbar", width = 0.3,
                                              colour=rep(tableau10medium[1:7], 2))
)
    ## So these genes are in both excit. & inhib. (& OPCs) - are they in certain sub-clusters?

sce.sub <- sce.dlpfc[ ,which(sce.dlpfc$cellType %in% c("Excit", "Inhib", "OPC"))]

sce.sub$prelimWithAnnot <- paste0(sce.sub$cellType, ".", sce.sub$prelimCluster)




pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/expression_UCHL1-CNTN4_10x-DLPFC_clusters_MNT.pdf",
    height=5, width=10)
# Collapsed clusters
print(
  plotExpression(sce.dlpfc, exprs_values = "logcounts", features=c("UCHL1", "CNTN4"),
                 x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,
                 add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                              geom = "crossbar", width = 0.3,
                                              colour=rep(tableau10medium[1:7], 2))
)
# Prelim clusters for nuclei that are neurons or OPCs
print(
  plotExpression(sce.sub, exprs_values = "logcounts", features=c("UCHL1", "CNTN4"),
                 x="prelimWithAnnot", colour_by="prelimWithAnnot", point_alpha=0.5, point_size=.7,
                 add_legend=F) + theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                 geom = "crossbar", width = 0.3, colour="black")
)
dev.off()




