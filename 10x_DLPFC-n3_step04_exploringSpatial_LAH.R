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
library(limma)
library(lattice)
library(RColorBrewer)
library(pheatmap)
# library(Seurat)


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
sce.dlpfc$PseudoSample = paste0(sce.dlpfc$donor, ":", sce.dlpfc$cellType)
table(sce.dlpfc$PseudoSample)
length(unique(sce.dlpfc$PseudoSample))
# [1] 53

cIndexes = splitit(sce.dlpfc$PseudoSample)
umiComb <- sapply(cIndexes, function(ii)
  rowSums(assays(sce.dlpfc)$counts[, ii, drop = FALSE]))

table(rowSums(umiComb) == 0)
# FALSE  TRUE 
# 29310  4228

phenoComb = colData(sce.dlpfc)[!duplicated(sce.dlpfc$PseudoSample), c(12:21)]
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
    # sce_raw <-
    #   SingleCellExperiment(
    #     list(counts = umiComb),
    #     colData = phenoComb,
    #     rowData = rowData(sce.dlpfc)
    #   )   # of mean(librarySizeFactors(sce_raw)) == 1, so this is ok

#save(sce_pseudobulk, file = "rda/dlpfc_snRNAseq_pseudobulked.Rdata")

## extract expression
mat <- assays(sce_pseudobulk)$logcounts

## Build a group model
#previoulsy prelimClusters
mod <- with(colData(sce_pseudobulk),
            model.matrix(~ 0 + cellType))

# colnames(mod) <- gsub('prelimCluster', '', colnames(mod))

## get duplicate correlation
## was 'sample' swap to 'donor'
corfit <- duplicateCorrelation(mat, mod,
                               block = sce_pseudobulk$donor)
corfit$consensus.correlation
        # 0.02392859
# [1] 0.09367473

#save(corfit, file = "rda/dlpfc_snRNAseq_pseudobulked_dupCor.Rdata")

## Next for each layer test that layer vs the rest
# cell_idx <- splitit(sce_pseudobulk$prelimCluster)
cell_idx <- splitit(sce_pseudobulk$cellType)

eb0_list_cell <- lapply(cell_idx, function(x) {
  res <- rep(0, ncol(sce_pseudobulk))
  res[x] <- 1
  m <- with(colData(sce_pseudobulk),
            model.matrix(~ res))
  eBayes(
    lmFit(
      mat,design = m,
      block = sce_pseudobulk$donor,
      correlation = corfit$consensus.correlation
    )
  )
})

    ## MNT addition: [\n] Warning messages:
     #1: Zero sample variances detected, have been offset away from zero
     #2: Zero sample variances detected, have been offset away from zero

    ## LAH Coefficents not estimabale: res

#save(eb0_list_cell, file = "rda/dlpfc_snRNAseq_pseudobulked_specific_Ts.Rdata")

## Extract the p-values
pvals0_contrasts_cell <- sapply(eb0_list_cell, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts_cell) = rowData(sce_pseudobulk)$gene_id

t0_contrasts_cell <- sapply(eb0_list_cell, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts_cell) = rowData(sce_pseudobulk)$gene_id
fdrs0_contrasts_cell = apply(pvals0_contrasts_cell, 2, p.adjust, 'fdr')

data.frame(
  'FDRsig' = colSums(fdrs0_contrasts_cell < 0.05 &
                       t0_contrasts_cell > 0),
  'Pval10-6sig' = colSums(pvals0_contrasts_cell < 1e-6 &
                            t0_contrasts_cell > 0),
  'Pval10-8sig' = colSums(pvals0_contrasts_cell < 1e-8 &
                            t0_contrasts_cell > 0)
)
#         FDRsig Pval10.6sig Pval10.8sig
# Astro     1296         306         189
# Excit_A    254          56          39
# Excit_B    288          47          28
# Excit_C    191          27          14
# Excit_D     44          12           6
# Excit_E     47           4           2
# Excit_F    694          87          23
# Inhib_A    210          49          24
# Inhib_B    138          36          17
# Inhib_C    105          17           8
# Inhib_D    150          27          19
# Inhib_E      3           0           0
# Inhib_F      2           2           2
# Micro      967         403         272
# Mural       54           7           6
# Oligo    20155        2767        1357
# OPC        502         111          66
# Tcell       89          13          10

############################
### correlate to layer?? ###
############################

## load modeling outputs
load("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/rda/eb_contrasts.Rdata", verbose=T)
    # eb_contrasts

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

## line up
all(rownames(pvals0_contrasts) %in% rownames(pvals0_contrasts_cell))
all(rownames(pvals0_contrasts_cell) %in% rownames(pvals0_contrasts))

common_genes <- intersect(rownames(pvals0_contrasts), rownames(pvals0_contrasts_cell))
length(common_genes)
# [1] 22331

# Layer
pvals0_contrasts = pvals0_contrasts[common_genes, ]
t0_contrasts = t0_contrasts[common_genes, ]
fdrs0_contrasts = fdrs0_contrasts[common_genes, ]

    # SN-PB'd
pvals0_contrasts_cell = pvals0_contrasts_cell[common_genes, ]
t0_contrasts_cell = t0_contrasts_cell[common_genes, ]
fdrs0_contrasts_cell = fdrs0_contrasts_cell[common_genes, ]

cor_t = cor(t0_contrasts_cell, t0_contrasts)
signif(cor_t, 2)

### just layer specific genes from ones left
layer_specific_indices = mapply(function(t) {
  oo = order(t, decreasing = TRUE)[1:100]},
  as.data.frame(t0_contrasts))

layer_ind = unique(as.numeric(layer_specific_indices))
length(layer_ind)

cor_t_layer = cor(t0_contrasts_cell[layer_ind, ],
                  t0_contrasts[layer_ind, ])
signif(cor_t_layer, 2)

### heatmap
theSeq = seq(-.81, .81, by = 0.01)
my.col <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq))

# ct = colData(sce_pseudobulk)
# ct = ct[!duplicated(sce_pseudobulk$collapsedCluster),]
# ct = ct[order(ct$cellType, ct$collapsedCluster),]
# ct$lab = paste0(ct$collapsedCluster, " (", ct$cellType,")")

cor_t_toPlot = cor_t[, c(1, 7:2)]
cor_t_layer_toPlot = cor_t_layer[, c(1, 7:2)]

pdf(here("pdfs/revision/exploration/overlap-ST-dlpfc_with_10x-snRNA-DLPFC-n3_LAH2021.pdf"),
    height=5)
print(
  levelplot(
    cor_t_toPlot,
    aspect = "fill",
    at = theSeq,
    col.regions = my.col,
    ylab = "",
    xlab = "",
    scales = list(x = list(rot = 90, cex = 0.9), y = list(cex = 1.5)),
    main = list("All Genes",side=1,line=0.5)
  )
)

print(
  levelplot(
    cor_t_layer_toPlot,
    aspect = "fill",
    at = theSeq,
    col.regions = my.col,
    ylab = "",
    xlab = "",
    scales = list(x = list(rot = 90, cex = 0.9), y = list(cex = 1.5)),
    main = list("Top 100 Layer Specific Genes",side=1,line=0.5)
  ) 
)


dev.off()

pdf("pdfs/revision/exploration/dlpfc_snRNAseq_overlap_pheatmap_LAH2021.pdf", width = 10)
print(
  pheatmap(
    cor_t_toPlot,
    main = "All Genes"
  )
)
print(
  pheatmap(
    cor_t_layer_toPlot,
    main = "Top 100 Layer Specific Genes"
  )
)
dev.off()

#### gene expression
g = c("SNAP25", "CAMK2A", "GAD2", "SOX11",
      "FOXP2", "PDGFRA", "MBP", "PLP1",
      "AQP4", "GFAP", "CD74")
g2 = rowData(sce.dlpfc)[g,]$gene_id

t0_contrasts_cell_markers = t0_contrasts_cell[g2,]

cc_cell_layer = cor(t(t0_contrasts_cell_markers), cor_t_layer)
rownames(cc_cell_layer) <- g
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


#### Compare with findMarkers ####
load(here("rdas/revision/markers-stats_DLPFC-n3_findMarkers-SN-LEVEL_LAHMay2021.rda"), verbose = TRUE)

## fix
logFC_fm <- sapply(markers.t.1vAll, function(x) {
  x <- x[, "std.logFC", drop = FALSE]
  x <- x[rownames(markers.t.1vAll[["Astro"]]),]
  return(x)
})

t0_fm_cell <- logFC_fm * ncol(sce.dlpfc)

rownames(t0_fm_cell) <- rowData(sce_pseudobulk[rownames(t0_fm_cell)])$gene_id
common_genes2 <- intersect(common_genes, rownames(t0_fm_cell))

t0_fm_cell <- t0_fm_cell[common_genes2, ]
t0_contrasts_cell2 <- t0_contrasts_cell[common_genes2,]

## Marker Finding Corelation 
cor_t_cells = cor(t0_contrasts_cell2, t0_fm_cell)
signif(cor_t_cells, 2)

t0_contrasts2 <- t0_contrasts[common_genes2,]
cor_t_fm = cor(t0_fm_cell, t0_contrasts2)
signif(cor_t_fm, 2)

layer_specific_indices2 <- mapply(function(t) {
  oo = order(t, decreasing = TRUE)[1:100]},
  as.data.frame(t0_contrasts2))

layer_ind2 = unique(as.numeric(layer_specific_indices2))
length(layer_ind2)
# [1] 691

cor_t_layer_fm = cor(t0_fm_cell[layer_ind2, ],
                  t0_contrasts[layer_ind2, ])
signif(cor_t_layer_fm, 2)

## Top100 cell Type genes
cellType_specific_indices <- mapply(function(t) {
  oo = order(t, decreasing = TRUE)[1:100]},
  as.data.frame(t0_fm_cell))

cellType_ind = unique(as.numeric(cellType_specific_indices))
length(cellType_ind)
# [1] 1402
cor_t_layer_ct_fm = cor(t0_fm_cell[cellType_ind, ],
                     t0_contrasts[cellType_ind, ])
signif(cor_t_layer_ct_fm, 2)


cor_fm <- list(cor_t_cells, cor_t_fm[, c(1, 7:2)], cor_t_layer_fm[, c(1, 7:2)], cor_t_layer_ct_fm[, c(1, 7:2)])
titles <- list("manual vs. findMarkers", "All Genes" , "Top 100 Layer Specific Genes", "Top 100 cell type specific genes")

theSeq3 = seq(-.12, .12, by = 0.005)
my.col3 <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq3))

pdf(here("pdfs/revision/exploration/overlap-ST-dlpfc_with_10x-snRNA-DLPFC-n3_findMarkers_LAH2021.pdf"),
    height=5)
for(i in 1:length(cor_fm)){
  print(
    levelplot(
      cor_fm[[i]],
      aspect = "fill",
      at = theSeq,
      col.regions = my.col,
      ylab = "",
      xlab = "",
      scales = list(x = list(rot = 90, cex = 0.9), y = list(cex = 1.5)),
      main = list(titles[[i]],side=1,line=0.5)
    )
  )
}

for(i in 3:length(cor_fm)){
  print(
    levelplot(
      cor_fm[[i]],
      aspect = "fill",
      at = theSeq3,
      col.regions = my.col3,
      ylab = "",
      xlab = "",
      scales = list(x = list(rot = 90, cex = 0.9), y = list(cex = 1.5)),
      main = list(titles[[i]],side=1,line=0.5)
    )
  )
}

dev.off()

