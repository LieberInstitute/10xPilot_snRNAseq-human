### MNT 10x snRNA-seq workflow: step 04 - downstream comparisons
###   **Region-specific analyses**
###     - (3x) NAc samples from: Br5161 & Br5212 & Br5287
###     - (2x) NeuN-sorted samples from: Br5207 & Br5182
###   * Comparison to Jeremy Day Lab's rat NAc samples (n=4)
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

### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

# ===


## load modeling outputs
# 10x-pilot human NAc
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_NAc_all-n5_manualContrasts_MNTApr2020.rda", verbose=T)
# eb_list.nac.all, sce.nac.all.PB, corfit.nac.all

# Day Lab Rat NAc
load("/dcl01/ajaffe/data/lab/singleCell/day_rat_snRNAseq/markers-stats_DayLab-ratNAc_manualContrasts_MNTApr2020.rda", verbose=T)
# eb_list.nac.rat, sce.nac.rat.PB, corfit.nac.rat


# Add EntrezID for human
hs.entrezIds <- mapIds(org.Hs.eg.db, keys=rowData(sce.nac.all.PB)$ID, 
                       column="ENTREZID", keytype="ENSEMBL")
# "'select()' returned 1:many mapping between keys and columns"
table(!is.na(hs.entrezIds))
    # 20,931 valid entries (remember this is already subsetted for those non-zero genes only)

# Add to rowData
rowData(sce.nac.all.PB) <- cbind(rowData(sce.nac.all.PB), hs.entrezIds)


## Bring in 'HomoloGene.ID' for human (already in rowData for rat SCE) ===
## JAX annotation info
hom = read.delim("http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt",
                 as.is=TRUE)

hom_hs <- hom[hom$Common.Organism.Name == "human", ]
    # of 19,124 entries
table(rowData(sce.nac.all.PB)$hs.entrezIds %in% hom_hs$EntrezGene.ID)
    # 17,491
table(rowData(sce.nac.all.PB)$Symbol %in% hom_hs$Symbol)
    # 17,150 - not a bad difference

        # So for mapping === == === ===
        # human.entrez > HomoloGene.ID < rat.Symbol
        #                ^ filter SCE's on this

# Human (by Entrez)
rowData(sce.nac.all.PB)$HomoloGene.ID <- hom_hs$HomoloGene.ID[match(rowData(sce.nac.all.PB)$hs.entrezIds,
                                                                 hom_hs$EntrezGene.ID)]


## Now set/match to shared homologous genes ===

length(intersect(rowData(sce.nac.all.PB)$HomoloGene.ID,
                 rowData(sce.nac.rat.PB)$HomoloGene.ID))  # 14,122 (well, 14,121)

sharedHomologs <- intersect(rowData(sce.nac.all.PB)$HomoloGene.ID,
                            rowData(sce.nac.rat.PB)$HomoloGene.ID)
# That first one is NA lol - get rid of it
sharedHomologs <- sharedHomologs[-1]

# Human not in rat
length(setdiff(rowData(sce.nac.all.PB)$HomoloGene.ID,
                 rowData(sce.nac.rat.PB)$HomoloGene.ID))  # 3187
# Rat not in human
length(setdiff(rowData(sce.nac.rat.PB)$HomoloGene.ID,
               rowData(sce.nac.all.PB)$HomoloGene.ID))  # 1533


# Subset for those
sce.rat.PBsub <- sce.nac.rat.PB[rowData(sce.nac.rat.PB)$HomoloGene.ID %in% sharedHomologs, ]   # 14247
sce.hsap.PBsub <- sce.nac.all.PB[rowData(sce.nac.all.PB)$HomoloGene.ID %in% sharedHomologs, ]  # 14178
    ## Many are duplicated...

rowData(sce.rat.PBsub)$Symbol[duplicated(rowData(sce.rat.PBsub)$HomoloGene.ID)]
    # shoot many genes are orthologs
rowData(sce.hsap.PBsub)$Symbol[duplicated(rowData(sce.hsap.PBsub)$HomoloGene.ID)]
    # same here, though less


### -> Take the higher-expressing of the duplicated - just mean across PB clusters:

    ## Rat ===
    duplicatedSet.rat <- which(duplicated(rowData(sce.rat.PBsub)$HomoloGene.ID))
    genes2compare.rat <- list()
    gene2keep.rat <- character()
    for(g in 1:length(duplicatedSet.rat)){
      genes2compare.rat[[g]] <- rownames(sce.rat.PBsub)[rowData(sce.rat.PBsub)$HomoloGene.ID ==
                                              rowData(sce.rat.PBsub)$HomoloGene.ID[duplicatedSet.rat[g]]]
      rowmeansmat <- rowMeans(assay(sce.rat.PBsub[genes2compare.rat[[g]], ], "logcounts"))
      gene2keep.rat[g] <- names(rowmeansmat[order(rowmeansmat, decreasing=TRUE)])[1]
    }
    
    # Now pull out those that not being compared, so can `c()`
    table(rownames(sce.rat.PBsub) %in% unlist(genes2compare.rat)) # 234   - why isn't this ==
    sum(lengths(genes2compare.rat))                               # 312 ????
    length(unique(unlist(genes2compare.rat))) # 234   - oh. also `length(unique(gene2keep.rat)) == 108`
    
    genesNoCompare.rat <- rownames(sce.rat.PBsub)[!(rownames(sce.rat.PBsub) %in% unlist(genes2compare.rat))]
    
    # Finally combine and subset
    sce.rat.PBsub <- sce.rat.PBsub[c(genesNoCompare.rat, unique(gene2keep.rat)), ]
    
    table(rowData(sce.rat.PBsub)$HomoloGene.ID %in% sharedHomologs) # 14121 TRUE
    table(duplicated(rowData(sce.rat.PBsub)$HomoloGene.ID)) # 14121 FALSE         dope.

    
    ## Human ===
    # First change rownames to EnsemblID
    rowData(sce.hsap.PBsub)$Symbol.unique <- rownames(sce.hsap.PBsub)
    rownames(sce.hsap.PBsub) <- rowData(sce.hsap.PBsub)$ID
        
    duplicatedSet.hsap <- which(duplicated(rowData(sce.hsap.PBsub)$HomoloGene.ID))
    genes2compare.hsap <- list()
    gene2keep.hsap <- character()
    for(g in 1:length(duplicatedSet.hsap)){
      genes2compare.hsap[[g]] <- rownames(sce.hsap.PBsub)[rowData(sce.hsap.PBsub)$HomoloGene.ID ==
                                                          rowData(sce.hsap.PBsub)$HomoloGene.ID[duplicatedSet.hsap[g]]]
      rowmeansmat <- rowMeans(assay(sce.hsap.PBsub[genes2compare.hsap[[g]], ], "logcounts"))
      gene2keep.hsap[g] <- names(rowmeansmat[order(rowmeansmat, decreasing=TRUE)])[1]
    }
    
    # Now pull out those that not being compared, so can `c()`
    table(rownames(sce.hsap.PBsub) %in% unlist(genes2compare.hsap)) # 112   - why isn't this ==
    sum(lengths(genes2compare.hsap))                               # 118 ????
    length(unique(unlist(genes2compare.hsap))) # 112   - oh. also `length(unique(gene2keep.hsap)) == 55`
    
    genesNoCompare.hsap <- rownames(sce.hsap.PBsub)[!(rownames(sce.hsap.PBsub) %in% unlist(genes2compare.hsap))]
        # of length 14066 (which + 55 == 14121)
    
    # Finally combine and subset
    sce.hsap.PBsub <- sce.hsap.PBsub[c(genesNoCompare.hsap, unique(gene2keep.hsap)), ]
    
    table(rowData(sce.hsap.PBsub)$HomoloGene.ID %in% sharedHomologs) # 14121 TRUE
    table(duplicated(rowData(sce.hsap.PBsub)$HomoloGene.ID)) # 14121 FALSE         dope.


    ## Match order and save
    sce.rat.PBsub <- sce.rat.PBsub[match(rowData(sce.hsap.PBsub)$HomoloGene.ID,
                                   rowData(sce.rat.PBsub)$HomoloGene.ID), ]

    table(rowData(sce.rat.PBsub)$HomoloGene.ID == rowData(sce.hsap.PBsub)$HomoloGene.ID)
        # all TRUE - good
    pheatmap(cor(assay(sce.rat.PBsub, "logcounts"), assay(sce.hsap.PBsub, "logcounts")), fontsize=5)
        # (ah but this is at the sample:cluster level)
    
    Readme <- "These two SCEs are subsetted and ordered for matching HomoloGene.ID in the rowData. This can be used to subset the nucleus-level SCEs in their respected Rdata files."
    save(sce.rat.PBsub, sce.hsap.PBsub, Readme, file="/dcl01/ajaffe/data/lab/singleCell/day_rat_snRNAseq//SCE_rat-NAc-PBd_w_matchingHsap-NAc-PBd_HomoloGene.IDs_MNT.rda")
    
    
    
### FINALLY resume comparisons === === === === ===

## Rat stats
pvals_rat <- sapply(eb_list.nac.rat, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals_rat) = rownames(sce.nac.rat.PB)

ts_rat <- sapply(eb_list.nac.rat, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(ts_rat) = rownames(sce.nac.rat.PB)



## Human stats
pvals_hsap <- sapply(eb_list.nac.all, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals_hsap) = rowData(sce.nac.all.PB)$ID

ts_hsap <- sapply(eb_list.nac.all, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(ts_hsap) =rowData(sce.nac.all.PB)$ID



### Subset and check matching 'HomoloGene.ID' === === === ===
pvals_rat <- pvals_rat[rownames(sce.rat.PBsub), ]
ts_rat <- ts_rat[rownames(sce.rat.PBsub), ]

pvals_hsap <- pvals_hsap[rowData(sce.hsap.PBsub)$ID, ]
ts_hsap <- ts_hsap[rowData(sce.hsap.PBsub)$ID, ]

rownames(ts_rat) <- rowData(sce.rat.PBsub)$HomoloGene.ID
rownames(pvals_rat) <- rowData(sce.rat.PBsub)$HomoloGene.ID

rownames(ts_hsap) <- rowData(sce.hsap.PBsub)$HomoloGene.ID
rownames(pvals_hsap) <- rowData(sce.hsap.PBsub)$HomoloGene.ID

table(rownames(ts_rat) == rownames(ts_hsap))
    ## all 14121 TRUE - good


## Now run correlation
cor_t = cor(ts_rat, ts_hsap)
signif(cor_t, 2)

## On just hsap cluster-specific homologous genes ===
hsap_specific_indices = mapply(function(t, p) {
    oo = order(t, decreasing = TRUE)[1:100]
  },
  as.data.frame(ts_hsap),
  as.data.frame(pvals_hsap)
)
hsap_ind = unique(as.numeric(hsap_specific_indices))

cor_t_hsap = cor(ts_rat[hsap_ind, ],
                  ts_hsap[hsap_ind, ])
signif(cor_t_hsap, 3)

## On just rat cluster-specific homologous genes ===
rat_specific_indices = mapply(function(t, p) {
    oo = order(t, decreasing = TRUE)[1:100]
  },
  as.data.frame(ts_rat),
  as.data.frame(pvals_rat)
)
rat_ind = unique(as.numeric(rat_specific_indices))

cor_t_rat = cor(ts_rat[rat_ind, ],
                 ts_hsap[rat_ind, ])
signif(cor_t_rat, 3)



### Heatmap
theSeq.all = seq(-.4, .4, by = 0.01)
my.col.all <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq))

ct = colData(sce.hsap.PBsub)
ct = ct[!duplicated(sce.hsap.PBsub$cellType.final), ]

cor_t_hsap_toPlot = cor_t_hsap
rownames(cor_t_hsap_toPlot) = paste0(rownames(cor_t_hsap_toPlot),"_","R.nor")
colnames(cor_t_hsap_toPlot) = paste0(colnames(cor_t_hsap_toPlot),"_","H.sap")

cor_t_rat_toPlot = cor_t_rat
rownames(cor_t_rat_toPlot) = paste0(rownames(cor_t_rat_toPlot),"_","R.nor")
colnames(cor_t_rat_toPlot) = paste0(colnames(cor_t_rat_toPlot),"_","H.sap")

cor_t_all_toPlot = cor_t
rownames(cor_t_all_toPlot) = paste0(rownames(cor_t_all_toPlot),"_","R.nor")
colnames(cor_t_all_toPlot) = paste0(colnames(cor_t_all_toPlot),"_","H.sap")


    ## MNT added 14Apr2020: Reorder to diagonal & threshold at 0.4 for all-gene correlation === === ===
        # Start from descending - easier to manually order
        #cor_t_all_toPlot <- cor_t_all_toPlot[ ,rev(1:ncol(cor_t_all_toPlot))]
    # This is useful:
    apply(cor_t_all_toPlot, 2, which.max)
        # If want to re-order human labels (but prefer re-ordering rat labels)
        #cor_t_all_toPlot <- cor_t_all_toPlot[ ,rev(c(14,5,3,4, 7,10,12,6, 9,8,2,1, 11,13))]
    cor_t_all_toPlot <- cor_t_all_toPlot[c(14,11:13,4,3, 2,8,7,9, 6,15,5,16,1,10), ]
    # Threshold at 0.4
    range(cor_t_all_toPlot)
    cor_t_all_toPlot <- ifelse(cor_t_all_toPlot >= 0.4, 0.4, cor_t_all_toPlot)
    
    
    ## Do for other gene subsets ===
    # Human
        #cor_t_hsap_toPlot <- cor_t_hsap_toPlot[ ,rev(1:ncol(cor_t_hsap_toPlot))]
        #cor_t_hsap_toPlot <- cor_t_hsap_toPlot[ ,rev(c(14,5,3,4, 7,10,12,6, 9,8,2,1, 11,13))]
    cor_t_hsap_toPlot <- cor_t_hsap_toPlot[c(14,11:13,4,3, 2,8,7,9, 6,15,5,16,1,10), ]
    # Threshold at 0.4
    range(cor_t_hsap_toPlot)
    cor_t_hsap_toPlot <- ifelse(cor_t_hsap_toPlot >= 0.4, 0.4, cor_t_hsap_toPlot)

    # Rat
        #cor_t_rat_toPlot <- cor_t_rat_toPlot[ ,rev(1:ncol(cor_t_rat_toPlot))]
        #cor_t_rat_toPlot <- cor_t_rat_toPlot[ ,rev(c(14,5,3,4, 7,10,12,6, 9,8,2,1, 11,13))]
    cor_t_rat_toPlot <- cor_t_rat_toPlot[c(14,11:13,4,3, 2,8,7,9, 6,15,5,16,1,10), ]
    # Threshold at 0.4
    range(cor_t_rat_toPlot)
    cor_t_rat_toPlot <- ifelse(cor_t_rat_toPlot >= 0.4, 0.4, cor_t_rat_toPlot)
    
    

#pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/DayLab-ratNAc/overlap-DayLab-ratNAc_with_LIBD-10x-NAc-n5_top100-or-all_Apr2020.pdf")
pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/DayLab-ratNAc/overlap-DayLab-ratNAc_with_LIBD-10x-NAc-n5_top100-or-all_v2_Apr2020.pdf")
# Most human-specific
print(
  levelplot(
    cor_t_hsap_toPlot,
    aspect = "fill",
    at = theSeq.all,
    col.regions = my.col.all,
    ylab = "",
    xlab = "",
    scales = list(x = list(rot = 90, cex = 1), y = list(cex = 1)),
    main="Correlation of cluster-specific t's \n (top 100 genes/human (LIBD) clusters)"
  )
)
# Most rat-specific
print(
  levelplot(
    cor_t_rat_toPlot,
    aspect = "fill",
    at = theSeq.all,
    col.regions = my.col.all,
    ylab = "",
    xlab = "",
    scales = list(x = list(rot = 90, cex = 1), y = list(cex = 1)),
    main="Correlation of cluster-specific t's \n (top 100 genes/rat (Day) clusters)"
  )
)
# All
print(
  levelplot(
    cor_t_all_toPlot,
    aspect = "fill",
    at = theSeq.all,
    col.regions = my.col.all,
    ylab = "",
    xlab = "",
    scales = list(x = list(rot = 90, cex = 1), y = list(cex = 1)),
    main="Correlation of cluster-specific t's \n (all shared 14,121 homologs)",
    fontsize = 20
  )
)
dev.off()





### Another comparison: Rat nuclei vs human subclusters (t stats) ====================================
# 10x-pilot human NAc stats
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_NAc_all-n5_manualContrasts_MNTApr2020.rda", verbose=T)
    # eb_list.nac.all, sce.nac.all.PB, corfit.nac.all

# Day Lab Rat NAc stats
load("/dcl01/ajaffe/data/lab/singleCell/day_rat_snRNAseq/markers-stats_DayLab-ratNAc_manualContrasts_MNTApr2020.rda", verbose=T)
    # eb_list.nac.rat, sce.nac.rat.PB, corfit.nac.rat

# Day Lab Rat NAc full nuclei-level SCE
load("/dcl01/ajaffe/data/lab/singleCell/day_rat_snRNAseq/SCE_rat-NAc_downstream-processing_MNT.rda", verbose=T)
    # sce.nac.rat, chosen.hvgs.nac.rat

# Already subsetted on shared homologous genes (14,121):
load("/dcl01/ajaffe/data/lab/singleCell/day_rat_snRNAseq/SCE_rat-NAc-PBd_w_matchingHsap-NAc-PBd_HomoloGene.IDs_MNT.rda", verbose=T)
    # sce.rat.PBsub, sce.hsap.PBsub, Readme


sce.nac.rat
    # class: SingleCellExperiment
    # dim: 32883 15631
    # metadata(1): Samples
    # assays(2): counts logcounts
    # rownames(32883): ENSRNOG00000046319 ENSRNOG00000047964 ...
    #   ENSRNOG00000053444 ENSRNOG00000060590
    # rowData names(4): ID Symbol Type HomoloGene.ID
    # colnames(15631): AAACCCATCGGACTGC_1 AAACCCATCTCCTGCA_1 ...
    #   TTTGGTTTCTCATAGG_4 TTTGTTGGTTCTCGTC_4
    # colData names(14): Sample Barcode ... Celltype Idents.All_Groups_log.
    # reducedDimNames(3): PCA TSNE UMAP
    # spikeNames(0):
    # altExpNames(0):


sce.nac.rat <- sce.nac.rat[rownames(sce.rat.PBsub), ]
rownames(sce.nac.rat) <- rowData(sce.nac.rat)$HomoloGene.ID


## Human setup (rownames already in EnsemblID)
# Specificity stats
pvals_hsap <- sapply(eb_list.nac.all, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals_hsap) <- rowData(sce.nac.all.PB)$ID

ts_hsap <- sapply(eb_list.nac.all, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(ts_hsap) <- rowData(sce.nac.all.PB)$ID


# Subset and check matching 'HomoloGene.ID' ===
pvals_hsap <- pvals_hsap[rownames(sce.hsap.PBsub), ]
ts_hsap <- ts_hsap[rownames(sce.hsap.PBsub), ]

rownames(ts_hsap) <- rowData(sce.hsap.PBsub)$HomoloGene.ID
rownames(pvals_hsap) <- rowData(sce.hsap.PBsub)$HomoloGene.ID


# Now both rat SCE and human stats have been subsetted to the row names & order of the matching 'sce.___.PBsub'
table(rownames(sce.nac.rat) == rownames(ts_hsap))
    # all TRUE - good.


### Now look at rat expression vs human t's ========
ratExprs <- as.matrix(assay(sce.nac.rat, "logcounts"))
colnames(ratExprs) <- sce.nac.rat$Celltype

cor_ratExprs.hsapTs <- cor(ratExprs, ts_hsap)

pdf("pdfs/exploration/DayLab-ratNAc/overlap-DayLab-ratNAc-nucleiExprs_with_LIBD-HsapNAc-ts_Apr2020.pdf", height=12)
# Break up into Rat 'Celltype''s
for(i in levels(sce.nac.rat$Celltype)){
  cor_temp <- cor_ratExprs.hsapTs[rownames(cor_ratExprs.hsapTs)==i, ]
  colnames(cor_temp) <- paste0(colnames(cor_temp),": mean.r=",round(apply(cor_temp,2,mean),3))
  corRange <- range(cor_temp)
  pheatmap(cor_temp, fontsize_row=2,main=paste0("Correlation of rat nuclei labelled '", i,
                                                "' (n=", nrow(cor_temp),
                                                ") to LIBD NAc cluster t's \n (with average Pearson's r)"),
           cluster_cols=FALSE, show_rownames=FALSE, breaks=seq(corRange[1],corRange[2],by=((corRange[2]-corRange[1])/99)))
}
dev.off()




### What if we used this to predict classification of rat nuclei?
sce.nac.rat$cellType.hsapPredxn <- colnames(cor_ratExprs.hsapTs)[apply(cor_ratExprs.hsapTs,1,which.max)]

# Originally
table(sce.nac.rat$Celltype)
    # Astrocyte          Drd1-MSN        Drd2-MSN-1        Drd2-MSN-2
    #      1118              2748              1983               173
    #  Drd3-MSN    GABA undefined     Glutamatergic          Grm8-MSN
    #       351              1058               113              1059
    # Microglia             Mural            Olig-1            Olig-2
    #       723               157              3147              1564
    #    Olig-3    Polydendrocyte Pvalb-Interneuron   Sst-Interneuron
    #       129               835               297               176

# Human prediction
table(sce.nac.rat$cellType.hsapPredxn)
    # Astro  Inhib.1  Inhib.4    Micro MSN.D1.1 MSN.D1.2 MSN.D1.3 MSN.D2.1
    #   766     3536        4      680     1086     2308      239     2131
    # Oligo      OPC
    #  4754      127

# Cross-tabulated
options(width=120)  # default = 80 (I think)
table(sce.nac.rat$Celltype, sce.nac.rat$cellType.hsapPredxn)
    #                   Astro Inhib.1 Inhib.4 Micro MSN.D1.1 MSN.D1.2 MSN.D1.3 MSN.D2.1 Oligo  OPC
    # Astrocyte           759     326       3     0       22        6        0        1     1    0
    # Drd1-MSN              0     194       0     0      171     1637       44      702     0    0
    # Drd2-MSN-1            0     111       0     0       51      619        3     1199     0    0
    # Drd2-MSN-2            0      24       0     0        9        0        0      140     0    0
    # Drd3-MSN              0     227       0     0       91       24        4        5     0    0
    # GABA undefined        0    1048       0     0        9        0        0        1     0    0
    # Glutamatergic         0     111       0     0        1        1        0        0     0    0
    # Grm8-MSN              0      85       0     0      713        4      185       72     0    0
    # Microglia             0      16       0   677        7        5        0        0    18    0
    # Mural                 7     107       1     3       11       12        3       11     1    1
    # Olig-1                0       2       0     0        0        0        0        0  3144    1
    # Olig-2                0       1       0     0        0        0        0        0  1563    0
    # Olig-3                0      93       0     0        1        0        0        0    27    8
    # Polydendrocyte        0     718       0     0        0        0        0        0     0  117
    # Pvalb-Interneuron     0     297       0     0        0        0        0        0     0    0
    # Sst-Interneuron       0     176       0     0        0        0        0        0     0    0


table(sce.nac.rat$cellType.hsapPredxn, sce.nac.rat$Sex_Stim)
    #          Female_Cocaine Female_Saline Male_Cocaine Male_Saline
    # Astro               312            69          116         269
    # Inhib.1             926           910          583        1117
    # Inhib.4               1             2            0           1
    # Micro               211           169          150         150
    # MSN.D1.1            340           215          276         255
    # MSN.D1.2            718           651          512         427
    # MSN.D1.3             69            41           72          57
    # MSN.D2.1            869           400          536         326
    # Oligo              1407          1190          903        1254
    # OPC                  39            13           33          42




### What if do with top 100 genes defining each human NAc cluster? === === === === ===
hsap_specific_indices = mapply(function(t, p) {
  oo = order(t, decreasing = TRUE)[1:100]
},
  t=as.data.frame(ts_hsap),
  p=as.data.frame(pvals_hsap)
)
hsap_ind = unique(as.integer(hsap_specific_indices))  # of length 1331

table(hsap_ind %in% rowData(sce.hsap.PBsub)$HomoloGene.ID)
    # only 596... OH this is a POSITION index, not the actual IDs...

table(rownames(ts_hsap) == rownames(ratExprs))
    # 14121 TRUE

### Now look at rat expression vs human t's
    #ratExprs <- as.matrix(assay(sce.nac.rat, "logcounts"))       # (previously done)
    #colnames(ratExprs) <- sce.nac.rat$Celltype

# Make subsets of these
ratExprs.100sub <- ratExprs[hsap_ind, ]
ts_hsap.100sub <- ts_hsap[hsap_ind, ]

cor_ratExprs.hsapTs.top100 <- cor(ratExprs.100sub, ts_hsap.100sub)


pdf("pdfs/exploration/DayLab-ratNAc/overlap-DayLab-ratNAc-nucleiExprs_with_LIBD-HsapNAc-ts_top100perClust_Apr2020.pdf", height=12)
# Break up into Rat 'Celltype''s
for(i in levels(sce.nac.rat$Celltype)){
  cor_temp <- cor_ratExprs.hsapTs.top100[rownames(cor_ratExprs.hsapTs.top100)==i, ]
  colnames(cor_temp) <- paste0(colnames(cor_temp),": mean.r=",round(apply(cor_temp,2,mean),3))
  corRange <- range(cor_temp)
  pheatmap(cor_temp, fontsize_row=2,main=paste0("Correlation of rat nuclei labelled '", i,
                                                "' (n=", nrow(cor_temp),
                                                ") to LIBD NAc cluster t's \n (top 100 homologous genes/cluster)"),
           cluster_cols=FALSE, show_rownames=FALSE, breaks=seq(corRange[1],corRange[2],by=((corRange[2]-corRange[1])/99)))
}
dev.off()


# Proportion of zeros in genes
sce.nac.rat$propZeros.top100hsap <- apply(ratExprs.100sub, 2, function(n){mean(n==0)})
sce.nac.rat$propZeros <- apply(assay(sce.nac.rat, "logcounts"), 2, function(n){mean(n==0)})


pdf("pdfs/exploration/DayLab-ratNAc/dayLab-RatNAc-n4_propZeros_top100hsapClustGenes-or-all_MNTApr2020.pdf")
par(mar=c(7,4,4,2))
boxplot(sce.nac.rat$propZeros.top100hsap ~ sce.nac.rat$Celltype, las=2, xlab=NULL, cex.axis=0.8,
        main="Proportion of zeros across top 1,331 genes defining human NAc clusters", cex.main=0.9)
boxplot(sce.nac.rat$propZeros ~ sce.nac.rat$Celltype, las=2, xlab=NULL, cex.axis=0.8,
        main="Proportion of zeros across all (14,121) homologous genes to H. sapiens", cex.main=0.9)
dev.off()



## Assignment by highest proportion of non-zero genes across top 100 human-specific genes === ===
asgnmtByPropNon0 <- sapply(colnames(hsap_specific_indices), function(c){
  apply(ratExprs, 2, function(n){
    mean(n[hsap_specific_indices[ ,c]] != 0)
  })
})

sce.nac.rat$asgnmtByPropNon0 <- colnames(asgnmtByPropNon0)[apply(asgnmtByPropNon0,1,which.max)]

options(width=100)
table(sce.nac.rat$Celltype, sce.nac.rat$asgnmtByPropNon0)
    #                   Astro Inhib.1 Inhib.2 Micro MSN.D1.1 MSN.D1.2 MSN.D1.3 MSN.D2.1 Oligo  OPC
    # Astrocyte           230      11     873     0        0        0        0        0     4    0
    # Drd1-MSN              0       6    2740     0        0        0        0        2     0    0
    # Drd2-MSN-1            0       2    1980     0        0        1        0        0     0    0
    # Drd2-MSN-2            0       0     173     0        0        0        0        0     0    0
    # Drd3-MSN              0       2     348     0        0        1        0        0     0    0
    # GABA undefined        0       2    1056     0        0        0        0        0     0    0
    # Glutamatergic         0       0     113     0        0        0        0        0     0    0
    # Grm8-MSN              0       0    1059     0        0        0        0        0     0    0
    # Microglia             0       0      79   639        0        0        0        0     5    0
    # Mural                 1      19     122     2        1        5        1        4     2    0
    # Olig-1                0       1     243     0        0        0        0        0  2903    0
    # Olig-2                0       0      74     0        0        0        0        0  1490    0
    # Olig-3                0       1     118     0        0        0        0        0    10    0
    # Polydendrocyte        0      39     780     0        0        0        0        0     0   16
    # Pvalb-Interneuron     0       2     295     0        0        0        0        0     0    0
    # Sst-Interneuron       0     103      73     0        0        0        0        0     0    0



## Write out .mtx & colData for Day Lab ========================================
library(Matrix)

load("rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda", verbose=T)
     # sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all, ref.sampleInfo

# Mtx
mat2write <- assay(sce.nac.all, "counts")
Matrix::writeMM(mat2write, file="pdfs/exploration/DayLab-ratNAc/10xCounts/libd_n3-hom_n2-NeuN_countMat.mtx")
    ## NULL     - uh what?  Lol

# Features
write.table(rowData(sce.nac.all), file="pdfs/exploration/DayLab-ratNAc/10xCounts/libd_n3-hom_n2-NeuN_features.tsv",
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Barcodes
write.table(sce.nac.all$Barcode, file="pdfs/exploration/DayLab-ratNAc/10xCounts/libd_n3-hom_n2-NeuN_barcodes.tsv",
            row.names=FALSE, col.names=FALSE, quote=FALSE)


## Now try reading back in and comparing === === ===
path.test <- file.path("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/DayLab-ratNAc/10xCounts")
sce.test <- read10xCounts(path.test, col.names=TRUE)
    ## issues such as "/genes.tsv': No such file or directory" (& mtx file & barcodes.tsv)...
     #      -> just re-named those in the test dir and re-did
     #      (Seurat will operate differently... going to keep those names (i.e. re-write) bc the Day
     #       Lab shared files with more specific prefixes..)

sce.test
    # class: SingleCellExperiment
    # dim: 33538 13241
    # metadata(1): Samples
    # assays(1): counts
    # rownames(33538): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
    #   ENSG00000268674
    # rowData names(3): ID Symbol NA
    # colnames(13241): AAACCCACATCGAACT-1 AAACCCATCCAACCAA-1 ...
    #   TTTGTTGGTACTGCGC-1 TTTGTTGGTTCCAGGC-1
    # colData names(2): Sample Barcode
    # reducedDimNames(0):
    # spikeNames(0):
    # altExpNames(0):

head(rownames(assay(sce.test, "counts"))) # ensemblID

table(rownames(sce.test) == rowData(sce.nac.all)$ID)  # all TRUE - so make rownames(sce.test)

rownames(sce.test) <- rownames(sce.nac.all)

all.equal(assay(sce.test, "counts"), assay(sce.nac.all, "counts"))
    ## [1] TRUE - dope


## Write out colData
colnames(colData(sce.nac.all))
    # [1] "Sample"                "Barcode"               "sum"
    # [4] "detected"              "percent_top_50"        "percent_top_100"
    # [7] "percent_top_200"       "percent_top_500"       "subsets_Mito_sum"
    # [10] "subsets_Mito_detected" "subsets_Mito_percent"  "high.mito"
    # [13] "sample"                "region"                "donor"
    # [16] "processDate"           "protocol"              "prelimCluster"
    # [19] "collapsedCluster"      "lessCollapsed"         "cellType"
    # [22] "cellType.split"        "prelimCluster.split"   "cellType.moreSplit"
    # [25] "cellType.final"

# Keep $Barcode so can rid of the rownames
    # (because `table(rownames(colData(sce.nac.all)) == sce.nac.all$Barcode)` == all TRUE)
pheno2write <- colData(sce.nac.all)[ ,c(2, 13, 3,4, 11,12, 14:18, 25)]

write.csv(pheno2write, row.names=FALSE, file="pdfs/exploration/DayLab-ratNAc/10xCounts/libd_n3-hom_n2-NeuN_metadata.csv")






### Comparison to Day Lab Rat with SN-LEVEL stats =============================================
  # Added MNT 11May2020

# Load rat SCE
load("/dcl01/ajaffe/data/lab/singleCell/day_rat_snRNAseq/SCE_rat-NAc_downstream-processing_MNT.rda", verbose=T)
    # sce.nac.rat, chosen.hvgs.nac.rat
  
load("/dcl01/ajaffe/data/lab/singleCell/day_rat_snRNAseq/markers-stats_DayLab-ratNAc_findMarkers-SN-LEVEL_MNTMay2020.rda",
     verbose=T)
    # markers.rat.t.1vAll

## Calculate and add t-statistic (= std.logFC * sqrt(N)) for rat clusters
 #      and fix row order to the first entry "Astrocyte"
fixTo <- rownames(markers.rat.t.1vAll[["Astrocyte"]])
for(x in names(markers.rat.t.1vAll)){
  markers.rat.t.1vAll[[x]]$t.stat <- markers.rat.t.1vAll[[x]]$std.logFC * sqrt(ncol(sce.nac.rat))
  markers.rat.t.1vAll[[x]] <- markers.rat.t.1vAll[[x]][fixTo, ]
}

# Pull out the t's
ts.rat <- sapply(markers.rat.t.1vAll, function(x){x$t.stat})
rownames(ts.rat) <- fixTo



## Human t stats subset/re-ordering ===
# Bring in human stats; create t's
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_NAc-n5_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.nac.t.design, markers.nac.t.1vAll
    rm(markers.nac.t.design)

# Need to add t's with N nuclei used in constrasts

load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda",
     verbose=T)
    # sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all,ref.sampleInfo
    rm(chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all,ref.sampleInfo)
    # First drop "ambig.lowNtrxts" (93 nuclei)
    sce.nac.all <- sce.nac.all[ ,sce.nac.all$cellType.final != "ambig.lowNtrxts"]
    sce.nac.all$cellType.final <- droplevels(sce.nac.all$cellType.final)

## As above, calculate and add t-statistic (= std.logFC * sqrt(N)) for rat clusters
 #      and fix row order to the first entry "Astrocyte"
fixTo <- rownames(markers.nac.t.1vAll[["Astro"]])

for(s in names(markers.nac.t.1vAll)){
  markers.nac.t.1vAll[[s]]$t.stat <- markers.nac.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.nac.all))
  markers.nac.t.1vAll[[s]] <- markers.nac.t.1vAll[[s]][fixTo, ]
}

# Pull out the t's
ts.nac <- sapply(markers.nac.t.1vAll, function(x){x$t.stat})
rownames(ts.nac) <- fixTo



## Bring in HomoloGene.ID info to subset/match order
load("/dcl01/ajaffe/data/lab/singleCell/day_rat_snRNAseq/SCE_rat-NAc-PBd_w_matchingHsap-NAc-PBd_HomoloGene.IDs_MNT.rda",
     verbose=T)
    # sce.rat.PBsub, sce.hsap.PBsub, Readme

table(rowData(sce.rat.PBsub)$HomoloGene.ID == rowData(sce.hsap.PBsub)$HomoloGene.ID)  # all TRUE - dope
    # (see above - these are the intersecting homologs)

# First give [human] ts.nac rownames their respective EnsemblID
#   (have to use the full sce bc rownames(sce.hsap.PBsub) is EnsemblID and we uniquified the $Symbol)
rownames(ts.nac) <- rowData(sce.nac.all)$ID[match(rownames(ts.nac), rownames(sce.nac.all))]

# How many human genes with rat homologs are in these?
table(rownames(sce.hsap.PBsub) %in% rownames(ts.nac)) # 14121 good

# Subset/re-order for these and set to HomoloGene.ID
ts.nac <- ts.nac[rownames(sce.hsap.PBsub), ]
rownames(ts.nac) <- rowData(sce.hsap.PBsub)$HomoloGene.ID


# Same for rat t's
table(rownames(sce.rat.PBsub) %in% rownames(ts.rat))

ts.rat <- ts.rat[rownames(sce.rat.PBsub), ]
rownames(ts.rat) <- rowData(sce.rat.PBsub)$HomoloGene.ID

table(rownames(ts.nac) == rownames(ts.rat)) # all 14121 TRUE (well duh)


    # Save the ts matrices to reduce work next time
    Readme <- "These t-statistic matrices are subsetted and matched for shared 'HomoloGene.ID', so `cor()` can simply be run or other gene subsets applied first."
    save(ts.nac, ts.rat, Readme, file="rdas/zTsMats_libd-NAc_and_DayLab-ratNAc_sharedGenes_29May2020.rda")



cor_t_nac <- cor(ts.nac, ts.rat)
rownames(cor_t_nac) = paste0(rownames(cor_t_nac),"_H")
colnames(cor_t_nac) = paste0(colnames(cor_t_nac),"_R")


### Heatmap - typically use levelplot (e.g. below), but will want pheatmap bc can cluster cols/rows
theSeq.all = seq(-.65, .65, by = 0.01)
my.col.all <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq.all)-1)


# or thresholded at .55
theSeq.th = seq(-.55, .55, by = 0.01)
my.col.th <- colorRampPalette(brewer.pal(7, "RdYlBu"))(length(theSeq.th)-1)

cor_t_nac.th <- cor_t_nac
cor_t_nac.th <- ifelse(cor_t_nac.th >= 0.55, 0.55, cor_t_nac.th)


#pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/DayLab-ratNAc/overlap-DayLab-ratNAc_with_LIBD-10x-NAc-n5_SN-LEVEL-stats_allGenes_May2020.pdf")
pheatmap(cor_t_nac,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=8.5, fontsize_col=8.5,
         main="Correlation of cluster-specific t's \n (all shared expressed genes)")
# or thresholded at .5
pheatmap(cor_t_nac.th,
         color=my.col.th,
         breaks=theSeq.th,
         fontsize_row=8.5, fontsize_col=8.5,
         main="Correlation of cluster-specific t's \n (all shared expressed genes, thresholded)")
#dev.off()


## Or with manual ordering



# Manually re-order rat labels
apply(cor_t_nac, 1, which.max)
cor_t_nac <- cor_t_nac[ ,c(1,16, 7,15,6, 9,10, 6,8, 2,5,4,3,11,12,14)]
cor_t_nac.th <- cor_t_nac.th[ ,c(1,16, 7,15,6, 9,10, 6,8, 2,5,4,3,11,12,14)]

pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/DayLab-ratNAc/overlap-DayLab-ratNAc_with_LIBD-10x-NAc-n5_SN-LEVEL-stats_allGenes_v2_May2020.pdf")
## no cutoff
# "BrBG"
my.col.all <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.all)-1)
pheatmap(cor_t_nac,
         cluster_cols=F, cluster_rows=F,
         angle_col=90,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=11.5, fontsize_col=11.5,
         main="Correlation of cluster-specific t's \n (all shared expressed genes)")
# Original "PRGn"
my.col.all <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq.all)-1)
pheatmap(cor_t_nac,
         cluster_cols=F, cluster_rows=F,
         angle_col=90,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=11.5, fontsize_col=11.5,
         main="Correlation of cluster-specific t's \n (all shared expressed genes)")

## or thresholded at .55
# "BrBG"
my.col.th <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.th)-1)
pheatmap(cor_t_nac.th,
         cluster_cols=F, cluster_rows=F,
         angle_col=90,
         color=my.col.th,
         breaks=theSeq.th,
         fontsize_row=11.5, fontsize_col=11.5,
         legend_breaks=round(c(seq(-0.55,0.55,by=0.275)),2),
         main="Correlation of cluster-specific t's \n (all shared expressed genes, thresholded)")
# Original "PRGn"
my.col.th <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq.th)-1)
pheatmap(cor_t_nac.th,
         cluster_cols=F, cluster_rows=F,
         angle_col=90,
         color=my.col.th,
         breaks=theSeq.th,
         fontsize_row=11.5, fontsize_col=11.5,
         legend_breaks=round(c(seq(-0.55,0.55,by=0.275)),2),
         main="Correlation of cluster-specific t's \n (all shared expressed genes, thresholded)")

dev.off()



