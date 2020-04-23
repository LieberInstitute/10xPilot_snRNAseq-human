### MNT 10x snRNA-seq workflow: step 04 - downstream comparisons
###   **Region-specific analyses**
###     - (3x) NAc samples from: Br5161 & Br5212 & Br5287
###   * Comparison to UCLA's Drop-seq on mouse medial amyg (MeA)
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


## load modeling outputs
# 10x-pilot human Amyg
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_Amyg-n2_manualContrasts_MNTMar2020.rda", verbose=T)
    # eb_contrasts.amy.broad, eb_list.amy.broad, sce.amy.PB

# UCLA mouse MeA Drop-seq
load("/dcl01/ajaffe/data/lab/singleCell/ucla_mouse-MeA/markers-stats_mouse-MeA-Drop-seq_manualContrasts_MNTApr2020.rda", verbose=T)
    # eb_list.amy.mm, corfit.amy.mm, sce.amy.mm.PB


# Add EntrezID for human
hs.entrezIds <- mapIds(org.Hs.eg.db, keys=rowData(sce.amy.PB)$ID, 
                       column="ENTREZID", keytype="ENSEMBL")
# "'select()' returned 1:many mapping between keys and columns"
table(!is.na(hs.entrezIds))
    # 20,578 valid entries (remember this is already subsetted for those non-zero genes only)

# Add to rowData
rowData(sce.amy.PB) <- cbind(rowData(sce.amy.PB), hs.entrezIds)


## Bring in 'HomoloGene.ID' for human (already in rowData for mm SCE) ===
## JAX annotation info
hom = read.delim("http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt",
                 as.is=TRUE)

hom_hs <- hom[hom$Common.Organism.Name == "human", ]
    # of 19,124 entries
table(rowData(sce.amy.PB)$hs.entrezIds %in% hom_hs$EntrezGene.ID)
    # 17,261
table(rowData(sce.amy.PB)$Symbol %in% hom_hs$Symbol)
    # 16,916 - not a bad difference

        # So for mapping === == === ===
        # human.entrez > HomoloGene.ID < mm.Symbol
        #                ^ filter SCE's on this

# Human (by Entrez)
rowData(sce.amy.PB)$HomoloGene.ID <- hom_hs$HomoloGene.ID[match(rowData(sce.amy.PB)$hs.entrezIds,
                                                                 hom_hs$EntrezGene.ID)]


## Now set/match to shared homologous genes ===

length(intersect(rowData(sce.amy.PB)$HomoloGene.ID,
                 rowData(sce.amy.mm.PB)$HomoloGene.ID))  # 13,444

sharedHomologs <- intersect(rowData(sce.amy.PB)$HomoloGene.ID,
                            rowData(sce.amy.mm.PB)$HomoloGene.ID)
# # That first one is NA - get rid of it
# sharedHomologs <- sharedHomologs[-1]

# Human not in mm
length(setdiff(rowData(sce.amy.PB)$HomoloGene.ID,
                 rowData(sce.amy.mm.PB)$HomoloGene.ID))  # 3657
# mm not in human
length(setdiff(rowData(sce.amy.mm.PB)$HomoloGene.ID,
               rowData(sce.amy.PB)$HomoloGene.ID))  # 928


# Subset for those
sce.mm.PBsub <- sce.amy.mm.PB[rowData(sce.amy.mm.PB)$HomoloGene.ID %in% sharedHomologs, ]   # 14247
sce.hsap.PBsub <- sce.amy.PB[rowData(sce.amy.PB)$HomoloGene.ID %in% sharedHomologs, ]  # 14178
    ## Many are duplicated...

rowData(sce.mm.PBsub)$Symbol[duplicated(rowData(sce.mm.PBsub)$HomoloGene.ID)]
    # shoot many genes are orthologs
rowData(sce.hsap.PBsub)$Symbol[duplicated(rowData(sce.hsap.PBsub)$HomoloGene.ID)]
    # same here, slightly less


### -> Take the higher-expressing of the duplicated - just mean across PB clusters:

    ## mm ===
    duplicatedSet.mm <- which(duplicated(rowData(sce.mm.PBsub)$HomoloGene.ID))
    genes2compare.mm <- list()
    gene2keep.mm <- character()
    for(g in 1:length(duplicatedSet.mm)){
      genes2compare.mm[[g]] <- rownames(sce.mm.PBsub)[rowData(sce.mm.PBsub)$HomoloGene.ID ==
                                              rowData(sce.mm.PBsub)$HomoloGene.ID[duplicatedSet.mm[g]]]
      rowmeansmat <- rowMeans(assay(sce.mm.PBsub[genes2compare.mm[[g]], ], "logcounts"))
      gene2keep.mm[g] <- names(rowmeansmat[order(rowmeansmat, decreasing=TRUE)])[1]
    }
    
    # Now pull out those that not being compared, so can `c()`
    table(rownames(sce.mm.PBsub) %in% unlist(genes2compare.mm)) # 133   - why isn't this ==
    sum(lengths(genes2compare.mm))                               # 328 ????
    length(unique(unlist(genes2compare.mm))) # 133   - oh. also `length(unique(gene2keep.mm)) == 52`
    
    genesNoCompare.mm <- rownames(sce.mm.PBsub)[!(rownames(sce.mm.PBsub) %in% unlist(genes2compare.mm))]
    
    # Finally combine and subset
    sce.mm.PBsub <- sce.mm.PBsub[c(genesNoCompare.mm, unique(gene2keep.mm)), ]
    
    table(rowData(sce.mm.PBsub)$HomoloGene.ID %in% sharedHomologs) # 13444 TRUE
    table(duplicated(rowData(sce.mm.PBsub)$HomoloGene.ID)) # 13444 FALSE         dope.

    
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
    table(rownames(sce.hsap.PBsub) %in% unlist(genes2compare.hsap)) # 109   - why isn't this ==
    sum(lengths(genes2compare.hsap))                               # 136 ????
    length(unique(unlist(genes2compare.hsap))) # 109   - oh. also `length(unique(gene2keep.hsap)) == 52`
    
    genesNoCompare.hsap <- rownames(sce.hsap.PBsub)[!(rownames(sce.hsap.PBsub) %in% unlist(genes2compare.hsap))]
        # of length 13392 (which + 52 == 13444)
    
    # Finally combine and subset
    sce.hsap.PBsub <- sce.hsap.PBsub[c(genesNoCompare.hsap, unique(gene2keep.hsap)), ]
    
    table(rowData(sce.hsap.PBsub)$HomoloGene.ID %in% sharedHomologs) # 13444 TRUE
    table(duplicated(rowData(sce.hsap.PBsub)$HomoloGene.ID)) # 13444 FALSE         dope.


    ## Match order and save
    sce.mm.PBsub <- sce.mm.PBsub[match(rowData(sce.hsap.PBsub)$HomoloGene.ID,
                                   rowData(sce.mm.PBsub)$HomoloGene.ID), ]

    table(rowData(sce.mm.PBsub)$HomoloGene.ID == rowData(sce.hsap.PBsub)$HomoloGene.ID)
        # all TRUE - good
    pheatmap(cor(assay(sce.mm.PBsub, "logcounts"), assay(sce.hsap.PBsub, "logcounts")), fontsize=5)
        # (ah but this is at the sample:cluster level)
    
    Readme <- "These two SCEs are subsetted and ordered for matching HomoloGene.ID in the rowData. This can be used to subset the nucleus-level SCEs in their respective Rdata files."
    save(sce.mm.PBsub, sce.hsap.PBsub, Readme, file="/dcl01/ajaffe/data/lab/singleCell/ucla_mouse-MeA/SCE_mm-MeA-PBd_w_matchingHsap-Amyg-PBd_HomoloGene.IDs_MNT.rda")
    
    
    
### FINALLY resume comparisons === === === === ===

## mm stats
pvals_mm <- sapply(eb_list.amy.mm, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals_mm) = rownames(sce.amy.mm.PB)

ts_mm <- sapply(eb_list.amy.mm, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(ts_mm) = rownames(sce.amy.mm.PB)



## Human stats
pvals_hsap <- sapply(eb_list.amy.broad, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals_hsap) = rowData(sce.amy.PB)$ID

ts_hsap <- sapply(eb_list.amy.broad, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(ts_hsap) = rowData(sce.amy.PB)$ID



### Subset and check matching 'HomoloGene.ID' === === === ===
pvals_mm <- pvals_mm[rownames(sce.mm.PBsub), ]
ts_mm <- ts_mm[rownames(sce.mm.PBsub), ]

pvals_hsap <- pvals_hsap[rowData(sce.hsap.PBsub)$ID, ]
ts_hsap <- ts_hsap[rowData(sce.hsap.PBsub)$ID, ]

rownames(ts_mm) <- rowData(sce.mm.PBsub)$HomoloGene.ID
rownames(pvals_mm) <- rowData(sce.mm.PBsub)$HomoloGene.ID

rownames(ts_hsap) <- rowData(sce.hsap.PBsub)$HomoloGene.ID
rownames(pvals_hsap) <- rowData(sce.hsap.PBsub)$HomoloGene.ID

table(rownames(ts_mm) == rownames(ts_hsap))
    ## all 14121 TRUE - good


## Now run correlation
cor_t = cor(ts_mm, ts_hsap)
signif(cor_t, 2)

## On just hsap cluster-specific homologous genes ===
hsap_specific_indices = mapply(function(t, p) {
    oo = order(t, decreasing = TRUE)[1:100]
  },
  as.data.frame(ts_hsap),
  as.data.frame(pvals_hsap)
)
hsap_ind = unique(as.numeric(hsap_specific_indices))

cor_t_hsap = cor(ts_mm[hsap_ind, ],
                  ts_hsap[hsap_ind, ])
signif(cor_t_hsap, 3)

## On just mouse cluster-specific homologous genes ===
mm_specific_indices = mapply(function(t, p) {
    oo = order(t, decreasing = TRUE)[1:100]
  },
  as.data.frame(ts_mm),
  as.data.frame(pvals_mm)
)
mm_ind = unique(as.numeric(mm_specific_indices))

cor_t_mm = cor(ts_mm[mm_ind, ],
                 ts_hsap[mm_ind, ])
signif(cor_t_mm, 3)



### Heatmap
theSeq.all = seq(-.5, .5, by = 0.01)
my.col.all <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq.all))

ct = colData(sce.hsap.PBsub)
ct = ct[!duplicated(sce.hsap.PBsub$cellType.final), ]

cor_t_hsap_toPlot = cor_t_hsap
rownames(cor_t_hsap_toPlot) = paste0(rownames(cor_t_hsap_toPlot),"_","M.mus")
colnames(cor_t_hsap_toPlot) = paste0(colnames(cor_t_hsap_toPlot),"_","H.sap")

cor_t_mm_toPlot = cor_t_mm
rownames(cor_t_mm_toPlot) = paste0(rownames(cor_t_mm_toPlot),"_","M.mus")
colnames(cor_t_mm_toPlot) = paste0(colnames(cor_t_mm_toPlot),"_","H.sap")

cor_t_all_toPlot = cor_t
rownames(cor_t_all_toPlot) = paste0(rownames(cor_t_all_toPlot),"_","M.mus")
colnames(cor_t_all_toPlot) = paste0(colnames(cor_t_all_toPlot),"_","H.sap")


    ## MNT added 14Apr2020: Reorder to diagonal & threshold at 0.4 for all-gene correlation === === ===
        # Start from descending - easier to manually order
        #cor_t_all_toPlot <- cor_t_all_toPlot[ ,rev(1:ncol(cor_t_all_toPlot))]
    # This is useful:
    apply(cor_t_all_toPlot, 2, which.max)
        # If want to re-order human labels (but prefer re-ordering mm labels)
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

    # mm
        #cor_t_mm_toPlot <- cor_t_mm_toPlot[ ,rev(1:ncol(cor_t_mm_toPlot))]
        #cor_t_mm_toPlot <- cor_t_mm_toPlot[ ,rev(c(14,5,3,4, 7,10,12,6, 9,8,2,1, 11,13))]
    cor_t_mm_toPlot <- cor_t_mm_toPlot[c(14,11:13,4,3, 2,8,7,9, 6,15,5,16,1,10), ]
    # Threshold at 0.4
    range(cor_t_mm_toPlot)
    cor_t_mm_toPlot <- ifelse(cor_t_mm_toPlot >= 0.4, 0.4, cor_t_mm_toPlot)
    
    

pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/overlap-mouse-MeA_with_LIBD-10x-Amyg_top100-or-all_Apr2020.pdf")
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
# Most mm-specific
print(
  levelplot(
    cor_t_mm_toPlot,
    aspect = "fill",
    at = theSeq.all,
    col.regions = my.col.all,
    ylab = "",
    xlab = "",
    scales = list(x = list(rot = 90, cex = 1), y = list(cex = 1)),
    main="Correlation of cluster-specific t's \n (top 100 genes/mouse MeA (UCLA) clusters)"
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
    main="Correlation of cluster-specific t's \n (all shared 13,444 homologs)",
    fontsize = 20
  )
)
dev.off()





### Another comparison: mm nuclei vs human subclusters (t stats) ====================================
# 10x-pilot human Amyg stats
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_NAc_all-n5_manualContrasts_MNTApr2020.rda", verbose=T)
    # eb_list.amy.broad, sce.amy.PB, corfit.amy.all

# Day Lab mm NAc stats
load("/dcl01/ajaffe/data/lab/singleCell/day_mm_snRNAseq/markers-stats_DayLab-mmNAc_manualContrasts_MNTApr2020.rda", verbose=T)
    # eb_list.amy.mm, sce.amy.mm.PB, corfit.amy.mm

# Day Lab mm NAc full nuclei-level SCE
load("/dcl01/ajaffe/data/lab/singleCell/day_mm_snRNAseq/SCE_mm-NAc_downstream-processing_MNT.rda", verbose=T)
    # sce.amy.mm, chosen.hvgs.amy.mm

# Already subsetted on shared homologous genes (14,121):
load("/dcl01/ajaffe/data/lab/singleCell/day_mm_snRNAseq/SCE_mm-NAc-PBd_w_matchingHsap-NAc-PBd_HomoloGene.IDs_MNT.rda", verbose=T)
    # sce.mm.PBsub, sce.hsap.PBsub, Readme


sce.amy.mm
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


sce.amy.mm <- sce.amy.mm[rownames(sce.mm.PBsub), ]
rownames(sce.amy.mm) <- rowData(sce.amy.mm)$HomoloGene.ID


## Human setup (rownames already in EnsemblID)
# Specificity stats
pvals_hsap <- sapply(eb_list.amy.broad, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals_hsap) <- rowData(sce.amy.PB)$ID

ts_hsap <- sapply(eb_list.amy.broad, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(ts_hsap) <- rowData(sce.amy.PB)$ID


# Subset and check matching 'HomoloGene.ID' ===
pvals_hsap <- pvals_hsap[rownames(sce.hsap.PBsub), ]
ts_hsap <- ts_hsap[rownames(sce.hsap.PBsub), ]

rownames(ts_hsap) <- rowData(sce.hsap.PBsub)$HomoloGene.ID
rownames(pvals_hsap) <- rowData(sce.hsap.PBsub)$HomoloGene.ID


# Now both mm SCE and human stats have been subsetted to the row names & order of the matching 'sce.___.PBsub'
table(rownames(sce.amy.mm) == rownames(ts_hsap))
    # all TRUE - good.


### Now look at mm expression vs human t's ========
mmExprs <- as.matrix(assay(sce.amy.mm, "logcounts"))
colnames(mmExprs) <- sce.amy.mm$Celltype

cor_mmExprs.hsapTs <- cor(mmExprs, ts_hsap)

pdf("pdfs/explommion/DayLab-mmNAc/overlap-DayLab-mmNAc-nucleiExprs_with_LIBD-HsapNAc-ts_Apr2020.pdf", height=12)
# Break up into mm 'Celltype''s
for(i in levels(sce.amy.mm$Celltype)){
  cor_temp <- cor_mmExprs.hsapTs[rownames(cor_mmExprs.hsapTs)==i, ]
  colnames(cor_temp) <- paste0(colnames(cor_temp),": mean.r=",round(apply(cor_temp,2,mean),3))
  corRange <- range(cor_temp)
  pheatmap(cor_temp, fontsize_row=2,main=paste0("Correlation of mm nuclei labelled '", i,
                                                "' (n=", nrow(cor_temp),
                                                ") to LIBD NAc cluster t's \n (with average Pearson's r)"),
           cluster_cols=FALSE, show_rownames=FALSE, breaks=seq(corRange[1],corRange[2],by=((corRange[2]-corRange[1])/99)))
}
dev.off()




### What if we used this to predict classification of mm nuclei?
sce.amy.mm$cellType.hsapPredxn <- colnames(cor_mmExprs.hsapTs)[apply(cor_mmExprs.hsapTs,1,which.max)]

# Originally
table(sce.amy.mm$Celltype)
    # Astrocyte          Drd1-MSN        Drd2-MSN-1        Drd2-MSN-2
    #      1118              2748              1983               173
    #  Drd3-MSN    GABA undefined     Glutamatergic          Grm8-MSN
    #       351              1058               113              1059
    # Microglia             Mural            Olig-1            Olig-2
    #       723               157              3147              1564
    #    Olig-3    Polydendrocyte Pvalb-Interneuron   Sst-Interneuron
    #       129               835               297               176

# Human prediction
table(sce.amy.mm$cellType.hsapPredxn)
    # Astro  Inhib.1  Inhib.4    Micro MSN.D1.1 MSN.D1.2 MSN.D1.3 MSN.D2.1
    #   766     3536        4      680     1086     2308      239     2131
    # Oligo      OPC
    #  4754      127

# Cross-tabulated
options(width=120)  # default = 80 (I think)
table(sce.amy.mm$Celltype, sce.amy.mm$cellType.hsapPredxn)
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


table(sce.amy.mm$cellType.hsapPredxn, sce.amy.mm$Sex_Stim)
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

table(rownames(ts_hsap) == rownames(mmExprs))
    # 14121 TRUE

### Now look at mm expression vs human t's
    #mmExprs <- as.matrix(assay(sce.amy.mm, "logcounts"))       # (previously done)
    #colnames(mmExprs) <- sce.amy.mm$Celltype

# Make subsets of these
mmExprs.100sub <- mmExprs[hsap_ind, ]
ts_hsap.100sub <- ts_hsap[hsap_ind, ]

cor_mmExprs.hsapTs.top100 <- cor(mmExprs.100sub, ts_hsap.100sub)


pdf("pdfs/explommion/DayLab-mmNAc/overlap-DayLab-mmNAc-nucleiExprs_with_LIBD-HsapNAc-ts_top100perClust_Apr2020.pdf", height=12)
# Break up into mm 'Celltype''s
for(i in levels(sce.amy.mm$Celltype)){
  cor_temp <- cor_mmExprs.hsapTs.top100[rownames(cor_mmExprs.hsapTs.top100)==i, ]
  colnames(cor_temp) <- paste0(colnames(cor_temp),": mean.r=",round(apply(cor_temp,2,mean),3))
  corRange <- range(cor_temp)
  pheatmap(cor_temp, fontsize_row=2,main=paste0("Correlation of mm nuclei labelled '", i,
                                                "' (n=", nrow(cor_temp),
                                                ") to LIBD NAc cluster t's \n (top 100 homologous genes/cluster)"),
           cluster_cols=FALSE, show_rownames=FALSE, breaks=seq(corRange[1],corRange[2],by=((corRange[2]-corRange[1])/99)))
}
dev.off()


# Proportion of zeros in genes
sce.amy.mm$propZeros.top100hsap <- apply(mmExprs.100sub, 2, function(n){mean(n==0)})
sce.amy.mm$propZeros <- apply(assay(sce.amy.mm, "logcounts"), 2, function(n){mean(n==0)})


pdf("pdfs/explommion/DayLab-mmNAc/dayLab-mmNAc-n4_propZeros_top100hsapClustGenes-or-all_MNTApr2020.pdf")
par(mar=c(7,4,4,2))
boxplot(sce.amy.mm$propZeros.top100hsap ~ sce.amy.mm$Celltype, las=2, xlab=NULL, cex.axis=0.8,
        main="Proportion of zeros across top 1,331 genes defining human NAc clusters", cex.main=0.9)
boxplot(sce.amy.mm$propZeros ~ sce.amy.mm$Celltype, las=2, xlab=NULL, cex.axis=0.8,
        main="Proportion of zeros across all (14,121) homologous genes to H. sapiens", cex.main=0.9)
dev.off()



## Assignment by highest proportion of non-zero genes across top 100 human-specific genes === ===
asgnmtByPropNon0 <- sapply(colnames(hsap_specific_indices), function(c){
  apply(mmExprs, 2, function(n){
    mean(n[hsap_specific_indices[ ,c]] != 0)
  })
})

sce.amy.mm$asgnmtByPropNon0 <- colnames(asgnmtByPropNon0)[apply(asgnmtByPropNon0,1,which.max)]

options(width=100)
table(sce.amy.mm$Celltype, sce.amy.mm$asgnmtByPropNon0)
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
     # sce.amy.all, chosen.hvgs.amy.all, pc.choice.amy.all, clusterRefTab.amy.all, ref.sampleInfo

# Mtx
mat2write <- assay(sce.amy.all, "counts")
Matrix::writeMM(mat2write, file="pdfs/explommion/DayLab-mmNAc/10xCounts/libd_n3-hom_n2-NeuN_countMat.mtx")
    ## NULL     - uh what?  Lol

# Features
write.table(rowData(sce.amy.all), file="pdfs/explommion/DayLab-mmNAc/10xCounts/libd_n3-hom_n2-NeuN_features.tsv",
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Barcodes
write.table(sce.amy.all$Barcode, file="pdfs/explommion/DayLab-mmNAc/10xCounts/libd_n3-hom_n2-NeuN_barcodes.tsv",
            row.names=FALSE, col.names=FALSE, quote=FALSE)


## Now try reading back in and comparing === === ===
path.test <- file.path("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/explommion/DayLab-mmNAc/10xCounts")
sce.test <- read10xCounts(path.test, col.names=TRUE)
    ## issues such as "/genes.tsv': No such file or directory" (& mtx file & barcodes.tsv)...
     #      -> just re-named those in the test dir and re-did
     #      (Seumm will opemme differently... going to keep those names (i.e. re-write) bc the Day
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

table(rownames(sce.test) == rowData(sce.amy.all)$ID)  # all TRUE - so make rownames(sce.test)

rownames(sce.test) <- rownames(sce.amy.all)

all.equal(assay(sce.test, "counts"), assay(sce.amy.all, "counts"))
    ## [1] TRUE - dope


## Write out colData
colnames(colData(sce.amy.all))
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
    # (because `table(rownames(colData(sce.amy.all)) == sce.amy.all$Barcode)` == all TRUE)
pheno2write <- colData(sce.amy.all)[ ,c(2, 13, 3,4, 11,12, 14:18, 25)]

write.csv(pheno2write, row.names=FALSE, file="pdfs/explommion/DayLab-mmNAc/10xCounts/libd_n3-hom_n2-NeuN_metadata.csv")



