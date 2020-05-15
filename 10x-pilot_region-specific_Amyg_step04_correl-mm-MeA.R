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







### Comparison to UCLA mouse MeA with SN-LEVEL stats ==================================
# Added MNT 14May2020

# Load mouse stats
load("/dcl01/ajaffe/data/lab/singleCell/ucla_mouse-MeA/2017Neuron/markers-stats_mouseMeA-2017-neuSubs_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.mmMeAneu.t.1vAll

# Load mouse SCE
load("/dcl01/ajaffe/data/lab/singleCell/ucla_mouse-MeA/2017Neuron/SCE_mouse-MeA-2017_neuronalSubclusters_HVGs_MNT.rda", verbose=T)
    # sce.amy.mm17hvgs

## Calculate and add t-statistic (= std.logFC * sqrt(N)) for mouse clusters
#      and fix row order to the first entry "Astrocyte"
fixTo <- rownames(markers.mmMeAneu.t.1vAll[[1]])
for(x in names(markers.mmMeAneu.t.1vAll)){
  markers.mmMeAneu.t.1vAll[[x]]$t.stat <- markers.mmMeAneu.t.1vAll[[x]]$std.logFC * sqrt(ncol(sce.amy.mm17hvgs))
  markers.mmMeAneu.t.1vAll[[x]] <- markers.mmMeAneu.t.1vAll[[x]][fixTo, ]
}

# Pull out the t's
ts.mmMeA <- sapply(markers.mmMeAneu.t.1vAll, function(x){x$t.stat})
rownames(ts.mmMeA) <- fixTo



## Human t stats subset/re-ordering ===
# Bring in human stats; create t's
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_Amyg-n2_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.amy.t.1vAll, markers.amy.t.design, markers.amy.wilcox.block
    rm(markers.amy.t.design, markers.amy.wilcox.block)

# Need to add t's with N nuclei used in constrasts
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    #sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo
    rm(chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy,ref.sampleInfo)

# First drop "ambig.lowNtrxts" (93 nuclei)
sce.amy <- sce.amy[ ,sce.amy$cellType.split != "Ambig.lowNtrxts"]
sce.amy$cellType.split <- droplevels(sce.amy$cellType.split)

## As above, calculate and add t-statistic (= std.logFC * sqrt(N)) from contrasts
#      and fix row order to the first entry "Astro"
fixTo <- rownames(markers.amy.t.1vAll[["Astro"]])

for(s in names(markers.amy.t.1vAll)){
  markers.amy.t.1vAll[[s]]$t.stat <- markers.amy.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.amy))
  markers.amy.t.1vAll[[s]] <- markers.amy.t.1vAll[[s]][fixTo, ]
}

# Pull out the t's
ts.amy <- sapply(markers.amy.t.1vAll, function(x){x$t.stat})
rownames(ts.amy) <- fixTo



## Bring in HomoloGene.ID info to subset/match order
load("/dcl01/ajaffe/data/lab/singleCell/ucla_mouse-MeA/2019Cell/SCE_mm-MeA-PBd_w_matchingHsap-Amyg-PBd_HomoloGene.IDs_MNT.rda",
     verbose=T)
    # sce.mm.PBsub, sce.hsap.PBsub, Readme

table(rowData(sce.mm.PBsub)$HomoloGene.ID == rowData(sce.hsap.PBsub)$HomoloGene.ID)  # all TRUE - dope
# (see above - these are the intersecting homologs)

## However!
table(rownames(ts.mmMeA) %in% rownames(sce.mm.PBsub)) # not all - so will need to get union
rm(sce.mm.PBsub, sce.hsap.PBsub, Readme)

## HomoloGene.ID for all human genes ====
    hom = read.delim("http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt",
                     as.is=TRUE)
    
    hom_hs <- hom[hom$Common.Organism.Name == "human", ]
    # of 19,124 entries
    
    # First Add EntrezID for human
    hs.entrezIds <- mapIds(org.Hs.eg.db, keys=rowData(sce.amy)$ID, 
                           column="ENTREZID", keytype="ENSEMBL")
    # "'select()' returned 1:many mapping between keys and columns"
    table(!is.na(hs.entrezIds))
        # 22,818 valid entries (remember this is already subsetted for those non-zero genes only)
    
    # Add to rowData
    rowData(sce.amy) <- cbind(rowData(sce.amy), hs.entrezIds)
    
    # Now how many in JAX db?
    table(rowData(sce.amy)$hs.entrezIds %in% hom_hs$EntrezGene.ID)
        # 18,865
    table(rowData(sce.amy)$Symbol %in% hom_hs$Symbol)
        # 18,479 - not a bad difference
    
          # So for mapping === === ===
          # human.entrez > HomoloGene.ID < mm.Symbol
          #                ^ filter SCE's on this
    
    # Human (by Entrez)
    rowData(sce.amy)$HomoloGene.ID <- hom_hs$HomoloGene.ID[match(rowData(sce.amy)$hs.entrezIds,
                                                                    hom_hs$EntrezGene.ID)]
    # end chunk ====


# Intersection?
table(rowData(sce.amy.mm17hvgs)$HomoloGene.ID %in% rowData(sce.amy)$HomoloGene.ID)
    # FALSE  TRUE
    #   135  2778


# First give [human] ts.amy rownames their respective EnsemblID
    #   (have to use the full sce bc rownames(sce.hsap.PBsub) is EnsemblID and we uniquified the $Symbol)
rownames(ts.amy) <- rowData(sce.amy)$ID[match(rownames(ts.amy), rownames(sce.amy))]
# Then to HomoloGene.ID
rownames(ts.amy) <- rowData(sce.amy)$HomoloGene.ID[match(rownames(ts.amy), rowData(sce.amy)$ID)]
    # Btw some are NA

# Subset for those with HomoloGene.ID
ts.amy <- ts.amy[!is.na(rownames(ts.amy)), ]



# Mouse - can just go to HomoloGene.ID
rownames(ts.mmMeA) <- rowData(sce.amy.mm17hvgs)$HomoloGene.ID[match(rownames(ts.mmMeA), rownames(sce.amy.mm17hvgs))]

# Intersecting?
table(rownames(ts.mmMeA) %in% rownames(ts.amy)) # all 14121 TRUE (well duh)
    # FALSE  TRUE
    #   193  2720

# Subset and match order
ts.mmMeA <- ts.mmMeA[rownames(ts.mmMeA) %in% rownames(ts.amy), ]
ts.amy <- ts.amy[rownames(ts.mmMeA), ]

cor_t_amyNeu <- cor(ts.amy, ts.mmMeA)
rownames(cor_t_amyNeu) = paste0(rownames(cor_t_amyNeu),"_","H")
colnames(cor_t_amyNeu) = paste0(colnames(cor_t_amyNeu),"_","M")
range(cor_t_amyNeu)
    #[1] -0.2557751  0.2577207      - kinda sucks lol...

### Heatmap - typically use levelplot (e.g. below), but will want pheatmap bc can cluster cols/rows
theSeq.all = seq(-.3, .3, by = 0.025)
my.col.all <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.all)-1)


# # or thresholded at .5    ( unneeded for now )
# theSeq.th = seq(-.5, .5, by = 0.025)
# my.col.th <- colorRampPalette(brewer.pal(7, "RdYlBu"))(length(theSeq.th)-1)

# cor_t_amyNeu.th <- cor_t_amyNeu
# cor_t_amyNeu.th <- ifelse(cor_t_amyNeu.th >= 0.5, 0.5, cor_t_amyNeu.th)


pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/HongLab-UCLA_mmMeA/overlap-mouseMeA-neuSubs_with_LIBD-10x-AMY_SN-LEVEL-stats_May2020.pdf")
pheatmap(cor_t_amyNeu,
         color=my.col.all,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.all,
         fontsize=9, fontsize_row=12, fontsize_col=12,
         main="Correlation of cluster-specific t's for mouse MeA neuronal subclusters \n (Wu et al., Neuron 2017)")
dev.off()


## Or with manual ordering (nvm not worth it with this comparison)
 #    - maybe for combined MeA stats...

# Manually re-order rat labels
apply(cor_t_amyNeu, 1, which.max)
cor_t_amyNeu <- cor_t_amyNeu[ ,c(1,16, 7,15,6, 9,10, 6,8, 2,5,4,3,11,12,14)]
cor_t_amyNeu.th <- cor_t_amyNeu.th[ ,c(1,16, 7,15,6, 9,10, 6,8, 2,5,4,3,11,12,14)]

#pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/HongLab-UCLA_mmMeA/overlap-mouseMeA-neuSubs_with_LIBD-10x-AMY_SN-LEVEL-stats_May2020.pdf")
## no cutoff
# "BrBG"
my.col.all <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.all)-1)
pheatmap(cor_t_amyNeu,
         cluster_cols=F, cluster_rows=F,
         angle_col=90,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=11.5, fontsize_col=11.5,
         main="Correlation of cluster-specific t's \n (all shared expressed genes)")

## or thresholded at .5
# "BrBG"
my.col.th <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.th)-1)
pheatmap(cor_t_amyNeu.th,
         cluster_cols=F, cluster_rows=F,
         angle_col=90,
         color=my.col.th,
         breaks=theSeq.th,
         fontsize_row=11.5, fontsize_col=11.5,
         legend_breaks=c(seq(-0.5,0.5,by=0.25)),
         main="Correlation of cluster-specific t's \n (all shared expressed genes, thresholded)")

#dev.off()




### Looking into some neuronal marker stats b/tw spp ===
markers.mmMeAneu.rowsFixed <- markers.mmMeAneu.t.1vAll

load("/dcl01/ajaffe/data/lab/singleCell/ucla_mouse-MeA/2017Neuron/markers-stats_mouseMeA-2017-neuSubs_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.mmMeAneu.t.1vAll
sapply(markers.mmMeAneu.t.1vAll, function(x){which(rownames(x)=="Npy")})
    # N4: 4; N5: 19       (oh but this is poorly expressed in mouse still...)


markers.amy.rowsFixed <- markers.amy.t.1vAll

load("rdas/markers-stats_Amyg-n2_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.amy.t.1vAll, markers.amy.t.design, markers.amy.wilcox.block
    rm(markers.amy.t.design, markers.amy.wilcox.block)

sapply(markers.amy.t.1vAll, function(x){which(rownames(x)=="NPY")}) #nothing
sapply(markers.amy.t.1vAll, function(x){which(rownames(x)=="TAC1")})  # Inhib.3 looks like

    

# Look at some notable genes from paper
plotExpression(sce.amy, exprs_values = "logcounts", features=c("NPY", "SST", "CCK", "TAC1", "CARTPT",
                                                               "SNAP25", "GAD1", "GAD2", "KCNQ2"),
               x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=5,
               add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                            geom = "crossbar", width = 0.3,
                                            colour=rep(tableau20[1:12], 9)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) 









