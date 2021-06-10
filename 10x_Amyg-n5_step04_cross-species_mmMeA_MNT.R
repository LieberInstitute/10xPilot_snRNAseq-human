### MNT 10x snRNA-seq workflow: step 04 - downstream comparisons
###   **Region-specific analyses**
###     - 5x AMY samples (incl'g revision samples)
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

### *** 09Jun 2021 update:
# JAX MGI database no longer reports a $HomoloGene.ID -> now use $DB.Class.Key
# (corresponded with David Shaw @ JAX/MGI)

# For mapping === == === ===
# human.entrez > DB.Class.Key < mouse.entrez < mm.Symbol
#                ^ filter SCE's on this - to be more descriptive, call 'JAX.geneID'


### Setting up homologous gene IDs, for mapping b/tw species =============

## 10x-pilot human Amyg SCE
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_Amyg-n5_cleaned-combined_SCE_MNT2021.rda",
     verbose=T)
#rm(chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo, annotationTab.amy, cell_colors.amy)

# First drop "drop.lowNTx_" (1138 nuclei)
sce.amy <- sce.amy[ ,-grep("drop.",sce.amy$cellType)]
sce.amy$cellType <- droplevels(sce.amy$cellType)

# Remove 0 genes across all nuclei
sce.amy <- sce.amy[!rowSums(assay(sce.amy, "counts"))==0, ]  # keeps same 29371 genes


## Add EntrezID for human
hs.entrezIds <- mapIds(org.Hs.eg.db, keys=rowData(sce.amy)$gene_id, 
                       column="ENTREZID", keytype="ENSEMBL")
# "'select()' returned 1:many mapping between keys and columns"
table(!is.na(hs.entrezIds))
    # 21,041 valid entries  - that's a lot lost (8330)
    withoutEntrez <- names(hs.entrezIds)[is.na(hs.entrezIds)]
    # Store those somewhere, maybe for later reference
    table(rowData(sce.amy)[rowData(sce.amy)$gene_id %in% withoutEntrez, ]$gene_id == withoutEntrez)
    names(withoutEntrez) <- rowData(sce.amy)[rowData(sce.amy)$gene_id %in% withoutEntrez, ]$gene_name
    
# Add to rowData
rowData(sce.amy) <- cbind(rowData(sce.amy), hs.entrezIds)


## Bring in 'DB.Class.Key' for human ===
# JAX annotation info:
hom = read.delim("http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt",
                 as.is=TRUE)
hom_hs <- hom[hom$Common.Organism.Name == "human", ]
    # of 22,522 entries
table(hs.entrezIds %in% hom_hs$EntrezGene.ID)
    # 16,986
table(rowData(sce.amy)$gene_name %in% hom_hs$Symbol)
    # 16,666 - not a bad difference (interestingly these numbers are lower than when done for preprint...)

# Human (by Entrez)
rowData(sce.amy)$JAX.geneID <- hom_hs$DB.Class.Key[match(rowData(sce.amy)$hs.entrezIds,
                                                                 hom_hs$EntrezGene.ID)]



## UCLA mouse MeA Drop-seq ===
load("/dcl01/ajaffe/data/lab/singleCell/ucla_mouse-MeA/2019Cell/SCE_mouse-MeA_downstream-processing_MNT.rda", verbose=T)
    # sce.amy.mm, chosen.hvgs.amy.mm

# This SCE has already been set up with the 'EntrezGene.ID' & [deprecated] 'HomoloGene.ID'
hom_mm <- hom[hom$Common.Organism.Name == "mouse, laboratory", ]

table(rowData(sce.amy.mm)$EntrezGene.ID %in% hom_mm$EntrezGene.ID)
    #FALSE  TRUE 
    #  107 14403

# Add this new 'DB.Class.Key'
rowData(sce.amy.mm)$JAX.geneID <- hom_mm$DB.Class.Key[match(rowData(sce.amy.mm)$EntrezGene.ID,
                                                            hom_mm$EntrezGene.ID)]




## Now set/match to shared homologous genes ===

length(intersect(rowData(sce.amy)$JAX.geneID,
                 rowData(sce.amy.mm)$JAX.geneID))  # 13,669

sharedHomologs <- intersect(rowData(sce.amy)$JAX.geneID,
                            rowData(sce.amy.mm)$JAX.geneID)
    # That first one is NA - get rid of it
    sharedHomologs <- sharedHomologs[-1]

# Human not in mm
length(setdiff(rowData(sce.amy)$JAX.geneID,
                 rowData(sce.amy.mm)$JAX.geneID))  # 2806
# mm not in human
length(setdiff(rowData(sce.amy.mm)$JAX.geneID,
               rowData(sce.amy)$JAX.geneID))  # 735


# Subset for those
sce.mm.sub <- sce.amy.mm[rowData(sce.amy.mm)$JAX.geneID %in% sharedHomologs, ]   # 13668
sce.hsap.sub <- sce.amy[rowData(sce.amy)$JAX.geneID %in% sharedHomologs, ]  # 13983
    ## Many are duplicated...

rowData(sce.mm.sub)$Symbol[duplicated(rowData(sce.mm.sub)$JAX.geneID)]
    # no orthologs in this space
rowData(sce.hsap.sub)$gene_name[duplicated(rowData(sce.hsap.sub)$JAX.geneID)]
    # up to 315 orthologs...


### -> Take the higher-expressing of the orthologs - just mean across PB clusters:

    # ## mm ===
    # duplicatedSet.mm <- which(duplicated(rowData(sce.mm.sub)$HomoloGene.ID))
    # genes2compare.mm <- list()
    # gene2keep.mm <- character()
    # for(g in 1:length(duplicatedSet.mm)){
    #   genes2compare.mm[[g]] <- rownames(sce.mm.sub)[rowData(sce.mm.sub)$HomoloGene.ID ==
    #                                           rowData(sce.mm.sub)$HomoloGene.ID[duplicatedSet.mm[g]]]
    #   rowmeansmat <- rowMeans(assay(sce.mm.sub[genes2compare.mm[[g]], ], "logcounts"))
    #   gene2keep.mm[g] <- names(rowmeansmat[order(rowmeansmat, decreasing=TRUE)])[1]
    # }
    # 
    # # Now pull out those that not being compared, so can `c()`
    # table(rownames(sce.mm.sub) %in% unlist(genes2compare.mm)) # 133   - why isn't this ==
    # sum(lengths(genes2compare.mm))                               # 328 ????
    # length(unique(unlist(genes2compare.mm))) # 133   - oh. also `length(unique(gene2keep.mm)) == 52`
    # 
    # genesNoCompare.mm <- rownames(sce.mm.sub)[!(rownames(sce.mm.sub) %in% unlist(genes2compare.mm))]
    # 
    # # Finally combine and subset
    # sce.mm.sub <- sce.mm.sub[c(genesNoCompare.mm, unique(gene2keep.mm)), ]
    # 
    # table(rowData(sce.mm.sub)$HomoloGene.ID %in% sharedHomologs) # 13444 TRUE
    # table(duplicated(rowData(sce.mm.sub)$HomoloGene.ID)) # 13444 FALSE         dope.

    
    ## Human ===
    # First change rownames to EnsemblID
    rowData(sce.hsap.sub)$Symbol.unique <- rownames(sce.hsap.sub)
    rownames(sce.hsap.sub) <- rowData(sce.hsap.sub)$ID
        
    duplicatedSet.hsap <- which(duplicated(rowData(sce.hsap.sub)$JAX.geneID))
    genes2compare.hsap <- list()
    gene2keep.hsap <- character()
    for(g in 1:length(duplicatedSet.hsap)){
      genes2compare.hsap[[g]] <- rownames(sce.hsap.sub)[rowData(sce.hsap.sub)$JAX.geneID ==
                                                          rowData(sce.hsap.sub)$JAX.geneID[duplicatedSet.hsap[g]]]
      rowmeansmat <- rowMeans(assay(sce.hsap.sub[genes2compare.hsap[[g]], ], "logcounts"))
      gene2keep.hsap[g] <- names(rowmeansmat[order(rowmeansmat, decreasing=TRUE)])[1]
    }
    
    # This is interesting...
    length(gene2keep.hsap)  #315
    length(unique(gene2keep.hsap))  #230
    
    # Ahh that's because many 'tested' might have been orthologous,
    #   b/tw themselves (i.e. 3+ orthologous genes):
    length(unique(rowData(sce.hsap.sub)$JAX.geneID[duplicatedSet.hsap]))
        # 230

    
    # Now pull out those that not being compared, so can `c()`
    genesNoCompare.hsap <- rownames(sce.hsap.sub)[!(rownames(sce.hsap.sub) %in% unlist(genes2compare.hsap))]

    # Finally combine and subset
    sce.hsap.sub <- sce.hsap.sub[c(genesNoCompare.hsap, unique(gene2keep.hsap)), ]
    
    table(rowData(sce.hsap.sub)$JAX.geneID %in% sharedHomologs) # 13668 TRUE
    table(duplicated(rowData(sce.hsap.sub)$JAX.geneID)) # 13668 FALSE         dope.

    table(rowData(sce.hsap.sub)$JAX.geneID %in% rowData(sce.mm.sub)$JAX.geneID) # 13668 TRUE

    
    ## Match order and save
    sce.mm.sub <- sce.mm.sub[match(rowData(sce.hsap.sub)$JAX.geneID,
                                   rowData(sce.mm.sub)$JAX.geneID), ]

    table(rowData(sce.mm.sub)$JAX.geneID == rowData(sce.hsap.sub)$JAX.geneID)
        # all TRUE - good
    
    Readme <- "These two SCEs are subsetted and ordered for matching 'JAX.geneID' in the rowData. This can be used to subset the nucleus-level SCEs in their respective Rdata files."
    save(sce.mm.sub, sce.hsap.sub, Readme, file="rdas/revision/SCE_mm-MeA_matched2_Hsap-Amyg_JAX.geneIDs_MNT2021.rda")
    
    
    

### Comparison to UCLA mouse MeA with SN-LEVEL stats ==================================
  # Added MNT 14May2020 - UPDATED 22May2020 to compare to 2019 dataset
  # (previously only 2017 neuronal subclusters), now with neuronal subcluster info

# Load mouse stats - can still use these (revision 2021)
load("/dcl01/ajaffe/data/lab/singleCell/ucla_mouse-MeA/2019Cell/markers-stats_mouseMeA-2019-with-16neuSubs_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.mmMeA.t.1vAll
    table(rownames(sce.mm.sub) %in% rownames(markers.mmMeA.t.1vAll[[1]]))
        # all 13668 (of 14510 entries) TRUE
        
# Load SCEs
load("rdas/revision/SCE_mm-MeA_matched2_Hsap-Amyg_JAX.geneIDs_MNT2021.rda", verbose=T)
    # sce.mm.sub, sce.hsap.sub, Readme

## Calculate and add t-statistic (= std.logFC * sqrt(N)) for mouse clusters
#      and fix row order to the first entry "AS"
fixTo <- rownames(sce.mm.sub)
for(x in names(markers.mmMeA.t.1vAll)){
  markers.mmMeA.t.1vAll[[x]]$t.stat <- markers.mmMeA.t.1vAll[[x]]$std.logFC * sqrt(ncol(sce.mm.sub))
  markers.mmMeA.t.1vAll[[x]] <- markers.mmMeA.t.1vAll[[x]][fixTo, ]
}

# Pull out the t's
ts.mmMeA <- sapply(markers.mmMeA.t.1vAll, function(x){x$t.stat})
rownames(ts.mmMeA) <- rowData(sce.mm.sub)$JAX.geneID[match(fixTo, rownames(sce.mm.sub))]



## Human t stats subset/re-ordering ===
# Bring in human stats; create t's
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/markers-stats_Amyg-n5_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.amy.t.pw, markers.amy.wilcox.block, markers.amy.t.1vAll, medianNon0.amy
    rm(markers.amy.t.pw, markers.amy.wilcox.block)

## As above, calculate and add t-statistic (= std.logFC * sqrt(N)) from contrasts
#      and fix row order to the first entry "Astro"
fixTo <- rownames(sce.hsap.sub)

## These stats have both an '_enriched' & '_depleted' result - take the '_enriched'
markers.amy.enriched <- lapply(markers.amy.t.1vAll, function(x){x[[2]]})

# Compute t-stat
for(s in names(markers.amy.enriched)){
  markers.amy.enriched[[s]]$t.stat <- markers.amy.enriched[[s]]$std.logFC * sqrt(ncol(sce.hsap.sub))
  markers.amy.enriched[[s]] <- markers.amy.enriched[[s]][fixTo, ]
}

# Pull out the t's
ts.amy <- sapply(markers.amy.enriched, function(x){x$t.stat})
rownames(ts.amy) <- rowData(sce.hsap.sub)$JAX.geneID[match(fixTo, rownames(sce.hsap.sub))]

# Do they all intersect?
table(rownames(ts.amy) == rownames(ts.mmMeA))
     # 13668 TRUE

cor_t_amy <- cor(ts.amy, ts.mmMeA)
rownames(cor_t_amy) = paste0(rownames(cor_t_amy),"_","H")
colnames(cor_t_amy) = paste0(colnames(cor_t_amy),"_","M")
range(cor_t_amy)
    # [1] -0.2857217  0.5872299      (previously {-0.2203968, 0.5023080 } w preprint pops)

### Heatmap - typically use levelplot (e.g. below), but will want pheatmap bc can cluster cols/rows
theSeq.all = seq(-.6, .6, by = 0.01)
my.col.all <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.all)-1)

# Re-order mouse labels - move MG/MU to after neuronal subclusters
cor_t_amy <- cor_t_amy[ ,c(1,2, 5,13:20,6:12, 3,4, 21:23)]

pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/revision/overlap-mouseMeA-2019_x_LIBD-10x-AMY_allGenes_MNT2021.pdf")
pheatmap(cor_t_amy,
         color=my.col.all,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.all,
         fontsize=11, fontsize_row=13, fontsize_col=12,
         main="Correlation of cluster-specific t's to mouse MeA \n subclusters (Chen-Hu-Wu et al., Cell 2019)")
dev.off()




## Iteration with top N spp:subcluster-specific genes: ========
 #   -> Basically just run through line 488, under ("Subset and match order")

# Save the ts matrices to reduce work next time
homologInfo <- data.frame(JAX.geneID = rownames(ts.amy),
                          H.sap.ID = rowData(sce.hsap.sub)$Symbol.uniq[
                            match(rownames(ts.amy), rowData(sce.hsap.sub)$JAX.geneID)],
                          M.mus.ID = rowData(sce.mm.sub)$Symbol[
                            match(rownames(ts.amy), rowData(sce.mm.sub)$JAX.geneID)]
                          )

Readme <- "These t-statistic matrices are subsetted and matched for shared 'JAX.geneID', so `cor()` can simply be run or other gene subsets applied first."
save(ts.amy, ts.mmMeA, homologInfo, Readme, file="rdas/revision/zTsMats_libd-AMY_and_ucla-mouseMeA-2019Cell_sharedGenes_MNT2021.rda")


    
    
## On just hsap cluster-specific homologous genes ===
hsap_specific_indices = mapply(function(t) {
  oo = order(t, decreasing = TRUE)[1:100]
  },
as.data.frame(ts.amy)
)
hsap_ind = unique(as.numeric(hsap_specific_indices))
length(hsap_ind)  # so of 1900 (100 x 19 cellType), 1399 unique

cor_t_hsap = cor(ts.amy[hsap_ind, ],
                 ts.mmMeA[hsap_ind, ])
rownames(cor_t_hsap) = paste0(rownames(cor_t_hsap),"_","H")
colnames(cor_t_hsap) = paste0(colnames(cor_t_hsap),"_","M")
range(cor_t_hsap)
    # [1] -0.3881377  0.7220732
        # (previously [1] -0.2738376  0.6612352)


## On just mouse cluster-specific homologous genes ===
mouse_specific_indices = mapply(function(t) {
  oo = order(t, decreasing = TRUE)[1:100]
},
as.data.frame(ts.mmMeA)
)
mouse_ind = unique(as.numeric(mouse_specific_indices))
length(mouse_ind)  # so of 2300 (100 x 23 subCluster), 1543 unique

cor_t_mouse = cor(ts.amy[mouse_ind, ],
                 ts.mmMeA[mouse_ind, ])
rownames(cor_t_mouse) = paste0(rownames(cor_t_mouse),"_","H")
colnames(cor_t_mouse) = paste0(colnames(cor_t_mouse),"_","M")
range(cor_t_mouse)
    # [1] -0.4078769  0.6786080
        # (previously [1] -0.2731605  0.6113445)

## Between these gene spaces:
length(intersect(rownames(ts.amy)[hsap_ind], rownames(ts.mmMeA)[mouse_ind]))
    # 480

toptop.genes <- intersect(rownames(ts.amy)[hsap_ind], rownames(ts.mmMeA)[mouse_ind])
cor_t_top <- cor(ts.amy[toptop.genes, ],
                 ts.mmMeA[toptop.genes, ])
rownames(cor_t_top) = paste0(rownames(cor_t_top),"_","H")
colnames(cor_t_top) = paste0(colnames(cor_t_top),"_","M")
range(cor_t_top)
    # [1] -0.4681437  0.7896791

## UPDATED heatmap:
theSeq.all = seq(-.70, .70, by = 0.01)
my.col.all <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.all)-1)

# Re-order mouse labels - move MG/MU to after neuronal subclusters
cor_t_hsap <- cor_t_hsap[ ,c(1,2, 5,13:20,6:12, 3,4, 21:23)]

cor_t_mouse <- cor_t_mouse[ ,c(1,2, 5,13:20,6:12, 3,4, 21:23)]

cor_t_top <- cor_t_top[ ,c(1,2, 5,13:20,6:12, 3,4, 21:23)]

# Print all four iterations
pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/revision/overlap-mouseMeA-2019_x_LIBD-10x-AMY_allGenes_wValues_MNT2021.pdf")
pheatmap(cor_t_amy,
         color=my.col.all,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.all,
         fontsize=11, fontsize_row=13, fontsize_col=12,
            display_numbers=T, number_format="%.2f", fontsize_number=6,
         legend_breaks=c(seq(-0.7,0.7,by=0.35)),
         main="Correlation of cluster-specific t's to mouse MeA \n subclusters (Chen-Hu-Wu et al., Cell 2019)")
# On human-specific genes
pheatmap(cor_t_hsap,
         color=my.col.all,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.all,
         fontsize=11, fontsize_row=13, fontsize_col=12,
            display_numbers=T, number_format="%.2f", fontsize_number=6,
         legend_breaks=c(seq(-0.7,0.7,by=0.35)),
         main="Correlation of top-100 cluster-specific t's (1399) to \n (Chen-Hu-Wu et al., Cell 2019) subclusters")
# On mm-MeA-specific genes
pheatmap(cor_t_mouse,
         color=my.col.all,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.all,
         fontsize=11, fontsize_row=13, fontsize_col=12,
            display_numbers=T, number_format="%.2f", fontsize_number=6,
         legend_breaks=c(seq(-0.7,0.7,by=0.35)),
         main="Correlation of LIBD-AMY subclusters to \n (Chen-Hu-Wu et al., Cell 2019) subcluster top-100 t's (1543)")
# On intersection between the top spp.-specific genes (480 genes)
theSeq.new = seq(-.80, .80, by = 0.01)
my.col.new <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.new)-1)
pheatmap(cor_t_top,
         color=my.col.new,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.new,
         fontsize=9.5, fontsize_row=13, fontsize_col=12,
         display_numbers=T, number_format="%.2f", fontsize_number=6,
         legend_breaks=c(seq(-0.8,0.8,by=0.4)),
         main="Correlation of LIBD-AMY subclusters to \n (Chen-Hu-Wu et al., Cell 2019) subcluster t's (shared top 100's, 480)")
dev.off()



## Intersecting some of the top markers =====================
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_Amyg-n2_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.amy.t.1vAll, markers.amy.t.design, markers.amy.wilcox.block
    rm(markers.amy.t.design, markers.amy.wilcox.block)

# Take top 100
markerList.t.hsap <- lapply(markers.amy.t.1vAll, function(x){
  rownames(x)[x$log.FDR < log10(1e-6)]
  }
)
genes.top100.hsap <- lapply(markerList.t.hsap, function(x){head(x, n=100)})


load("/dcl01/ajaffe/data/lab/singleCell/ucla_mouse-MeA/2019Cell/markers-stats_mouseMeA-2019-with-16neuSubs_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.mmMeA.t.1vAll

# Just `toupper()` it
markerList.t.mm <- lapply(markers.mmMeA.t.1vAll, function(x){
  rownames(x)[x$log.FDR < log10(1e-6)]
  }
)
genes.top100.mm <- lapply(markerList.t.mm, function(x){toupper(head(x, n=100))})
genes.top100.mm <- sapply(genes.top100.mm, cbind)

## sapply
sapply(genes.top100.hsap, function(x){
  apply(genes.top100.mm,2,function(y){length(intersect(x,y))})
})
    #        Astro Excit.1 Excit.2 Excit.3 Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5 Micro Oligo OPC
    # AS        20       0       1       0       0       1       0       0       1     0     0   1
    # EN         2       1       1       0       1       2       1       0       0     0     0   2
    # MG         0       0       0       1       0       0       0       0       0    19     0   0
    # MU         2       1       0       1       0       0       1       1       0     0     0   0
    # N.1        1       4       2       0      14       8       3       4       9     0     1   1
    # N.10       0       6       1       4       7       7       2       0       6     0     0   0
    # N.11       1      10       5       2       8       3       4       6       8     0     0   4
    # N.12       2       7       4       3       7       5       2       3       5     0     2   2
    # N.13       1       2       1       3       1       1       0       0       5     1     0   1
    # N.14       0       7       2       4       9       6       0       4       7     1     1   2
    # N.15       0       7       1       6       0       1       1       0       1     0     0   1
    # N.16       1       3       4       1       7       3       3       6       4     0     0   4
    # N.2        2       6       2       1       9       5       2       3       6     0     0   3
    # N.3        2       3       1       4       0       3       0       0       2     0     0   0
    # N.4        2       5       3       1      10       7       3      10       6     1     1   3
    # N.5        0       4       3       2       4       4       1       2       5     0     0   2
    # N.6        1       2       3       0      13      10       6       8       9     0     3   2
    # N.7        0       4      10       1       1       3       1       2       2     0     0   1
    # N.8        1       7       4       4       6       6       2       3      19     1     1   3
    # N.9        0       3       1       1      10       5       2       5       4     0     0   1
    # OL         0       0       2       0       0       0       0       0       0     0    19   0
    # OPC        0       0       0       0       0       1       0       0       0     0     0  26
    # OPC.OL     0       0       0       1       0       0       0       0       1     0     5   7

  
  
## Amonst top 40 ===
genes.top40.hsap <- lapply(markerList.t.hsap, function(x){head(x, n=40)})

genes.top40.mm <- lapply(markerList.t.mm, function(x){toupper(head(x, n=40))})
genes.top40.mm <- sapply(genes.top40.mm, cbind)

sapply(genes.top40.hsap, function(x){
  apply(genes.top40.mm,2,function(y){length(intersect(x,y))})
})
    #       Astro Excit.1 Excit.2 Excit.3 Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5 Micro Oligo OPC
    # AS         7       0       0       0       0       0       0       0       0     0     0   0
    # EN         1       0       0       0       0       0       0       0       0     0     0   0
    # MG         0       0       0       0       0       0       0       0       0     4     0   0
    # MU         0       0       0       0       0       0       0       0       0     0     0   0
    # N.1        0       0       1       0       1       0       0       0       0     0     0   0
    # N.10       0       0       0       0       1       2       0       0       2     0     0   0
    # N.11       0       4       0       0       2       0       0       2       0     0     0   1
    # N.12       1       2       2       0       0       1       0       0       2     0     0   1
    # N.13       0       0       0       1       0       0       0       0       1     0     0   0
    # N.14       0       2       0       1       0       1       0       0       1     1     0   0
    # N.15       0       3       0       0       0       0       0       0       1     0     0   0
    # N.16       0       1       1       0       0       1       0       0       1     0     0   2
    # N.2        0       1       1       0       1       0       0       0       1     0     0   0
    # N.3        0       1       0       0       0       2       0       0       0     0     0   0
    # N.4        0       1       1       0       3       3       0       1       0     0     0   1
    # N.5        0       0       0       0       1       0       0       0       0     0     0   1
    # N.6        0       1       0       0       2       2       0       0       0     0     0   0
    # N.7        0       0       2       0       0       1       0       1       0     0     0   0
    # N.8        0       1       0       0       1       1       1       2       1     0     0   0
    # N.9        0       0       0       0       0       1       0       1       1     0     0   0
    # OL         0       0       1       0       0       0       0       0       0     0     7   0
    # OPC        0       0       0       0       0       0       0       0       0     0     0  10
    # OPC.OL     0       0       0       0       0       0       0       0       0     0     1   1

# Inhib.5 : N.8 genes ==
  intersect(genes.top40.hsap[["Inhib.5"]], genes.top100.mm[ ,"N.8"])
      # [1] "NPFFR2" "SV2C"   "OTOF"   "GRM8"   "OLFM3"  "FOXP2"

  # round(ts.mmMeA["49202", ],3)  # (Tll1 - looking because a highlighted gene in text)
  #     # AS      EN      MG      MU     N.1    N.10    N.11    N.12    N.13    N.14    N.15    N.16
  #     # -5.939  -5.932  -6.699   1.698   8.835   2.691 107.521  -5.323  20.345  86.122  -5.484  -5.423
  #     # N.2     N.3     N.4     N.5     N.6     N.7     N.8     N.9      OL     OPC  OPC.OL
  #     # 13.117  -5.297  33.339  16.283  -6.203  -5.520 108.310  22.783  -5.886  -4.273  -5.318
  # 
  plotExpression(sce.amy.mm, exprs_values="logcounts", x="subCluster", colour_by="subCluster", features="Tll1")
  #     # ahh nothing but a few outliers
  
  sce.amy.mm.sub <- sce.amy.mm[ ,grep("N.", sce.amy.mm$subCluster)]
  sce.amy.mm.sub$subCluster <- droplevels(sce.amy.mm.sub$subCluster)
  plotExpression(sce.amy.mm.sub, exprs_values="logcounts", x="subCluster", colour_by="subCluster",
                 features=c("Npffr2","Sv2c","Otof","Grm8","Olfm3","Foxp2"))
      # Actually nothing suuuper convicing - mostly outlier.  These just happen to have _more_ lol
  
  # N.8 top genes include Pcdh8 & Lamp5
  plotExpression(sce.amy.mm.sub, exprs_values="logcounts", x="subCluster", colour_by="subCluster",
                 features=c("Pcdh8","Lamp5"))

  # N.12 reported marker genes (reported in supplementals "mmc2.xlsx" with paper)
  plotExpression(sce.amy.mm.sub, exprs_values="logcounts", x="subCluster", colour_by="subCluster",
                 features=c("Eomes","Dsp","Nhlh2","Samd3","Trpc3","Cdhr1","Lhx1"))
      # Oh six of these were of the top 10 from my own test and plotted lol.  Well good.
  
  
# (and btw) ===
table(sce.amy$cellType.split, sce.amy$donor)


# Glucocorticoid receptors? (in relation to TLL1, as per https://doi.org/10.1016/j.molbrainres.2005.09.016)
plotExpression(sce.amy, exprs_values="logcounts", x="cellType.split", colour_by="cellType.split",
               features=c("NR3C1","NR3C2")) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                          geom = "crossbar", width = 0.3,
                                                          colour=rep(tableau20[1:12], 2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25))
    # No particular high/specific expression in Inhib.5



### FINAL GRANT VERSION ===
  # Remove EN, MU, OPC.OL, N.12 & N.15 (the latter two bc < 50 cells)
load("rdas/zTsMats_libd-AMY_and_ucla-mouseMeA-2019Cell_sharedGenes_25May2020.rda", verbose=T)
    # ts.amy, ts.mmMeA, Readme


cor_t_amy <- cor(ts.amy, ts.mmMeA)
rownames(cor_t_amy) = paste0(rownames(cor_t_amy),"_","H")
colnames(cor_t_amy) = paste0(colnames(cor_t_amy),"_","M")

# Re-order mouse labels - move EN/MG/MU to after neuronal subclusters
cor_t_amy <- cor_t_amy[ ,c(1, 5,13:20,6:12, 3,2,4, 21,23,22)]
# Remove those selected
cor_t_amy_sub <- cor_t_amy[ ,-which(colnames(cor_t_amy) %in% c("EN_M", "MU_M", "OPC.OL_M",
                                                               "N.12_M", "N.15_M"))]
range(cor_t_amy_sub)
    #[1] -0.2203968  0.5023080  --> Threshold to 0.5

cor_t_amy_sub <- ifelse(cor_t_amy_sub >= 0.5, 0.5, cor_t_amy_sub)

### Heatmap - typically use levelplot (e.g. below), but will want pheatmap bc can cluster cols/rows
theSeq.all = seq(-.5, .5, by = 0.01)
my.col.all <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.all)-1)

pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/HongLab-UCLA_mmMeA/overlap-mouseMeA-2019-fullSubclusters_with_LIBD-10x-AMY_SN-LEVEL-stats_FINAL_May2020.pdf",width=8)
pheatmap(cor_t_amy_sub,
         color=my.col.all,
         cluster_cols=F, cluster_rows=F,
         angle_col=90,
         breaks=theSeq.all,
         fontsize=11.5, fontsize_row=17, fontsize_col=15,
         legend_breaks=c(seq(-0.5,0.5,by=0.25)),
         main="Correlation of cluster-specific t's to mouse MeA \n subclusters (Chen-Hu-Wu et al., Cell 2019)")
dev.off()





## For supplement: Print top markers for 'Inhib.5' & corresponding in MeA 'N.8' === ===
# (load AMY SCE - already done in session)

# Prep mouse MeA
load("/dcl01/ajaffe/data/lab/singleCell/ucla_mouse-MeA/2019Cell/SCE_mouse-MeA_downstream-processing_MNT.rda", verbose=T)
    # sce.amy.mm, chosen.hvgs.amy.mm

sce.amy.mm.sub <- sce.amy.mm[ ,grep("N.", sce.amy.mm$subCluster)]
sce.amy.mm.sub$subCluster <- droplevels(sce.amy.mm.sub$subCluster)

genes2print <- c("Npffr2", "Tll1", "Grm8", "Foxp2")

pdf("pdfs/pubFigures/suppFig_AMY-vs-MeA_topInhib.5markers_MNTSep2020.pdf", height=2.5, width=5)
# Human AMY
print(
  plotExpression(sce.amy, exprs_values = "logcounts", features=toupper(genes2print),
                 x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=1.0, ncol=4,
                 add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                              geom = "crossbar", width = 0.3,
                                              colour=rep(tableau20[1:12], length(genes2print))) +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9.5),
          axis.title.y = element_text(angle = 90, size = 10),
          panel.grid.major=element_line(colour="grey95", size=0.8),
          panel.grid.minor=element_line(colour="grey95", size=0.4))
)

# mouse MeA
print(
  plotExpression(sce.amy.mm.sub, exprs_values = "logcounts", features=genes2print,
                 x="subCluster", colour_by="subCluster", point_alpha=0.5, point_size=1.0, ncol=4,
                 add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                              geom = "crossbar", width = 0.3,
                                              colour=rep(tableau20[1:16], length(genes2print))) +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9.5),
          axis.title.y = element_text(angle = 90, size = 10),
          panel.grid.major=element_line(colour="grey95", size=0.8),
          panel.grid.minor=element_line(colour="grey95", size=0.4))
)
dev.off()


## Heatmap version ===
# Take more overlapping, from above exploration
genes2print <- c("Npffr2", "Tll1", "Grm8", "Foxp2", "Sv2c", "Olfm3")

pdf("pdfs/pubFigures/suppFig_AMY-vs-MeA_topInhib.5markers_heatmap_MNTSep2020.pdf", width=5, height=5)
dat <- assay(sce.amy, "logcounts")
cell.idx <- splitit(sce.amy$cellType.split)
current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMeans(dat[toupper(genes2print), ii])))
pheatmap(current_dat, cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(0.02, 4.0, length.out = 101),
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "BuGn"))(100),
         fontsize_row = 18, fontsize_col=16)

dat <- assay(sce.amy.mm.sub, "logcounts")
cell.idx <- splitit(sce.amy.mm.sub$subCluster)
current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMeans(dat[genes2print, ii])))
pheatmap(current_dat, cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(0.02, 1, length.out = 101),
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "BuGn"))(100),
         fontsize_row = 16, fontsize_col=16)
dev.off()



