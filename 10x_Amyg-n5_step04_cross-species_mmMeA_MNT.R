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

source('plotExpressionCustom.R')

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


### Unmapped [-across species] markers ======
  # I.e. 'What are significantly enriched but weren't kept in the homology mapping?'

## Human vs. shared homology space ===
markerList.1vAll <- lapply(markers.amy.enriched, function(x){rownames(x)[x$log.FDR < log(0.05) & x$non0median==TRUE]})
length(unique(unlist(markerList.1vAll)))
    # 8483
length(setdiff(unique(unlist(markerList.1vAll)), rownames(sce.hsap.sub)))
    # 1154 - wow so ~1/8th of these markers were 'missing' in homology

    # of PW markers?
    markerList.pw <- lapply(markers.amy.t.pw, function(x){rownames(x)[x$FDR < 0.05 & x$non0median==TRUE]})
    length(unique(unlist(markerList.pw)))
        # 2461
    length(setdiff(unique(unlist(markerList.pw)), rownames(sce.hsap.sub)))
        # 388

## By cell type?
sapply(markerList.1vAll, function(x){length(setdiff(x, rownames(sce.hsap.sub)))})
    # Astro_A Astro_B    Endo Excit_A Excit_B Excit_C Inhib_A Inhib_B Inhib_C Inhib_D 
    #     169       8      38     662     131     317     517     330     266     298 
    # Inhib_E Inhib_F Inhib_G Inhib_H   Micro   Mural   Oligo     OPC   Tcell 
    #     238     316      86     129      59      51      78     121      30 

# % 'not conserved'?
round(sapply(markerList.1vAll, function(x){length(setdiff(x, rownames(sce.hsap.sub)))}) / lengths(markerList.1vAll) * 100, 3)
    # Astro_A Astro_B    Endo Excit_A Excit_B Excit_C Inhib_A Inhib_B Inhib_C Inhib_D 
    #  11.388   9.756   7.170  13.749  11.245  16.374  11.471   9.169  10.666   8.851 
    # Inhib_E Inhib_F Inhib_G Inhib_H   Micro   Mural   Oligo     OPC   Tcell 
    #  14.050  10.701  14.651   9.656   6.191   8.472   8.657   7.511   7.371
    # - mean 10.4%

        # PW results:
        sapply(markerList.pw, function(x){length(setdiff(x, rownames(sce.hsap.sub)))})
            # Astro_A Astro_B    Endo Excit_A Excit_B Excit_C Inhib_A Inhib_B Inhib_C Inhib_D 
            #      50       0      11      15      32      63       8      14       2       9 
            # Inhib_E Inhib_F Inhib_G Inhib_H   Micro   Mural   Oligo     OPC   Tcell 
            #      37      14       1       7      30      16      38      24      17
        
        # % 'not conserved'?
        round(sapply(markerList.pw, function(x){length(setdiff(x, rownames(sce.hsap.sub)))}) / lengths(markerList.pw) * 100, 3)
            #  Astro_A Astro_B    Endo Excit_A Excit_B Excit_C Inhib_A Inhib_B Inhib_C Inhib_D 
            #   15.337   0.000   6.111  25.424  18.497  44.366  30.769  15.909  14.286   6.870 
            #  Inhib_E Inhib_F Inhib_G Inhib_H   Micro   Mural   Oligo     OPC   Tcell 
            #   25.170  32.558  25.000   8.861  10.563   9.302  12.418  13.953  15.179 


## What are some of those genes?
sapply(markerList.1vAll, function(x){head(setdiff(x, rownames(sce.hsap.sub)))})



### Mouse MeA vs shared homology space
markerList.mm <- lapply(markers.mmMeA.t.1vAll, function(x){rownames(x)[x$log.FDR < log(0.05)]})
length(unique(unlist(markerList.mm)))
    # 14478
length(setdiff(unique(unlist(markerList.mm)), rownames(sce.mm.sub)))  # 838; 5.8% of genes

sapply(markerList.mm, head)

# Checking out a handful of these ('N.11' markers)
plotExpressionCustom(sce.mm.sub, anno_name="subCluster", features=c("Neurod6", "Olfm1", "Synpr", "Slc30a3", "Snap25", "Slc17a7"),
                     ncol=3, features_name="mmMeA custom-selected")

  # Observation - the 'markers' of these mouse clusters have quite sparse expression, that
  #               a lot of them don't exhibit optimal distribution patterns (Slc7a7)...
  #               In other words the results from the 'cluster-vs-all-others' marker
  #               test is heavily noise-driven (see Slc30a3, Neurod6)



## Intersecting some of the top ('1vAll') markers =====================
load("rdas/revision/markers-stats_Amyg-n5_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.amy.t.pw, markers.amy.wilcox.block, markers.amy.t.1vAll, medianNon0.amy
    rm(markers.amy.t.pw, markers.amy.wilcox.block, medianNon0.amy)

# Take top 100
markerList.t.hsap <- lapply(markers.amy.t.1vAll, function(x){
  rownames(x[[2]])[x[[2]]$log.FDR < log(1e-6) & x[[2]]$non0median==TRUE]
  }
)
genes.top100.hsap <- lapply(markerList.t.hsap, function(x){head(x, n=100)})


## mm MeA ===
load("/dcl01/ajaffe/data/lab/singleCell/ucla_mouse-MeA/2019Cell/markers-stats_mouseMeA-2019-with-16neuSubs_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.mmMeA.t.1vAll

## MNT 15Jun2021 - let's add the same marker restriction as applying to revision human data:
 #                 i.e. for any given subCluster, a marker must have non-0 median expression of that gene
load("/dcl01/ajaffe/data/lab/singleCell/ucla_mouse-MeA/2019Cell/SCE_mouse-MeA_downstream-processing_MNT.rda", verbose=T)
    # sce.amy.mm, chosen.hvgs.amy.mm

cellSubtype.idx <- splitit(sce.amy.mm$subCluster)
medianNon0.amy.mm <- lapply(cellSubtype.idx, function(x){
  apply(as.matrix(assay(sce.amy.mm, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

sapply(medianNon0.amy.mm, table)
    #          AS    EN    MG    MU   N.1   N.2   N.3   N.4   N.5   N.6   N.7   N.8
    # FALSE 14417 14400 14382 14420 14406 14401 14408 14373 14339 14371 14077 14330
    # TRUE     93   110   128    90   104   109   102   137   171   139   433   180
    
    #         N.9  N.10  N.11  N.12  N.13  N.14  N.15  N.16    OL   OPC OPC.OL
    # FALSE 14400 14493 14254 14234 14495 14485 14486 14392 14088 14250  14110
    # TRUE    110    17   256   276    15    25    24   118   422   260    400

# Add this to the markers stats object
for(i in names(markers.mmMeA.t.1vAll)){
  markers.mmMeA.t.1vAll[[i]] <- cbind(markers.mmMeA.t.1vAll[[i]],
                                           medianNon0.amy.mm[[i]][match(rownames(markers.mmMeA.t.1vAll[[i]]),
                                                                     names(medianNon0.amy.mm[[i]]))])
  colnames(markers.mmMeA.t.1vAll[[i]])[5] <- "non0median"
}

# Go ahead and save that into that markers .rda
save(markers.mmMeA.t.1vAll, medianNon0.amy.mm,
     file="/dcl01/ajaffe/data/lab/singleCell/ucla_mouse-MeA/2019Cell/markers-stats_mouseMeA-2019-with-16neuSubs_findMarkers-SN-LEVEL_MNTMay2020.rda")



# Map through the 'homologInfo' [if it exists in the shared space]
markerList.t.mm <- lapply(markers.mmMeA.t.1vAll, function(x){
  rownames(x)[x$log.FDR < log(1e-6) & x$non0median==TRUE]
  }
)
    # MNT comment: this may be too strict, given the technology and that looking at
    #              mean expression looks more somewhat convincing
    #            - Will leave the below numbers from the iteration without `non0median==TRUE`,
    #              unless otherwise specified
genes.top100.mm <- lapply(markerList.t.mm, function(x){head(x, n=100)})

# Convert to H. sap homolog [if applicable]
for(i in names(genes.top100.mm)){ 
  homolog.space.idx <- which(genes.top100.mm[[i]] %in% homologInfo$M.mus.ID)
  genes.top100.mm[[i]] <- c(genes.top100.mm[[i]][setdiff(1:100, homolog.space.idx)],
                            homologInfo$H.sap.ID[match(genes.top100.mm[[i]][homolog.space.idx],
                                                       homologInfo$M.mus.ID)])
}
  

genes.top100.mm <- sapply(genes.top100.mm, cbind)

## sapply
sapply(genes.top100.hsap, function(x){
  apply(genes.top100.mm,2,function(y){length(intersect(x,y))})
})
    #        Astro_A Astro_B Endo Excit_A Excit_B Excit_C Inhib_A Inhib_B Inhib_C
    # AS          19       8    0       0       1       0       1       0       0
    # EN           3       1   34       1       1       0       2       1       0
    # MG           0       0    1       0       0       0       0       0       0
    # MU           1       0    3       0       0       0       0       2       0
    # N.1          1       1    0       3       3       2       3       5       7
    # N.10         0       0    0       5       1       3       3       3       8
    # N.11         1       1    0      14       8       6       5       6       8
    # N.12         2       2    0      11       5       5       2       4       6
    # N.13         1       1    0       3       1       2       1       1       7
    # N.14         0       0    0       8       2       4       3       6      10
    # N.15         1       0    1       6       1       5       3       0       4
    # N.16         1       1    0       5       6       3       2       6       1
    # N.2          2       2    0       7       4       4       3       4       5
    # N.3          2       1    0       2       6       1       1       0       5
    # N.4          2       2    0       9       5       6       3      12       6
    # N.5          0       0    0       1       2       1       3       2       6
    # N.6          1       0    1       1       2       0       0       9       5
    # N.7          0       0    0       6       8       2       2       3       3
    # N.8          1       1    0       5       6       4       4       4      17
    # N.9          0       0    0       1       2       1       1       4       4
    # OL           0       0    0       0       0       0       0       0       0
    # OPC          0       0    2       0       0       0       2       0       1
    # OPC.OL       0       0    2       0       0       0       1       0       1

    #        Inhib_D Inhib_E Inhib_F Inhib_G Inhib_H Micro Mural Oligo OPC Tcell
    # AS           2       0       0       0       0     0     1     0   1     0
    # EN           2       1       0       0       2     0     4     0   1     3
    # MG           0       0       0       0       1    19     2     0   0     4
    # MU           1       0       1       0       0     0    28     0   0     0
    # N.1          3       1       5       3       5     0     0     1   1     1
    # N.10         5       0       4       4       2     0     0     0   0     0
    # N.11         4       2       5       4       5     0     0     0   4     0
    # N.12         5       3       4       4       2     0     2     2   1     1
    # N.13         0       3       1       4       3     1     0     0   1     0
    # N.14         5       3       4       7       2     1     0     1   2     0
    # N.15         0       3       1       2       4     0     0     0   1     1
    # N.16         2       2       5       4       5     0     0     0   4     0
    # N.2          2       2       2       2       5     0     0     0   3     1
    # N.3          2       1       0       2       1     0     1     0   0     0
    # N.4          9       2       2       2       2     1     0     1   2     0
    # N.5          1       2       1       4       6     0     0     0   2     0
    # N.6          8       1       9       3       4     0     0     3   2     0
    # N.7          1       4       2       3       1     0     2     0   0     1
    # N.8          6       4       2       5       4     0     0     1   3     0
    # N.9          1       2       2       3       2     0     1     0   0     1
    # OL           0       0       2       0       0     0     0    23   0     0
    # OPC          1       0       2       1       0     0     0     0  26     0
    # OPC.OL       0       0       0       0       0     0     0     5   7     0



# Inhib_C : N.8 genes ==
  intersect(genes.top100.hsap[["Inhib_C"]], genes.top100.mm[ ,"N.8"])
    # [1] "NPFFR2"  "FOXP2"   "OLFM3"   "OTOF"    "EPHA6"   "SV2C"    "ANO3"   
    # [8] "RASGRP1" "SYT1"    "VSTM2A"  "KCNJ3"   "CACNA1E" "GRIA1"   "CPNE5"  
    # [15] "PLPPR4"  "KCNH5"   "GABRB3"
  
      # # (also closely correlating across t's)
      # printThese <- intersect(genes.top100.hsap[["Inhib_C"]], genes.top100.mm[ ,"N.5"])
      #     # "NELL2"  "SYT1"   "GRIA1"  "PLPPR4" "GABRB3"
      #     # Oh well these are just a subset of what's shared b/tw 'N.8' & Inhib_C...
      # 
      # plotExpressionCustom(sce.amy.mm, anno_name="subCluster", features_name="subCluster", ncol=3,
      #                      features=homologInfo$M.mus.ID[match(printThese, homologInfo$H.sap.ID)])
      #     # These aren't 'N.5'-specific/most-expressing...

  plotExpressionCustom(sce.mm.sub, anno_name="subCluster", features_name="subCluster", ncol=3,
                       features=c('Npffr2','Foxp2','Olfm3','Otof','Epha6','Sv2c'))
      # Actually nothing suuuper convicing - mostly outlier.  These just happen to have _more_ lol
  
  # N.8 top genes include Pcdh8 & Lamp5 - these are way more convincing
  plotExpressionCustom(sce.mm.sub, anno_name="subCluster", features_name="subCluster", ncol=3,
                 features=c("Pcdh8","Lamp5"))

      # Employing a non0median Boolean like with the human data:
      printThese <- intersect(genes.top100.hsap[["Inhib_C"]], genes.top100.mm[ ,"N.8"])
          # [1] "NELL2"  "SYT1"   "VSTM2A" "GRIA1"  "PLPPR4" "KCNH5"  "GABRB3"
      plotExpressionCustom(sce.amy.mm, anno_name="subCluster", features_name="subCluster", ncol=3,
                           features=homologInfo$M.mus.ID[match(printThese, homologInfo$H.sap.ID)])
          # These look good.  Kcnh5 is the most 'N.8' defining
          which(genes.top100.hsap[["Inhib_C"]] == "KCNH5")  # 90, nice
      
  

# (and btw) ===
# Human AMY (n=5)
table(sce.amy$cellType, sce.amy$donor)
    #         br5161 br5212 br5276 br5400 br5701
    # Astro_A    484    350    230    111    380
    # Astro_B      7     10     49     12      5
    # Endo         0      0     21      7      3
    # Excit_A    106    203      5     14     16
    # Excit_B      0     39      0      0      5
    # Excit_C      5     43      0      7      0
    # Inhib_A      0      0    366    362      0
    # Inhib_B     36    115     71    245     74
    # Inhib_C    128     17    284     85     11
    # Inhib_D     36     75    124    271     49
    # Inhib_E      0      0    405      7      2
    # Inhib_F     24     68     36     81      7
    # Inhib_G      0      0     76      9      1
    # Inhib_H      0      0     50      2      0
    # Micro      411    304     14    117    355
    # Mural        2      0     24      7      6
    # Oligo     1688   1736    304    309   2043
    # OPC        340    290    199     93    537
    # Tcell        3      7      3      3     15

# mm MeA (n=13)
table(sce.amy.mm$subCluster)
    #    AS     EN     MG     MU    N.1    N.2    N.3    N.4    N.5    N.6    N.7 
    # 12640   2467   6366   1114   2376   1586    694    594    234    291    235 
    #   N.8    N.9   N.10   N.11   N.12   N.13   N.14   N.15   N.16     OL    OPC 
    #   196    176   1768    176     43   1440    195     38     81   4770   5226 
    #OPC.OL 
    #   639


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

dat <- assay(sce.mm.sub, "logcounts")
cell.idx <- splitit(sce.mm.sub$subCluster)
current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMeans(dat[genes2print, ii])))
pheatmap(current_dat, cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(0.02, 1, length.out = 101),
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "BuGn"))(100),
         fontsize_row = 16, fontsize_col=16)
dev.off()


### Session info for 16Jun2021 ===============
sessionInfo()
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
#   [1] pheatmap_1.0.12             RColorBrewer_1.1-2          lattice_0.20-41            
# [4] limma_3.46.0                jaffelab_0.99.30            rafalib_1.0.0              
# [7] DropletUtils_1.10.3         batchelor_1.6.2             scran_1.18.5               
# [10] scater_1.18.6               ggplot2_3.3.3               org.Hs.eg.db_3.12.0        
# [13] EnsDb.Hsapiens.v86_2.99.0   ensembldb_2.14.1            AnnotationFilter_1.14.0    
# [16] GenomicFeatures_1.42.3      AnnotationDbi_1.52.0        SingleCellExperiment_1.12.0
# [19] SummarizedExperiment_1.20.0 Biobase_2.50.0              GenomicRanges_1.42.0       
# [22] GenomeInfoDb_1.26.7         IRanges_2.24.1              S4Vectors_0.28.1           
# [25] BiocGenerics_0.36.1         MatrixGenerics_1.2.1        matrixStats_0.58.0         
# 
# loaded via a namespace (and not attached):
#   [1] googledrive_1.0.1         ggbeeswarm_0.6.0          colorspace_2.0-0         
# [4] ellipsis_0.3.2            scuttle_1.0.4             bluster_1.0.0            
# [7] XVector_0.30.0            BiocNeighbors_1.8.2       rstudioapi_0.13          
# [10] farver_2.1.0              bit64_4.0.5               fansi_0.4.2              
# [13] xml2_1.3.2                splines_4.0.4             R.methodsS3_1.8.1        
# [16] sparseMatrixStats_1.2.1   cachem_1.0.4              Rsamtools_2.6.0          
# [19] ResidualMatrix_1.0.0      dbplyr_2.1.1              R.oo_1.24.0              
# [22] HDF5Array_1.18.1          compiler_4.0.4            httr_1.4.2               
# [25] dqrng_0.2.1               assertthat_0.2.1          Matrix_1.3-2             
# [28] fastmap_1.1.0             lazyeval_0.2.2            BiocSingular_1.6.0       
# [31] prettyunits_1.1.1         tools_4.0.4               rsvd_1.0.3               
# [34] igraph_1.2.6              gtable_0.3.0              glue_1.4.2               
# [37] GenomeInfoDbData_1.2.4    dplyr_1.0.5               rappdirs_0.3.3           
# [40] Rcpp_1.0.6                vctrs_0.3.6               Biostrings_2.58.0        
# [43] rhdf5filters_1.2.0        rtracklayer_1.50.0        DelayedMatrixStats_1.12.3
# [46] stringr_1.4.0             beachmat_2.6.4            lifecycle_1.0.0          
# [49] irlba_2.3.3               statmod_1.4.35            XML_3.99-0.6             
# [52] edgeR_3.32.1              zlibbioc_1.36.0           scales_1.1.1             
# [55] hms_1.0.0                 ProtGenerics_1.22.0       rhdf5_2.34.0             
# [58] curl_4.3                  memoise_2.0.0             gridExtra_2.3            
# [61] segmented_1.3-3           biomaRt_2.46.3            stringi_1.5.3            
# [64] RSQLite_2.2.7             BiocParallel_1.24.1       rlang_0.4.10             
# [67] pkgconfig_2.0.3           bitops_1.0-7              purrr_0.3.4              
# [70] Rhdf5lib_1.12.1           labeling_0.4.2            GenomicAlignments_1.26.0 
# [73] cowplot_1.1.1             bit_4.0.4                 tidyselect_1.1.1         
# [76] magrittr_2.0.1            R6_2.5.0                  generics_0.1.0           
# [79] DelayedArray_0.16.3       DBI_1.1.1                 pillar_1.6.0             
# [82] withr_2.4.2               RCurl_1.98-1.3            tibble_3.1.1             
# [85] crayon_1.4.1              utf8_1.2.1                BiocFileCache_1.14.0     
# [88] viridis_0.6.0             progress_1.2.2            locfit_1.5-9.4           
# [91] grid_4.0.4                blob_1.2.1                digest_0.6.27            
# [94] R.utils_2.10.1            openssl_1.4.3             munsell_0.5.0            
# [97] beeswarm_0.3.1            viridisLite_0.4.0         vipor_0.4.5              
# [100] askpass_1.1 
