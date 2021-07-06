### MNT 10x snRNA-seq workflow: step 04 - downstream comparisons
###   **Region-specific analyses**
###     - (5x) NAc samples from Oct2020
###     - (3x) revision samples, incl'g female donors
###   * Comparison to Jeremy Day Lab's rat NAc samples (n=4)
#####################################################################

library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
#BiocManager::install("org.Rn.eg.db")
library(org.Rn.eg.db)
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
# JAX MGI database no longer reports a $JAX.geneID -> now use $DB.Class.Key
# (corresponded with David Shaw @ JAX/MGI)

# For mapping === == === ===
# human.entrez > DB.Class.Key < mouse.entrez < rn.Symbol
#                ^ filter SCE's on this - to be more descriptive, call 'JAX.geneID'


### Setting up homologous gene IDs, for mapping b/tw species =============

## Load SCE with new info
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda",
     verbose=T)
    # sce.nac, chosen.hvgs.nac, pc.choice.nac, ref.sampleInfo, annotationTab.nac, cell_colors.nac

# First drop the "drop." clusters (669 nuclei, total)
sce.nac <- sce.nac[ ,-grep("drop.", sce.nac$cellType)]
sce.nac$cellType <- droplevels(sce.nac$cellType)

# Drop genes with all 0's
sce.nac <- sce.nac[!rowSums(assay(sce.nac, "counts"))==0, ]


# Add EntrezID for human
hs.entrezIds <- mapIds(org.Hs.eg.db, keys=rowData(sce.nac)$gene_id, 
                       column="ENTREZID", keytype="ENSEMBL")
    # "'select()' returned 1:many mapping between keys and columns"
table(!is.na(hs.entrezIds))
    # 21,191 valid entries (remember this is already subsetted for those non-zero genes only)
    withoutEntrez <- names(hs.entrezIds)[is.na(hs.entrezIds)]
    # Store those somewhere, maybe for later reference
    table(rowData(sce.nac)[rowData(sce.nac)$gene_id %in% withoutEntrez, ]$gene_id == withoutEntrez)
    names(withoutEntrez) <- rowData(sce.nac)[rowData(sce.nac)$gene_id %in% withoutEntrez, ]$gene_name


# Add to rowData
rowData(sce.nac) <- cbind(rowData(sce.nac), hs.entrezIds)


## Bring in 'DB.Class.Key' for human ===
# JAX annotation info
hom = read.delim("http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt",
                 as.is=TRUE)

hom_hs <- hom[hom$Common.Organism.Name == "human", ]
    # of 22,514 entries
table(rowData(sce.nac)$hs.entrezIds %in% hom_hs$EntrezGene.ID)
    # 17,063
table(rowData(sce.nac)$gene_name %in% hom_hs$Symbol)
    # 16,739 - not a bad difference

# Human JAX.geneID (by Entrez)
rowData(sce.nac)$JAX.geneID <- hom_hs$DB.Class.Key[match(rowData(sce.nac)$hs.entrezIds,
                                                       hom_hs$EntrezGene.ID)]


## Do the same for rat NAc ===
load("/dcl01/ajaffe/data/lab/singleCell/day_rat_snRNAseq/SCE_rat-NAc_downstream-processing_MNT.rda", verbose=T)
    # sce.nac.rat, chosen.hvgs.nac.rat
    # (This object was already processed to contain the now-deprecated $HomoloGene.ID)

# Find/add Entrez ID
rn.entrezIds <- mapIds(org.Rn.eg.db, keys=rowData(sce.nac.rat)$ID, 
                       column="ENTREZID", keytype="ENSEMBL")

table(!is.na(hs.entrezIds)) #21866

# Add to rowData
rowData(sce.nac.rat) <- cbind(rowData(sce.nac.rat), rn.entrezIds)

hom_rn <- hom[hom$Common.Organism.Name == "rat", ]
    # of 17,655 entries
table(rowData(sce.nac.rat)$rn.entrezIds %in% hom_rn$EntrezGene.ID)
    # 16,201
table(rowData(sce.nac.rat)$Symbol %in% hom_rn$Symbol)
    # 15,893

# Now add rat JAX.geneID (by Entrez)
rowData(sce.nac.rat)$JAX.geneID <- hom_rn$DB.Class.Key[match(rowData(sce.nac.rat)$rn.entrezIds,
                                                         hom_rn$EntrezGene.ID)]


## Now set/match to shared homologous genes ===
length(intersect(rowData(sce.nac)$JAX.geneID,
                 rowData(sce.nac.rat)$JAX.geneID))  # 14,008

sharedHomologs <- intersect(rowData(sce.nac)$JAX.geneID,
                            rowData(sce.nac.rat)$JAX.geneID)
    # That first one is NA - rm
sharedHomologs <- sharedHomologs[-1]

# Human not in rat
length(setdiff(rowData(sce.nac)$JAX.geneID,
                 rowData(sce.nac.rat)$JAX.geneID))  # 2524
# Rat not in human
length(setdiff(rowData(sce.nac.rat)$JAX.geneID,
               rowData(sce.nac)$JAX.geneID))  # 1741


# Subset for those
sce.rn.sub <- sce.nac.rat[rowData(sce.nac.rat)$JAX.geneID %in% sharedHomologs, ]   # 14213
sce.hsap.sub <- sce.nac[rowData(sce.nac)$JAX.geneID %in% sharedHomologs, ]  # 14408
    ## Many are duplicated...

rowData(sce.rn.sub)$Symbol[duplicated(rowData(sce.rn.sub)$JAX.geneID)]
    # many orthologs
rowData(sce.hsap.sub)$gene_name[duplicated(rowData(sce.hsap.sub)$JAX.geneID)]
    # many, almost double


### -> Take the higher-expressing of the duplicated

    ## Rat ===
    duplicatedSet.rat <- which(duplicated(rowData(sce.rn.sub)$JAX.geneID))
    genes2compare.rat <- list()
    gene2keep.rat <- character()
    for(g in 1:length(duplicatedSet.rat)){
      genes2compare.rat[[g]] <- rownames(sce.rn.sub)[rowData(sce.rn.sub)$JAX.geneID ==
                                              rowData(sce.rn.sub)$JAX.geneID[duplicatedSet.rat[g]]]
      rowmeansmat <- rowMeans(assay(sce.rn.sub[genes2compare.rat[[g]], ], "logcounts"))
      gene2keep.rat[g] <- names(rowmeansmat[order(rowmeansmat, decreasing=TRUE)])[1]
    }
    
    length(genes2compare.rat) # 206
    length(unique(gene2keep.rat)) # 181
        # This is because many 'tested' might have been orthologous,
        #   b/tw themselves (i.e. 3+ orthologous genes):
    length(unique(rowData(sce.rn.sub)$JAX.geneID[duplicatedSet.rat])) # 181 - good
    
    genesNoCompare.rat <- rownames(sce.rn.sub)[!(rownames(sce.rn.sub) %in% unlist(genes2compare.rat))]
    
    # Finally combine and subset
    sce.rn.sub <- sce.rn.sub[c(genesNoCompare.rat, unique(gene2keep.rat)), ]
    
    table(rowData(sce.rn.sub)$JAX.geneID %in% sharedHomologs) # 14007 TRUE
    table(duplicated(rowData(sce.rn.sub)$JAX.geneID)) # 14007 FALSE         dope.

    
    ## Human ===
    # First change rownames to EnsemblID
    rownames(sce.hsap.sub) <- rowData(sce.hsap.sub)$gene_id
        
    duplicatedSet.hsap <- which(duplicated(rowData(sce.hsap.sub)$JAX.geneID))
    genes2compare.hsap <- list()
    gene2keep.hsap <- character()
    for(g in 1:length(duplicatedSet.hsap)){
      genes2compare.hsap[[g]] <- rownames(sce.hsap.sub)[rowData(sce.hsap.sub)$JAX.geneID ==
                                                          rowData(sce.hsap.sub)$JAX.geneID[duplicatedSet.hsap[g]]]
      rowmeansmat <- rowMeans(assay(sce.hsap.sub[genes2compare.hsap[[g]], ], "logcounts"))
      gene2keep.hsap[g] <- names(rowmeansmat[order(rowmeansmat, decreasing=TRUE)])[1]
    }
    
    length(gene2keep.hsap)  # 401
    length(unique(gene2keep.hsap))  # 274
        # This is because many 'tested' might have been orthologous,
        #   b/tw themselves (i.e. 3+ orthologous genes):
    length(unique(rowData(sce.hsap.sub)$JAX.geneID[duplicatedSet.hsap]))  # 274 - good
    
    genesNoCompare.hsap <- rownames(sce.hsap.sub)[!(rownames(sce.hsap.sub) %in% unlist(genes2compare.hsap))]
        # of length 13733 (which + 274 == 14007)
    
    # Finally combine and subset
    sce.hsap.sub <- sce.hsap.sub[c(genesNoCompare.hsap, unique(gene2keep.hsap)), ]
    
    table(rowData(sce.hsap.sub)$JAX.geneID %in% sharedHomologs) # 14007 TRUE
    table(duplicated(rowData(sce.hsap.sub)$JAX.geneID)) # 14007 FALSE         dope.


    ## Match order and save
    sce.rn.sub <- sce.rn.sub[match(rowData(sce.hsap.sub)$JAX.geneID,
                                   rowData(sce.rn.sub)$JAX.geneID), ]

    table(rowData(sce.rn.sub)$JAX.geneID == rowData(sce.hsap.sub)$JAX.geneID)
        # all TRUE - good

    Readme <- "These two SCEs are subsetted and ordered for matching 'JAX.geneID' in the rowData. This can be used to subset the nucleus-level SCEs in their respected Rdata files."
    save(sce.rn.sub, sce.hsap.sub, Readme, file="/dcl01/ajaffe/data/lab/singleCell/day_rat_snRNAseq/SCE_rat-NAc_matched2_Hsap-NAc_JAX.geneIDs_MNT2021.rda")
    # Also save in the main project dir
    save(sce.rn.sub, sce.hsap.sub, Readme, file="rdas/revision/SCE_rat-NAc_matched2_Hsap-NAc_JAX.geneIDs_MNT2021.rda")
    
    
    ## MNT 06Jul2021: oops - the rat SCE wasn't first subsetted for 0-capture genes ===
    load("rdas/revision/SCE_rat-NAc_matched2_Hsap-NAc_JAX.geneIDs_MNT2021.rda", verbose=T)
        # sce.rn.sub, sce.hsap.sub, Readme
    
    sce.rn.sub <- sce.rn.sub[!rowSums(assay(sce.rn.sub, "counts"))==0, ]
    # Subset the corresponding human-NAc SCE
    sce.hsap.sub <- sce.hsap.sub[rowData(sce.hsap.sub)$JAX.geneID %in% rowData(sce.rn.sub)$JAX.geneID, ]
    # re-save
    save(sce.rn.sub, sce.hsap.sub, Readme, file="/dcl01/ajaffe/data/lab/singleCell/day_rat_snRNAseq/SCE_rat-NAc_matched2_Hsap-NAc_JAX.geneIDs_MNT2021.rda")
    # Also save in the main project dir
    save(sce.rn.sub, sce.hsap.sub, Readme, file="rdas/revision/SCE_rat-NAc_matched2_Hsap-NAc_JAX.geneIDs_MNT2021.rda")
    
    
    
### FINALLY resume comparisons === === === === ===




## Write out .mtx & colData for Day Lab ========================================
library(Matrix)

load("rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda", verbose=T)
     # sce.nac, chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac, ref.sampleInfo

# Mtx
mat2write <- assay(sce.nac, "counts")
Matrix::writeMM(mat2write, file="pdfs/exploration/DayLab-ratNAc/10xCounts/libd_n3-hom_n2-NeuN_countMat.mtx")
    ## NULL     - uh what?  Lol

# Features
write.table(rowData(sce.nac), file="pdfs/exploration/DayLab-ratNAc/10xCounts/libd_n3-hom_n2-NeuN_features.tsv",
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Barcodes
write.table(sce.nac$Barcode, file="pdfs/exploration/DayLab-ratNAc/10xCounts/libd_n3-hom_n2-NeuN_barcodes.tsv",
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

table(rownames(sce.test) == rowData(sce.nac)$ID)  # all TRUE - so make rownames(sce.test)

rownames(sce.test) <- rownames(sce.nac)

all.equal(assay(sce.test, "counts"), assay(sce.nac, "counts"))
    ## [1] TRUE - dope


## Write out colData
colnames(colData(sce.nac))
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
    # (because `table(rownames(colData(sce.nac)) == sce.nac$Barcode)` == all TRUE)
pheno2write <- colData(sce.nac)[ ,c(2, 13, 3,4, 11,12, 14:18, 25)]

write.csv(pheno2write, row.names=FALSE, file="pdfs/exploration/DayLab-ratNAc/10xCounts/libd_n3-hom_n2-NeuN_metadata.csv")






### Comparison to Day Lab Rat with SN-LEVEL stats =============================================

# Load rat SCE
load("/dcl01/ajaffe/data/lab/singleCell/day_rat_snRNAseq/SCE_rat-NAc_downstream-processing_MNT.rda", verbose=T)
    # sce.nac.rat, chosen.hvgs.nac.rat
  
load("/dcl01/ajaffe/data/lab/singleCell/day_rat_snRNAseq/markers-stats_DayLab-ratNAc_findMarkers-SN-LEVEL_MNTMay2020.rda",
     verbose=T)
    # markers.rat.t.1vAll

table(rownames(sce.rn.sub) %in% rownames(markers.rat.t.1vAll[[1]]))
    #FALSE  TRUE 
    #  562 13445    # Why is this... (should be all 14007 TRUE)

    # Ohhh this is bc not all genes in that SCE are expressed:
    table(rowSums(assay(sce.rn.sub, "counts"))==0)  # 562

## Calculate and add t-statistic (= std.logFC * sqrt(N)) for rat clusters
 #      and fix row order to the first entry "Astrocyte"
fixTo <- rownames(markers.rat.t.1vAll[["Astrocyte"]])
for(x in names(markers.rat.t.1vAll)){
  markers.rat.t.1vAll[[x]]$t.stat <- markers.rat.t.1vAll[[x]]$std.logFC * sqrt(ncol(sce.nac.rat))
  markers.rat.t.1vAll[[x]] <- markers.rat.t.1vAll[[x]][fixTo, ]
}

# Pull out the t's
ts.rat <- sapply(markers.rat.t.1vAll, function(x){x$t.stat})
rownames(ts.rat) <- rowData(sce.rn.sub)$JAX.geneID[match(fixTo, rownames(sce.rn.sub))]



## Human t stats subset/re-ordering ===
# Bring in human stats; create t's
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_NAc-n5_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.nac.t.design, markers.nac.t.1vAll
    rm(markers.nac.t.design)

# Need to add t's with N nuclei used in constrasts

load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda",
     verbose=T)
    # sce.nac, chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac,ref.sampleInfo
    rm(chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac,ref.sampleInfo)
    # First drop "ambig.lowNtrxts" (93 nuclei)
    sce.nac <- sce.nac[ ,sce.nac$cellType.final != "ambig.lowNtrxts"]
    sce.nac$cellType.final <- droplevels(sce.nac$cellType.final)

## As above, calculate and add t-statistic (= std.logFC * sqrt(N)) for rat clusters
 #      and fix row order to the first entry "Astrocyte"
fixTo <- rownames(markers.nac.t.1vAll[["Astro"]])

for(s in names(markers.nac.t.1vAll)){
  markers.nac.t.1vAll[[s]]$t.stat <- markers.nac.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.nac))
  markers.nac.t.1vAll[[s]] <- markers.nac.t.1vAll[[s]][fixTo, ]
}

# Pull out the t's
ts.nac <- sapply(markers.nac.t.1vAll, function(x){x$t.stat})
rownames(ts.nac) <- fixTo



## Bring in HomoloGene.ID info to subset/match order
load("/dcl01/ajaffe/data/lab/singleCell/day_rat_snRNAseq/SCE_rat-NAc-PBd_w_matchingHsap-NAc-PBd_HomoloGene.IDs_MNT.rda",
     verbose=T)
    # sce.rn.sub, sce.hsap.sub, Readme

table(rowData(sce.rn.sub)$JAX.geneID == rowData(sce.hsap.sub)$JAX.geneID)  # all TRUE - dope
    # (see above - these are the intersecting homologs)

# First give [human] ts.nac rownames their respective EnsemblID
#   (have to use the full sce bc rownames(sce.hsap.sub) is EnsemblID and we uniquified the $Symbol)
rownames(ts.nac) <- rowData(sce.nac)$ID[match(rownames(ts.nac), rownames(sce.nac))]

# How many human genes with rat homologs are in these?
table(rownames(sce.hsap.sub) %in% rownames(ts.nac)) # 14121 good

# Subset/re-order for these and set to HomoloGene.ID
ts.nac <- ts.nac[rownames(sce.hsap.sub), ]
rownames(ts.nac) <- rowData(sce.hsap.sub)$JAX.geneID


# Same for rat t's
table(rownames(sce.rn.sub) %in% rownames(ts.rat))

ts.rat <- ts.rat[rownames(sce.rn.sub), ]
rownames(ts.rat) <- rowData(sce.rn.sub)$JAX.geneID

table(rownames(ts.nac) == rownames(ts.rat)) # all 14121 TRUE (well duh)


    # Save the ts matrices to reduce work next time
    Readme <- "These t-statistic matrices are subsetted and matched for shared 'HomoloGene.ID', so `cor()` can simply be run or other gene subsets applied first."
    save(ts.nac, ts.rat, Readme, file="rdas/zTsMats_libd-NAc_and_DayLab-ratNAc_sharedGenes_29May2020.rda")



cor_t_nac <- cor(ts.nac, ts.rat)
rownames(cor_t_nac) = paste0(rownames(cor_t_nac),"_H")
colnames(cor_t_nac) = paste0(colnames(cor_t_nac),"_R")
range(cor_t_nac)
    # [1] -0.3791969  0.6132932

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
#cor_t_nac <- cor_t_nac[ ,c(1,16, 7,15,6, 9,10, 6,8, 2,5,4,3,11,12,14)] # oops - duplicated col (6)
cor_t_nac <- cor_t_nac[ ,c(1,16, 7,15,6, 9,10, 8, 2,5,4,3,11:14)]
cor_t_nac.th <- cor_t_nac.th[ ,c(1,16, 7,15,6, 9,10, 8, 2,5,4,3,11:14)]

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



### A few more iterations (These looked great with the mm MeA revision analyses) ======
## On just hsap cluster-specific homologous genes ===
hsap_specific_indices = mapply(function(t) {
  oo = order(t, decreasing = TRUE)[1:100]
  },
as.data.frame(ts.nac)
)
hsap_ind = unique(as.numeric(hsap_specific_indices))
length(hsap_ind)  # so of 1400 (100 x 14 cellType), 1077 unique

cor_t_hsap = cor(ts.nac[hsap_ind, ],
                 ts.rat[hsap_ind, ])
rownames(cor_t_hsap) = paste0(rownames(cor_t_hsap),"_","H")
colnames(cor_t_hsap) = paste0(colnames(cor_t_hsap),"_","R")
range(cor_t_hsap)
    # [1] -0.3872332  0.7588219


## On just rat cluster-specific homologous genes ===
rat_specific_indices = mapply(function(t) {
  oo = order(t, decreasing = TRUE)[1:100]
},
as.data.frame(ts.rat)
)
rat_ind = unique(as.numeric(rat_specific_indices))
length(rat_ind)  # so of 1600 (100 x 16 subCluster), 1210 unique

cor_t_rat = cor(ts.nac[rat_ind, ],
                  ts.rat[rat_ind, ])
rownames(cor_t_rat) = paste0(rownames(cor_t_rat),"_","H")
colnames(cor_t_rat) = paste0(colnames(cor_t_rat),"_","R")
range(cor_t_rat)
    # [1] -0.4604021  0.7284520

## Between these gene spaces:
length(intersect(rownames(ts.nac)[hsap_ind], rownames(ts.rat)[rat_ind]))
    # 476

toptop.genes <- intersect(rownames(ts.nac)[hsap_ind], rownames(ts.rat)[rat_ind])
cor_t_top <- cor(ts.nac[toptop.genes, ],
                 ts.rat[toptop.genes, ])
rownames(cor_t_top) = paste0(rownames(cor_t_top),"_","H")
colnames(cor_t_top) = paste0(colnames(cor_t_top),"_","R")
range(cor_t_top)
    # [1] -0.4404794  0.8510041

# Print these
theSeq.all = seq(-.75, .75, by = 0.01)
my.col.all <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.all)-1)

# Reorder, as before
#cor_t_nac <- cor_t_nac[ ,c(1,16, 7,15,6, 9,10, 6,8, 2,5,4,3,11,12,14)] OOPS
cor_t_nac <- cor_t_nac[ ,c(1,16, 7,15,6, 9,10, 8, 2,5,4,3,11:14)]
cor_t_hsap <- cor_t_hsap[ ,c(1,16, 7,15,6, 9,10, 8, 2,5,4,3,11:14)]
cor_t_rat <- cor_t_rat[ ,c(1,16, 7,15,6, 9,10, 8, 2,5,4,3,11:14)]
cor_t_top <- cor_t_top[ ,c(1,16, 7,15,6, 9,10, 8, 2,5,4,3,11:14)]

pdf("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/DayLab-ratNAc/overlap-DayLab-ratNAc_with_LIBD-10x-NAc-n5_top-X-iterations_MNT2021.pdf")
pheatmap(cor_t_nac,
         color=my.col.all,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.all,
         fontsize=11, fontsize_row=13, fontsize_col=12,
         display_numbers=T, number_format="%.2f", fontsize_number=6,
         legend_breaks=c(seq(-0.75,0.75,by=0.375)),
         main="Correlation of cluster-specific t's to rat NAc \n subclusters (Savell et al., Sci Adv 2020)")
# On human-specific genes
pheatmap(cor_t_hsap,
         color=my.col.all,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.all,
         fontsize=11, fontsize_row=13, fontsize_col=12,
         display_numbers=T, number_format="%.2f", fontsize_number=6,
         legend_breaks=c(seq(-0.75,0.75,by=0.375)),
         main="Correlation of top-100 cluster-specific t's (1077) to \n (Savell et al., Sci Adv 2020) subclusters")
# On rat-NAc-specific genes
pheatmap(cor_t_rat,
         color=my.col.all,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.all,
         fontsize=11, fontsize_row=13, fontsize_col=12,
         display_numbers=T, number_format="%.2f", fontsize_number=6,
         legend_breaks=c(seq(-0.75,0.75,by=0.375)),
         main="Correlation of LIBD-NAc subclusters to \n (Savell et al., Sci Adv 2020) subcluster top-100 t's (1210)")
# On intersection between the top spp.-specific genes (476 genes)
theSeq.new = seq(-.85, .85, by = 0.01)
my.col.new <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.new)-1)
pheatmap(cor_t_top,
         color=my.col.new,
         cluster_cols=F, cluster_rows=F,
         breaks=theSeq.new,
         fontsize=9.5, fontsize_row=13, fontsize_col=12,
         display_numbers=T, number_format="%.2f", fontsize_number=6,
         legend_breaks=c(seq(-0.85,0.85,by=0.425)),
         main="Correlation of LIBD-NAc subclusters to \n (Savell et al., Sci Adv 2020) subcluster t's (shared top 100's, 476)")
dev.off()




## Final small exploration: Rat Grm8-MSNs ~ MSN.D1.3? ================
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_NAc-n5_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.nac.t.design, markers.nac.t.1vAll
load("/dcl01/ajaffe/data/lab/singleCell/day_rat_snRNAseq/markers-stats_DayLab-ratNAc_findMarkers-SN-LEVEL_MNTMay2020.rda",
     verbose=T)
    # markers.rat.t.1vAll

# Rat SCE for gene IDs
load("/dcl01/ajaffe/data/lab/singleCell/day_rat_snRNAseq/SCE_rat-NAc_downstream-processing_MNT.rda", verbose=T)
    # sce.nac.rat, chosen.hvgs.nac.rat

for(i in names(markers.rat.t.1vAll)){
  rownames(markers.rat.t.1vAll[[i]]) <- rowData(sce.nac.rat)$Symbol[match(rownames(markers.rat.t.1vAll[[i]]),
                                                                          rowData(sce.nac.rat)$ID)]
}


# Intersecting top 40?  Just use brute-force `toupper()`
intersect(toupper(head(rownames(markers.rat.t.1vAll[["Grm8-MSN"]]), n=40)),
          head(rownames(markers.nac.t.1vAll[["MSN.D1.3"]]), n=40))
    # "FOXP2"   "NTNG1"   "CNTN5"   "SLC35F1"
        # *More w/ MSN.D1.1: "EYA2"  "VWC2L" "FOXP2" "KCNJ3" "GNG4"  "SEZ6L" "PPM1E"
        # Three with MSN.D1.2: "THSD7B" "KCNH5"  "LYPD1"

# Plot it
sce.nac <- sce.nac[ ,sce.nac$cellType.final != "ambig.lowNtrxts"]
sce.nac$cellType.final <- droplevels(sce.nac$cellType.final)

plotExpression(sce.nac, x="cellType.final", colour_by="cellType.final",
               exprs_values="logcounts", features="GRM8") +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.3,
               colour=rep(tableau20[1:14], 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12))
    ## D1.1 expresses the most, followed by D1.2
  
which(rownames(markers.nac.t.1vAll[["MSN.D1.3"]]) == "GRM8")
    # 29084 - interesting     * for 'MSN.D1.1' it's the 157th-top marker



# What does Grm8 look like in the rat data?
rowData(sce.nac.rat[which(rowData(sce.nac.rat)$Symbol %in% c("Grm8","Drd1","Drd2","Drd3")), ])
    # "ENSRNOG00000021468", "ENSRNOG00000023688", "ENSRNOG00000008428", "ENSRNOG00000060806"
plotExpression(sce.nac.rat, x="Celltype", colour_by="Celltype",
               exprs_values="logcounts", point_alpha=0.6, point_size=0.8,
               features=c("ENSRNOG00000023688","ENSRNOG00000008428","ENSRNOG00000060806","ENSRNOG00000021468")) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.3,
               colour=rep(tableau20[1:16], 4)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  ggtitle("Drd1, Drd2, Drd3, and Grm8 expression, respectively (Savell et al. 2020 Rat NAc)")
    # (screencapped this)

# Shared genes b/tw Grm8-MSNs : 'MSN.D1.1' vs. 'MSN.D1.3' ? Taking top 100 for each
shared.w.d1.1 <- intersect(toupper(head(rownames(markers.rat.t.1vAll[["Grm8-MSN"]]), n=100)),
                           head(rownames(markers.nac.t.1vAll[["MSN.D1.1"]]), n=100))
    # 20 genes
shared.w.d1.3 <- intersect(toupper(head(rownames(markers.rat.t.1vAll[["Grm8-MSN"]]), n=100)),
                           head(rownames(markers.nac.t.1vAll[["MSN.D1.3"]]), n=100))
    # 19 genes

intersect(shared.w.d1.1, shared.w.d1.3)
    # only 9 of the above overlaps are the same


# Correlation coefficients ===
load("rdas/zTsMats_libd-NAc_and_DayLab-ratNAc_sharedGenes_29May2020.rda", verbose=T)
    # ts.nac, ts.rat, Readme
cor_t_nac <- cor(ts.nac, ts.rat)
rownames(cor_t_nac) = paste0(rownames(cor_t_nac),"_H")
colnames(cor_t_nac) = paste0(colnames(cor_t_nac),"_R")
cor_t_nac[ ,"Grm8-MSN_R"]
    #      Astro_H    Inhib.1_H    Inhib.2_H    Inhib.3_H    Inhib.4_H      Micro_H
    #  0.005444515  0.037596773  0.025707518 -0.018774007 -0.035129239 -0.083253870
    #   MSN.D1.1_H   MSN.D1.2_H   MSN.D1.3_H   MSN.D1.4_H   MSN.D2.1_H   MSN.D2.2_H
    #  0.192243866  0.079574932  0.262900534  0.234844145  0.134149319  0.059178530
    #      Oligo_H        OPC_H
    # -0.345118419  0.149465693





### Test MNT Apr2021 ======================================
  # What about trying `fastMNN` as an integration approach?
library(batchelor)

# Load human NAc .rda (preprint)
load("rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda", verbose=T)
    # sce.nac, chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac, ref.sampleInfo

# Load rat SCE
load("/dcl01/ajaffe/data/lab/singleCell/day_rat_snRNAseq/SCE_rat-NAc_downstream-processing_MNT.rda", verbose=T)
    # sce.nac.rat, chosen.hvgs.nac.rat

# Load previous object with 'pseudo-bulked' SCEs subsetted for aligned $JAX.geneID
load("/dcl01/ajaffe/data/lab/singleCell/day_rat_snRNAseq/SCE_rat-NAc-PBd_w_matchingHsap-NAc-PBd_HomoloGene.IDs_MNT.rda",
     verbose=T)
    # sce.rn.sub, sce.hsap.sub, Readme

table(rowData(sce.rn.sub)$JAX.geneID == rowData(sce.hsap.sub)$JAX.geneID)  # all TRUE - dope


## Subset/match order for those
sce.nac <- sce.nac[rowData(sce.nac)$ID %in% rowData(sce.hsap.sub)$ID, ]
sce.nac.rat <- sce.nac.rat[rowData(sce.nac.rat)$ID %in% rowData(sce.rn.sub)$ID, ]

# The rat SCE already has its $JAX.geneID; add it to 'sce.nac
rowData(sce.nac)$JAX.geneID <- rowData(sce.hsap.sub)$JAX.geneID[match(
  rowData(sce.nac)$ID, rowData(sce.hsap.sub)$ID
)]

# Match the order
sce.nac.rat <- sce.nac.rat[match(rowData(sce.nac)$JAX.geneID, rowData(sce.nac.rat)$JAX.geneID), ]

table(rowData(sce.nac)$JAX.geneID == rowData(sce.nac.rat)$JAX.geneID)
    # good

# Duplicate, in case mess anything up
sce.rat <- sce.nac.rat
sce.hsap <- sce.nac

# Clean up and get 'intersecting' colData cols
sce.hsap$Sample <- ss(sce.hsap$Sample,"/",9)
sce.hsap$cellType <- sce.hsap$cellType.final
sce.hsap$batch.a <- sce.hsap$processDate
sce.hsap$batch.b <- sce.hsap$protocol
sce.hsap$species <- "H.sap"

sce.rat$sum <- sce.rat$nCount_RNA
sce.rat$detected <- sce.rat$nFeature_RNA
sce.rat$cellType <- sce.rat$Celltype
sce.rat$batch.a <- sce.rat$Stim
sce.rat$batch.b <- sce.rat$Sex
sce.rat$species <- "R.nor"
# Get rid of that weird $Barcode duplicate (it's col 3 that's the same as colnames)
colData(sce.rat) <- colData(sce.rat)[ ,-2]

# Keep:
keepInfo <- c("Sample", "Barcode", "species", "sum", "detected", "batch.a", "batch.b", "cellType")

colData(sce.rat) <- colData(sce.rat)[ ,keepInfo]
colData(sce.hsap) <- colData(sce.hsap)[ ,keepInfo]

# Remove all reducedDims
reducedDim(sce.rat) <- NULL   # had to do this 3x to remove each of the three
reducedDim(sce.hsap) <- NULL  # 4x

# Remove logcounts
assay(sce.rat, "logcounts") <- NULL
assay(sce.hsap, "logcounts") <- NULL

table(rowData(sce.rat)$JAX.geneID == rowData(sce.hsap)$JAX.geneID)

# Make rownames the $JAX.geneID
rownames(sce.rat) <- rowData(sce.rat)$JAX.geneID
rownames(sce.hsap) <- rowData(sce.hsap)$JAX.geneID 


## Requiring that the rowData basically match exactly - have to drop other info
rowData(sce.rat) <- rowData(sce.rat)$JAX.geneID
rowData(sce.hsap) <- rowData(sce.hsap)$JAX.geneID

sce.comb <- cbind(sce.hsap, sce.rat)


## Save this for future work (along with original SCEs for their rowData)
table(rownames(sce.comb) == rowData(sce.nac)$JAX.geneID) # good

Readme <- "These SCEs are all subsetted for matching 'HomoloGene.ID'; the original spp. SCEs are saved for convenience of their gene info (rowData)"
save(sce.comb, sce.nac, sce.nac.rat, Readme,
     file="rdas/zPiloting_human-rat-combinedSCEs_for-fastMNN-integration_MNT.rda")

sce.comb
    # class: SingleCellExperiment 
    # dim: 14121 28872 
    # metadata(6): Samples Samples ... Samples Samples
    # assays(1): counts
    # rownames(14121): 34983 6980 ... 55587 40764
    # rowData names(1): value
    # colnames(28872): AAACCCACATCGAACT-1 AAACCCATCCAACCAA-1 ...
    #   TTTGGTTTCTCATAGG_4 TTTGTTGGTTCTCGTC_4
    # colData names(8): Sample Barcode ... batch.b cellType
    # reducedDimNames(0):
    # altExpNames(0):



### Re-normalizing counts; dimensionality reduction === === ===

# Distribution across spp./Samples??
sample.idx <- splitit(sce.comb$Sample)
sapply(sample.idx, function(x){quantile(sce.comb$sum[x])})
    #      Br5161_NAc Br5182_NAc_NeuN Br5207_NAc_NeuN Br5212_NAc Br5287_NAc cocaineF
    # 0%          108           131.0          105.00     101.00      102.0    509.0
    # 25%        4893         17699.5        18366.25    3025.25     7855.5   2135.0
    # 50%        6702         22951.0        22404.00    4963.50    10226.0   4493.5
    # 75%        9475         28709.0        26982.50   10292.25    14735.5   7431.5
    # 100%      69559         85460.0        84761.00   76583.00    78013.0  26180.0
    #      cocaineM  salineF  salineM
    # 0%        524   540.00   596.00
    # 25%      2501  2438.25  2373.75
    # 50%      5228  5386.00  4179.50
    # 75%      9400  9318.75  7631.75
    # 100%    38400 47393.00 37781.00


# Use `multiBatchNorm()` to compute log-normalized counts, matching the scaling across samples
#   (since it's so variable, above)
sce.comb <- multiBatchNorm(sce.comb, batch=sce.comb$Sample)

# Use the simple `modelGeneVar` - this makes more sense over `combineVar`, since the
#   cell composition is already known to be quite different (with NeuN selection)
geneVar.comb <- modelGeneVar(sce.comb)
chosen.hvgs.comb <- geneVar.comb$bio > 0
sum(chosen.hvgs.comb)
    # [1] 7087

    # Alternatively, take top 1000
    chosen.hvgs.1000 <- getTopHVGs(geneVar.comb, n=1000)
        # (instead of a Boolean, this is a character vector of n defined length)
    
    # Top 1000, fitting (design=) ~ species
    mod <- model.matrix(~ species, data=colData(sce.comb))
    geneVar.spp <- modelGeneVar(sce.comb, design=mod)
    chosen.hvgs.1000spp <- getTopHVGs(geneVar.spp, n=1000)
        length(intersect(chosen.hvgs.1000, chosen.hvgs.1000spp))
            # 860 - nah, just take top 500 lol

    chosen.hvgs.500 <- getTopHVGs(geneVar.comb, n=500)
        
      
### Run `fastMNN` (internally uses `multiBatchPCA`), taking 50 PCs ===
  # Do two iterations: R.nor merged to H.sap (merge indiv. samples first); then vice versa

# 'to H.sap' ===
set.seed(109)
mnn.hold <-  fastMNN(sce.comb, batch=sce.comb$Sample,
                     merge.order=list(list("Br5161_NAc","Br5212_NAc","Br5287_NAc",
                                           "Br5207_NAc_NeuN","Br5182_NAc_NeuN"),
                                      list("salineF","salineM","cocaineF","cocaineM")
                                      ),
                     #subset.row=chosen.hvgs.comb, d=50,
                     #subset.row=chosen.hvgs.1000, d=50,
                     subset.row=chosen.hvgs.500, d=50,
                     correct.all=TRUE, get.variance=TRUE,
                     BSPARAM=BiocSingular::IrlbaParam())
    # This temp file just used for getting batch-corrected components (drops a variety of entries)

table(colnames(mnn.hold) == colnames(sce.comb))  # all TRUE
table(mnn.hold$batch == sce.comb$Sample) # all TRUE

# Add them to the SCE, as well as the metadata (though the latter might not be so usefl)
#reducedDim(sce.comb, "PCA_corrected_2H") <- reducedDim(mnn.hold, "corrected")
#reducedDim(sce.comb, "PCA_corrected_2H.1000") <- reducedDim(mnn.hold, "corrected")
reducedDim(sce.comb, "PCA_corrected_2H.500") <- reducedDim(mnn.hold, "corrected")


    # # 'to R.nor' ===
    # set.seed(109)
    # mnn.hold <-  fastMNN(sce.comb, batch=sce.comb$species,
    #                      merge.order=c("R.nor", "H.sap"),
    #                      subset.row=chosen.hvgs.comb, d=50,
    #                      correct.all=TRUE, get.variance=TRUE,
    #                      BSPARAM=BiocSingular::IrlbaParam())
    
    # Add them to the SCE, as well as the metadata (though the latter might not be so usefl)
    #reducedDim(sce.comb, "PCA_corrected_2R") <- reducedDim(mnn.hold, "corrected") # 100 components


# # Save with these new dims
# save(sce.comb, sce.nac, sce.nac.rat, Readme, chosen.hvgs.comb,
#      file="rdas/zPiloting_human-rat-combinedSCEs_for-fastMNN-integration_MNT.rda")


    # # Are those coordinates similar?
    # sapply(c(1:50), function(x){round(cor(reducedDim(sce.comb, "PCA_corrected_2H")[ ,x],
    #                                 reducedDim(sce.comb, "PCA_corrected_2R")[ ,x]), 3)})
    #     # [1]  0.999  0.987  0.832 -0.418  0.918  0.986  0.740  0.852 -0.309  0.600
    #     # [11]  0.874  0.789  0.195  0.808 -0.067  0.548  0.768  0.738 -0.033  0.623
    #     # [21]  0.770  0.506  0.475  0.617  0.424  0.768  0.796  0.719  0.765  0.690
    #     # [31]  0.590  0.742  0.456  0.773  0.607  0.811  0.658  0.753  0.790  0.807
    #     # [41]  0.550  0.521  0.614  0.690  0.681  0.647  0.670  0.643  0.762  0.807
    #     # Interesting.

## t-SNE
reducedDim(sce.comb, "TSNE") <- NULL
set.seed(109)
#sce.comb <- runTSNE(sce.comb, dimred="PCA_corrected_2H")
#sce.comb <- runTSNE(sce.comb, dimred="PCA_corrected_2H.1000")
sce.comb <- runTSNE(sce.comb, dimred="PCA_corrected_2H.500")


# How do these look? 
#pdf("pdfs/exploration/DayLab-ratNAc/zPiloting_human-rat-fastMNN-integration_posBioComp_MNT.pdf")
#pdf("pdfs/exploration/DayLab-ratNAc/zPiloting_human-rat-fastMNN-integration_top1000hvgs_MNT.pdf")
pdf("pdfs/exploration/DayLab-ratNAc/zPiloting_human-rat-fastMNN-integration_top500hvgs_MNT.pdf")
plotReducedDim(sce.comb, dimred="TSNE", colour_by="species", point_alpha=0.1, point_size=2.0)
plotReducedDim(sce.comb, dimred="TSNE", colour_by="Sample", point_alpha=0.25, point_size=2.0)
    # Observation: these are quite still sample (esp. human)-batch-y
    #     -> Will probs want to merge those samples, first, then across species
    #reducedDim(sce.comb) <- NULL  # 3x
plotReducedDim(sce.comb, dimred="TSNE", colour_by="cellType", point_alpha=0.25, point_size=2.0,
               text_by="cellType", text_size=3)
dev.off()


# Save progress
save(sce.comb, sce.nac, sce.nac.rat, Readme,
     file="rdas/zPiloting_human-rat-combinedSCEs_for-fastMNN-integration_MNT.rda")




# ## UMAP
# set.seed(109)
# sce.comb <- runUMAP(sce.comb, dimred="PCA_corrected_2H")
# 
# plotReducedDim(sce.comb, dimred="UMAP", colour_by="species")





