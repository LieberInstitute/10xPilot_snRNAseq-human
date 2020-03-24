################################################################################
### LIBD pilot (2x) NeuN-sorted NAc snRNA-seq samples
### STEP 01: QC datasets, separately (integration in step 2 script(s))
### Initiated: MNT 04Mar2020   
################################################################################

library(SingleCellExperiment)
library(scRNAseq)
library(batchelor)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(uwot)
library(DropletUtils)
library(jaffelab)
library(Rtsne)
library(gridExtra)


### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

### =======

## Read in raw UMI x barcode matrix - **use pre-mRNA-aligned reads
path.5207.nacNeuN <- file.path("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/Br5207_NAc_NeuN/outs/raw_feature_bc_matrix")
nac.neun.5207 <- read10xCounts(path.5207.nacNeuN, col.names=TRUE)

path.5182.nacNeuN <- file.path("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/Br5182_NAc_NeuN/outs/raw_feature_bc_matrix")
nac.neun.5182 <- read10xCounts(path.5182.nacNeuN, col.names=TRUE)

rm(list=ls(pattern="path."))


## Gene annotation (from scater)
pilot.neun <- list(nac.neun.5207, nac.neun.5182)
names(pilot.neun) <- c("nac.neun.5207", "nac.neun.5182")

# Save memory
rm(nac.neun.5207, nac.neun.5182)


for(i in 1:length(pilot.neun)){
  rownames(pilot.neun[[i]]) <- uniquifyFeatureNames(rowData(pilot.neun[[i]])$ID, rowData(pilot.neun[[i]])$Symbol)
}

# In case mess anything up
save(pilot.neun, file="rdas/NeuN-sortedNAc_n2_processing-QC_MNTMar2020.rda")

### Quality control ============================================================
## - Going to ignore the adaptive NMAD-approach to outlier detection for UMI/feature count
#    because this hasn't been as straightforward in past experience (might throw away neurons)
## - Vignette for the 10x PBMC dataset (OSCA Ch.24) only does mito & droplet QC anyhow
#       - (mention that for a sample with very heterogeneous cell comp., don't want
#          to drop potential cells with low RNA content)


## Cell detection (droplet exclusion, rather)
# Can use UMI count vs barcode rank (knee/inflection plot) to decide threshold, but
#      "this unnecessarily discards libraries derived from cell types with low RNA content" (OSCA, Ch. 6)
#      -> Instead should prefer this Monte Carlo-simulation-based empty droplet test:
# Additionally:
# For any Sig==FALSE & Limited==TRUE, may need to increase n iterations (default = 10000) with 'niters='
#   - this field = whether "the computed p-value for a...barcode is bounded by the number of iterations"


e.out <- list()
for(i in 1:length(pilot.neun)){
  cat(paste0("Simulating empty drops for: ",names(pilot.neun)[i],"... \n"))
  
  set.seed(109)
  e.out[[i]] <- emptyDrops(counts(pilot.neun[[i]]), niters=20000)
      # Bc in exploratory phase, not all samples passed with niters=15000 (default 10,000)
}

save(pilot.neun, e.out, file="rdas/NeuN-sortedNAc_n2_processing-QC_MNTMar2020.rda")

# As needed:
#load("rdas/NeuN-sortedNAc_n2_processing-QC_MNTMar2020.rda")

names(e.out) <- names(pilot.neun)

for(i in 1:length(e.out)){
  print(names(e.out)[[i]])
  print(table(Signif = e.out[[i]]$FDR <= 0.001, Limited = e.out[[i]]$Limited))
}
    ##[1] "nac.neun.5207"
    #      Limited
    #Signif  FALSE  TRUE
    #  FALSE 18112     0
    #  TRUE     10  4672
    
    #[1] "nac.neun.5182"
    #      Limited
    #Signif  FALSE TRUE
    #  FALSE  3552    0
    #  TRUE     23 4579              - all are good and not lower-p-value-bound-limited

# Subset in for-loop:
sapply(pilot.neun, dim)
    #nac.neun.5207 nac.neun.5182
    #[1,]         33538         33538
    #[2,]       6794880       6794880

for(i in 1:length(pilot.neun)){
  pilot.neun[[i]] <- pilot.neun[[i]][ ,which(e.out[[i]]$FDR <= 0.001)]
}
# Check
sapply(pilot.neun, dim)
    #     nac.neun.5207 nac.neun.5182
    #[1,]         33538         33538
    #[2,]          4682          4602

# Temp save 
save(pilot.neun, e.out, file="rdas/NeuN-sortedNAc_n2_processing-QC_MNTMar2020.rda")




### Mito rate QC ===
table(rownames(pilot.neun[[1]])==rownames(pilot.neun[[1]]))  # and checked various other pairs

location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(pilot.neun[[1]])$ID, 
                   column="SEQNAME", keytype="GENEID")
    ## Warning message: [\n] Unable to map 312 of 33538 requested IDs... is this an issue??
     #                  - some of them, ENSG00000284194 (SCO2), are not just scaffold/contig genes... (chr22)

# ID those mito genes
stats <- list()
for(i in 1:length(pilot.neun)){
  stats[[i]] <- perCellQCMetrics(pilot.neun[[i]], subsets=list(Mito=which(location=="MT")))
}
names(stats) <- names(pilot.neun)


## Lapply: MAD approach for mito rate thresholding
high.mito <- lapply(stats, function(x) isOutlier(x$subsets_Mito_percent, nmads=3, type="higher"))
high.mito.table <- lapply(high.mito, table)
lapply(high.mito.table, function(x) x[2]/sum(x))  # 5-7%
sapply(high.mito, attributes)
    ## $nac.neun.5207.thresholds
    #   lower    higher
    #    -Inf 0.1116235
    
    #  $nac.neun.5182.thresholds
    #   lower    higher
    #    -Inf 0.3087009


# Bind stats to colData
for(i in 1:length(pilot.neun)){
  colData(pilot.neun[[i]]) <- cbind(colData(pilot.neun[[i]]), stats[[i]])
}
for(i in 1:length(pilot.neun)){
  colData(pilot.neun[[i]]) <- cbind(colData(pilot.neun[[i]]), high.mito[[i]])
}
for(i in 1:length(pilot.neun)){
  colnames(colData(pilot.neun[[i]]))[13] <- "high.mito"
}

# $sum vs. $total ??
for(i in 1:length(pilot.neun)){
  print(table(pilot.neun[[i]]$sum == pilot.neun[[i]]$total))
}
    ## all TRUE so can remove this second column:
    for(i in 1:length(pilot.neun)){
      pilot.neun[[i]]$total <- NULL
    }


# Store original for comparison/plotting
pilot.neun.unfiltered <- pilot.neun

## Subset - remove those indexed as high.mito
for(i in 1:length(pilot.neun)){
  pilot.neun[[i]] <- pilot.neun[[i]][ ,!high.mito[[i]]]
}
sapply(pilot.neun, dim)

## Plot metrics

mitoCutoffs <- unlist(lapply(high.mito, function(x){attributes(x)$thresholds["higher"]}))
mean(mitoCutoffs)
    #[1] 0.3903217
median(mitoCutoffs)
    #[1] 0.138046
mitoCutoffs <- round(mitoCutoffs, 2)

pdf("pdfs/NeuN-sortedNAc-n2_QCmetrics_high-mitoColored_MNTMar2020.pdf", height=5)
for(i in 1:length(pilot.neun.unfiltered)){
  grid.arrange(
    plotColData(pilot.neun.unfiltered[[i]], y="sum", colour_by="high.mito") +
      scale_y_log10() + ggtitle(paste0("Total count: ", names(pilot.neun.unfiltered)[[i]])),
    plotColData(pilot.neun.unfiltered[[i]], y="detected", colour_by="high.mito") +
      scale_y_log10() + ggtitle("Detected features"),
    plotColData(pilot.neun.unfiltered[[i]], y="subsets_Mito_percent",
                colour_by="high.mito") + ggtitle(paste0("Mito % (cutoff = ", mitoCutoffs[i],")")),
    ncol=3
  )
}
dev.off()


## Save!
save(pilot.neun, pilot.neun.unfiltered, e.out, file="rdas/NeuN-sortedNAc_n2_processing-QC_MNTMar2020.rda")



    ### To add in, manually, for each region (but will put here for the pipeline) === === === ===
    ## Aside - testing add of rowRanges data to SCE's - from [ST_project_dir]/Layer_Notebook.R
    #         (see https://github.com/LieberInstitute/HumanPilot/blob/c8a3a31b991081d656ededee59da45aa0494b334/Analysis/Layer_Notebook.R#L78-L87)
    library(rtracklayer)      
    
    ## get annotation
    map = read.delim("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/151675/151675_raw_feature_bc_matrix__features.tsv.gz",
                     as.is=TRUE, header=FALSE,col.names=c("EnsemblID", "Symbol", "Type"))
    ## get GTF, this seems like what they used
    gtf = import("/dcl01/ajaffe/data/lab/singleCell/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf")
        ## of length 2565061
    gtf = gtf[gtf$type    == "gene"]
        ## of length 33538
    names(gtf) = gtf$gene_id
    table(names(gtf) == map$EnsemblID)
        #TRUE
        #33538      - (so the reordering below isn't actually necessary, but safe)
    gtf = gtf[map$EnsemblID]
    seqlevels(gtf)[1:25] = paste0("chr", seqlevels(gtf)[1:25])
    mcols(gtf) = mcols(gtf)[,c(5:9)]
    
        table(gtf$gene_biotype)
            ##      antisense       IG_C_gene IG_C_pseudogene       IG_D_gene       IG_J_gene
            #           5497              14               9              37              18
            #IG_J_pseudogene       IG_V_gene IG_V_pseudogene         lincRNA  protein_coding
            #              3             144             188            7484           19912
            #      TR_C_gene       TR_D_gene       TR_J_gene TR_J_pseudogene       TR_V_gene
            #              6               4              79               4             106
            #TR_V_pseudogene
            #             33
    
    save(gtf,file="rdas/zref_genes-GTF-fromGRCh38-3.0.0_33538.rda")


    
    # === === === === === === === === === === ===
    # And end here -> proceed to 'step02' scripts
    # === === === === === === === === === === ===




