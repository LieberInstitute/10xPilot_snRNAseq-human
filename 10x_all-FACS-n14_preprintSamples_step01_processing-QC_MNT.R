################################################################################
### LIBD pilot 10x-Frankenstein (n=12) snRNA-seq samples
### STEP 01: QC datasets, separately, and integrate, for PAN-BRAIN analyses
### Initiated: MNT 29Jan2020   
### Intention: To generate/have a streamlined, easy-to-follow pipeline
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
path.5161amy <- file.path("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/Br5161_Amy/outs/raw_feature_bc_matrix")
amy.5161 <- read10xCounts(path.5161amy, col.names=TRUE)

path.5212amy <- file.path("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/Br5212_Amy/outs/raw_feature_bc_matrix")
amy.5212 <- read10xCounts(path.5212amy, col.names=TRUE)

path.5161dlpfc <- file.path("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/Br5161_DLPFC/outs/raw_feature_bc_matrix")
dlpfc.5161 <- read10xCounts(path.5161dlpfc, col.names=TRUE)

path.5212dlpfc <- file.path("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/Br5212_DLPFC/outs/raw_feature_bc_matrix")
dlpfc.5212 <- read10xCounts(path.5212dlpfc, col.names=TRUE)

    #path.dlpfcFACS <- file.path("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/Br5287_DLPFC/outs/raw_feature_bc_matrix")
    #dlpfc.5287 <- read10xCounts(path.dlpfcFACS, col.names=TRUE)

path.hpcFACS <- file.path("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/Br5161_HPC/outs/raw_feature_bc_matrix")
hpc.5161 <- read10xCounts(path.hpcFACS, col.names=TRUE)

path.5212hpc <- file.path("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/Br5212_HPC/outs/raw_feature_bc_matrix")
hpc.5212 <- read10xCounts(path.5212hpc, col.names=TRUE)

path.5287hpc <- file.path("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/Br5287_HPC/outs/raw_feature_bc_matrix")
hpc.5287 <- read10xCounts(path.5287hpc, col.names=TRUE)
##
path.5161nac <- file.path("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/Br5161_NAc/outs/raw_feature_bc_matrix")
nac.5161 <- read10xCounts(path.5161nac, col.names=TRUE)

path.5212nac <- file.path("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/Br5212_NAc/outs/raw_feature_bc_matrix")
nac.5212 <- read10xCounts(path.5212nac, col.names=TRUE)

path.5287nac <- file.path("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/Br5287_NAc/outs/raw_feature_bc_matrix")
nac.5287 <- read10xCounts(path.5287nac, col.names=TRUE)

path.5161sacc <- file.path("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/Br5161_sACC/outs/raw_feature_bc_matrix")
sacc.5161 <- read10xCounts(path.5161sacc, col.names=TRUE)

path.5212sacc <- file.path("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/Br5212_sACC/outs/raw_feature_bc_matrix")
sacc.5212 <- read10xCounts(path.5212sacc, col.names=TRUE)

rm(list=ls(pattern="path."))


## Gene annotation (from scater)
pilot.data <- list(amy.5161, amy.5212,
                   dlpfc.5161, dlpfc.5212, #dlpfc.5287,
                   hpc.5161, hpc.5212, hpc.5287,
                   nac.5161, nac.5212, nac.5287,
                   sacc.5161, sacc.5212)

# Save memory
rm(amy.5161, amy.5212,
   dlpfc.5161, dlpfc.5212, #dlpfc.5287,
   hpc.5161, hpc.5212, hpc.5287,
   nac.5161, nac.5212, nac.5287,
   sacc.5161, sacc.5212)

names(pilot.data) <- c("amy.5161", "amy.5212",
                       "dlpfc.5161", "dlpfc.5212", #"dlpfc.5287",
                       "hpc.5161", "hpc.5212", "hpc.5287",
                       "nac.5161", "nac.5212", "nac.5287",
                       "sacc.5161", "sacc.5212")

#lapply(pilot.data,
#       function(x) {rownames(x) <- uniquifyFeatureNames(rowData(x)$ID, rowData(x)$Symbol)}
#  )
## didn't work - just printed _all_ the gene symbols...

for(i in 1:length(pilot.data)){
  rownames(pilot.data[[i]]) <- uniquifyFeatureNames(rowData(pilot.data[[i]])$ID, rowData(pilot.data[[i]])$Symbol)
}

# In case mess anything up
#dir.create("rdas/")
save(pilot.data, file="rdas/all-FACS-homogenates_n12_processing-QC_MNTJan2020.rda")

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
for(i in 1:length(pilot.data)){
  cat(paste0("Simulating empty drops for: ",names(pilot.data)[i],"... \n"))
  
  set.seed(109)
  e.out[[i]] <- emptyDrops(counts(pilot.data[[i]]), niters=20000)
      # Bc in exploratory phase, not all samples passed with niters=15000 (default 10,000)
}

save(pilot.data, e.out, file="rdas/all-FACS-homogenates_n12_processing-QC_MNTJan2020.rda")

# As needed:
#load("rdas/all-FACS-homogenates_n12_processing-QC_MNTJan2020.rda")

names(e.out) <- names(pilot.data)

for(i in 1:length(e.out)){
  print(names(e.out)[[i]])
  print(table(Signif = e.out[[i]]$FDR <= 0.001, Limited = e.out[[i]]$Limited))
}
    ##[1] "amy.5161"
    #Limited
    #Signif  FALSE TRUE
    #FALSE  1020    0
    #TRUE    149 3537
    
    #[1] "amy.5212"
    #Limited
    #Signif  FALSE  TRUE
    #FALSE 16780     0
    #TRUE     30  3631
    
    #[1] "dlpfc.5161"
    #Limited
    #Signif  FALSE TRUE
    #FALSE  1229    0
    #TRUE     75 4713
    
    #[1] "dlpfc.5212"
    #Limited
    #Signif  FALSE TRUE
    #FALSE   755    0
    #TRUE     35 1971
    
    #[1] "hpc.5161"
    #Limited
    #Signif  FALSE TRUE
    #FALSE   735    0
    #TRUE     91 4951
    
    #[1] "hpc.5212"
    #Limited
    #Signif  FALSE TRUE
    #FALSE   358    0
    #TRUE     59 4479
    
    #[1] "hpc.5287"
    #Limited
    #Signif  FALSE TRUE
    #FALSE   325    0
    #TRUE     35 2154
    
    #[1] "nac.5161"
    #Limited
    #Signif  FALSE TRUE
    #FALSE    41    0
    #TRUE     14 2280
    
    #[1] "nac.5212"
    #Limited
    #Signif  FALSE TRUE
    #FALSE   456    0
    #TRUE     37 1864
    
    #[1] "nac.5287"
    #Limited
    #Signif  FALSE TRUE
    #FALSE   116    0
    #TRUE     18  752
    
    #[1] "sacc.5161"
    #Limited
    #Signif  FALSE TRUE
    #FALSE  3762    0
    #TRUE     93 3454
    
    #[1] "sacc.5212"
    #Limited
    #Signif  FALSE  TRUE
    #FALSE 55255     0
    #TRUE      0  4398              - all are good and not lower-p-value-bound-limited

# Subset in for-loop:
for(i in 1:length(pilot.data)){
  pilot.data[[i]] <- pilot.data[[i]][ ,which(e.out[[i]]$FDR <= 0.001)]
}
# Check
sapply(pilot.data, dim)


# Save 
save(pilot.data, e.out, file="rdas/all-FACS-homogenates_n12_processing-QC_MNTJan2020.rda")



### Mito rate QC ===
table(rownames(pilot.data[[1]])==rownames(pilot.data[[6]]))  # and checked various other pairs

location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(pilot.data[[1]])$ID, 
                   column="SEQNAME", keytype="GENEID")

# ID those mito genes
stats <- list()
for(i in 1:length(pilot.data)){
  stats[[i]] <- perCellQCMetrics(pilot.data[[i]], subsets=list(Mito=which(location=="MT")))
}
names(stats) <- names(pilot.data)


## Lapply: MAD approach for mito rate thresholding
high.mito <- lapply(stats, function(x) isOutlier(x$subsets_Mito_percent, nmads=3, type="higher"))
high.mito.table <- lapply(high.mito, table)
lapply(high.mito.table, function(x) x[2]/sum(x))  # about 10% across all (besides weird 5212dlpfc, ~50%)
sapply(high.mito, attributes)



# Bind stats to colData
for(i in 1:length(pilot.data)){
  colData(pilot.data[[i]]) <- cbind(colData(pilot.data[[i]]), stats[[i]])
}
for(i in 1:length(pilot.data)){
  colData(pilot.data[[i]]) <- cbind(colData(pilot.data[[i]]), high.mito[[i]])
}
for(i in 1:length(pilot.data)){
  colnames(colData(pilot.data[[i]]))[13] <- "high.mito"
}

# $sum vs. $total ??
for(i in 1:length(pilot.data)){
  print(table(pilot.data[[i]]$sum == pilot.data[[i]]$total))
}
    ## all TRUE so can remove this second column:
    for(i in 1:length(pilot.data)){
      pilot.data[[i]]$total <- NULL
    }


# Store original for comparison/plotting
pilot.data.unfiltered <- pilot.data

## Subset - remove those indexed as high.mito
for(i in 1:length(pilot.data)){
  pilot.data[[i]] <- pilot.data[[i]][ ,!high.mito[[i]]]
}
sapply(pilot.data, dim)

## Plot metrics

mitoCutoffs <- unlist(lapply(high.mito, function(x){attributes(x)$thresholds["higher"]}))
mean(mitoCutoffs)
    #[1] 0.3903217
median(mitoCutoffs)
    #[1] 0.138046
mitoCutoffs <- round(mitoCutoffs, 2)

#dir.create("pdfs/")
pdf("pdfs/all-FACS-homogenates_n12_QCmetrics_high-mitoColored_MNTJan2020.pdf", height=5)
for(i in 1:length(pilot.data.unfiltered)){
  grid.arrange(
    plotColData(pilot.data.unfiltered[[i]], y="sum", colour_by="high.mito") +
      scale_y_log10() + ggtitle(paste0("Total count: ", names(pilot.data.unfiltered)[[i]])),
    plotColData(pilot.data.unfiltered[[i]], y="detected", colour_by="high.mito") +
      scale_y_log10() + ggtitle("Detected features"),
    plotColData(pilot.data.unfiltered[[i]], y="subsets_Mito_percent",
                colour_by="high.mito") + ggtitle(paste0("Mito % (cutoff = ", mitoCutoffs[i],")")),
    ncol=3
  )
}
dev.off()


## Save!
save(pilot.data, pilot.data.unfiltered, e.out, file="rdas/all-FACS-homogenates_n12_processing-QC_MNTJan2020.rda")



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


### DISREGARD CHUNK: Merging all n=12 samples into one SCE =================================
  # ** Because not using multiBatchNorm() anymore -- too conservative for scaling down counts
load("rdas/all-FACS-homogenates_n12_processing-QC_MNTJan2020.rda", verbose=T)

# First normalize across a scaled library size factor
pilot.data <- lapply(pilot.data, function(x) logNormCounts(x))

    # Re-save, so that can easily pull from, for region-specific analyses
    save(pilot.data, pilot.data.unfiltered, e.out, file="rdas/all-FACS-homogenates_n12_processing-QC_MNTJan2020.rda")

# Model gene var across each sample
geneVar <- lapply(pilot.data, function(x) modelGeneVar(x))

# Combine var to identify hvgs
combined.var <- combineVar(geneVar[[1]], geneVar[[2]], geneVar[[3]], geneVar[[4]],
                           geneVar[[5]], geneVar[[6]], geneVar[[7]], geneVar[[8]],
                           geneVar[[9]], geneVar[[10]], geneVar[[11]], geneVar[[12]])

chosen.hvgs.n12 <- combined.var$bio > 0
sum(chosen.hvgs.n12)
    ## [1] 9225

## Rescale to account for differences in sequencing depth
rescaled.pilot <- multiBatchNorm(pilot.data[[1]], pilot.data[[2]], pilot.data[[3]], pilot.data[[4]],
                                 pilot.data[[5]], pilot.data[[6]], pilot.data[[7]], pilot.data[[8]],
                                 pilot.data[[9]], pilot.data[[10]], pilot.data[[11]], pilot.data[[12]])
names(rescaled.pilot) <- names(pilot.data)


# Add region.donor info into $sample variable
for(i in 1:length(rescaled.pilot)){
  rescaled.pilot[[i]]$sample <- names(rescaled.pilot)[i]
}

# Cbind to make single SCE
sce.pilot.n12 <- cbind(rescaled.pilot[[1]], rescaled.pilot[[2]], rescaled.pilot[[3]], rescaled.pilot[[4]],
                       rescaled.pilot[[5]], rescaled.pilot[[6]], rescaled.pilot[[7]], rescaled.pilot[[8]],
                       rescaled.pilot[[9]], rescaled.pilot[[10]], rescaled.pilot[[11]], rescaled.pilot[[12]])


colnames(colData(sce.pilot.n12))
table(sce.pilot.n12$sample)

# Save
save(sce.pilot.n12, chosen.hvgs.n12, file="rdas/all-FACS-homogenates_n12_cleaned-combined_SCE_MNTJan2020.rda")

# Wrap up this portion
rm(list=ls())

### Dimensionality reduction ================================================================
load("rdas/all-FACS-homogenates_n12_cleaned-combined_SCE_MNTJan2020.rda", verbose=T)

# Starting afresh: run PCA, taking top 250 (instead of default 50 PCs)
set.seed(109)
sce.pilot.n12 <- runPCA(sce.pilot.n12, subset_row=chosen.hvgs.n12, ncomponents=250,
                        BSPARAM=BiocSingular::RandomParam())

# Save into a new data file, which will dedicate for pan-brain-analyses
save(sce.pilot.n12, chosen.hvgs.n12, file="rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTJan2020.rda")






    ## Identify optimal working PC space with 'getClustered - will do interactively.





