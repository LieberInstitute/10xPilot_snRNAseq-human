################################################################################
### LIBD 10x snRNA-seq pilot (n=14) re-processing (Bioc v3.12)
### STEP 01.batchJob: Read in SCEs and perform nuclei calling (`emptyDrops()`)
### Initiated: MNT 03Mar2021
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
library(rtracklayer)


#### AS SUBMITTED JOB ====

### Read in preprint 'samples.manifest.full' for streamlining
samples.prepr <- read.table("/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-10-15_Tran2021_published/samples.manifest.full",
                            sep="\t", header=F)$V5

    # Drop Br5287-DLPFC (poor quality sample; dropped for preprint) and the test sucrose samples
    samples.prepr <- samples.prepr[-c(grep("Br5287_DLPFC", samples.prepr),
                                      grep("_suc", samples.prepr))]
    # Add '_NeuN' suffix to the two NeuN sorts
    samples.prepr[c(6,14)] <- paste0(samples.prepr[c(6,14)],"_NeuN")

# Make list of paths
paths.rawCounts <- c(paste0("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/",
                            samples.prepr,"/outs/raw_feature_bc_matrix"))

# Simpler names for individual SCEs
names(paths.rawCounts) <- gsub("_",".", tolower(samples.prepr))


## Read in raw UMI x barcode matrix - **use pre-mRNA-aligned reads
pilot.data <- lapply(paths.rawCounts, function(x){ read10xCounts(x, col.names=TRUE) })
names(pilot.data) <- names(paths.rawCounts)



### Gene annotation (from scater) ===
# Pull in GTF information
gtf = import("/dcl01/ajaffe/data/lab/singleCell/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf")
gtf = gtf[gtf$type == "gene"]
length(gtf)
names(gtf) = gtf$gene_id
table(names(gtf) == rowData(pilot.data[[1]])$ID)

seqlevels(gtf)[1:25] = paste0("chr", seqlevels(gtf)[1:25])
mcols(gtf) = mcols(gtf)[,c(5:9)]


for(i in 1:length(pilot.data)){
  rowRanges(pilot.data[[i]]) <- gtf
  # Because some gene names are the same:
  rownames(pilot.data[[i]]) <- uniquifyFeatureNames(rowData(pilot.data[[i]])$gene_id, rowData(pilot.data[[i]])$gene_name)
  rowData(pilot.data[[i]])$Symbol.uniq <- rownames(pilot.data[[i]])
}


# Preliminary save
save(pilot.data, file="rdas/revision/all-FACS-n14_preprint_SCEs_processing-QC_MNTMar2021.rda")


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

# -> In exploratory phase (preprint), not all samples passed with niters=15000 (default 10,000), so use 20,000

e.out <- list()
for(i in 1:length(pilot.data)){
  cat(paste0("Simulating empty drops for: ",names(pilot.data)[i],"... \n"))
  
  set.seed(109)
  e.out[[i]] <- emptyDrops(counts(pilot.data[[i]]), niters=20000)
  cat(paste0("\n\t...Simulations complete. \n\t", date(), "\n\n\n"))
  date()
}

names(e.out) <- names(pilot.data)

save(pilot.data, e.out, file="rdas/revision/all-FACS-n14_preprint_SCEs_processing-QC_MNTMar2021.rda")


rm(list=ls())
sessionInfo()
