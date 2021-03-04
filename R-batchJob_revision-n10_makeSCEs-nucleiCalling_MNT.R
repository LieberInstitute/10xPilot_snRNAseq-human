################################################################################
### LIBD 10x snRNA-seq [pilot] revision (n=10)
### STEP 01.batchJob: Read in SCEs and perform nuclei calling (`emptyDrops()`)
### Initiated: MNT 25Feb2021
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

### Read in (2021) 'samples.manifest' for streamlining
samples.revision <- read.table("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/FASTQ/Feb2021/samples.manifest",
                               sep="\t", header=F)$V1

# Make list of paths
paths.rawCounts <- c(paste0("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/Feb2021/",
                            samples.revision,"/outs/raw_feature_bc_matrix"))

# Simpler names for individual SCEs
names(paths.rawCounts) <- gsub("_",".", tolower(samples.revision))


## Read in raw UMI x barcode matrix - **use pre-mRNA-aligned reads
pilot.data.2 <- lapply(paths.rawCounts, function(x){ read10xCounts(x, col.names=TRUE) })
names(pilot.data.2) <- names(paths.rawCounts)



### Gene annotation (from scater) ===
# Pull in GTF information
gtf = import("/dcl01/ajaffe/data/lab/singleCell/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf")
gtf = gtf[gtf$type == "gene"]
length(gtf)
names(gtf) = gtf$gene_id
table(names(gtf) == rowData(pilot.data.2[[1]])$ID)

seqlevels(gtf)[1:25] = paste0("chr", seqlevels(gtf)[1:25])
mcols(gtf) = mcols(gtf)[,c(5:9)]


for(i in 1:length(pilot.data.2)){
  rowRanges(pilot.data.2[[i]]) <- gtf
  # Because some gene names are the same:
  rownames(pilot.data.2[[i]]) <- uniquifyFeatureNames(rowData(pilot.data.2[[i]])$gene_id, rowData(pilot.data.2[[i]])$gene_name)
  rowData(pilot.data.2[[i]])$Symbol.uniq <- rownames(pilot.data.2[[i]])
}


# Preliminary save
save(pilot.data.2, file="rdas/revision/all-FACS-n10_2021rev_SCEs_processing-QC_MNTMar2021.rda")

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

e.out.2 <- list()
for(i in 1:length(pilot.data.2)){
  cat(paste0("Simulating empty drops for: ",names(pilot.data.2)[i],"... \n"))
  
  set.seed(109)
  e.out.2[[i]] <- emptyDrops(counts(pilot.data.2[[i]]), niters=20000)
  cat(paste0("\n\t...Simulations complete. \n\t", date(), "\n\n\n"))
  date()
}

names(e.out.2) <- names(pilot.data.2)

save(pilot.data.2, e.out.2, file="rdas/revision/all-FACS-n10_2021rev_SCEs_processing-QC_MNTMar2021.rda")


rm(list=ls())
sessionInfo()
