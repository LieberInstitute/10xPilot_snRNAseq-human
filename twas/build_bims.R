## Load plink before starting R:
# module load plink/1.90b6.6
## Also load twas fusion code
# module load fusion_twas/github

library("SummarizedExperiment")
library("jaffelab")
library("data.table")
library("sessioninfo")
library("getopt")
library("BiocParallel")
library("tidyr")
library("here")

## For styling this script
# library("styler")
# styler::style_file("build_bims.R", transformers = biocthis::bioc_style())

## Without this, the memory use blows up
## getDTthreads() will detect 64 threads in some cases here
setDTthreads(threads = 1)

## Flags that are supplied with RScript
spec <- matrix(c(
    "cores", "c", 1, "integer", "Number of cores to use. Use a small number",
    "help", "h", 0, "logical", "Display help",
    "test", "t", 2, "logical", "Test run? TRUE/FALSE"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

## For an interactive test
if(FALSE) {
    opt <- list("cores" = 1, "test" = TRUE)
}

opt$region <- "NAc"
opt$feature <- "gene"

# create the NAc_gene dir
# dir.create(opt$region, showWarnings = FALSE)
# dir.create(file.path(opt$region, opt$feature), showWarnings = FALSE)
dir.create(paste0(opt$region, "_", opt$feature), showWarnings = FALSE)

# default arguments for flags
if (is.null(opt$test)) {
    opt$test <- FALSE
}
if (is.null(opt$cores)) {
    ## Default to 1 core
    opt$cores <- 1
}

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}

## Show the options used
message(paste(Sys.time(), "options used"))
print(opt)


load_rse <- function(feat, reg) {
    message(paste(Sys.time(), "loading expression data"))
    
    # expmnt data
    load(here("twas", "filter_data", "rda", "NAc_Nicotine_hg38_rseGene_rawCounts_allSamples_n205.Rdata"), verbose = TRUE)
    ## Could be more complicated later on for exon, jxn, tx
    rse <- rse_gene
    assays(rse)$raw_expr <- assays(rse_gene)$RPKM
    
    ## Define the main model with effects to remove from the expression
    load(here("twas", "filter_data", "rda", "genePCs.Rdata"), verbose = TRUE)
    mod <- model.matrix(~ Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + genePCs, data = colData(rse)
    )
    
    ## Regress out effects. If we had a diagnosis variable (Dx), we would use it
    ## first, then use P = 2 in cleaningY()
    message(paste(Sys.time(), "cleaning expression"))
    assays(rse)$clean_expr <- cleaningY(log2(assays(rse)$raw_expr + 1), mod, P = 1)
        
    message(paste(Sys.time(), "switch column names to BrNum"))
    stopifnot(!any(duplicated(colnames(rse))))
    colnames(rse) <- colData(rse)$BrNum
    
    return(rse)
}

rse_file <- file.path("NAc_gene/working_rse.Rdata")

if (!file.exists(rse_file) == TRUE) {
    message(paste(Sys.time(), "rse file does not already exist, generating now", rse_file))
    rse <- load_rse(opt$feature, opt$region)
    message(paste(Sys.time(), "saving the rse file for later at", rse_file))
    save(rse, file = rse_file)
} else if (file.exists(rse_file) == TRUE) {
    message(paste(Sys.time(), "loading previous rse file", rse_file))
    load(rse_file, verbose = TRUE)
}
message(Sys.time(), " working RSE dimensions")
print(dim(rse))

## Using the files where the SNV names have been made unique
## using make.names(unique = TRUE)
## Details at filter_data/filter_snps.R
bim_file <- "/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/twas/filter_data/unique_snps_bim/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38_filtered_NAc_Nicotine_uniqueSNPs"

message(paste(Sys.time(), "reading the bim file", bim_file))
bim <- fread(
    paste0(bim_file, ".bim"),
    col.names = c("chr", "snp", "position", "basepair", "allele1", "allele2")
)

# convert 23 to X, as is std in plink
bim$chr <- as.character(bim$chr)
bim[chr == "23", ]$chr <- "X"

bim_gr <- GRanges(
    paste0("chr", bim$chr),
    IRanges(bim$basepair, width = 1)
)
mcols(bim_gr) <- bim[, -c("chr", "basepair")]

## Based on http://gusevlab.org/projects/fusion/#computing-your-own-functional-weights
## they restrict to 500kb on each side
rse_window <- resize(rowRanges(rse), width(rowRanges(rse)) + 500000 * 2, fix = "center")
mcols(rse_window) <- NULL

## Number of SNPs per feature window
# table(countOverlaps(rse_window, bim_gr))

## Keep only those feature windows with some SNPs nearby
keep_feat <- which(countOverlaps(rse_window, bim_gr) > 0)
message(paste(Sys.time(), "number of features kept:", length(keep_feat)))
rse <- rse[keep_feat, ]
rse_window <- rse_window[keep_feat, ]
stopifnot(nrow(rse) == length(rse_window))

message(Sys.time(), " number of genes per chr")
table(seqnames(rse))

print("Final RSE feature dimensions:")
print(dim(rse))

rse_file_subset <- file.path("NAc_gene/subsetted_rse.Rdata")

if (!file.exists(rse_file_subset)){
    message(paste(Sys.time(), "saving the subsetted rse file for later at", rse_file_subset))
    save(rse, file = rse_file_subset)
} else {
    message(paste(Sys.time(), "subsetted rse file already exists:", rse_file))
}

## Simplify RSE file to reduce mem
assays(rse) <- list("clean_expr" = assays(rse)$clean_expr)
mcols(rowRanges(rse)) <- NULL

## Subset to features
setwd(file.path("NAc_gene"))
dir.create("snp_files", showWarnings = FALSE)
dir.create("bim_files", showWarnings = FALSE)
message(Sys.time(), " current files:")
print(dir())

## For testing
if (opt$test == TRUE) {
    rse <- head(rse, n = 20)
}

## For checking if the triplet of bim files exists
check_bim <- function(x) {
    all(file.exists(paste0(x, c(".fam", ".bed", ".bim"))))
}

## For matching, due to complicated BrNums
## this function was used in filter_data/filter_snps.R
brnumerical <- function(x) {
    as.integer(gsub("Br|_.*", "", x))
}

## Create the bim files
i <- seq_len(nrow(rse))
i.names <- rownames(rse)

if (file.exists((file.path("i_info.Rdata")))) {
    message(Sys.time(), " loading pre-existing i_info.Rdata")
    load(file.path("i_info.Rdata"), verbose = TRUE)
} else {
    message(Sys.time(), " creating the i_info.Rdata file")
    save(i, i.names, file = "i_info.Rdata")
}

message(paste(Sys.time(), "length of i and i.names"))
print(length(i))
print(length(i.names))

## Save the ids
write.table(data.frame(i, i.names), file = "input_ids.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
# rm(b, q, feat_id, base, j, filt_fam, m)
# setwd("NAc_gene/")
# load("pre-bimFilesFxn.RData")
# save.image(file = "testing_bim_files.RData")

## Manual testing:
# i <- 1
# feat_id <- i.names[1]
# clean <- TRUE

bim_files <- bpmapply(function(i, feat_id, clean = TRUE) {
    if (i == 1 || i %% 1000 == 0) {
        message("*******************************************************************************")
        message(paste(Sys.time(), "building bim file for i =", i, "corresponding to feature", feat_id))
    }
    
    base <- paste0(opt$region, "_", opt$feature, "_", i)
    filt_snp <- paste0("LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38_filtered_NAc_Nicotine_", base, ".txt")
    
    # change this file
    filt_bim <- gsub(".txt", "", filt_snp)
    filt_snp <- file.path("snp_files", filt_snp)
    dir.create(file.path("bim_files", base), showWarnings = FALSE)
    filt_bim <- file.path("bim_files", base, filt_bim)
    
    ## Re-use the bim file if it exists already
    if (check_bim(filt_bim)) {
        return(TRUE)
    }
    
    j <- subjectHits(findOverlaps(rse_window[i], bim_gr))
    
    fwrite(
        bim[j, "snp"], # %>% unique()
        file = filt_snp,
        sep = "\t", col.names = FALSE
    )
    
    system(paste(
        "plink --bfile", bim_file, "--extract", filt_snp,
        "--make-bed --out", filt_bim, "--memory 5000 --threads 1"
    ))

    ## Edit the "phenotype" column of the fam file
    filt_fam <- fread(paste0(filt_bim, ".fam"),
        col.names = c("famid", "w_famid", "w_famid_fa", "w_famid_mo", "sex_code", "phenotype")
    )

    ## Note BrNums might be duplicated, hence the use of match()
    m <- match(brnumerical(filt_fam$famid), brnumerical(colData(rse)$BrNum))
    stopifnot(all(!is.na(m)))

    ## Use cleaned expression for now. Could be an argument for the code.
    filt_fam$phenotype <- assays(rse)$clean_expr[i, m]

    ## Ovewrite fam file (for the phenotype info)
    fwrite(filt_fam, file = paste0(filt_bim, ".fam"), sep = " ", col.names = FALSE)
    ## Clean up
    if (clean) unlink(filt_snp)

    return(check_bim(filt_bim))
}, i, i.names, BPPARAM = MulticoreParam(workers = opt$cores), SIMPLIFY = FALSE)

message(paste(Sys.time(), "finished creating the bim files"))
stopifnot(all(unlist(bim_files)))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

message("-- plink version information --")
system("plink --version")
