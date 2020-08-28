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
library("sva")
library("recount")
library("tidyr")
library("styler")

setDTthreads(threads = 1)

## Flags that are supplied with RScript
spec <- matrix(c(
    "cores", "c", 1, "integer", "Number of cores to use. Use a small number",
    "help", "h", 0, "logical", "Display help",
    "degradation", "d", 2, "logical", "degradation data present? T/F",
    "test", "t", 2, "logical", "Test run? T/F"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

opt$region <- "NAc"
opt$feature <- "gene"

# create the NAc_gene dir
dir.create(paste0(opt$region, "_", opt$feature), showWarnings = F)

# default arguments for flags
if ( is.null(opt$degradation ) ) { opt$degradation = FALSE }
if ( is.null(opt$test ) ) { opt$test = FALSE }
if ( is.null(opt$cores ) ) { opt$cores = 12 }

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}

# dir.create(opt$region, showWarnings = FALSE)
# dir.create(file.path(opt$region, opt$feature), showWarnings = FALSE)

load_rse <- function(feat, reg) {
    message(paste(Sys.time(), "loading expression data"))

    # expmnt data
    load("/dcl01/lieber/ajaffe/lab/Nicotine/NAc/RNAseq/paired_end_n239/count_data/NAc_Nicotine_hg38_rseGene_rawCounts_allSamples_n239.rda", verbose = T)

    rse <- rse_gene
    assays(rse)$raw_expr <- getRPKM(rse_gene, "Length")

    # genotype file, contains mds object
    load("/dcl01/lieber/ajaffe/lab/Nicotine/NAc/RNAseq/paired_end_n239/genotype_data/Nicotine_NAc_Genotypes_n206.rda", verbose=T)

    # load("/dcl01/lieber/ajaffe/Brain/Imputation/Merged/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_Quads_maf005_hwe6_geno10_updatedMap.rda", verbose = TRUE)

    table(rse$BrNum %in% rownames(mds)) # matching samples == 195
    # T = 195

    mds <- na.omit(mds) # omits nas, inits test var

    # subset rows in mds that are present as samples in rse
    mds <- mds[rownames(mds) %in% rse$BrNum,] %>%
        na.omit()

    dim(mds)

    dim(rse)

    # only keep columns in rse that are present as rows in mds
    rse <- rse[,rse$BrNum %in% rownames(mds)]

    stopifnot(identical(nrow(mds), ncol(rse)))

    mod = model.matrix(~ Sex + as.matrix(mds [,1:5]), data = colData(rse))

    pcaGene = prcomp(t(log2(assays(rse)$raw_expr + 1)))
    kGene = num.sv(log2(assays(rse)$raw_expr + 1), mod)

    genePCs = pcaGene$x[,1:kGene]

    pcs <- genePCs

    colData(rse) <- cbind(colData(rse), pcs, mds) #cbind mds[,1:]

    mod <- model.matrix(~ Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + pcs,
        data = colData(rse))

    message(paste(Sys.time(), "cleaning expression"))
    assays(rse) <- List(
        "raw_expr" = assays(rse)$raw_expr,
        "clean_expr" = cleaningY(log2(assays(rse)$raw_expr + 1), mod, P = 2)
    )

    message(paste(Sys.time(), "switch column names to BrNum"))
    stopifnot(!any(duplicated(colnames(rse))))
    colnames(rse) <- colData(rse)$BrNum

    return(rse)
}

# requires degredation data - saving for later
if (opt$degradation == TRUE) {

    library("recount")
    load_rse <- function(feat, reg) {
        message(paste(Sys.time(), "loading expression data"))
        load("/dcl01/lieber/ajaffe/lab/Nicotine/NAc/RNAseq/paired_end_n239/count_data/NAc_Nicotine_hg38_rseGene_rawCounts_allSamples_n239.rda", verbose = TRUE)
        if (feat == "gene") {
            rse <- rse_gene_joint
            assays(rse)$raw_expr <- recount::getRPKM(rse, "Length")
        } else if (feat == "exon") {
            rse <- rse_exon_joint
            assays(rse)$raw_expr <- recount::getRPKM(rse, "Length")
        } else if (feat == "jxn") {
            rse <- rse_jxn_joint
            rowRanges(rse)$Length <- 100
            assays(rse)$raw_expr <- recount::getRPKM(rse, "Length")
        } else if (feat == "tx") {
            rse <- rse_tx_joint
            assays(rse)$raw_expr <- assays(rse)$tpm
        }

        message(paste(Sys.time(), "subsetting to age and region data"))
        keepInd <- which(colData(rse)$Age > 13 & colData(rse)$Region == reg)
        rse <- rse[, keepInd]

        message(paste(Sys.time(), "loading model pieces"))
        if (reg == "HIPPO") {
            ## NOTE:
            ## As written right now, this won't ever run since HIPPO will be processed with the BSP2 files
            load("/dcl01/ajaffe/data/lab/dg_hippo/eQTL_all_SNPs_n596/rdas/pcs_4features_hippo_333.rda", verbose = TRUE)
            load("/dcl01/ajaffe/data/lab/dg_hippo/genotype_data/merged_dg_hippo_allSamples_n596.rda", verbose = TRUE)
            mds <- mds[rse$BrNum, ]
        } else if (reg == "DentateGyrus") {
            load("/dcl01/ajaffe/data/lab/dg_hippo/eQTL_all_SNPs_n596/rdas/pcs_4features_dg_263.rda", verbose = TRUE)
            load("/dcl01/ajaffe/data/lab/dg_hippo/genotype_data/astellas_dg_genotype_data_n263.rda", verbose = TRUE)
            stopifnot(identical(rse$BrNum, rownames(mds)))
        }

        message(paste(Sys.time(), "building model"))
        if (feat == "gene") {
            pcs <- genePCs
        } else if (feat == "exon") {
            pcs <- exonPCs
        } else if (feat == "jxn") {
            pcs <- jxnPCs
        } else if (feat == "tx") {
            pcs <- txPCs
        }
        colData(rse) <- cbind(colData(rse), pcs, mds)
        mod <- model.matrix(~ Dx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + pcs,
            data = colData(rse)
        )

        message(paste(Sys.time(), "cleaning expression"))
        assays(rse) <- List(
            "raw_expr" = assays(rse)$raw_expr,
            "clean_expr" = cleaningY(log2(assays(rse)$raw_expr + 1), mod, P = 2)
        )

        message(paste(Sys.time(), "switch column names to BrNum"))
        stopifnot(!any(duplicated(colnames(rse))))
        colnames(rse) <- colData(rse)$BrNum

        return(rse)
    }
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

bim_file <- "/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/twas/filter_data/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38_filtered_NAc_Nicotine"

message(paste(Sys.time(), "reading the bim file", bim_file))
bim <- fread(
    paste0(bim_file, ".bim"),
    col.names = c("chr", "snp", "position", "basepair", "allele1", "allele2")
)

# newsnp <- make.names(bim$snp, unique = T) # make snp names unique, see 188 - 190
# bim$snp <- newsnp

# convert 23 to X, as is std in plink
bim$chr <- as.character(bim$chr)
bim[chr == "23",]$chr <- "X"

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

# drops non-overlapping chromosomes
rse <- rse[!seqnames(rowRanges(rse)) %in% c("chrM", "chrY"),]

# > seqnames(rowRanges(rse))
# factor-Rle of length 56888 with 23 runs
#   Lengths:  5178  3971  2911  2505  2868 ...  2926  1372   795  1281  2332
#   Values :  chr1  chr2  chr3  chr4  chr5 ... chr19 chr20 chr21 chr22  chrX
# Levels(25): chr1 chr2 chr3 chr4 chr5 chr6 ... chr20 chr21 chr22 chrX chrY chrM <--- Important?


print("Final RSE feature dimensions:")
print(dim(rse))

rse_file_subset <- file.path("NAc_genes/subsetted_rse.Rdata")
message(paste(Sys.time(), "saving the subsetted rse file for later at", rse_file))
save(rse, file = rse_file)

## Simplify RSE file to reduce mem
assays(rse) <- list("clean_expr" = assays(rse)$clean_expr)
mcols(rowRanges(rse)) <- NULL

## Subset to features
setwd(file.path("NAc_gene"))
dir.create("snp_files", showWarnings = FALSE)
dir.create("bim_files", showWarnings = FALSE)

## For testing
if(opt$test == T) {
    rse <- head(rse, n = 20)
}

## For checking if the triplet of bim files exists
check_bim <- function(x) {
    all(file.exists(paste0(x, c(".fam", ".bed", ".bim"))))
}

## Create the bim files
i <- seq_len(nrow(rse))
i.names <- rownames(rse)

if (file.exists((file.path("i_info.Rdata")))) {
    load(file.path("i_info.Rdata"), verbose = T)
} else {
    save(i, i.names, file = "i_info.Rdata")
}

## Save the ids
write.table(data.frame(i, i.names), file = "input_ids.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

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
    dir.create(file.path("bim_files", base), showWarnings = F)
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
        "--make-bed --out", filt_bim, "--memory 8000 --threads 1"
    ))

    ## Edit the "phenotype" column of the fam file
    filt_fam <- fread(paste0(filt_bim, ".fam"),
        col.names = c("famid", "w_famid", "w_famid_fa", "w_famid_mo", "sex_code", "phenotype")
    )

    ## Note BrNums might be duplicated, hence the use of match()
    m <- match(filt_fam$famid, colData(rse)$BrNum)

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
