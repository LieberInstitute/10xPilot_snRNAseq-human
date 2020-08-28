## Load plink before starting R:
# module load plink/1.90b6.6
## Also load twas fusion code
# module load fusion_twas/github

library("SummarizedExperiment")
library("sessioninfo")
library("getopt")
library("BiocParallel")

## Flags that are supplied with RScript
spec <- matrix(c(
    "cores", "c", 1, "integer", "Number of cores to use. Use a small number",
    "pgconly", "p", 1, "logical", "Subset to only PGC loci?",
    "help", "h", 0, "logical", "Display help"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

opt$region <- "NAc"
opt$feat <- "gene"

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}
## Not sure what if(FALSE) really does
if (FALSE) {
    # feat = opt$feature; reg = opt$reg
    opt <- list(region = "NAc", feature = "gene", cores = 3, "pgconly" = FALSE)
}

stopifnot(opt$region %in% c("NAc"))
stopifnot(opt$feature %in% c("gene"))

## Use the rse file from build_bims.R
rse_file <- file.path(opt$region, paste0(opt$feature, ifelse(opt$pgconly, "_pgconly", "")), "subsetted_rse.Rdata") # can't make .Rdata
stopifnot(file.exists(rse_file))


## Change to the wd
setwd(file.path(opt$region, paste0(opt$feature, ifelse(opt$pgconly, "_pgconly", ""))))


## Create the bim files with half the cores used for the later stage
## to reduce memory blow up
load("i_info.Rdata", verbose = TRUE)

## Compute TWAS weights in parallel
output_status <- bpmapply(function(i, feat_id, clean_bim = FALSE) {

    ## Check that the output file exists
    out_file <- file.path("out_files", paste0(opt$feature, "_", i))

    ## Clean up
    if (clean_bim) unlink(dirname(filt_bim), recursive = TRUE)

    return(file.exists(paste0(out_file, ".wgt.RDat")))
}, i, i.names, BPPARAM = MulticoreParam(workers = opt$cores), SIMPLIFY = FALSE)
output_status <- unlist(output_status)

message("*******************************************************************************")
message(paste(Sys.time(), "summary output status (TRUE means that there is a file)"))
table(output_status)

## Load RSE
message(paste(Sys.time(), "loading the rse file"))
load(basename(rse_file), verbose = TRUE)

message(paste(Sys.time(), "saving the output_status.Rdata file"))
names(output_status) <- paste0(opt$feature, "_", i)
save(output_status, file = "output_status.Rdata")

message(paste(Sys.time(), "creating the wgt profile files"))
rdat_files <- dir("out_files", ".wgt.RDat", full.names = TRUE)
stopifnot(length(rdat_files) == sum(output_status))

wglist <- paste0(opt$region, "_", opt$feature, ".list")
write.table(rdat_files, file = wglist, row.names = FALSE, col.names = FALSE, quote = FALSE)
system(paste0("Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/utils/FUSION.profile_wgt.R ", wglist, " > ", opt$region, "_", opt$feature, ".profile.err 2>&1"))
## Rename the wglist_summary.txt file to keep the naming convention consistent
system(paste0("mv wglist_summary.txt ", opt$region, "_", opt$feature, ".profile"))

message(paste(Sys.time(), "creating the .pos file"))
if (opt$feature %in% c("gene", "exon")) {
    var <- "gencodeID"
    # vars <- 'Symbol'
} else if (opt$feature == "jxn") {
    var <- "gencodeGeneID"
    # vars <- 'Symbol'
} else {
    var <- "gene_id"
    # vars <- 'gene_name'
}

pos_match <- match(gsub("out_files/|\\.wgt\\.RDat", "", rdat_files), names(output_status))
stopifnot(all(!is.na(pos_match)))

stopifnot(identical(gsub("out_files/|\\.wgt\\.RDat", "", rdat_files), names(output_status)[pos_match]))


pos_info <- data.frame(
    "WGT" = rdat_files,
    ## use the feature id
    "ID" = names(rowRanges(rse))[pos_match],
    "CHR" = gsub("chr", "", seqnames(rowRanges(rse)[pos_match])),
    "P0" = start(rowRanges(rse)[pos_match]),
    "P1" = end(rowRanges(rse)[pos_match]),
    "geneID" = mcols(rowRanges(rse))[, var][pos_match],
    stringsAsFactors = FALSE
)
sapply(pos_info, function(x) sum(is.na(x)))
pos_info$ID[pos_info$ID == ""] <- NA
pos_info$geneID[pos_info$geneID == ""] <- NA
write.table(pos_info[, -which(colnames(pos_info) == "geneID")],
    file = paste0(opt$region, "_", opt$feature, ".pos"),
    row.names = FALSE, col.names = TRUE, quote = FALSE
)
save(pos_info, file = "pos_info.Rdata")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

message("-- plink version information --")
system("plink --version")
