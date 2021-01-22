## Based on:
# https://github.com/LieberInstitute/brainseq_phase2/blob/master/twas/read_twas.R
# https://github.com/LieberInstitute/brainseq_phase2/blob/master/twas/explore_twas.R
# https://github.com/LieberInstitute/brainseq_phase2/blob/master/twas/explore_twas_psycm.R


library("readr")
library("purrr")
library("dplyr")
library("stringr")
library("sessioninfo")
library("SummarizedExperiment")
library("ggplot2")
library("gplots")
library("VennDiagram")
library("RColorBrewer")
library("readr")
library("here")
library("getopt")
library("GenomicRanges")
library("dplyr")
library("data.table")

spec <- matrix(
    c('region', 'r', 1, 'character', 'Either Amygdala or SACC'),
    byrow = TRUE,
    ncol = 5
)

opt <- getopt(spec)

# Redundant to do this but just keeping the flag options
# in here in case it makes a difference which brain
# subregion we source our rse from
opt$region <- "NAc"

opt$region <- tolower(opt$region)

## Main options to work through
regions <- c("NAc")
types <- c("aoi", "si", "cpd", "dpw", "sc")
features <- "gene"
file_types <- c("all", "included", "dropped")

## Function for locating the files
locate_files <- function(ftype, path, type) {
    ## Determine the file pattern to use
    fpatt <- case_when(
        ftype == "all" ~ "",
        ftype == "included" ~ "\\.analysis\\.joint_included",
        ftype == "dropped" ~ "\\.analysis\\.joint_dropped"
    )
    patt <- paste0(type,
        "\\.[[:digit:]]*",
        fpatt,
        "\\.dat$")

    ## Find the files
    dat_files <- dir(path = path,
        pattern = patt,
        full.names = TRUE)

    ## Extract the chromosome
    names(dat_files) <- str_extract(basename(dat_files),
        "[:digit:]+")

    ## Done
    return(dat_files)
}


## Now read in the data for all file types
## Note that each file type has different numbers of columns
## when why I'm keeping them as different elements of the
## 'twas' list
twas <- map(file_types, function(ftype) {
    ## data.frame with all the combinations of arguments
    arg_grid <- expand.grid(
        region = regions,
        type = types,
        feature = features,
        stringsAsFactors = FALSE
    )


    pmap_dfr(arg_grid, function(region, type, feature) {
        ## Construct the path to the files
        path <- file.path(
            "/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/twas",
            paste0(region, "_", feature),
            type
        )
        
        ## Now locate the files
        dat_files <- locate_files(ftype = ftype,
            path = path,
            type = type)

        ## Next read the files and add the chromosome info
        result <- map2_dfr(dat_files,
            names(dat_files),
            function(f, chr) {
                res <- suppressMessages(read_tsv(f))
                res$chr <- chr
                return(res)
            })

        ## Next add the region, feature and type information
        result$region <- region
        result$feature <- feature
        result$type <- type

        ## Done
        return(result)
    })
})
names(twas) <- file_types

## Explore the resulting dimensions
map_dfr(twas, dim)
# # A tibble: 2 x 3
#     all included dropped
#   <int>    <int>   <int>
# 1 31609       79   16136
# 2    24       12      12

## Save the data for later use
# dir.create("rda", showWarnings = FALSE)
# save(twas, file = "rda/twas.Rdata")

## Read in the RSE info
## adapted from https://github.com/LieberInstitute/brainseq_phase2/blob/master/development/load_funs.R
load_rse <- function(type) {
    # expmnt data
    load_file <-
        Sys.glob(here(
            "twas",
            "filter_data",
            paste0(opt$region, "_rda"),
            "NAc_Nicotine_hg38_rseGene_rawCounts_allSamples_n205.Rdata")
        )

    stopifnot(file.exists(load_file))
    load(load_file)

    ## Get the appropriate object
    rse <- rse_gene

    return(rse)
}

rse <- map(unique(twas$all$feature), load_rse)
names(rse) <- unique(twas$all$feature)

## Add the gene id
twas_exp <- map(twas, function(tw) {
    ## For testing the inner part of the function
    # tw <- twas$included

    by_feat <- split(tw, tw$feature)
    ## Make sure it's in the right order
    if (!identical(names(by_feat), names(rse))) {
        message(paste(Sys.time(), "fixing order"))
        by_feat <- by_feat[names(rse)]
    }

    ## Now add the gene gencode ID and symbol
    result <- pmap_dfr(list(by_feat, rse, names(rse)),
        function(info, rs, feature) {
            ## Find the appropriate variables
            gene_var <- case_when(
                feature == "gene" ~ "gencodeID",
                feature == "exon" ~ "gencodeID",
                feature == "jxn" ~ "newGeneID",
                feature == "tx" ~ "gene_id"
            )
            symbol_var <- case_when(
                feature == "gene" ~ "Symbol",
                feature == "exon" ~ "Symbol",
                feature == "jxn" ~ "newGeneSymbol",
                feature == "tx" ~ "gene_name"
            )

            ## Match by id
            m <- match(info$ID, names(rowRanges(rs)))
            stopifnot(!is.na(m))

            ## Add the gene id/symbol
            info$geneid <- mcols(rowRanges(rs))[[gene_var]][m]
            info$genesymbol <- mcols(rowRanges(rs))[[symbol_var]][m]

            ## Done
            return(info)
        })

    return(result)
})
names(twas_exp) <- names(twas)

# get rse ranges
rse_ranges <- rowRanges(rse$gene) %>%
    ranges() %>%
    as.data.table()

# Set std. col names, specifically ID
colnames(rse_ranges) <- c("start", "end", "width", "ID")

# convert to dt
twas_exp_dt <- as.data.table(twas_exp$all)

# get ranges from rse
twas_exp_ranges <- merge(twas_exp_dt, rse_ranges, by = "ID")

# find midpoint for each gene
twas_mean_dist <- round(twas_exp_ranges[, c(width / 2) + start])

# add midpoint col
twas_exp_fin <- cbind(twas_exp_ranges, twas_mean_dist)

## Save the data for later use
dir.create("rda", showWarnings = FALSE)
save(twas_exp_fin, file = "rda/twas_exp_ranges.Rdata")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ??? Session info ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#  setting  value
#  version  R version 4.0.2 Patched (2020-06-24 r78746)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2020-10-26
#
# ??? Packages ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#  package              * version  date       lib source
#  assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.0.0)
#  backports              1.1.10   2020-09-15 [2] CRAN (R 4.0.2)
#  Biobase              * 2.48.0   2020-04-27 [2] Bioconductor
#  BiocGenerics         * 0.34.0   2020-04-27 [2] Bioconductor
#  bitops                 1.0-6    2013-08-17 [2] CRAN (R 4.0.0)
#  caTools                1.18.0   2020-01-17 [2] CRAN (R 4.0.0)
#  cli                    2.1.0    2020-10-12 [2] CRAN (R 4.0.2)
#  colorout             * 1.2-2    2020-09-24 [1] Github (jalvesaq/colorout@726d681)
#  colorspace             1.4-1    2019-03-18 [2] CRAN (R 4.0.0)
#  crayon                 1.3.4    2017-09-16 [2] CRAN (R 4.0.0)
#  data.table           * 1.13.2   2020-10-19 [2] CRAN (R 4.0.2)
#  DelayedArray         * 0.14.1   2020-07-14 [2] Bioconductor
#  dplyr                * 1.0.2    2020-08-18 [2] CRAN (R 4.0.2)
#  ellipsis               0.3.1    2020-05-15 [2] CRAN (R 4.0.0)
#  fansi                  0.4.1    2020-01-08 [2] CRAN (R 4.0.0)
#  formatR                1.7      2019-06-11 [2] CRAN (R 4.0.0)
#  futile.logger        * 1.4.3    2016-07-10 [2] CRAN (R 4.0.0)
#  futile.options         1.0.1    2018-04-20 [2] CRAN (R 4.0.0)
#  generics               0.0.2    2018-11-29 [2] CRAN (R 4.0.0)
#  GenomeInfoDb         * 1.24.2   2020-06-15 [2] Bioconductor
#  GenomeInfoDbData       1.2.3    2020-05-18 [2] Bioconductor
#  GenomicRanges        * 1.40.0   2020-04-27 [2] Bioconductor
#  getopt               * 1.20.3   2019-03-22 [2] CRAN (R 4.0.0)
#  ggplot2              * 3.3.2    2020-06-19 [1] CRAN (R 4.0.2)
#  glue                   1.4.2    2020-08-27 [2] CRAN (R 4.0.2)
#  gplots               * 3.1.0    2020-09-18 [2] CRAN (R 4.0.2)
#  gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.0.0)
#  gtools                 3.8.2    2020-03-31 [2] CRAN (R 4.0.0)
#  here                 * 0.1      2017-05-28 [1] CRAN (R 4.0.2)
#  hms                    0.5.3    2020-01-08 [2] CRAN (R 4.0.0)
#  IRanges              * 2.22.2   2020-05-21 [2] Bioconductor
#  KernSmooth             2.23-17  2020-04-26 [3] CRAN (R 4.0.2)
#  lambda.r               1.2.4    2019-09-18 [2] CRAN (R 4.0.0)
#  lattice                0.20-41  2020-04-02 [3] CRAN (R 4.0.2)
#  lifecycle              0.2.0    2020-03-06 [2] CRAN (R 4.0.0)
#  magrittr               1.5      2014-11-22 [2] CRAN (R 4.0.0)
#  Matrix                 1.2-18   2019-11-27 [3] CRAN (R 4.0.2)
#  matrixStats          * 0.57.0   2020-09-25 [2] CRAN (R 4.0.2)
#  munsell                0.5.0    2018-06-12 [2] CRAN (R 4.0.0)
#  pillar                 1.4.6    2020-07-10 [2] CRAN (R 4.0.2)
#  pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.0.0)
#  ps                     1.4.0    2020-10-07 [2] CRAN (R 4.0.2)
#  purrr                * 0.3.4    2020-04-17 [2] CRAN (R 4.0.0)
#  R6                     2.4.1    2019-11-12 [2] CRAN (R 4.0.0)
#  RColorBrewer         * 1.1-2    2014-12-07 [1] CRAN (R 4.0.2)
#  RCurl                  1.98-1.2 2020-04-18 [2] CRAN (R 4.0.0)
#  readr                * 1.4.0    2020-10-05 [2] CRAN (R 4.0.2)
#  rlang                  0.4.8    2020-10-08 [1] CRAN (R 4.0.2)
#  rprojroot              1.3-2    2018-01-03 [2] CRAN (R 4.0.0)
#  rstudioapi             0.11     2020-02-07 [2] CRAN (R 4.0.0)
#  S4Vectors            * 0.26.1   2020-05-16 [2] Bioconductor
#  scales                 1.1.1    2020-05-11 [2] CRAN (R 4.0.0)
#  sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 4.0.0)
#  stringi                1.5.3    2020-09-09 [2] CRAN (R 4.0.2)
#  stringr              * 1.4.0    2019-02-10 [1] CRAN (R 4.0.2)
#  SummarizedExperiment * 1.18.2   2020-07-09 [2] Bioconductor
#  tibble                 3.0.4    2020-10-12 [2] CRAN (R 4.0.2)
#  tidyselect             1.1.0    2020-05-11 [2] CRAN (R 4.0.0)
#  utf8                   1.1.4    2018-05-24 [2] CRAN (R 4.0.0)
#  vctrs                  0.3.4    2020-08-29 [2] CRAN (R 4.0.2)
#  VennDiagram          * 1.6.20   2018-03-28 [2] CRAN (R 4.0.0)
#  withr                  2.3.0    2020-09-22 [2] CRAN (R 4.0.2)
#  XVector                0.28.0   2020-04-27 [2] Bioconductor
#  zlibbioc               1.34.0   2020-04-27 [2] Bioconductor
#
# [1] /users/aseyedia/R/4.0
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0/R/4.0/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0/R/4.0/lib64/R/library
