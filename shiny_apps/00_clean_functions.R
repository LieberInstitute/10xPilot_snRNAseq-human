explore_sce_original <- function(sce_original) {
    print("Dimensions:")
    print(dim(sce_original))
    print("Number of unique cell names:")
    print(length(unique(colnames(sce_original))))
    print("Repeated cell names:")
    col_tab <- table(colnames(sce_original))
    print(col_tab[col_tab > 1])
    print("Number of unique genes names:")
    print(length(unique(rownames(sce_original))))
}

create_small_sce <- function(sce_original, cell_var = "cellType.final") {
    library("SingleCellExperiment")
    library("rafalib")

    sce_original <- sce_original[, tolower(colData(sce_original)[[cell_var]]) != tolower("ambig.lowNtrxts")]
    colData(sce_original)[[cell_var]] <- factor(colData(sce_original)[[cell_var]])
    colData(sce_original)[[cell_var]] <- factor(colData(sce_original)[[cell_var]])

    message(Sys.time(), " reducing the sce object")
    sce_small <- sce_original
    assays(sce_small) <- assays(sce_small)["logcounts"]
    rowData(sce_small)
    sce_small$sample <- factor(sce_small$sample)
    sce_small$region <- factor(sce_small$region)
    sce_small$donor <- factor(sce_small$donor)
    sce_small$processDate <- factor(sce_small$processDate)
    sce_small$protocol <- factor(sce_small$protocol)
    colData(sce_small) <- colData(sce_small)[, !colnames(colData(sce_small)) %in% c("Sample", "prelimCluster", "collapsedCluster", "lessCollapsed", "cellType", "prelimCluster.split", "cellType.moreSplit", "cellType.split", cell_var)]
    sce_small$cell_type <- factor(colData(sce_original)[[cell_var]])

    ## Make the rows more browsable
    colnames(sce_small) <- paste0(sce_small$sample, '_', sce_small$Barcode)
    sce_small$Barcode <- make.names(sce_small$Barcode, unique = TRUE)
    metadata(sce_small) <- list()

    ## It's all "Gene Expression", so we can remove it
    rowData(sce_small)$Type <- NULL

    message(Sys.time(), " computing propNucleiExprs")

    rowData(sce_small)$propNucleiExprs <- apply(
        assay(sce_original, "counts"),
        1,
        function(x) {
            mean(x != 0)
        }
    )
    # The above, by cell type ===
    cellType.idx <- splitit(colData(sce_original)[[cell_var]])
    rowdat.sce <- rowData(sce_small)
    for(i in names(cellType.idx)){
        message(Sys.time(), " computing propNucleiExprs for ", i)
        rowdat.sce[, paste0("propExprsIn.", i)] <- apply(
            assay(sce_original, "counts")[, cellType.idx[[i]]],
            1,
            function(x){
                mean(x != 0)
            }
        )
    }
    rowData(sce_small) <- rowdat.sce

    print(pryr::object_size(sce_original))
    print(pryr::object_size(sce_small))

    return(sce_small)
}

## Donor map info from Matt
donor_map <- paste0("donor", seq_len(8))
names(donor_map) <- c("br5161", "br5212", "br5287", "br5400", "br5276", "br5207", "br5182", "br5701")

## 2021 data version
create_small_sce_2021 <- function(sce_original) {
    library("SingleCellExperiment")
    library("rafalib")

    ## Drop some nuclei to begin with
    sce_original <- sce_original[ , !grepl("drop", sce_original$cellType)]
    sce_original$cellType <- droplevels(sce_original$cellType)
    stopifnot(all(unique(sce_original$donor) %in% names(donor_map)))
    sce_original$donor <- unname(donor_map[sce_original$donor])

    message(Sys.time(), " reducing the sce object")
    sce_small <- sce_original
    assays(sce_small) <- assays(sce_small)["logcounts"]
    sce_small$donor <- factor(sce_small$donor)
    sce_small$sex <- factor(sce_small$sex)
    sce_small$processBatch <- factor(sce_small$processBatch)
    sce_small$protocol <- factor(sce_small$protocol)
    sce_small$sequencer <- factor(sce_small$sequencer)
    sce_small$cell_type <- factor(sce_small$cellType)
    colData(sce_small) <- colData(sce_small)[, !colnames(colData(sce_small)) %in% c("Sample", "prelimCluster", "collapsedCluster", "lessCollapsed", "cellType", "prelimCluster.split", "cellType.moreSplit", "cellType.split", "sampleID", "high.mito", "region")]


    ## Make the rows more browsable
    sce_small$Barcode <- make.names(sce_small$Barcode, unique = TRUE)
    colnames(sce_small) <- paste0(sce_small$donor, '_', sce_small$Barcode)
    metadata(sce_small) <- list()

    ## Fix the "Symbol" and "ID" variables
    colnames(rowData(sce_small))[colnames(rowData(sce_small)) == "gene_id"] <- "ID"
    colnames(rowData(sce_small))[colnames(rowData(sce_small)) == "Symbol.uniq"] <- "Symbol"

    ## Drop things we don't need
    rowData(sce_small)$gene_version <- NULL
    rowData(sce_small)$gene_name <- NULL
    rowData(sce_small)$gene_source <- NULL

    ## I guess that this one could be of use
    rowData(sce_small)$gene_biotype <- factor(rowData(sce_small)$gene_biotype)

    message(Sys.time(), " computing propNucleiExprs")

    rowData(sce_small)$propNucleiExprs <- apply(
        assay(sce_original, "counts"),
        1,
        function(x) {
            mean(x != 0)
        }
    )
    # The above, by cell type ===
    cellType.idx <- splitit(sce_small$cell_type)
    rowdat.sce <- rowData(sce_small)
    for(i in names(cellType.idx)){
        message(Sys.time(), " computing propNucleiExprs for ", i)
        rowdat.sce[, paste0("propExprsIn.", i)] <- apply(
            assay(sce_original, "counts")[, cellType.idx[[i]]],
            1,
            function(x){
                mean(x != 0)
            }
        )
    }
    rowData(sce_small) <- rowdat.sce

    print(c(lobstr::obj_size(sce_original), lobstr::obj_size(sce_small)) )

    return(sce_small)
}


save_sce_small <- function(sce_small, region, prefix = "tran2020_") {
    region_dir <- paste0(prefix, region)
    dir.create(here::here("shiny_apps", region_dir), showWarnings = FALSE)
    saveRDS(sce_small, file = here("shiny_apps", region_dir, paste0("sce_", tolower(region), "_small.rds")))
}

save_cell_colors <- function(cell_colors, region, prefix = "tran2021_") {
    region_dir <- paste0(prefix, region)
    dir.create(here::here("shiny_apps", region_dir), showWarnings = FALSE)
    saveRDS(cell_colors, file = here("shiny_apps", region_dir, paste0("cell_colors_", tolower(region), ".rds")))
}

cellmarkers_fig_s7_2021 <- paste0(c('SNAP25','SLC17A6','SLC17A7','SLC17A8','GAD1','GAD2','DRD1','DRD2','AQP4','GFAP','CLDN5','FLT1',
           'CD163','SIGLEC1','C3','CD74','COL1A2','PDGFRB','MBP','PDGFRA','VCAN','SKAP1','CD247'), collapse = "\n")

create_app <- function(sce_small, region, cellmarkers = "SNAP25\nMBP\nPCP4", prefix = "tran2020_") {
    library("whisker")
    library("usethis")
    library("here")
    data <- list(
        region = region,
        regionlower = tolower(region),
        region_n = length(unique(sce_small$donor)),
        cellone = colnames(sce_small)[1],
        cellmarkers = cellmarkers
    )
    region_dir <- paste0(prefix, region)
    new_template <- whisker::whisker.render(
        usethis:::read_utf8(here::here("shiny_apps", "templates", "app.R")),
        data
    )
    writeLines(new_template, here::here("shiny_apps", region_dir, "app.R"))
    new_template <- whisker::whisker.render(
        usethis:::read_utf8(here::here("shiny_apps", "templates", "deploy.R")),
        data
    )
    writeLines(new_template, here::here("shiny_apps", region_dir, "deploy.R"))
}
