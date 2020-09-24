
library("SingleCellExperiment")
library("iSEE")
library("shiny")

sce <- readRDS("sce_nac_small.rds")

packageVersion("iSEE")


initial <- list()

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot", Type = "PCA", XAxis = 1L, YAxis = 2L,
    ColorByColumnData = "cell_type", ColorByFeatureNameAssay = "logcounts",
    ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "donor",
    SizeByColumnData = "sum", FacetByRow = "---", FacetByColumn = "---",
    ColorBy = "Column data", ColorByDefaultColor = "#000000",
    ColorByFeatureName = "MIR1302-2HG", ColorByFeatureSource = "---",
    ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "nac.5161_AAACCCACATCGAACT-1",
    ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
    ShapeBy = "Column data", SizeBy = "None", SelectionEffect = "Transparent",
    SelectionColor = "#FF0000", SelectionAlpha = 0.1, ZoomData = numeric(0),
    BrushData = list(), VisualBoxOpen = FALSE, VisualChoices = c("Color",
        "Shape"), ContourAdd = FALSE, ContourColor = "#0000FF", PointSize = 1,
    PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
    FontSize = 1, LegendPosition = "Bottom", PanelId = c(ReducedDimensionPlot = 1L),
    PanelHeight = 600L, PanelWidth = 4L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, RowSelectionType = "Active",
    RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE,
    ColumnSelectionType = "Active", ColumnSelectionSaved = 0L,
    SelectionHistory = list())

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new("RowDataTable", Selected = "MOBP", Search = "", SearchColumns = c("",
    ""), PanelId = c(RowDataTable = 1L), PanelHeight = 600L, PanelWidth = 3L,
    SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, RowSelectionType = "Active",
    RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE,
    ColumnSelectionType = "Active", ColumnSelectionSaved = 0L,
    SelectionHistory = list())

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data",
    XAxisColumnData = "cell_type", XAxisFeatureName = "MIR1302-2HG",
    XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE,
    YAxisFeatureName = "MOBP", YAxisFeatureSource = "RowDataTable1",
    YAxisFeatureDynamicSource = TRUE, ColorByColumnData = "cell_type",
    ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000",
    ShapeByColumnData = "donor", SizeByColumnData = "sum", FacetByRow = "donor",
    FacetByColumn = "---", ColorBy = "Column data", ColorByDefaultColor = "#000000",
    ColorByFeatureName = "MIR1302-2HG", ColorByFeatureSource = "---",
    ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "nac.5161_AAACCCACATCGAACT-1",
    ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None", SizeBy = "None", SelectionEffect = "Transparent",
    SelectionColor = "#FF0000", SelectionAlpha = 0.1, ZoomData = numeric(0),
    BrushData = list(), VisualBoxOpen = FALSE, VisualChoices = c("Color",
        "Facet"), ContourAdd = FALSE, ContourColor = "#0000FF", PointSize = 1,
    PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
    FontSize = 1, LegendPosition = "Bottom", PanelId = c(FeatureAssayPlot = 1L),
    PanelHeight = 600L, PanelWidth = 5L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, RowSelectionType = "Active",
    RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE,
    ColumnSelectionType = "Active", ColumnSelectionSaved = 0L,
    SelectionHistory = list())

################################################################################
# Settings for Column data plot 1
################################################################################

initial[["ColumnDataPlot1"]] <- new("ColumnDataPlot", XAxis = "Column data", YAxis = "detected",
    XAxisColumnData = "cell_type", ColorByColumnData = "cell_type",
    ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000",
    ShapeByColumnData = "high.mito", SizeByColumnData = "sum",
    FacetByRow = "donor", FacetByColumn = "---", ColorBy = "Column data",
    ColorByDefaultColor = "#000000", ColorByFeatureName = "MIR1302-2HG",
    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
    ColorBySampleName = "nac.5161_AAACCCACATCGAACT-1", ColorBySampleSource = "---",
    ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
    SelectionEffect = "Transparent", SelectionColor = "#FF0000",
    SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(),
    VisualBoxOpen = FALSE, VisualChoices = c("Color", "Facet"
    ), ContourAdd = FALSE, ContourColor = "#0000FF", PointSize = 1,
    PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
    FontSize = 1, LegendPosition = "Bottom", PanelId = c(ColumnDataPlot = 1L),
    PanelHeight = 600L, PanelWidth = 5L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, RowSelectionType = "Active",
    RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE,
    ColumnSelectionType = "Active", ColumnSelectionSaved = 0L,
    SelectionHistory = list())

################################################################################
# Settings for Sample assay plot 1
################################################################################

initial[["SampleAssayPlot1"]] <- new("SampleAssayPlot", Assay = "logcounts", XAxis = "None", XAxisRowData = "ID",
    XAxisSampleName = "nac.5161_AAACCCACATCGAACT-1", XAxisSampleSource = "---",
    XAxisSampleDynamicSource = FALSE, YAxisSampleName = "nac.5161_AAACCCACATCGAACT-1",
    YAxisSampleSource = "ColumnDataTable1", YAxisSampleDynamicSource = TRUE,
    ColorByRowData = "ID", ColorBySampleNameAssay = "logcounts",
    ColorByFeatureNameColor = "#FF0000", ShapeByRowData = NA_character_,
    SizeByRowData = "", FacetByRow = "---", FacetByColumn = "---",
    ColorBy = "None", ColorByDefaultColor = "#000000", ColorByFeatureName = "MIR1302-2HG",
    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
    ColorBySampleName = "nac.5161_AAACCCACATCGAACT-1", ColorBySampleSource = "---",
    ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
    SelectionEffect = "Transparent", SelectionColor = "#FF0000",
    SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(),
    VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE,
    ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1,
    Downsample = FALSE, DownsampleResolution = 200, FontSize = 1,
    LegendPosition = "Bottom", PanelId = c(SampleAssayPlot = 1L),
    PanelHeight = 600L, PanelWidth = 3L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, RowSelectionType = "Active",
    RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE,
    ColumnSelectionType = "Active", ColumnSelectionSaved = 0L,
    SelectionHistory = list())

################################################################################
# Settings for Column data table 1
################################################################################

initial[["ColumnDataTable1"]] <- new("ColumnDataTable", Selected = "nac.5161_AAACCCACATCGAACT-1",
    Search = "", SearchColumns = c("", "", "", "", "", "", "",
        "", "", "", "", "", "", "", "", "", ""), PanelId = c(ColumnDataTable = 1L),
    PanelHeight = 600L, PanelWidth = 4L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, RowSelectionType = "Active",
    RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE,
    ColumnSelectionType = "Active", ColumnSelectionSaved = 0L,
    SelectionHistory = list())

iSEE(sce, appTitle = "M.N. Tran et al 2020, NAc, n = 5", initial = initial)
