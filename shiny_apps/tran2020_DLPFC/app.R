
library("SingleCellExperiment")
library("iSEE")
library("shiny")

sce <- readRDS("sce_dlpfc_small.rds")
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
    ColorByFeatureName = "MOBP", ColorByFeatureSource = "---",
    ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "dlpfc.5161_AAACCCACACCGTCGA-1",
    ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None", SizeBy = "None", SelectionEffect = "Transparent",
    SelectionColor = "#FF0000", SelectionAlpha = 0.1, ZoomData = numeric(0),
    BrushData = list(), VisualBoxOpen = FALSE, VisualChoices = c("Color",
        "Shape"), ContourAdd = FALSE, ContourColor = "#0000FF", PointSize = 1,
    PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
    FontSize = 1, LegendPosition = "Bottom", PanelId = c(ReducedDimensionPlot = 1L),
    PanelHeight = 600L, PanelWidth = 6L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, RowSelectionType = "Active",
    RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE,
    ColumnSelectionType = "Active", ColumnSelectionSaved = 0L,
    SelectionHistory = list())

################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot", Assay = "logcounts", CustomRows = TRUE,
    CustomRowsText = "SNAP25
MBP
PCP4", ClusterRows = FALSE,
    ClusterRowsDistance = "spearman", ClusterRowsMethod = "ward.D2",
    DataBoxOpen = FALSE, VisualChoices = "Annotations",
    ColumnData = c("cell_type", "donor"),
    RowData = character(0), CustomBounds = FALSE, LowerBound = NA_real_,
    UpperBound = NA_real_, AssayCenterRows = FALSE, AssayScaleRows = FALSE,
    DivergentColormap = "purple < black < yellow", ShowDimNames = "Rows",
    LegendPosition = "Bottom", LegendDirection = "Horizontal",
    VisualBoxOpen = FALSE, SelectionEffect = "Color", SelectionColor = "#FF0000",
    PanelId = 1L, PanelHeight = 600L, PanelWidth = 6L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    RowSelectionDynamicSource = FALSE, RowSelectionType = "Active",
    RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE,
    ColumnSelectionType = "Active", ColumnSelectionSaved = 0L,
    SelectionHistory = list())

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new("RowDataTable", Selected = "MOBP", Search = "", SearchColumns = c("",
    "", "", "", "", "", "", "", "", "", "", "", "", "", ""), PanelId = c(RowDataTable = 1L),
    PanelHeight = 600L, PanelWidth = 6L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, RowSelectionType = "Active",
    RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE,
    ColumnSelectionType = "Active", ColumnSelectionSaved = 0L,
    SelectionHistory = list())

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data",
    XAxisColumnData = "cell_type", XAxisFeatureName = "MOBP",
    XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE,
    YAxisFeatureName = "MOBP", YAxisFeatureSource = "RowDataTable1",
    YAxisFeatureDynamicSource = TRUE, ColorByColumnData = "cell_type",
    ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000",
    ShapeByColumnData = "donor", SizeByColumnData = "sum", FacetByRow = "---",
    FacetByColumn = "---", ColorBy = "Column data", ColorByDefaultColor = "#000000",
    ColorByFeatureName = "MOBP", ColorByFeatureSource = "---",
    ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "dlpfc.5161_AAACCCACACCGTCGA-1",
    ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None", SizeBy = "None", SelectionEffect = "Transparent",
    SelectionColor = "#FF0000", SelectionAlpha = 0.1, ZoomData = numeric(0),
    BrushData = list(), VisualBoxOpen = FALSE, VisualChoices = "Color",
    ContourAdd = FALSE, ContourColor = "#0000FF", PointSize = 1,
    PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
    FontSize = 1, LegendPosition = "Bottom", PanelId = c(FeatureAssayPlot = 1L),
    PanelHeight = 600L, PanelWidth = 6L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, RowSelectionType = "Active",
    RowSelectionSaved = 0L, ColumnSelectionDynamicSource = FALSE,
    ColumnSelectionType = "Active", ColumnSelectionSaved = 0L,
    SelectionHistory = list())

iSEE(sce, appTitle = "M.N. Tran et al 2020, DLPFC region https://bit.ly/LIBD10xHuman", initial = initial)
