library(ggplot2)
library(ggrepel)
library(dplyr)
library(data.table)
library(plotly)
library(htmlwidgets)
library(sessioninfo)
library(ggpubr)
library(tools)
library(GGally)
library(tidyverse)
library(patchwork)

data.table::setDTthreads(threads = 1)

# Sourcing Data/Inst. Vars. ####
load("rda/twas_exp_ranges.Rdata")
# load("twas_exp_ranges.Rdata")

dir.create(file.path("analysis/plots"),
           showWarnings = FALSE,
           recursive = TRUE)
dir.create(file.path("analysis/tables"),
           showWarnings = FALSE,
           recursive = TRUE)

# Filter N/A Z scores
twas_z <- twas_exp_fin %>% filter(!is.na(TWAS.Z))

# twas_z_NAc <- twas_z[twas_z$region == "NAc",]

don <- list()

axisdf <- list()

don_key <- list()

p <- list()

intctv_plot <- list()

fin_plot <- list()

for (i in 1:5) {
    # Preprocessing Data ####
    if (i == 1) {
        twas_var <- twas_z[type == "aoi"]
    } else if (i == 2) {
        twas_var <- twas_z[type == "cpd"]
    } else if (i == 3) {
        twas_var <- twas_z[type == "dpw"]
    } else if (i == 4) {
        twas_var <- twas_z[type == "sc"]
    } else if (i == 5) {
        twas_var <- twas_z[type == "si"]
    }
    
    don[[i]] <-
        twas_var %>%
        # Compute chromosome size
        group_by(CHR) %>%
        summarise(chr_len = max(end)) %>%
        
        # Calculate cumulative position of each chromosome
        mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
        select(-chr_len) %>%
        
        # Add this info to the initial dataset
        left_join(twas_var,
                  .,
                  by = c("CHR" = "CHR")) %>%
        
        # Add a cumulative position of each SNP
        arrange(CHR, twas_mean_dist) %>%
        mutate(BPcum = twas_mean_dist + tot)
    
    
    axisdf[[i]] = don[[i]] %>% group_by(CHR) %>% summarise(center = (max(BPcum) + min(BPcum)) / 2)
    
    # Prepare text description for each SNP:
    don[[i]]$text <-
        paste0(
            "Gene Symbol: ",
            don[[i]]$genesymbol,
            "\nENSEMBL Gene ID: ",
            don[[i]]$geneid,
            "\nBrain Subregion: ",
            don[[i]]$region,
            "\nChromosome: ",
            don[[i]]$CHR,
            "\nStart Position: ",
            don[[i]]$start,
            "\nEnd Position: ",
            don[[i]]$end,
            "\nZ score: ",
            don[[i]]$TWAS.Z %>% round(2)
        )
    
    don_key[[i]] <-
        highlight_key(don[[i]], ~ genesymbol, group = "Gene Symbol")
}

# TWAS Z Manhattan Plot ####
pdf(file = "analysis/plots/10x_NAc_TWAS_ManhattanPlot.pdf")
# storing ggplot as an object3

sig <- qnorm(1 - 0.025 / table(twas_exp_fin$type))

for (i in 1:5) {
    # Bonferroni Correction
    sig_bonf <- sig[[i]]
    
    p[[i]] <-
        ggplot(don_key[[i]], aes(x = BPcum, y = TWAS.Z, text = text)) +
        
        ggtitle(paste0(
            "Gene Windows of ",
            case_when(
                i == 1 ~ "Age of Initiation",
                i == 2 ~ "Cigarettes Per Day",
                i == 3 ~ "Drinks Per Week",
                i == 4 ~ "Smoking Cessation",
                i == 5 ~ "Smoking Initiation"
            ) ,
            " TWAS"
        )) +
        # Show all points
        geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
        scale_color_manual(values = rep(c("#861657", "#D56AA0"), 22)) +
        geom_hline(
            yintercept = c(sig_bonf, -sig_bonf),
            color = "grey40",
            linetype = "dashed"
        ) +
        
        # custom X axis:
        scale_x_continuous(labels = axisdf[[i]]$CHR, breaks = axisdf[[i]]$center) +
        scale_y_continuous(expand = c(0, 0)) +     # remove space between plot area and x axis
        
        # Custom the theme:
        theme_bw() +
        theme(
            legend.position = "none",
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()
        )
    
    print(p[[i]])
}
dev.off()

# Z scores threshold
twas_z_NAc_threshold <- list()

for (i in 1:5) {
    twas_z_NAc_threshold[[i]] <-
        rbind(twas_z[TWAS.Z > sig[[i]],], twas_z[TWAS.Z < -sig[[i]],])
}

# Interactive TWAS Z Manhattan Plots ####
for (i in 1:5) {
    ##### Plotly
    intctv_plot[[i]] <- ggplotly(p[[i]], tooltip = "text")
    
    fin_plot[[i]] <- highlight(
        intctv_plot[[i]],
        on = "plotly_click",
        off = "plotly_doubleclick",
        color = "#60D394",
        selectize = TRUE
    )
    
    saveWidget(fin_plot[[i]],
               file.path(paste0(
                   # "analysis/plots/",
                   "NAc_TWAS_",
                   case_when(i == 1 ~ "aoi",
                             i == 2 ~ "cpd",
                             i == 3 ~ "dpw",
                             i == 4 ~ "sc",
                             i == 5 ~ "si"),
                   "_ManhattanPlotly.html"
               )))
}

# Plotly cannot save directly to relative path for whatever reason
system("mv *_ManhattanPlotly.html analysis/plots/")

## Correlations of TWAS Z scores across different GWAS sets ####
pdf(file = "analysis/plots/10x_NAc_TWAS_GWAS_Correlations.pdf")

final_2 <- bind_rows(don) %>%
    select(-c(FILE, text)) %>% 
    select(geneid, BPcum, genesymbol, type, TWAS.Z, TWAS.P) %>%
    as.data.table()

final_wide <-
    dcast(final_2,
          geneid + genesymbol ~ type,
          value.var = c("TWAS.Z"))

# sample_test <- sample_n(final_2, 10000)

make_plot <- function(df, var1, var2) {
    df$var1 <- var1
    df$var2 <- var2
    xlabel <- if (var2 == "si") var1 else NULL
    ylabel <- if (var1 == "aoi") var2 else NULL
    color_label <- if (var1 == "aoi") "type" else NULL
    xaxis <- if (var2 != "si") theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) else NULL
    yaxis <- if (var1 != "aoi") theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) else NULL
    ggplot(df, aes(x = .data[[var1]], y = .data[[var2]], color = var1)) +
        geom_point(alpha = 0.15) +
        scale_color_manual(values = cols, labels = var1) +
        labs(x = xlabel, y = ylabel, color = color_label) +
        xaxis +
        yaxis
}

var1 <- names(final_wide)[-c(1:2)]
var2 <- lapply(1:4, function(x) var1[(x + 1):5])

cols <- scales::hue_pal()(5)
cols <- setNames(cols, var1)
labs <- setNames(var1, var1)

plot_list <- map2(var1[1:4], var2, function(var1, var2) {
    blank <- rep(list(plot_spacer()), 4 - length(var2))
    map(var2, ~ make_plot(final_wide, var1, .x)) %>% c(blank, .)
})
plot_list <- reduce(plot_list, c)

wrap_plots(plot_list, nrow = 4, byrow = FALSE) +
    plot_layout(guides = "collect") &
    theme(legend.spacing.y = unit(0, "pt"),
          legend.margin = margin(1, 1, 1, 1 , "pt"),
          legend.title = element_text(margin = margin(0, 0, 3, 0, "pt")))

ggpairs(final_wide %>% select(-c(geneid, genesymbol)))

ggcorr(final_wide %>% select(-c(geneid, genesymbol)), method = c("pairwise", "pearson"))

dev.off()
dev.off()

# # Scatter plots ####
# pdf(
#     'analysis/plots/NAc_TWAS_ScatterPlots.pdf',
#     useDingbats = FALSE,
#     width = 10,
#     height = 10
# )
# 

# 
# # Render Z scores and P values horizontally by region
# twas_z_wide <-
#     dcast(twas_z_select,
#           geneid + genesymbol ~ region,
#           value.var = c("TWAS.Z", "TWAS.P"))
# 
# # FDR calculation per subregion
# twas_z_wide$amygdala.fdr.p <-
#     p.adjust(twas_z_wide$TWAS.P_amygdala, 'fdr')
# twas_z_wide$sacc.fdr.p <- p.adjust(twas_z_wide$TWAS.P_sacc, 'fdr')
# 
# # Indicate in both
# twas_z_wide$in_both <-
#     ifelse(!is.na(twas_z_wide$TWAS.Z_amygdala &
#                       twas_z_wide$TWAS.Z_sacc),
#            TRUE,
#            FALSE)
# 
# # FDR cutoffs
# twas_z_wide$FDR.5perc <- 'None'
# twas_z_wide$FDR.5perc[twas_z_wide$amygdala.fdr.p < 0.05] <-
#     'amygdala'
# twas_z_wide$FDR.5perc[twas_z_wide$sacc.fdr.p < 0.05] <- 'sACC'
# twas_z_wide$FDR.5perc[twas_z_wide$amygdala.fdr.p < 0.05 &
#                           twas_z_wide$sacc.fdr.p < 0.05] <- 'Both'
# 
# # Remove NAs
# twas_z_wide[is.na(twas_z_wide)] <- 0
# 
# twas_z_wide$FDR.5perc <-
#     factor(twas_z_wide$FDR.5perc,
#            levels = c('None', 'amygdala', 'sACC', 'Both'))
# 
# ggplot(twas_z_wide,
#        aes(x = TWAS.Z_amygdala,
#            y = TWAS.Z_sacc,
#            color = FDR.5perc)) +
#     xlab("Amygdala Z-Score") +
#     ylab("sACC Z-Score") +
#     labs(color = "FDR < 5%") +
#     geom_point() +
#     coord_fixed() +
#     theme_bw(base_size = 20) +
#     scale_color_manual(values = c('grey80', 'dark orange', 'skyblue3', 'purple')) # you can define names
# 
# dev.off()
# 
# ## Z-Score Correlation Test ####
# 
# both_genes_sacc <-
#     twas_z_sacc[twas_z_sacc$geneid %in% twas_z_amyg$geneid]
# both_genes_amyg <-
#     twas_z_amyg[twas_z_amyg$geneid %in% twas_z_sacc$geneid]
# 
# z_score_cor <-
#     cor.test(both_genes_sacc$TWAS.Z, both_genes_amyg$TWAS.Z, method = "pearson")
# z_score_cor
# 
# both_z_scores <- both_genes_sacc %>%
#     select(TWAS.Z) %>%
#     mutate(
#         sacc_TWAS_Z_both = TWAS.Z,
#         amyg_TWAS_Z_both = both_genes_amyg$TWAS.Z,
#         .keep = "unused"
#     )
# 
# pdf(
#     'analysis/plots/BD_TWAS_Z_Correlation.pdf',
#     useDingbats = FALSE,
#     width = 10,
#     height = 10
# )
# 
# 
# ggscatter(
#     both_z_scores,
#     x = "amyg_TWAS_Z_both",
#     y = "sacc_TWAS_Z_both",
#     add = "reg.line",
#     add.params = list(color = "red", fill = "orangered1"),
#     conf.int = TRUE,
#     cor.coef = TRUE,
#     cor.method = "pearson",
#     alpha = 0.5,
#     xlab = toTitleCase("TWAS Z-Scores for amygdala genes"),
#     ylab = toTitleCase("TWAS Z-Scores for sACC genes")
# ) +
#     grids(size = 1) + bgcolor("gray90") + border("gray90")
# 
# dev.off()
# 
# # Differential Expression ####
# load(
#     "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/case_control/bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation.rda",
#     verbose = TRUE
# )
# 
# statOutGene <- statOut[statOut$Type == "Gene",] %>%
#     as.data.table(keep.rownames = "geneid") %>%
#     select(geneid, t_Amyg, t_sACC)
# 
# 
# merged_t <- merge(statOutGene, twas_z_wide, by = "geneid")
# 
# pdf(
#     'analysis/plots/BD_TWAS_t-stat.pdf',
#     useDingbats = FALSE,
#     width = 10,
#     height = 10
# )
# 
# amyg_rho <-
#     format(
#         cor(merged_t$TWAS.Z_amygdala, merged_t$t_Amyg),
#         digits = 3,
#         scientific = TRUE
#     )
# sacc_rho <-
#     format(
#         cor(merged_t$TWAS.Z_sacc, merged_t$t_sACC),
#         digits = 3,
#         scientific = TRUE
#     )
# 
# ggplot(merged_t,
#        aes(x = TWAS.Z_amygdala,
#            y = t_Amyg)) + geom_point() + labs(
#                title = toTitleCase("TWAS vs BD differential expression in Amygdala"),
#                x = "TWAS Z-Score",
#                y = "BD vs Control t-Statistic"
#            ) + annotate(
#                "text",
#                x = -5.5,
#                y = 6,
#                label = paste0("rho == ", formatC(amyg_rho, format = "e")),
#                parse = TRUE
#            ) + scale_y_continuous(breaks = c(-6, -3, 0, 3, 6)) + xlim(-6, 6) +
#     theme_bw(base_size = 20)
# 
# ggplot(merged_t,
#        aes(x = TWAS.Z_sacc,
#            y = t_sACC)) + geom_point() + labs(
#                title = toTitleCase("TWAS vs BD differential expression in sACC"),
#                x = "TWAS Z-Score",
#                y = "BD vs Control t-Statistic"
#            ) + annotate(
#                "text",
#                x = -5.5,
#                y = 6,
#                label = paste0("rho == ", formatC(sacc_rho, format = "e")),
#                parse = TRUE
#            ) + scale_y_continuous(breaks = c(-6, -3, 0, 3, 6)) +
#     scale_x_continuous(breaks = waiver()) +
#     theme_bw(base_size = 20)
# 
# dev.off()
# 
# # XLSX Output ####
# 
# # xlsx does not work with conda_R/4.0, needs 4.0.x
# save(
#     twas_z_amyg_threshold,
#     twas_z_sacc_threshold,
#     twas_z_wide,
#     merged_t,
#     "rda/xlsx_output.RData"
# )
# 
# # for enrichment test
# save.image("rda/generate_plots_data.RData")

## Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# > print("Reproducibility information:")
# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2020-11-13 17:04:12 EST"
# > proc.time()
#     user   system  elapsed
#    53.14     9.85 25140.28
# > options(width = 120)
# > session_info()
# - Session info -------------------------------------------------------------------------------------------------------
#  setting  value
#  version  R version 4.0.3 (2020-10-10)
#  os       Windows 10 x64
#  system   x86_64, mingw32
#  ui       RStudio
#  language (EN)
#  collate  English_United States.1252
#  ctype    English_United States.1252
#  tz       America/New_York
#  date     2020-11-13
#
# - Packages -----------------------------------------------------------------------------------------------------------
#  package     * version date       lib source
#  assertthat    0.2.1   2019-03-21 [1] CRAN (R 4.0.2)
#  backports     1.1.10  2020-09-15 [1] CRAN (R 4.0.2)
#  callr         3.5.1   2020-10-13 [1] CRAN (R 4.0.3)
#  cli           2.1.0   2020-10-12 [1] CRAN (R 4.0.3)
#  colorspace    1.4-1   2019-03-18 [1] CRAN (R 4.0.2)
#  crayon        1.3.4   2017-09-16 [1] CRAN (R 4.0.2)
#  crosstalk     1.1.0.1 2020-03-13 [1] CRAN (R 4.0.3)
#  data.table  * 1.13.2  2020-10-19 [1] CRAN (R 4.0.3)
#  desc          1.2.0   2018-05-01 [1] CRAN (R 4.0.2)
#  devtools    * 2.3.2   2020-09-18 [1] CRAN (R 4.0.3)
#  digest        0.6.25  2020-02-23 [1] CRAN (R 4.0.2)
#  dplyr       * 1.0.2   2020-08-18 [1] CRAN (R 4.0.2)
#  ellipsis      0.3.1   2020-05-15 [1] CRAN (R 4.0.2)
#  fansi         0.4.1   2020-01-08 [1] CRAN (R 4.0.2)
#  farver        2.0.3   2020-01-16 [1] CRAN (R 4.0.2)
#  fastmap       1.0.1   2019-10-08 [1] CRAN (R 4.0.3)
#  fs            1.5.0   2020-07-31 [1] CRAN (R 4.0.2)
#  generics      0.1.0   2020-10-31 [1] CRAN (R 4.0.3)
#  ggplot2     * 3.3.2   2020-06-19 [1] CRAN (R 4.0.3)
#  ggrepel     * 0.8.2   2020-03-08 [1] CRAN (R 4.0.2)
#  glue          1.4.2   2020-08-27 [1] CRAN (R 4.0.2)
#  gtable        0.3.0   2019-03-25 [1] CRAN (R 4.0.2)
#  htmltools     0.5.0   2020-06-16 [1] CRAN (R 4.0.2)
#  htmlwidgets * 1.5.2   2020-10-03 [1] CRAN (R 4.0.3)
#  httpuv        1.5.4   2020-06-06 [1] CRAN (R 4.0.3)
#  httr          1.4.2   2020-07-20 [1] CRAN (R 4.0.2)
#  jsonlite      1.7.1   2020-09-07 [1] CRAN (R 4.0.2)
#  labeling      0.4.2   2020-10-20 [1] CRAN (R 4.0.3)
#  later         1.1.0.1 2020-06-05 [1] CRAN (R 4.0.2)
#  lazyeval      0.2.2   2019-03-15 [1] CRAN (R 4.0.2)
#  lifecycle     0.2.0   2020-03-06 [1] CRAN (R 4.0.2)
#  magrittr      1.5     2014-11-22 [1] CRAN (R 4.0.2)
#  memoise       1.1.0   2017-04-21 [1] CRAN (R 4.0.2)
#  mime          0.9     2020-02-04 [1] CRAN (R 4.0.0)
#  munsell       0.5.0   2018-06-12 [1] CRAN (R 4.0.2)
#  pillar        1.4.6   2020-07-10 [1] CRAN (R 4.0.2)
#  pkgbuild      1.1.0   2020-07-13 [1] CRAN (R 4.0.2)
#  pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.0.2)
#  pkgload       1.1.0   2020-05-29 [1] CRAN (R 4.0.2)
#  plotly      * 4.9.2.1 2020-04-04 [1] CRAN (R 4.0.3)
#  prettyunits   1.1.1   2020-01-24 [1] CRAN (R 4.0.2)
#  processx      3.4.4   2020-09-03 [1] CRAN (R 4.0.2)
#  promises      1.1.1   2020-06-09 [1] CRAN (R 4.0.2)
#  ps            1.3.4   2020-08-11 [1] CRAN (R 4.0.2)
#  purrr         0.3.4   2020-04-17 [1] CRAN (R 4.0.2)
#  R6            2.5.0   2020-10-28 [1] CRAN (R 4.0.3)
#  Rcpp          1.0.5   2020-07-06 [1] CRAN (R 4.0.2)
#  remotes       2.2.0   2020-07-21 [1] CRAN (R 4.0.2)
#  rlang         0.4.7   2020-07-09 [1] CRAN (R 4.0.2)
#  rprojroot     1.3-2   2018-01-03 [1] CRAN (R 4.0.2)
#  rstudioapi    0.11    2020-02-07 [1] CRAN (R 4.0.2)
#  scales        1.1.1   2020-05-11 [1] CRAN (R 4.0.2)
#  sessioninfo * 1.1.1   2018-11-05 [1] CRAN (R 4.0.2)
#  shiny         1.5.0   2020-06-23 [1] CRAN (R 4.0.3)
#  testthat    * 2.3.2   2020-03-02 [1] CRAN (R 4.0.2)
#  tibble        3.0.4   2020-10-12 [1] CRAN (R 4.0.3)
#  tidyr         1.1.2   2020-08-27 [1] CRAN (R 4.0.2)
#  tidyselect    1.1.0   2020-05-11 [1] CRAN (R 4.0.2)
#  usethis     * 1.6.3   2020-09-17 [1] CRAN (R 4.0.2)
#  utf8          1.1.4   2018-05-24 [1] CRAN (R 4.0.2)
#  vctrs         0.3.4   2020-08-29 [1] CRAN (R 4.0.2)
#  viridisLite   0.3.0   2018-02-01 [1] CRAN (R 4.0.2)
#  withr         2.3.0   2020-09-22 [1] CRAN (R 4.0.3)
#  xtable        1.8-4   2019-04-21 [1] CRAN (R 4.0.3)
#  yaml          2.2.1   2020-02-01 [1] CRAN (R 4.0.0)
#
# [1] C:/Users/artas/Documents/R/win-library/4.0
# [2] C:/Program Files/R/R-4.0.3/library
