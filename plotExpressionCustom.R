
## TEST
# test_markers <- list(A = paste0("Gene_000", 1:4),
#                      B = paste0("Gene_000", 5:8))
# 
# plotExpressionCustom(sce = example_sce, features = test_markers$A, "A", anno_name = "Cell_Cycle")


plotExpressionCustom <- function(sce, features, features_name, anno_name = "cellType",
                                 point_alpha=0.2, point_size=0.7, ncol=2){
  scater::plotExpression(sce, 
                         exprs_values = "logcounts", 
                         features = features,
                         x = anno_name, 
                         colour_by = anno_name,
                         ncol = ncol,
                         point_alpha = point_alpha, 
                         point_size = point_size,
                         add_legend = F) +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(face = "italic")) +  
  ggtitle(label=paste0(features_name, " markers"))
}

# MNT note: if want to add a custom title, just call this function and can 'overwrite' that with
#           another `+ ggtitle(...)`
