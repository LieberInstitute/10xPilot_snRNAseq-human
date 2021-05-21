
## TEST
# test_markers <- list(A = paste0("Gene_000", 1:4),
#                      B = paste0("Gene_000", 5:8))
# 
# plotExpressionCustom(sce = example_sce, features = test_markers$A, "A", anno_name = "Cell_Cycle")


plotExpressionCustom <- function(sce, features, features_name ,anno_name = "cellType"){
  scater::plotExpression(sce, 
                         exprs_values = "logcounts", 
                         features = features,
                         x = anno_name, 
                         colour_by = anno_name, 
                         point_alpha = 0.2, 
                         point_size = .7,
                         add_legend = F) +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(face = "italic")) +  
  ggtitle(label=paste0(features_name, " markers"))
}
