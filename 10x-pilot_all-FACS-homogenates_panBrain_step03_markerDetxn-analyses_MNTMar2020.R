### MNT 10x snRNA-seq workflow: step 03 - marker detection
###   **Pan-brain analyses**
###     - n=12 samples from 5 regions, up to three donors
### MNT Mar2020
#####################################################################

library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(limma)
library(scater)
library(scran)
library(batchelor)
library(DropletUtils)
library(jaffelab)
library(dendextend)
library(dynamicTreeCut)

### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

# ===


### Cell type marker gene detection =======================================
#   ** Approach - pseudo-bulk on sample:cellType stratification, then treat as SCE, so that
#                 can use 'findMarkers()' function that utilizes different tests

load("/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda",
     verbose=T)
    ## sce.all.n12, chosen.hvgs.all.n12, pc.choice.n12, ref.sampleInfo, clusterRefTab.all.n12

    # **MNT 16Apr: Deciding to remove the clusters won't focus on:
    #              'Ambig.hiVCAN' & 'Excit.4' & those .RS-annot'd as 'Ambig.lowNtrxts'

sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType.RS != "Ambig.lowNtrxts"] # 445
sce.all.n12$cellType.RS <- droplevels(sce.all.n12$cellType.RS)

sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType != "Ambig.hiVCAN"] # 32 nuclei
sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType != "Excit.4"]  # 33 nuclei
sce.all.n12$cellType <- droplevels(sce.all.n12$cellType)


# For reference
table(sce.all.n12$sample, sce.all.n12$cellType)
    ## 12 x 17

# Make the pseudo-bulked SCE
sce.all.n12.PB <- aggregateAcrossCells(sce.all.n12, ids=paste0(sce.all.n12$sample,":",sce.all.n12$cellType),
                                   use_exprs_values="counts")

# Drop genes with all 0's
sce.all.n12.PB <- sce.all.n12.PB[!rowSums(assay(sce.all.n12.PB, "counts"))==0, ]
    ## keeps 30,031 genes (x 149 sample:cellType cols)

# Clean up colData
colData(sce.all.n12.PB) <- colData(sce.all.n12.PB)[ ,c(13:17,19:20)]

# Remove stored `sizeFactors()` because this will mess you up
#     * Also, to be safe, can always provide manually-computed SFs:
sizeFactors(sce.all.n12.PB) <- NULL
LSFvec <- librarySizeFactors(sce.all.n12.PB)
sce.all.n12.PB <- logNormCounts(sce.all.n12.PB, size_factors=LSFvec)


### `findMarkers()` BioC approach (skip as of Apr2020) ====================================
## Find markers using stringent [max-p-value-of-all-pw-comparisons] test ('pval.type="all"')
 #    - first try no 'design=' argument
markers.all.n12.t <- findMarkers(sce.all.n12.PB, groups=sce.all.n12.PB$cellType,
                             assay.type="logcounts",
                             direction="up", pval.type="all", full.stats=T)

sapply(markers.all.n12.t, function(x){table(x$FDR<0.05)})
    # actually micros has (224); oligos (37); OPC (8); Astro.1 (8); Inhib.3 (1)



## With 'design=' - first try raw [sum of] "counts"
design.PB <- model.matrix(~region + processDate + donor, data=colData(sce.all.n12.PB))
design.PB <- design.PB[ , -1, drop=F] # 'drop=F' to keep as matrix - otherwise turns into numeric vec

markers.all.n12.t.design.counts <- findMarkers(sce.all.n12.PB, groups=sce.all.n12.PB$cellType,
                                           assay.type="counts", design=design.PB,
                                           direction="up", pval.type="all", full.stats=T)

sapply(markers.all.n12.t.design.counts, function(x){table(x$FDR<0.05)[2]})
    ##  Ambig.hiVCAN.NA    Astro.1.TRUE      Astro.2.NA      Excit.1.NA    Excit.2.TRUE
    #                NA             771              NA              NA              98
    #      Excit.3.TRUE      Excit.4.NA    Excit.5.TRUE      Excit.6.NA    Excit.7.TRUE
    #               125              NA              43              NA              26
    #      Excit.8.TRUE    Inhib.1.TRUE      Inhib.2.NA    Inhib.3.TRUE      Inhib.4.NA
    #                23              18              NA              27              NA
    #        Inhib.5.NA      Micro.TRUE      Oligo.TRUE        OPC.TRUE
    #                NA             864            9788             241


## normalized and transformed "logcounts"
markers.all.n12.t.design.log <- findMarkers(sce.all.n12.PB, groups=sce.all.n12.PB$cellType,
                                        assay.type="logcounts", design=design.PB,
                                        direction="up", pval.type="all", full.stats=T)

sapply(markers.all.n12.t.design.log, function(x){table(x$FDR<0.05)[2]})
    ##  Ambig.hiVCAN.NA    Astro.1.TRUE      Astro.2.NA      Excit.1.NA    Excit.2.TRUE
    #                NA             317              NA              NA              28
    #      Excit.3.TRUE      Excit.4.NA    Excit.5.TRUE    Excit.6.TRUE    Excit.7.TRUE
    #                 6              NA              16               1              46
    #      Excit.8.TRUE    Inhib.1.TRUE    Inhib.2.TRUE    Inhib.3.TRUE    Inhib.4.TRUE
    #               201              16               6              14               8
    #        Inhib.5.NA      Micro.TRUE      Oligo.TRUE        OPC.TRUE
    #                NA            1000             315             147




## Normalized $counts, but not log-transformed
#assay(sce.all.n12.PB, "countsNormd") <- t(apply(assay(sce.all.n12.PB, "counts"), 1, function(x) {x/LSFvec}))

markers.all.n12.t.design.countsN <- findMarkers(sce.all.n12.PB, groups=sce.all.n12.PB$cellType,
                                            assay.type="countsNormd", design=design.PB,
                                            direction="up", pval.type="all", full.stats=T)

sapply(markers.all.n12.t.design.countsN, function(x){table(x$FDR<0.05)[2]})
    ##  Ambig.hiVCAN.TRUE      Astro.1.TRUE      Astro.2.TRUE      Excit.1.TRUE
    #                   2               291                53                 7
    #        Excit.2.TRUE      Excit.3.TRUE      Excit.4.TRUE      Excit.5.TRUE
    #                  66                13                 1                40
    #        Excit.6.TRUE      Excit.7.TRUE      Excit.8.TRUE      Inhib.1.TRUE
    #                   7               127               523                43
    #        Inhib.2.TRUE      Inhib.3.TRUE      Inhib.4.TRUE        Inhib.5.NA
    #                  11                85                24                NA
    #          Micro.TRUE        Oligo.TRUE          OPC.TRUE
    #                1537               696               210


      ## Using 'block=' is N/A in this PB'd setting:
      #    -> warnings such as "no within-block comparison between ____ and ____"


# What if add donor into the mix? (not expected to be estimable within these other covars)
#design.PB <- model.matrix(~region + processDate, data=colData(sce.all.n12.PB))
#design.PB <- design.PB[ , -1, drop=F]
    ## Run ANOVA first (below chunk)


# Save prelim results
#save(markers.all.n12.t.design.log, markers.all.n12.t.design.countsN,
#     file="/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_panBrain_findMarkers_MNTMar2020.rda")


### ANOVA to find the most variable covars ================================================
  # Following `/dcl01/ajaffe/data/lab/singleCell/10x_pilot/check_varExpl.R`
library(edgeR)
library(doMC)
registerDoMC(cores=16)

# Pull pseudo-bulked, normalized expression
lc = assays(sce.all.n12.PB)$logcounts

## do regression
varCompAnalysis = foreach(i = 1:nrow(lc)) %dopar% {
# or
#varCompAnalysis.don1st = foreach(i = 1:nrow(lc)) %dopar% {
  if(i %% 1000 == 0) cat(".")
  fit = lm(as.numeric(lc[i,]) ~ cellType + region + processDate + donor,
  # or
  #fit = lm(as.numeric(lc[i,]) ~ cellType + region + donor + processDate,
#  fit = lm(as.numeric(lc[i,]) ~ region + processDate + donor + cellType,
#  fit = lm(as.numeric(lc[i,]) ~ cellType + processDate + region + donor,
#  fit = lm(as.numeric(lc[i,]) ~ cellType + region + processDate + donor + protocol,
           data=colData(sce.all.n12.PB))
  
  full = anova(fit)
  fullSS = full$"Sum Sq"
  signif(cbind(full, PctExp = fullSS/sum(fullSS)*100), 3)
}

names(varCompAnalysis) = rownames(lc)
names(varCompAnalysis.don1st) = rownames(lc)

save(sce.all.n12.PB, varCompAnalysis, varCompAnalysis.don1st,
     file="/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/varCompAnalysis_MNTApr2020.rda")


## make boxplots
varExpl = t(sapply(varCompAnalysis.don1st, function(x) x[,"PctExp"]))
colnames(varExpl) = rownames(varCompAnalysis.don1st[[1]])
boxplot(varExpl)


apply(varExpl, 2, function(x){quantile(x)})
    ##     cellType   region processDate    donor Residuals
    # 0%       3.65  0.00132    2.98e-04 4.68e-06     0.151
    # 25%     18.40  1.12000    6.61e-01 2.51e-01    49.200
    # 50%     31.40  1.99000    1.39e+00 6.70e-01    62.800
    # 75%     46.10  3.20000    2.56e+00 1.50e+00    74.800
    # 100%    99.80 39.50000    3.98e+01 2.51e+01    93.000

    # or if listed donor before processDate
    ##     cellType   region    donor processDate Residuals
    # 0%       3.65  0.00132 1.53e-04    5.05e-04     0.151
    # 25%     18.40  1.12000 3.22e-01    6.19e-01    49.200
    # 50%     31.40  1.99000 7.93e-01    1.29e+00    62.800
    # 75%     46.10  3.20000 1.62e+00    2.42e+00    74.800
    # 100%    99.80 39.50000 6.22e+01    3.75e+01    93.000

        ## (Looks like the 'protocol' effect couldn't get estimated when added to the model)



#varCompAnalysis$SNAP25




###### Direct limma approach ####
#################################

## Extract the count data
mat <- assays(sce.all.n12.PB)$logcounts

## Build a group model
mod <- with(colData(sce.all.n12.PB), model.matrix(~ 0 + cellType + region + processDate))
colnames(mod) <- gsub('cellType', '', colnames(mod))

corfit <- duplicateCorrelation(mat, mod, block = sce.all.n12.PB$donor)
    # "Coefficients not estimable: protocolpseudoSort"
    #       MNT comment: this is ok because it doesn't appear this technical flaw did much in the end
corfit$consensus.correlation
    # [1] -0.004489989


fit <-
  lmFit(
    mat,
    design = mod,
    block = sce.all.n12.PB$donor,
    correlation = corfit$consensus.correlation
  )

eb <- eBayes(fit)


## Contrasts for pairwise comparison
cellType_combs <- combn(colnames(mod), 2)
    # will be `choose(ncol(x), 2)` columns long, of course
cellType_contrasts <- apply(cellType_combs, 2, function(x) {
  z <- paste(x, collapse = '-')
  makeContrasts(contrasts = z, levels = mod)
})
rownames(cellType_contrasts) <- colnames(mod)
colnames(cellType_contrasts) <- apply(cellType_combs, 2, paste, collapse = '-')

eb_contrasts <- eBayes(contrasts.fit(fit, cellType_contrasts))

# ## Tabulating significant hits    - too big
# pvals_contrasts <- eb_contrasts$p.value
# 
# data.frame(
#   'FDRsig' = colSums(apply(pvals_contrasts, 2, p.adjust, 'fdr') < 0.05),
#   'Pval10-6sig' = colSums(pvals_contrasts < 1e-6),
#   'Pval10-8sig' = colSums(pvals_contrasts < 1e-8)
# )




## Then each cellType vs the rest (leave 'protocol' out, bc not estimable)
cellType_idx <- splitit(sce.all.n12.PB$cellType)

eb0_list <- lapply(cellType_idx, function(x) {
  res <- rep(0, ncol(sce.all.n12.PB))
  res[x] <- 1
  m <- model.matrix(~ res + region + processDate, data=colData(sce.all.n12.PB))
  eBayes(
    lmFit(
      mat,
      design = m,
      block = sce.all.n12.PB$donor,
      correlation = corfit$consensus.correlation
    )
  )
})

## Extract the p-values
pvals0_contrasts <- sapply(eb0_list, function(x) {
  x$p.value[, 2, drop = FALSE]
})


## Extract the tstats
t0_contrasts_cell <- sapply(eb0_list, function(x) {
  x$t[, 2, drop = FALSE]
})


data.frame(
  'FDRsig.05' = colSums(apply(pvals0_contrasts, 2, p.adjust, 'fdr') < 0.05),
  'FDRsig.01' = colSums(apply(pvals0_contrasts, 2, p.adjust, 'fdr') < 0.01),
  'Pval10-6sig' = colSums(pvals0_contrasts < 1e-6),
  'Pval10-8sig' = colSums(pvals0_contrasts < 1e-8)
)

# For only (+) t-stats
data.frame(
  'FDRsig.05' = colSums(apply(pvals0_contrasts, 2, p.adjust, 'fdr') < 0.05 &
                       t0_contrasts_cell > 0),
  'FDRsig.01' = colSums(apply(pvals0_contrasts, 2, p.adjust, 'fdr') < 0.01 &
                          t0_contrasts_cell > 0),
  'Pval10-6sig' = colSums(pvals0_contrasts < 1e-6 &
                            t0_contrasts_cell > 0),
  'Pval10-8sig' = colSums(pvals0_contrasts < 1e-8 &
                            t0_contrasts_cell > 0)
)

    ## With t > 0
    #        FDRsig.05 FDRsig.01 Pval10.6sig Pval10.8sig
    # Astro.1      4531      3209        1385         915
    # Astro.2       510       350          81          37
    # Excit.1       478       269         101          38
    # Excit.2       398       243          96          52
    # Excit.3       198       127          52          24
    # Excit.5       326       220         100          48
    # Excit.6       300       154          12           5
    # Excit.7       461       318         149          76
    # Excit.8       452       324         190         124
    # Inhib.1      1364       732         241         137
    # Inhib.2       984       479         143          67
    # Inhib.3       740       405         173         101
    # Inhib.4       253       120          62          30
    # Inhib.5       358       196          83          27
    # Micro        5368      3838        1967        1526
    # Oligo        3590      2449        1064         726
    # OPC          2478      1662         753         483




## Save for later
eb_contrasts.all.n12.broad <- eb_contrasts
eb_list.all.n12.broad <- eb0_list

save(sce.all.n12.PB, eb_list.all.n12.broad, eb_contrasts.all.n12.broad, 
     file = 'rdas/markers-stats_panBrain_manualContrasts-w-pairwise_MNTApr2020.rda')





### Top markers to look more into =========
load("rdas/markers-stats_panBrain_manualContrasts-w-pairwise_MNTApr2020.rda", verbose=T)
    # sce.all.n12.PB, eb_list.all.n12.broad, eb_contrasts.all.n12.broad


## Extract the p-values
pvals0_contrasts <- sapply(eb_list.all.n12.broad, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts) <- rownames(sce.all.n12.PB)

fdrs0_contrasts = apply(pvals0_contrasts, 2, p.adjust, "fdr")

## Extract the tstats
t0_contrasts_cell <- sapply(eb_list.all.n12.broad, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts_cell) <- rownames(sce.all.n12.PB)


# Make list of top markers at FDR < 0.01
markerList.PB.manual <- lapply(colnames(fdrs0_contrasts), function(x){
  rownames(fdrs0_contrasts)[fdrs0_contrasts[ ,x] < 0.01 & t0_contrasts_cell[ ,x] > 0]
})
names(markerList.PB.manual) <- colnames(fdrs0_contrasts)
lengths(markerList.PB.manual)
    # Astro.1 Astro.2 Excit.1 Excit.2 Excit.3 Excit.5 Excit.6 Excit.7 Excit.8 Inhib.1
    #    3209     350     269     243     127     220     154     318     324     732
    # Inhib.2 Inhib.3 Inhib.4 Inhib.5   Micro   Oligo     OPC
    #     479     405     120     196    3838    2449    1662

markerTs.fdr.01 <- lapply(colnames(fdrs0_contrasts), function(x){
  as.matrix(t0_contrasts_cell[fdrs0_contrasts[ ,x] < 0.01 & t0_contrasts_cell[ ,x] > 0, x])
})

names(markerTs.fdr.01) <- colnames(fdrs0_contrasts)

markerList.sorted <- lapply(markerTs.fdr.01, function(x){
  x[,1][order(x, decreasing=TRUE)]
})

# Check, from broad cell type marker expression plots (SLC17A6 seemed very specific to "Excit.8")
"SLC17A6" %in% names(markerList.sorted[["Excit.8"]])  # TRUE

# Let's print top 40 and to PNGs instead of pdfs
genes2plot <- lapply(markerList.sorted, function(x){head(x, n=40)})


## Let's plot some expression of these to see how much are 'real' (not driven by outliers)
load("/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda",
     verbose=T)
## sce.all.n12, chosen.hvgs.all.n12, pc.choice.n12, ref.sampleInfo, clusterRefTab.all.n12

# As before, remove the clusters that won't focus on:
#     'Ambig.hiVCAN' & 'Excit.4' & those .RS-annot'd as 'Ambig.lowNtrxts'

sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType.RS != "Ambig.lowNtrxts"] # 445
sce.all.n12$cellType.RS <- droplevels(sce.all.n12$cellType.RS)

sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType != "Ambig.hiVCAN"] # 32 nuclei
sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType != "Excit.4"]  # 33 nuclei
sce.all.n12$cellType <- droplevels(sce.all.n12$cellType)

#dir.create("pdfs/exploration/panBrainMarkers/")

library(grDevices)
for(i in names(genes2plot)){
  png(paste0("pdfs/exploration/panBrainMarkers/panBrain_top40markers-",i,"_logExprs_Apr2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.all.n12, exprs_values = "logcounts", features=c(names(genes2plot[[i]])),
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:17], length(genes2plot[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(i, " top 40 markers"))
  )
  dev.off()
}





### Brief test with `findMarkers()`, WMW at single-nucleus-level ======================
  # MNT 23Apr2020

## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.all.n12), model.matrix(~ region + processDate + donor))
mod <- mod[ ,-1]
    # This matrix treats the 'amy' region as the baseline:
    #design.full <- with(colData(sce.all.n12), model.matrix(~region + processDate))
markers.n12.t.design <- findMarkers(sce.all.n12, groups=sce.all.n12$cellType,
                                        assay.type="logcounts", design=mod, test="t",
                                        direction="up", pval.type="all", full.stats=T)

sapply(markers.n12.t.design, function(x){table(x$FDR<0.01)})
    #       Astro.1 Astro.2 Excit.1 Excit.2 Excit.3 Excit.5 Excit.6 Excit.7 Excit.8
    # FALSE   32852   33458   33200   33306   33392   33248   33382   32939   32930
    # TRUE      686      80     338     232     146     290     156     599     608
    #       Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5 Micro Oligo   OPC
    # FALSE   33325   33302   33301   33383   33222 32241 32926 33183
    # TRUE      213     236     237     155     316  1297   612   355


## WMW: Blocking on sample (this test doesn't take 'design=' argument) ===
markers.n12.wilcox.block <- findMarkers(sce.all.n12, groups=sce.all.n12$cellType,
                                          assay.type="logcounts", block=sce.all.n12$donor, test="wilcox",
                                          direction="up", pval.type="all", full.stats=T)
    # Warning messages:
    #   1: In .pairwise_blocked_template(x, group.vals, nblocks = length(block),  :
    #                                      no within-block comparison between Excit.8 and Excit.5
    #                                    2: In .pairwise_blocked_template(x, group.vals, nblocks = length(block),  :
    #                                                                       no within-block comparison between Excit.8 and Excit.7


sapply(markers.n12.wilcox.block, function(x){table(x$FDR<0.05)})
    ## if 'block=' on $sample, very poor... all excitatory subtypes will have 0 markers
     
     # but blocking on $donor:

    #       Astro.1 Astro.2 Excit.1 Excit.2 Excit.3 Excit.5 Excit.6 Excit.7 Excit.8
    # FALSE   33147   33537   33359   33423   33498   33378   33529   33131   33425
    # TRUE      391       1     179     115      40     160       9     407     113
    #       Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5 Micro Oligo   OPC
    # FALSE   33441   33326   33380   33455   33339 32899 33076 33301
    # TRUE       97     212     158      83     199   639   462   237

## Binomial ===
markers.n12.binom.block <- findMarkers(sce.all.n12, groups=sce.all.n12$cellType,
                                        assay.type="logcounts", block=sce.all.n12$donor, test="binom",
                                        direction="up", pval.type="all", full.stats=T)
    # [48x] Warning messages:
    # [...]
    # 10: In pbinom(host.nzero, size, p, log.p = TRUE) :
    #   pbeta(*, log.p=TRUE) -> bpser(a=1645, b=39, x=0.577685,...) underflow to -Inf
    # 11: In pbinom(host.nzero - 1, size, p, lower.tail = FALSE,  ... :
    #   pbeta(*, log.p=TRUE) -> bpser(a=1586, b=35, x=0.577685,...) underflow to -Inf
    # 12: In pbinom(host.nzero - 1, size, p, lower.tail = FALSE,  ... :
    #   pbeta(*, log.p=TRUE) -> bpser(a=1720, b=37, x=0.577685,...) underflow to -Inf
    # 13: In pbinom(host.nzero, size, p, log.p = TRUE) :
    #   pbeta(*, log.p=TRUE) -> bpser(a=5102, b=39, x=0.858593,...) underflow to -Inf
    # 14: In pbinom(host.nzero, size, p, log.p = TRUE) :
    #   pbeta(*, log.p=TRUE) -> bpser(a=3670, b=39, x=0.785647,...) underflow to -Inf
    # 15: In pbinom(host.nzero, size, p, log.p = TRUE) :
    #   pbeta(*, log.p=TRUE) -> bpser(a=3334, b=31, x=0.785647,...) underflow to -Inf
    # 16: In pbinom(host.nzero - 1, size, p, lower.tail = FALSE,  ... :
    # [...]

sapply(markers.n12.binom.block, function(x){table(x$FDR<0.05)})
    # ~ dozen to 100+ across most; none for Astro.2 or Excit.6

## Save all these for future reference
save(markers.n12.t.design, markers.n12.wilcox.block, markers.n12.binom.block,
     file="rdas/markers-stats_panBrain_allTests-single-nuc-lvl_MNTApr2020.rda")



## For funsies, print those from the WMW result: === === ===
 # - these are already ordered from best stats to worst

    ## btw, good sign that:
    which(rownames(markers.n12.wilcox.block[["Excit.8"]])=="SLC17A6")
    # [1] 34

markerList.wmw <- lapply(markers.n12.wilcox.block, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)

genes2plot.wmw <- sapply(markerList.wmw, function(x){head(x, n=40)})

## Plot WMW results ======
for(i in setdiff(names(genes2plot.wmw), c("Astro.2","Excit.6"))){
  png(paste0("pdfs/exploration/panBrainMarkers/panBrain_wmw-top40markers-",i,"_logExprs_Apr2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.all.n12, exprs_values = "logcounts", features=genes2plot.wmw[[i]],
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:17], length(genes2plot.wmw[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
      ggtitle(label=paste0(i, " top 40 markers"))
  )
  dev.off()
}

# Then for Astro.2
png("pdfs/exploration/panBrainMarkers/panBrain_wmw-top40markers-Astro.2_logExprs_Apr2020.png", height=1900/8, width=1200/6)
print(
  plotExpression(sce.all.n12, exprs_values = "logcounts", features=genes2plot.wmw[["Astro.2"]],
                 x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, #ncol=5,
                 add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                              geom = "crossbar", width = 0.3,
                                              colour=rep(tableau20[1:17], length(genes2plot.wmw[["Astro.2"]]))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
    ggtitle(label="Astro.2 top 40 markers")
)
dev.off()

# And Excit.6
png("pdfs/exploration/panBrainMarkers/panBrain_wmw-top40markers-Excit.6_logExprs_Apr2020.png", height=1900/2, width=1200)
print(
  plotExpression(sce.all.n12, exprs_values = "logcounts", features=genes2plot.wmw[["Excit.6"]],
                 x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                 add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                              geom = "crossbar", width = 0.3,
                                              colour=rep(tableau20[1:17], length(genes2plot.wmw[["Excit.6"]]))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
    ggtitle(label="Excit.6 top 40 markers")
)
dev.off()

    ## ====


# Compare to single-nucleus-level t-test? ===
markerList.t <- lapply(markers.n12.t.design, function(x){
  rownames(x)[x$FDR < 0.01]
  }
)

# FDR < 0.05 for WMW ===
lengths(markerList.wmw)

# FDR < 0.05 for t-test ===
lengths(markerList.t)

# Intersection ===
sapply(names(markerList.t), function(g){
  length(intersect(markerList.t[[g]],
                   markerList.wmw[[g]]))
})

## What about top 40 vs. top 40?
genes.top40.t <- lapply(markerList.t, function(x){head(x, n=40)})

# Intersection
sapply(names(genes.top40.t), function(g){
  length(intersect(genes2plot.wmw[[g]],
                   genes.top40.t[[g]]))
})
    # Astro.1 Astro.2 Excit.1 Excit.2 Excit.3 Excit.5 Excit.6 Excit.7 Excit.8 Inhib.1
    #      20       1      11      23      23      18       3      17      15      22
    # Inhib.2 Inhib.3 Inhib.4 Inhib.5   Micro   Oligo     OPC
    #      24      22      23      22      21      36      24     - so about half


## Let's plot the t-test results too ===
for(i in names(genes.top40.t)){
  png(paste0("pdfs/exploration/panBrainMarkers/panBrain_t-sn-level-top40markers-",i,"_logExprs_Apr2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.all.n12, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:17], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers: single-nucleus-level p.w. t-tests"))
  )
  dev.off()
}


pdf("pdfs/exploration/zForBoG-poster-SOX11-in-panBrain_Excit.8-marker.pdf", height=2.5, width=4)
plotExpression(sce.all.n12, exprs_values = "logcounts", features=c("SOX11"),
               x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,
               add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                            geom = "crossbar", width = 0.3,
                                            colour=rep(tableau20[1:17],1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()






### Cluster-vs-all single-nucleus-level iteration ====================================
  # Added 20May2020 - because this was done at all region-specific-level analyses,
  #                   once we accepted the sn-level tests to yield more robust results

# Load SCE
load("/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda",
     verbose=T)
    # sce.all.n12, chosen.hvgs.all.n12, pc.choice.n12, ref.sampleInfo, clusterRefTab.all.n12

sce.all.n12

# As before, remove the clusters that won't focus on:
#     'Ambig.hiVCAN' & 'Excit.4' & those .RS-annot'd as 'Ambig.lowNtrxts'

sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType.RS != "Ambig.lowNtrxts"] # 445
sce.all.n12$cellType.RS <- droplevels(sce.all.n12$cellType.RS)

sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType != "Ambig.hiVCAN"] # 32 nuclei
sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType != "Excit.4"]  # 33 nuclei
sce.all.n12$cellType <- droplevels(sce.all.n12$cellType)

# Remove 0 genes across all nuclei
sce.all.n12 <- sce.all.n12[!rowSums(assay(sce.all.n12, "counts"))==0, ]  # keeps same 30031

sce.all.n12



## Traditional t-test with design as above, for pairwise result
mod <- with(colData(sce.all.n12), model.matrix(~ region + processDate + donor))
mod <- mod[ ,-1]    # intercept otherwise automatically dropped by `findMarkers()`

markers.n12.t.1vAll <- list()
for(i in levels(sce.all.n12$cellType)){
  # Make temporary contrast
  sce.all.n12$contrast <- ifelse(sce.all.n12$cellType==i, 1, 0)
  # Test cluster vs. all
  markers.n12.t.1vAll[[i]] <- findMarkers(sce.all.n12, groups=sce.all.n12$contrast,
                                           assay.type="logcounts", design=mod, test="t",
                                           direction="up", pval.type="all", full.stats=T)
}

## Then, temp set of stats to get the standardized logFC
temp.1vAll <- list()
for(i in levels(sce.all.n12$cellType)){
  # Make temporary contrast
  sce.all.n12$contrast <- ifelse(sce.all.n12$cellType==i, 1, 0)
  # Test cluster vs. all
  temp.1vAll[[i]] <- findMarkers(sce.all.n12, groups=sce.all.n12$contrast,
                                 assay.type="logcounts", design=mod, test="t",
                                 std.lfc=TRUE,
                                 direction="up", pval.type="all", full.stats=T)
}


## For some reason all the results are in the second List entry (first is always empty)

# Replace that empty slot with the entry with the actul stats
markers.n12.t.1vAll <- lapply(markers.n12.t.1vAll, function(x){ x[[2]] })
# Same for that with std.lfc
temp.1vAll <- lapply(temp.1vAll, function(x){ x[[2]] })

# Now just pull from the 'stats.0' DataFrame column
markers.n12.t.1vAll <- lapply(markers.n12.t.1vAll, function(x){ x$stats.0 })
temp.1vAll <- lapply(temp.1vAll, function(x){ x$stats.0 })

# Re-name std.lfc column and add to the first result
for(i in names(temp.1vAll)){
  colnames(temp.1vAll[[i]])[1] <- "std.logFC"
  markers.n12.t.1vAll[[i]] <- cbind(markers.n12.t.1vAll[[i]], temp.1vAll[[i]]$std.logFC)
  # Oh the colname is kept weird
  colnames(markers.n12.t.1vAll[[i]])[4] <- "std.logFC"
  # Then re-organize
  markers.n12.t.1vAll[[i]] <- markers.n12.t.1vAll[[i]][ ,c("logFC","std.logFC","log.p.value","log.FDR")]
}


## Let's save this along with the previous pairwise results
load("rdas/markers-stats_panBrain_allTests-single-nuc-lvl_MNTApr2020.rda", verbose=T)
save(markers.n12.t.design, markers.n12.wilcox.block, markers.n12.binom.block,
     markers.n12.t.1vAll,
     file="rdas/markers-stats_panBrain_allTests-single-nuc-lvl_MNTApr2020.rda")


## Print these to pngs
markerList.t.1vAll <- lapply(markers.n12.t.1vAll, function(x){
  rownames(x)[x[ ,"log.FDR"] < log10(0.000001)]
  }
)
genes.top40.t.1vAll <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t.1vAll)){
  png(paste0("/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/panBrainMarkers/panBrain_t-sn-level_1vALL_top40markers-",i,"_logExprs_May2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.all.n12, exprs_values = "logcounts", features=genes.top40.t.1vAll[[i]],
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:17], length(genes.top40.t.1vAll[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers: single-nucleus-level t-tests, cluster-vs-all"))
  )
  dev.off()
}



## Write these top 40 lists to a csv ===
markerList.t.pw <- lapply(markers.n12.t.design, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)

names(markerList.t.pw) <- paste0(names(markerList.t.pw),"_pw")
names(markerList.t.1vAll) <- paste0(names(markerList.t.1vAll),"_1vAll")

# PW result for "Inhib.1" doesn't have 40 markers:
#markerList.t.pw[["Inhib.1_pw"]] <- c(markerList.t.pw[["Inhib.1_pw"]], rep("",9))

top40genes <- cbind(sapply(markerList.t.pw, function(x) head(x, n=40)),
                    sapply(markerList.t.1vAll, function(y) head(y, n=40)))
top40genes <- top40genes[ ,sort(colnames(top40genes))]

write.csv(top40genes, file="tables/top40genesLists_panBrain-n12_cellType_SN-LEVEL-tests_May2020.csv",
          row.names=FALSE)





















### Miscellaneous - For KeMa/AnJa 20May2020 =============================
# Plot separately bc doesn't look like you can plot two genes at once if facetting by region
pdf("pdfs/exploration/zCRHBP-OXTR_panBrain-clusters-by-region_MNT.pdf", height=3.5, width=10)
plotExpression(sce.all.n12, exprs_values = "logcounts", features="CRHBP",
               x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,
               add_legend=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) + facet_grid(~ sce.all.n12$region)

plotExpression(sce.all.n12, exprs_values = "logcounts", features="OXTR",
               x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,
               add_legend=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) + facet_grid(~ sce.all.n12$region)
dev.off()

pdf("pdfs/exploration/zNPY-SST_panBrain-clusters-by-region_MNT.pdf", height=3.5, width=10)
plotExpression(sce.all.n12, exprs_values = "logcounts", features="NPY",
               x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,
               add_legend=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) + facet_grid(~ sce.all.n12$region)

plotExpression(sce.all.n12, exprs_values = "logcounts", features="C",
               x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,
               add_legend=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) + facet_grid(~ sce.all.n12$region)
dev.off()

# CACNA1C
pdf("pdfs/exploration/zCACNA1C_panBrain-clusters-by-region_MNT.pdf", height=3.5, width=10)
plotExpression(sce.all.n12, exprs_values = "logcounts", features="CACNA1C",
               x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7,
               add_legend=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) + facet_grid(~ sce.all.n12$region)
dev.off()



dlpfc.inhib4 <- which(sce.all.n12$region=="dlpfc" & sce.all.n12$cellType=="Inhib.4")  # 133

table(sce.all.n12$region=="dlpfc" & sce.all.n12$cellType=="Inhib.4" & assay(sce.all.n12,"logcounts")[])

# DLPFC:Inhib.4
table(assay(sce.all.n12,"logcounts")["OXTR",dlpfc.inhib4] > 0 &      # 14
        assay(sce.all.n12,"logcounts")["CRHBP",dlpfc.inhib4] > 0 &   # 7
        assay(sce.all.n12,"logcounts")["SST",dlpfc.inhib4] > 0 &     # 6
        assay(sce.all.n12,"logcounts")["NPY",dlpfc.inhib4] > 0)      # 4

# All nuclei...
table(assay(sce.all.n12,"logcounts")["OXTR", ] > 0 &      # 1038
        assay(sce.all.n12,"logcounts")["CRHBP", ] > 0 &   # 103
        assay(sce.all.n12,"logcounts")["SST", ] > 0 &     # 32
        assay(sce.all.n12,"logcounts")["NPY", ] > 0)      # 15

allFourPos <- which(assay(sce.all.n12,"logcounts")["OXTR", ] > 0 &      # 1038
        assay(sce.all.n12,"logcounts")["CRHBP", ] > 0 &   # 103
        assay(sce.all.n12,"logcounts")["SST", ] > 0 &     # 32
        assay(sce.all.n12,"logcounts")["NPY", ] > 0)      # 15

table(sce.all.n12$cellType[allFourPos], sce.all.n12$region[allFourPos])

assay(sce.all.n12,"logcounts")["CACNA1C", allFourPos]



