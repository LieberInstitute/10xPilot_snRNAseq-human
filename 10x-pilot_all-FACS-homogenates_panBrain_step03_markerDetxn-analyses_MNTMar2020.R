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

load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda",
     verbose=T)
    ## sce.all.n12, chosen.hvgs.all.n12, pc.choice.n12, ref.sampleInfo, clusterRefTab.all.n12

    ## ** First drop ambiguous clusters if needed

# For reference
table(sce.all.n12$sample, sce.all.n12$cellType)
    ## 12 x 19 table; but only 164 of 228 cells with at least 1 cell

# Make the pseudo-bulked SCE
sce.all.n12.PB <- aggregateAcrossCells(sce.all.n12, ids=paste0(sce.all.n12$sample,":",sce.all.n12$cellType),
                                   use_exprs_values="counts")

# Drop genes with all 0's
sce.all.n12.PB <- sce.all.n12.PB[!rowSums(assay(sce.all.n12.PB, "counts"))==0, ]
    ## keeps 30038 genes (28312 `>=5`)

# Clean up colData
colData(sce.all.n12.PB) <- colData(sce.all.n12.PB)[ ,c(13:17,19:20)]

# Remove stored `sizeFactors()` because this will mess you up
#     * Also, to be safe, can always provide manually-computed SFs:
sizeFactors(sce.all.n12.PB) <- NULL
LSFvec <- librarySizeFactors(sce.all.n12.PB)
sce.all.n12.PB <- logNormCounts(sce.all.n12.PB, size_factors=LSFvec)


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
#     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_panBrain_findMarkers_MNTMar2020.rda")


### ANOVA to find the most variable covars =========
  # Following `/dcl01/ajaffe/data/lab/singleCell/10x_pilot/check_varExpl.R`
library(edgeR)
library(doMC)
registerDoMC(cores=16)

# Pull pseudo-bulked, normalized expression
lc = assays(sce.all.n12.PB)$logcounts

## do regression
varCompAnalysis = foreach(i = 1:nrow(lc)) %dopar% {
  if(i %% 1000 == 0) cat(".")
#  fit = lm(as.numeric(lc[i,]) ~ region + processDate + donor + cellType,
#  fit = lm(as.numeric(lc[i,]) ~ cellType + region + processDate + donor,
#  fit = lm(as.numeric(lc[i,]) ~ cellType + processDate + region + donor,
#  fit = lm(as.numeric(lc[i,]) ~ region + processDate + donor + protocol + cellType,
  fit = lm(as.numeric(lc[i,]) ~ cellType + region + processDate + donor + protocol,
           data=colData(sce.all.n12.PB))
  
  full = anova(fit)
  fullSS = full$"Sum Sq"
  signif(cbind(full, PctExp = fullSS/sum(fullSS)*100), 3)
}

names(varCompAnalysis) = rownames(lc)

save(sce.all.n12.PB, varCompAnalysis,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/varCompAnalysis_MNTMar2020.rda")


## make boxplots
varExpl = t(sapply(varCompAnalysis, function(x) x[,"PctExp"]))
colnames(varExpl) = rownames(varCompAnalysis[[1]])
#boxplot(varExpl)


apply(varExpl, 2, function(x){quantile(x)})
    ##      region processDate    donor cellType Residuals
    #0%    0.00786    1.09e-34 4.91e-35     1.92     0.141
    #25%   1.12000    6.25e-01 2.84e-01    19.70    49.200
    #50%   1.94000    1.32e+00 7.15e-01    33.20    61.300
    #75%   3.06000    2.42e+00 1.50e+00    46.00    73.800
    #100% 39.60000    3.76e+01 2.11e+01    98.10    92.000


## What if listed in different order?
 #    ~ cellType + region + processDate + donor:
    ##    cellType   region processDate    donor Residuals
    #0%       4.34  0.00107    6.06e-04 2.81e-06     0.141
    #25%     19.80  1.07000    5.84e-01 2.57e-01    49.200
    #50%     33.50  1.84000    1.23e+00 6.55e-01    61.300
    #75%     46.40  2.92000    2.28e+00 1.41e+00    73.800
    #100%    99.90 38.50000    3.63e+01 2.14e+01    92.000

 #    ~ cellType + processDate + region + donor:
    ##    cellType processDate   region    donor Residuals
    #0%       4.34     0.00172  0.00102 2.81e-06     0.141
    #25%     19.80     0.70300  0.97100 2.57e-01    49.200
    #50%     33.50     1.35000  1.73000 6.55e-01    61.300
    #75%     46.40     2.34000  2.83000 1.41e+00    73.800
    #100%    99.90    28.30000 62.90000 2.14e+01    92.000

## Add protocol?
 # ~ region + processDate + donor + protocol + cellType
    # Looks like didn't get estimated...






#varCompAnalysis$SNAP25

###### Direct limma approach ####
#################################


## Load in sce.dlpfc and cluster:sample bulk
# Clean up colData
colData(sce.all.n12.PB) <- colData(sce.all.n12.PB)[ ,c(13:17,19:20)]
    # (already done up top)



## Extract the count data
mat <- assays(sce.all.n12.PB)$logcounts

## Build a group model
mod <- with(colData(sce.all.n12.PB), model.matrix(~ 0 + cellType + region + processDate + protocol))
colnames(mod) <- gsub('cellType', '', colnames(mod))

corfit <- duplicateCorrelation(mat, mod, block = sce.all.n12.PB$donor)
    # "Coefficients not estimable: protocolpseudoSort"
corfit$consensus.correlation
    # [1] -0.0031468

    ## What about mod without processDate?
    #mod.noProcessDate <- with(colData(sce.all.n12.PB), model.matrix(~ 0 + cellType + region))# + protocol))
    #corfit.noProcessDate <- duplicateCorrelation(mat, mod.noProcessDate, block = sce.all.n12.PB$donor)
    #corfit.noProcessDate$consensus.correlation  # -0.000193579
        # (And without protocol: -0.0002819187)
    
fit <-
  lmFit(
    mat,
    design = mod,
    block = sce.all.n12.PB$donor,
    correlation = corfit$consensus.correlation
  )

eb <- eBayes(fit)


## Skip this chunk for now - too many 'cellType's ====
## Contrasts for pairwise comparison
cellType_combs <- combn(colnames(mod), 2) # will be `choose(ncol(x), 2)` columns long, of course
cellType_contrasts <- apply(cellType_combs, 2, function(x) {
  z <- paste(x, collapse = '-')
  makeContrasts(contrasts = z, levels = mod)
})
rownames(cellType_contrasts) <- colnames(mod)
colnames(cellType_contrasts) <- apply(cellType_combs, 2, paste, collapse = '-')

eb_contrasts <- eBayes(contrasts.fit(fit, cellType_contrasts))

## Tabulating significant hits
pvals_contrasts <- eb_contrasts$p.value

data.frame(
  'FDRsig' = colSums(apply(pvals_contrasts, 2, p.adjust, 'fdr') < 0.05),
  'Pval10-6sig' = colSums(pvals_contrasts < 1e-6),
  'Pval10-8sig' = colSums(pvals_contrasts < 1e-8)
)
# end skip ====



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
  'FDRsig' = colSums(apply(pvals0_contrasts, 2, p.adjust, 'fdr') < 0.05),
  'Pval10-6sig' = colSums(pvals0_contrasts < 1e-6),
  'Pval10-8sig' = colSums(pvals0_contrasts < 1e-8)
)

# For only (+) t-stats
data.frame(
  'FDRsig' = colSums(apply(pvals0_contrasts, 2, p.adjust, 'fdr') < 0.01 &
                       t0_contrasts_cell > 0),
  'Pval10-6sig' = colSums(pvals0_contrasts < 1e-6 &
                            t0_contrasts_cell > 0),
  'Pval10-8sig' = colSums(pvals0_contrasts < 1e-8 &
                            t0_contrasts_cell > 0)
)

    ## With t > 0
    #             FDRsig Pval10.6sig Pval10.8sig
    #Ambig.hiVCAN    104           5           2
    #Astro.1        4900        1471         990
    #Astro.2         537          83          39
    #Excit.1         509         114          53
    #Excit.2         454         116          60
    #Excit.3         216          65          28
    #Excit.4         122           6           3
    #Excit.5         383         117          56
    #Excit.6         243          26           5
    #Excit.7         532         172          92
    #Excit.8         508         211         144
    #Inhib.1        1675         272         152
    #Inhib.2        1153         177          80
    #Inhib.3         889         198         108
    #Inhib.4         339          69          41
    #Inhib.5         358          77          21
    #Micro          5630        2031        1587
    #Oligo          3827        1087         743
    #OPC            2752         784         519




## Save for later
#eb_contrasts.all.n12.broad <- eb_contrasts
eb_list.all.n12.broad <- eb0_list

save(sce.all.n12.PB, eb_list.all.n12.broad, #eb_contrasts.all.n12.broad, 
     file = 'rdas/zmarkers-stats_panBrain_manualContrasts_MNTMar2020.rda')



