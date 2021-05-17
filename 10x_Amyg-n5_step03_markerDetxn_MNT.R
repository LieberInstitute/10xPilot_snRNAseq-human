### MNT 10x snRNA-seq workflow: step 03 - marker detection
###   **Region-specific analyses**
###     - (2x) amygdala samples from: Br5161 & Br5212
### Initiated MNT 12Feb2020
#####################################################################

library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(batchelor)
library(DropletUtils)
library(jaffelab)
library(limma)

### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

# ===




## Load SCE with new info
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_Amyg-n5_cleaned-combined_SCE_MNT2021.rda",
     verbose=T)
    # sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo, annotationTab.amy

table(sce.amy$cellType)
    #   ambig.glial       Astro_A       Astro_B drop.lowNTx_A drop.lowNTx_B 
    #            31          1555            83          1067            71 
    #          Endo       Excit_A       Excit_B       Inhib_A       Inhib_B 
    #            70           399            44           780           541 
    #       Inhib_C       Inhib_D       Inhib_E       Inhib_F         Micro 
    #           525           500           555           216          1201 
    #         Oligo           OPC 
    #          6080          1459


# First drop "drop.lowNTx_" (1138 nuclei)
sce.amy <- sce.amy[ ,-grep("drop.",sce.amy$cellType)]
sce.amy$cellType <- droplevels(sce.amy$cellType)

# Remove 0 genes across all nuclei
sce.amy <- sce.amy[!rowSums(assay(sce.amy, "counts"))==0, ]  # keeps same 29371 genes


### Run ANOVA real quick ====
library(edgeR)
library(doMC)
registerDoMC(cores=16)

mat = assay(sce.amy, "counts")

## do regression
varCompAnalysis.amy = foreach(i = 1:nrow(mat)) %dopar% {
  if(i %% 100 == 0) cat(".")
  fit = lm(as.numeric(mat[i,]) ~ cellType + donor + processBatch + sequencer + sex + protocol,
           data=colData(sce.amy))
  full = anova(fit)
  fullSS = full$"Sum Sq"
  signif(cbind(full, PctExp = fullSS/sum(fullSS)*100), 3)
}

names(varCompAnalysis.amy) = rownames(mat)

## make boxplot
varExpl = t(sapply(varCompAnalysis.amy, function(x) x[,"PctExp"]))
colnames(varExpl) = rownames(varCompAnalysis.amy[[1]])



#pdf("pdfs/revision/anova_Amyg_MNT2021.pdf")
boxplot(varExpl, main="ANOVA on human Amyg 10x snRNA-seq \n (sn-level, n=5)",
        ylab="Percent Var explained (%)")
#dev.off()
# ok so DEF keep 'individualID'

save(varCompAnalysis.amy, file="rdas/revision/anova_Amyg-n5_MNT2021.rda")

apply(varExpl, 2, function(x){quantile(x, na.rm=T)})





## Traditional t-test === ===
 #    - this has been the most reliable (i.e. consistently yielding results) method, for every population

mod <- with(colData(sce.amy), model.matrix(~ donor))
    #If try to `+ processBatch`, get Error in .ranksafeQR(design, error = rank.error) :\n  design matrix is not of full rank
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`
    ## e.g. Warning message:
    #       In .fit_lm_internal(x, subset.row, groups, design = design, direction = direction,  :
    #       automatically removed intercept column

# Run pairwise t-tests
markers.amy.t.pw <- findMarkers(sce.amy, groups=sce.amy$cellType,
                                assay.type="counts", design=mod, test="t",
                                direction="up", pval.type="all", full.stats=T)

sapply(markers.amy.t.pw, function(x){table(x$FDR<0.05)})
    #       ambig.glial Astro_A Astro_B  Endo Excit_A Excit_B Inhib_A Inhib_B Inhib_C
    # FALSE       29059   28890   29368 28961   26281   29146   25588   29042   29292
    # TRUE          312     481       3   410    3090     225    3783     329      79
    
    #       Inhib_D Inhib_E Inhib_F Micro Oligo   OPC
    # FALSE   29218   28911   27881 28969 29198 29154
    # TRUE      153     460    1490   402   173   217

        ## modeling with 'logcounts' - made interactively:
        sapply(markers.amy.t.pw.logcounts, function(x){table(x$FDR<0.05)})
            # ambig.glial Astro_A Astro_B  Endo Excit_A Excit_B Inhib_A Inhib_B Inhib_C
            # FALSE       28554   28712   29255 28537   28946   28658   29052   29167   29309
            # TRUE          817     659     116   834     425     713     319     204      62
            # Inhib_D Inhib_E Inhib_F Micro Oligo   OPC
            # FALSE   28936   29126   29221 28550 28860 29037
            # TRUE      435     245     150   821   511   334


## What is that 'ambig.glial'?
head(rownames(markers.amy.t.pw[["ambig.glial"]]), n=20)
    #   [1] "SKAP1"     "SLFN12L"   "CD2"       "IKZF3"     "SLAMF6"    "RUNX3"    
    #   [7] "P2RY8"     "ITK"       "GRAP2"     "IL7R"      "SLAMF7"    "LINC00861"
    # [13] "THEMIS"    "PRF1"      "SLAMF1"    "GPR174"    "LCK"       "IL32"     
    #  [19] "TRBC2"     "CD247" 


## WMW: Blocking on donor (this test doesn't take 'design=' argument) ===
markers.amy.wilcox.block <- findMarkers(sce.amy, groups=sce.amy$cellType,
                                        assay.type="counts", block=sce.amy$donor, test="wilcox",
                                        direction="up", pval.type="all", full.stats=T)

# WMW FDR<0.05
sapply(markers.amy.wilcox.block, function(x){table(x$FDR<0.05)["TRUE"]})
    #ambig.glial.TRUE     Astro_A.TRUE       Astro_B.NA          Endo.NA 
    #               8              258               NA               NA 
    #      Excit_A.NA       Excit_B.NA     Inhib_A.TRUE     Inhib_B.TRUE 
    #              NA               NA             2315               60 
    #    Inhib_C.TRUE       Inhib_D.NA     Inhib_E.TRUE     Inhib_F.TRUE 
    #               4               NA              100               29 
    #      Micro.TRUE       Oligo.TRUE         OPC.TRUE 
    #             102              122              124


## Binomial ===
markers.amy.binom.block <- findMarkers(sce.amy, groups=sce.amy$cellType,
                                       assay.type="counts", block=sce.amy$donor, test="binom",
                                       direction="up", pval.type="all", full.stats=T)

sapply(markers.amy.binom.block, function(x){table(x$FDR<0.05)["TRUE"]})
    # only a couple dozen hits for a few - disregard these
    #ambig.glial.NA   Astro_A.TRUE     Astro_B.NA        Endo.NA     Excit_A.NA 
    #            NA             51             NA             NA             NA 
    #    Excit_B.NA     Inhib_A.NA     Inhib_B.NA     Inhib_C.NA     Inhib_D.NA 
    #            NA             NA             NA             NA             NA 
    #    Inhib_E.NA     Inhib_F.NA       Micro.NA     Oligo.TRUE       OPC.TRUE 
    #            NA             NA             NA              8             24

## Save all these for future reference
save(markers.amy.t.pw, markers.amy.wilcox.block, #markers.amy.binom.block,
     file="rdas/revision/markers-stats_Amyg-n5_findMarkers-SN-LEVEL_MNT2021.rda")


        # Btw - some have 0 p.value's and FDR's
        head(markers.amy.t.pw[["Excit_A"]][ ,1:4])
        head(markers.amy.t.pw[["Endo"]][ ,1:4])
        


# Print these to pngs
markerList.t <- lapply(markers.amy.t.pw, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)
genes.top40.t <- lapply(markerList.t, function(x){head(x, n=40)})

#dir.create("pdfs/revision/Amyg/")
for(i in names(genes.top40.t)){
  png(paste0("pdfs/revision/Amyg/Amyg_t_pairwise_top40markers-", i, "_logExprs_MNT2021.png"), height=1900, width=1200)
  print(
    plotExpression(sce.amy, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType", colour_by="cellType", point_alpha=0.3, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun = median, fun.min = median, fun.max = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:15], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers: single-nucleus-level p.w. t-tests"))
  )
  dev.off()
}






### Cluster-vs-all single-nucleus-level iteration ======

## Load SCE with new info
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo

table(sce.amy$cellType)

# First drop "Ambig.lowNtrxts" (50 nuclei)
sce.amy <- sce.amy[ ,sce.amy$cellType != "Ambig.lowNtrxts"]
sce.amy$cellType <- droplevels(sce.amy$cellType)

# Remove 0 genes across all nuclei
sce.amy <- sce.amy[!rowSums(assay(sce.amy, "counts"))==0, ]  # keeps same 28464 genes


## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.amy), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

markers.amy.t.1vAll <- list()
for(i in levels(sce.amy$cellType)){
  # Make temporary contrast
  sce.amy$contrast <- ifelse(sce.amy$cellType==i, 1, 0)
  # Test cluster vs. all
  markers.amy.t.1vAll[[i]] <- findMarkers(sce.amy, groups=sce.amy$contrast,
                                            assay.type="logcounts", design=mod, test="t",
                                            direction="up", pval.type="all", full.stats=T)
}

    ## Then, temp set of stats to get the standardized logFC
    temp.1vAll <- list()
    for(i in levels(sce.amy$cellType)){
      # Make temporary contrast
      sce.amy$contrast <- ifelse(sce.amy$cellType==i, 1, 0)
      # Test cluster vs. all
      temp.1vAll[[i]] <- findMarkers(sce.amy, groups=sce.amy$contrast,
                                     assay.type="logcounts", design=mod, test="t",
                                     std.lfc=TRUE,
                                     direction="up", pval.type="all", full.stats=T)
    }


## For some reason all the results are in the second List entry (first is always empty)

# Replace that empty slot with the entry with the actul stats
markers.amy.t.1vAll <- lapply(markers.amy.t.1vAll, function(x){ x[[2]] })
# Same for that with std.lfc
temp.1vAll <- lapply(temp.1vAll, function(x){ x[[2]] })

# Now just pull from the 'stats.0' DataFrame column
markers.amy.t.1vAll <- lapply(markers.amy.t.1vAll, function(x){ x$stats.0 })
temp.1vAll <- lapply(temp.1vAll, function(x){ x$stats.0 })

# Re-name std.lfc column and add to the first result
for(i in names(temp.1vAll)){
  colnames(temp.1vAll[[i]])[1] <- "std.logFC"
  markers.amy.t.1vAll[[i]] <- cbind(markers.amy.t.1vAll[[i]], temp.1vAll[[i]]$std.logFC)
  # Oh the colname is kept weird
  colnames(markers.amy.t.1vAll[[i]])[4] <- "std.logFC"
  # Then re-organize
  markers.amy.t.1vAll[[i]] <- markers.amy.t.1vAll[[i]][ ,c("logFC","std.logFC","log.p.value","log.FDR")]
}


## Let's save this along with the previous pairwise results
save(markers.amy.t.1vAll, markers.amy.t.pw, markers.amy.wilcox.block,
     file="rdas/markers-stats_Amyg-n2_findMarkers-SN-LEVEL_MNTMay2020.rda")




## Print these to pngs
markerList.t.1vAll <- lapply(markers.amy.t.1vAll, function(x){
  rownames(x)[x[ ,"log.FDR"] < log10(0.001)]
  }
)
genes.top40.t.1vAll <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t.1vAll)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/Amyg/Amyg_t_1vALL_top40markers-",gsub(":",".",i),"_logExprs_May2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.amy, exprs_values = "logcounts", features=genes.top40.t.1vAll[[i]],
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun = median, fun.min = median, fun.max = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:12], length(genes.top40.t.1vAll[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers: single-nucleus-level p.w. t-tests, cluster-vs-all"))
  )
  dev.off()
}






## How do these top 40 intersect? ===
sapply(names(genes.top40.t), function(c){
  length(intersect(genes.top40.t[[c]],
                   genes.top40.t.1vAll[[c]]))
})
    #  Astro Excit.1 Excit.2 Excit.3 Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5   Micro
    #     35       5      37      32       3      21      32      15      24      37
    #  Oligo     OPC
    #     30      27



## Write these top 40 lists to a csv
names(markerList.t) <- paste0(names(markerList.t),"_pw")
names(markerList.t.1vAll) <- paste0(names(markerList.t.1vAll),"_1vAll")

# PW result for "Inhib.1" doesn't have 40 markers:
markerList.t[["Inhib.1_pw"]] <- c(markerList.t[["Inhib.1_pw"]], rep("",9))

top40genes <- cbind(sapply(markerList.t, function(x) head(x, n=40)),
                    sapply(markerList.t.1vAll, function(y) head(y, n=40)))
top40genes <- top40genes[ ,sort(colnames(top40genes))]

write.csv(top40genes, file="tables/top40genesLists_Amyg-n2_cellType_SN-LEVEL-tests_May2020.csv",
          row.names=FALSE)



## 15May2020 for AnJa - print t's and FDRs ===
# Bring in human stats; create t's
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_Amyg-n2_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.amy.t.1vAll, markers.amy.t.pw, markers.amy.wilcox.block
    rm(markers.amy.t.pw, markers.amy.wilcox.block)
    
# Need to add t's with N nuclei used in constrasts
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    #sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo
    rm(chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy,ref.sampleInfo)

# First drop "ambig.lowNtrxts" (93 nuclei)
sce.amy <- sce.amy[ ,sce.amy$cellType != "Ambig.lowNtrxts"]
sce.amy$cellType <- droplevels(sce.amy$cellType)

## As above, calculate and add t-statistic (= std.logFC * sqrt(N)) from contrasts
#      and fix row order to the first entry "Astro"
fixTo <- rownames(markers.amy.t.1vAll[["Astro"]])
for(s in names(markers.amy.t.1vAll)){
  markers.amy.t.1vAll[[s]]$t.stat <- markers.amy.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.amy))
  markers.amy.t.1vAll[[s]] <- markers.amy.t.1vAll[[s]][fixTo, ]
}
# Pull out the t's
ts.amy <- sapply(markers.amy.t.1vAll, function(x){x$t.stat})
rownames(ts.amy) <- fixTo
# Change to EnsemblID
rownames(ts.amy) <- rowData(sce.amy)$ID[match(rownames(ts.amy), rownames(sce.amy))]

# logFDRs
logFDRs.amy <- sapply(markers.amy.t.1vAll, function(x){x$log.FDR})
rownames(logFDRs.amy) <- fixTo
# EnsemblID
rownames(logFDRs.amy) <- rownames(ts.amy)

save(ts.amy, logFDRs.amy, file="rdas/zForAnJa_AMY-sn-level-markerStats_MNT.rda")

apply(logFDRs.amy, 2, function(x) {table(x<log10(1e-6))})
    # Astro Excit.1 Excit.2 Excit.3 Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5
    # FALSE 23194   19604   25841   25201   21979   23391   26181   26787   24329
    # TRUE   5270    8860    2623    3263    6485    5073    2283    1677    4135
    # Micro Oligo   OPC
    # FALSE 24520 25698 24702
    # TRUE   3944  2766  3762



### For fig: Plot some top markers in vlnplot array (12Jun2020) === === === ===
  #     * had been done on local machine, but realized chunk isn't in script
  #       - updating 31Aug2020

load("rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
load("rdas/markers-stats_Amyg-n2_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.amy.t.1vAll, markers.amy.t.pw, markers.amy.wilcox.block
    # Focus on the pairwise result (".pw") bc more specific
    rm(markers.amy.t.1vAll, markers.amy.wilcox.block)

# First drop "ambig.lowNtrxts" (50 nuclei)
sce.amy <- sce.amy[ ,sce.amy$cellType != "Ambig.lowNtrxts"]
sce.amy$cellType <- droplevels(sce.amy$cellType)


# Take top two for broad glia
topToPrint <- as.data.frame(sapply(markers.amy.t.pw, function(x) {head(rownames(x),n=2)}))
topToPrint <- c(topToPrint$Astro, c("NPTX1", "SLC30A3", # Excit.1
                                    "SLC17A6", "SOX4", "SOX11", #Excit.2
                                    "MCHR2", "CDH22", # Excit.3
                                    "CCK", "CALB2", "KIT", # Inhib.1/2/4
                                    "CNTNAP3", "CNTNAP3B", "CALB1", # Inhib.3
                                    "NPFFR2", "TLL1"), # Inhib.5
                topToPrint$Micro, topToPrint$Oligo, topToPrint$OPC)

table(topToPrint %in% rownames(sce.amy)) # good

# With top 2 per glial
pdf("pdfs/pubFigures/Amyg_topMarkers-ARRAY_logExprs_Jun2020_v1.pdf", height=17, width=4)
print(
  plotExpression(sce.amy, exprs_values = "logcounts", features=topToPrint,
                 x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=1,
                 add_legend=F) + stat_summary(fun = median, fun.min = median, fun.max = median,
                                              geom = "crossbar", width = 0.3,
                                              colour=rep(tableau20[1:12], length(topToPrint))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12), plot.title = element_text(size = 25)) +  
    ggtitle(label="AMY top marker array (with glial)")
)
dev.off()

# Neuronal markers only (highlighted in paper)
topToPrint <- c("NRN1", "NPTX1", "SLC30A3", # Excit.1
  "SLC17A6", "VCAN", #Excit.2
  "MCHR2", "RBFOX3", # Excit.3
  "CCK", "CALB2", "KIT", # Inhib.1/2/4
  "CRH", "CALB1", # Inhib.3
  "NPFFR2", "TLL1")

pdf("pdfs/pubFigures/Amyg_topMarkers-ARRAY_logExprs_Jun2020_v2.pdf", height=11, width=4.5)
print(
  plotExpression(sce.amy, exprs_values = "logcounts", features=topToPrint,
                 x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=1,
                 add_legend=F, scales="free_y") + stat_summary(fun = median, fun.min = median, fun.max = median,
                                              geom = "crossbar", width = 0.3,
                                              colour=rep(tableau20[1:12], length(topToPrint))) +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 13),
          axis.text.y = element_text(size = 6.5),
          axis.title.y = element_text(angle = 90, size = 16),
          plot.title = element_text(size = 20),
          panel.grid.major=element_line(colour="grey95", size=0.8),
          panel.grid.minor=element_line(colour="grey95", size=0.4))# +  
    #ggtitle(label="AMY top marker array (neuronal only)")
)
dev.off()





### MNT add 18Nov2020 =================================
# -> What if add param/requirement that for any given subcluster, median expression has to > 0?
load("rdas/markers-stats_Amyg-n2_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.amy.t.1vAll, markers.amy.t.pw, markers.amy.wilcox.block

## Load SCE 
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)
    # sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo

table(sce.amy$cellType)

# First drop "Ambig.lowNtrxts" (50 nuclei)
sce.amy <- sce.amy[ ,sce.amy$cellType != "Ambig.lowNtrxts"]
sce.amy$cellType <- droplevels(sce.amy$cellType)

# Remove 0 genes across all nuclei
sce.amy <- sce.amy[!rowSums(assay(sce.amy, "counts"))==0, ]


## Make list of Boolean param / cell subtype ===
cellSubtype.idx <- splitit(sce.amy$cellType)
medianNon0.idx <- lapply(cellSubtype.idx, function(x){
  apply(as.matrix(assay(sce.amy, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

lengths(medianNon0.idx)
sapply(medianNon0.idx, head)

# Add these to the stats for each set of markers
for(i in names(markers.amy.t.1vAll)){
  markers.amy.t.1vAll[[i]] <- cbind(markers.amy.t.1vAll[[i]],
                                     medianNon0.idx[[i]][match(rownames(markers.amy.t.1vAll[[i]]),
                                                               names(medianNon0.idx[[i]]))])
  colnames(markers.amy.t.1vAll[[i]])[5] <- "non0median"
}


## Use these restrictions to print (to png) a refined top 40, as before ===
markerList.t.1vAll <- lapply(markers.amy.t.1vAll, function(x){
  rownames(x)[x$log.FDR < log10(0.000001) & x$non0median==TRUE]
  }
)
    # lengths(markerList.t.1vAll)     # ( **without $non0median==TRUE restriction )
        #  Astro Excit.1 Excit.2 Excit.3 Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5   Micro
        #   5270    8860    2623    3263    6485    5073    2283    1677    4135    3944
        #  Oligo     OPC
        #   2766    3762

lengths(markerList.t.1vAll)
    #  Astro Excit.1 Excit.2 Excit.3 Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5   Micro
    #   1326    4344     948    1795    3424    3031    1647    1121    2500     823
    #  Oligo     OPC
    #    976    1361

genes.top40.t <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/Amyg/Amyg_t_1vALL_top40markers-REFINED-",gsub(":",".",i),"_logExprs_Nov2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.amy, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun = median, fun.min = median, fun.max = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:12], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers, refined: single-nucleus-level p.w. t-tests, cluster-vs-all"))
  )
  dev.off()
}



## Do the same with the pairwise result ('markers.amy.t.pw') === === ===
# Add these to the stats for each set of markers
for(i in names(markers.amy.t.pw)){
  markers.amy.t.pw[[i]] <- cbind(markers.amy.t.pw[[i]],
                                      medianNon0.idx[[i]][match(rownames(markers.amy.t.pw[[i]]),
                                                                names(medianNon0.idx[[i]]))])
  colnames(markers.amy.t.pw[[i]])[14] <- "non0median"
}

markerList.t <- lapply(markers.amy.t.pw, function(x){
  rownames(x)[x$FDR < 0.05 & x$non0median==TRUE]
  }
)
    # lengths(markerList.t)     # ( **without $non0median==TRUE restriction )
        #  Astro Excit.1 Excit.2 Excit.3 Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5   Micro
        #    713      98     750     653      31     223     330     143     249     989
        #  Oligo     OPC
        #    434     294

lengths(markerList.t)
    #  Astro Excit.1 Excit.2 Excit.3 Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5   Micro
    #    403      44     226     350       7     142     182      34      86     405
    #  Oligo     OPC
    #    371     186


genes.top40.t <- lapply(markerList.t, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/pdfs/exploration/Amyg/Amyg_t_pairwise_top40markers-REFINED-", i, "_logExprs_Nov2020.png"), height=1900, width=1200)
  print(
    plotExpression(sce.amy, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType", colour_by="cellType", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun = median, fun.min = median, fun.max = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:12], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers, refined: single-nucleus-level p.w. t-tests"))
  )
  dev.off()
}

## Then write a new CSV of these refined top 40 genes ===
names(markerList.t) <- paste0(names(markerList.t),"_pw")
names(markerList.t.1vAll) <- paste0(names(markerList.t.1vAll),"_1vAll")

# Many of the PW results don't have at least 40 markers:
extend.idx <- names(which(lengths(markerList.t) < 40))
for(i in extend.idx){
  markerList.t[[i]] <- c(markerList.t[[i]], rep("", 40-length(markerList.t[[i]])))
}

top40genes <- cbind(sapply(markerList.t, function(x) head(x, n=40)),
                    sapply(markerList.t.1vAll, function(y) head(y, n=40)))
top40genes <- top40genes[ ,sort(colnames(top40genes))]

write.csv(top40genes, file="tables/top40genesLists-REFINED_Amyg-n2_cellType_Nov2020.csv",
          row.names=FALSE)





## Aside: add in 't.stat' as in 'step04' analyses to save for LoHu/LeCo ===
for(s in names(markers.amy.t.1vAll)){
  markers.amy.t.1vAll[[s]]$t.stat <- markers.amy.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.amy))
}

save(markers.amy.t.1vAll, markers.amy.t.pw, sce.amy,
     file="rdas/markerStats-and-SCE_AMY-n2_sn-level_cleaned_MNTNov2020.rda")








### Aside - difference b/tw dropping 'Ambig.lowNtrxts' before or after PB'ing ===========
  # MNT 05Mar2020
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)


# First drop "Ambig.lowNtrxts" (50 nuclei)
sce.amy.noAmbig <- sce.amy[ ,sce.amy$cellType != "Ambig.lowNtrxts"]
sce.amy.noAmbig$cellType <- droplevels(sce.amy.noAmbig$cellType)
# Then make the pseudo-bulked SCE
sce.amy.PBafterDrop <- aggregateAcrossCells(sce.amy.noAmbig, ids=paste0(sce.amy.noAmbig$sample,":",sce.amy.noAmbig$cellType),
                                              use_exprs_values="counts")
# Drop genes with all 0's
sce.amy.PBafterDrop <- sce.amy.PBafterDrop[!rowSums(assay(sce.amy.PBafterDrop, "counts"))==0, ]
    ## 28464 remaining genes

## OR

# Just make the pseudo-bulked SCE, without dropping that cluster
sce.amy.PB <- aggregateAcrossCells(sce.amy, ids=paste0(sce.amy$sample,":",sce.amy$cellType),
                                     use_exprs_values="counts")
# Drop genes with all 0's
sce.amy.PB <- sce.amy.PB[!rowSums(assay(sce.amy.PB, "counts"))==0, ]
    ## 28470

## Then
genesOfInterest <- setdiff(rownames(sce.amy.PB), rownames(sce.amy.PBafterDrop))
#  [1] "FOXC1"      "CSAG1"      "AP000879.1" "TBX3"       "FOXL1"        "FOXS1"

rowSums(assay(sce.amy.PB, "counts")[genesOfInterest, ])
## all 1's# 166 zeros

    ##    - so this is definitely just a poor, lowly-captured cluster; explains the poor
    ##      intra-cluster correlation, compared to other clusters



## And old...
# For BoG abstract ===
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_Amyg_n2_MNTFeb2020.rda", verbose=T)
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_DLPFC_n2_MNTFeb2020.rda", verbose=T)

intersectingMarkers <- list()
for(i in names(markerList.PB.tDesign.amy)){
  intersectingMarkers[[i]] <- intersect(markerList.PB.tDesign.amy[[i]],
                                        markerList.PB.tDesign.dlpfc[[i]])
}
# CEACAM21 a microglial marker - turns out to have been implicated in SCZ
#     (Jewish-Israeli familial study - I don't think in large GWAS though)

specificMarkers.amy <- list()
for(i in names(markerList.PB.tDesign.amy)){
  specificMarkers.amy[[i]] <- setdiff(markerList.PB.tDesign.amy[[i]],
                                      markerList.PB.tDesign.dlpfc[[i]])  
}

specificMarkers.dlpfc <- list()
for(i in names(markerList.PB.tDesign.amy)){
  specificMarkers.dlpfc[[i]] <- setdiff(markerList.PB.tDesign.dlpfc[[i]],
                                        markerList.PB.tDesign.amy[[i]])
}