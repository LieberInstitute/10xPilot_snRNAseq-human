################################################################################
### LIBD 10x snRNA-seq [pilot] revision (n=10)
### STEP 01: Read in SCEs and perform nuclei calling and QC
### Initiated: MNT 25Feb2021
################################################################################

library(SingleCellExperiment)
library(scRNAseq)
library(batchelor)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(uwot)
library(DropletUtils)
library(jaffelab)
library(Rtsne)
library(gridExtra)
library(rtracklayer)

### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

### =======

#### THIS CHUNK DONE IN SUBMITTED JOB ====
        
        ### Read in (2021) 'samples.manifest' for streamlining
        samples.revision <- read.table("/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-10-15_Tran2021_published/Feb2021/samples.manifest",
                                       sep="\t", header=F)$V1
        
        # Make list of paths
        paths.rawCounts <- c(paste0("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/Feb2021/",
                                    samples.revision,"/outs/raw_feature_bc_matrix"))
        # Make sure works
        sapply(paths.rawCounts, list.files) # good
        
        # Make names for individual SCEs
        names(paths.rawCounts) <- gsub("_",".", tolower(samples.revision))
        
        
        ## Read in raw UMI x barcode matrix - **use pre-mRNA-aligned reads
        pilot.data.2 <- lapply(paths.rawCounts, function(x){ read10xCounts(x, col.names=TRUE) })
        names(pilot.data.2) <- names(paths.rawCounts)
        
        
        
        ### Gene annotation (from scater) ===
        # Pull in GTF information
        gtf = import("/dcl01/ajaffe/data/lab/singleCell/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf")
        gtf = gtf[gtf$type == "gene"]
        length(gtf)
        names(gtf) = gtf$gene_id
        table(names(gtf) == rowData(pilot.data.2[[1]])$ID)
        
        seqlevels(gtf)[1:25] = paste0("chr", seqlevels(gtf)[1:25])
        mcols(gtf) = mcols(gtf)[,c(5:9)]
        
        
        for(i in 1:length(pilot.data.2)){
          rowRanges(pilot.data.2[[i]]) <- gtf
          # Because some gene names are the same:
          rownames(pilot.data.2[[i]]) <- uniquifyFeatureNames(rowData(pilot.data.2[[i]])$gene_id, rowData(pilot.data.2[[i]])$gene_name)
          rowData(pilot.data.2[[i]])$Symbol.uniq <- rownames(pilot.data.2[[i]])
        }
        
        
        # Preliminary save
        save(pilot.data.2, file="rdas/revision/all-FACS-n10_2021rev_SCEs_processing-QC_MNTFeb2021.rda")
        
        ### Quality control ============================================================
        ## - Going to ignore the adaptive NMAD-approach to outlier detection for UMI/feature count
        #    because this hasn't been as straightforward in past experience (might throw away neurons)
        ## - Vignette for the 10x PBMC dataset (OSCA Ch.24) only does mito & droplet QC anyhow
        #       - (mention that for a sample with very heterogeneous cell comp., don't want
        #          to drop potential cells with low RNA content)
        
        
        ## Cell detection (droplet exclusion, rather)
        # Can use UMI count vs barcode rank (knee/inflection plot) to decide threshold, but
        #      "this unnecessarily discards libraries derived from cell types with low RNA content" (OSCA, Ch. 6)
        #      -> Instead should prefer this Monte Carlo-simulation-based empty droplet test:
        # Additionally:
        # For any Sig==FALSE & Limited==TRUE, may need to increase n iterations (default = 10000) with 'niters='
        #   - this field = whether "the computed p-value for a...barcode is bounded by the number of iterations"
        
            # -> In exploratory phase (preprint), not all samples passed with niters=15000 (default 10,000), so use 20,000
        
        e.out.2 <- list()
        for(i in 1:length(pilot.data.2)){
          cat(paste0("Simulating empty drops for: ",names(pilot.data.2)[i],"... \n"))
          
          set.seed(109)
          e.out.2[[i]] <- emptyDrops(counts(pilot.data.2[[i]]), niters=20000)
          cat(paste0("\n\t...Simulations complete. \n\t", date(), "\n\n\n"))
          date()
        }
        
        names(e.out.2) <- names(pilot.data.2)
        
        save(pilot.data.2, e.out.2, file="rdas/revision/all-FACS-n10_2021rev_SCEs_processing-QC_MNTFeb2021.rda")
        
        #### ** END JOB - pick up interactive assessment, below ** ====

        
        

### (Interactive:) Read in data with `emptyDrops` stats =====
load("rdas/revision/all-FACS-n10_2021rev_SCEs_processing-QC_MNTMar2021.rda", verbose=T)
    # pilot.data.2, e.out.2

for(i in 1:length(e.out.2)){
  print(names(e.out.2)[[i]])
  print(table(Signif = e.out.2[[i]]$FDR <= 0.001, Limited = e.out.2[[i]]$Limited))
}
        # [1] "br5276.sacc.neun"
        #       Limited
        # Signif  FALSE TRUE
        #   FALSE  6649    0
        #   TRUE     26  900
        
        # [1] "br5400.nac"
        #       Limited
        # Signif  FALSE  TRUE
        #   FALSE 66248     0
        #   TRUE      0  4859
        
        # [1] "br5276.nac"
        #       Limited
        # Signif  FALSE  TRUE
        #   FALSE 42164     0
        #   TRUE      0  3012
        
        # [1] "br5701.nac.neun"
        #       Limited
        # Signif  FALSE TRUE
        #   FALSE  8280    0
        #   TRUE      0  814
        
        # [1] "br5701.sacc.neun"
        #       Limited
        # Signif  FALSE  TRUE
        #   FALSE 73689     0
        #   TRUE      0  4518
        
        # [1] "br5207.dlpfc"
        #       Limited
        # Signif  FALSE  TRUE
        #   FALSE 74753     0
        #   TRUE      0  6453
        
        # [1] "br5276.amy.neun"
        #       Limited
        # Signif  FALSE  TRUE
        #   FALSE 68568  1918     <- needs more iters, just do interactively (below)
        #   TRUE      1   816
        
        # [1] "br5400.amy.neun"
        #       Limited
        # Signif  FALSE  TRUE
        #   FALSE 55624     0
        #   TRUE      0  3185
        
        # [1] "br5400.sacc"
        #       Limited
        # Signif  FALSE  TRUE
        #   FALSE 67139     0
        #   TRUE      0  4364
        
        # [1] "br5701.amy"
        #       Limited
        # Signif  FALSE  TRUE
        #   FALSE 66271     0
        #   TRUE     25  4174
        #         - all are good and not lower-p-value-bound-limited?  All but 'br5276.amy.neun'
        
        # Increase 'niters' for br5276.amy.neun
        e.out.5276amy <- emptyDrops(counts(pilot.data.2[["br5276.amy.neun"]]), niters=30000)
        table(Signif = e.out.5276amy$FDR <= 0.001, Limited = e.out.5276amy$Limited)
                #       Limited
                #Signif  FALSE  TRUE
                # FALSE 68565     0
                # TRUE      1  2737     ok good.
                #         
        
        # Replace that entry in e.out.2
        e.out.2[["br5276.amy.neun"]] <- e.out.5276amy


# Subset in for-loop:
for(i in 1:length(pilot.data.2)){
  pilot.data.2[[i]] <- pilot.data.2[[i]][ ,which(e.out.2[[i]]$FDR <= 0.001)]
}
# Check
sapply(pilot.data.2, dim)


## Save this
save(pilot.data.2, e.out.2, file="rdas/revision/all-FACS-n10_2021rev_SCEs_processing-QC_MNTMar2021.rda")



### Mito rate QC ==================
table(rownames(pilot.data.2[[7]])==rownames(pilot.data.2[[10]]))  # and checked various other pairs

location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(pilot.data.2[[1]])$gene_id, 
                   column="SEQNAME", keytype="GENEID")
    #Warning message:
    #Unable to map 312 of 33538 requested IDs. - ok bc none of these are MT genes (13 pt-coding; `table(location)`)

# ID those mito genes
stats <- list()
for(i in 1:length(pilot.data.2)){
  stats[[i]] <- perCellQCMetrics(pilot.data.2[[i]], subsets=list(Mito=which(location=="MT")))
}
names(stats) <- names(pilot.data.2)


### Trick: Add a pseudo-count==1 for a 'MT transcript' ===
  # Note: This was implemented because we realized samples with mito rate distributions that
  #       were 'clean' and tightly distributed about 0 would yield a 3x MAD = 0, thus over-penalizing
  #       nuclei even if they had a single MT transcript (throwing out upwards of 50% of the sample)

# First check computation of mito percent:
table(stats[[8]]$subsets_Mito_percent == (stats[[8]]$subsets_Mito_sum/stats[[8]]$sum)*100)
    # All TRUE

test.stats <- stats
    
    for(i in 1:length(test.stats)){
      test.stats[[i]]$pseudo_subsets_Mito_sum <- test.stats[[i]]$subsets_Mito_sum + 1
      test.stats[[i]]$pseudo_subsets_Mito_percent <- test.stats[[i]]$pseudo_subsets_Mito_sum / (test.stats[[i]]$sum+1) * 100
    }

    ## Lapply: MAD approach for mito rate thresholding
    pseudo.high.mito <- lapply(test.stats, function(x) isOutlier(x$pseudo_subsets_Mito_percent, nmads=3, type="higher"))
    pseudo.high.mito.table <- lapply(pseudo.high.mito, table)
    # Percept dropped
    sapply(pseudo.high.mito.table, function(x) round(x[2]/sum(x), 3))
        #br5276.sacc.neun.TRUE       br5400.nac.TRUE       br5276.nac.TRUE 
        #                0.081                 0.155                 0.128 
        # br5701.nac.neun.TRUE br5701.sacc.neun.TRUE     br5207.dlpfc.TRUE 
        #                0.205                 0.158                 0.180 
        # br5276.amy.neun.TRUE  br5400.amy.neun.TRUE      br5400.sacc.TRUE 
        #                0.100                 0.173                 0.093 
        #      br5701.amy.TRUE 
        #                0.161
    
    # Thresholds
    sapply(pseudo.high.mito, function(x){round(attributes(x)[["thresholds"]]["higher"], 4)})
        #br5276.sacc.neun.higher       br5400.nac.higher       br5276.nac.higher 
        #                 5.6623                  0.0921                  1.3381 
        # br5701.nac.neun.higher br5701.sacc.neun.higher     br5207.dlpfc.higher 
        #                 3.5073                  0.0981                  0.0591 
        # br5276.amy.neun.higher  br5400.amy.neun.higher      br5400.sacc.higher 
        #                 5.2953                  0.6188                  0.0712 
        #      br5701.amy.higher 
        #                 0.2177 

    
## Bind [true] stats to colData
for(i in 1:length(pilot.data.2)){
  colData(pilot.data.2[[i]]) <- cbind(colData(pilot.data.2[[i]]), stats[[i]],
                                      #high.mito[[i]]
                                      pseudo.high.mito[[i]]
                                      )
  colnames(colData(pilot.data.2[[i]]))[9] <- "high.mito"
}

# $sum vs. $total ??
for(i in 1:length(pilot.data.2)){
  print(table(pilot.data.2[[i]]$sum == pilot.data.2[[i]]$total))
}
    ## all TRUE so can remove this 'duplicate' column:
    for(i in 1:length(pilot.data.2)){
      pilot.data.2[[i]]$total <- NULL
    }

# Store original for comparison/plotting
pilot.data.2.unfiltered <- pilot.data.2

## Subset - remove those indexed as high.mito
for(i in 1:length(pilot.data.2)){
  pilot.data.2[[i]] <- pilot.data.2[[i]][ ,!pilot.data.2[[i]]$high.mito]
}
sapply(pilot.data.2, dim)


## Plot metrics === ===

mitoCutoffs <- unlist(lapply(high.mito, function(x){attributes(x)$thresholds["higher"]}))
#mitoCutoffs <- unlist(lapply(pseudo.high.mito, function(x){attributes(x)$thresholds["higher"]}))
mean(mitoCutoffs)
    # [1] 1.453033;;     0.3903217 for first batch (n=12)
    ## with pseudo-MT count:
    # [1] 1.696016
median(mitoCutoffs)
    # [1] 0.1657229;;    0.138046 for first batch (n=12)
    ## with pseudo-MT count:
    # [1] 0.4182892
mitoCutoffs <- round(mitoCutoffs, 3)

#dir.create("pdfs/revision")
pdf("pdfs/revision/all-FACS-n10_2021rev_QCmetrics_high-mitoColored_MNT.pdf", height=4)
#pdf("pdfs/revision/all-FACS-n10_2021rev_QCmetrics_high-mitoColored_wPseudoMTcount_MNT.pdf", height=4)
for(i in 1:length(pilot.data.2.unfiltered)){
  grid.arrange(
    plotColData(pilot.data.2.unfiltered[[i]], y="sum", colour_by="high.mito") +
      scale_y_log10() + ggtitle(paste0("Total count: ", names(pilot.data.2.unfiltered)[[i]])),
    plotColData(pilot.data.2.unfiltered[[i]], y="detected", colour_by="high.mito") +
      scale_y_log10() + ggtitle("Detected features"),
    plotColData(pilot.data.2.unfiltered[[i]], y="subsets_Mito_percent",
                colour_by="high.mito") + ggtitle(paste0("Mito % (cutoff = ", mitoCutoffs[i],")")),
    ncol=3
  )
  # Mito rate vs n detected features
  print(
    plotColData(pilot.data.2.unfiltered[[i]], x="detected", y="subsets_Mito_percent",
                colour_by="high.mito", point_size=2.5, point_alpha=0.5) +
      ggtitle(paste0("Sample: ", names(pilot.data.2.unfiltered)[[i]],
                     ";   pre-QC nNuclei: ", ncol(pilot.data.2.unfiltered[[i]]),";    ",
                     "nNuclei kept: ", ncol(pilot.data.2[[i]])," (",
                     round(ncol(pilot.data.2[[i]]) / ncol(pilot.data.2.unfiltered[[i]]), 2), "%)"
      ))
  )
  # Detected features vs total count
  print(
    plotColData(pilot.data.2.unfiltered[[i]], x="sum", y="detected",
                colour_by="high.mito", point_size=2.5, point_alpha=0.5) +
      ggtitle(paste0("Sample: ", names(pilot.data.2.unfiltered)[[i]],
                     ";   pre-QC nNuclei: ", ncol(pilot.data.2.unfiltered[[i]]),";    ",
                     "nNuclei kept: ", ncol(pilot.data.2[[i]])," (",
                     round(ncol(pilot.data.2[[i]]) / ncol(pilot.data.2.unfiltered[[i]]), 2), "%)"
      ))
  )
}
dev.off()


## Save!
save(pilot.data.2, pilot.data.2.unfiltered, e.out.2, file="rdas/revision/all-FACS-n10_2021rev_SCEs_processing-QC_MNTMar2021.rda")



### Doublet detection / removal ==============================
  # Use default params, because this is at the single-sample-level
  # (within-region normalization, PCA, etc. will be performed with corresponding samples)

library(scDblFinder)

## To speed up, run on sample-level top-HVGs - just take top 1000 ===
pilot.data.normd <- lapply(pilot.data.2, logNormCounts)
geneVar.samples <- lapply(pilot.data.normd, modelGeneVar)
topHVGs <- lapply(geneVar.samples, function(x) {getTopHVGs(x, n=1000)})

# Generate doublet density scores
set.seed(109)
dbl.dens.focused <- lapply(names(pilot.data.normd), function(x) {
  computeDoubletDensity(pilot.data.normd[[x]], subset.row=topHVGs[[x]])})
names(dbl.dens.focused) <- names(pilot.data.normd)

sapply(dbl.dens.focused, function(x) round(quantile(x, probs=seq(0,1,by=0.05)),3))
    #      br5276.sacc.neun br5400.nac br5276.nac br5701.nac.neun br5701.sacc.neun
    # 0%              0.015      0.000      0.005           0.041            0.000
    # 5%              0.043      0.016      0.032           0.092            0.023
    # 10%             0.071      0.025      0.053           0.135            0.046
    # 15%             0.124      0.033      0.074           0.185            0.053
    # 20%             0.204      0.049      0.095           0.239            0.068
    # 25%             0.397      0.066      0.116           0.313            0.076
    # 30%             0.645      0.090      0.142           0.399            0.084
    # 35%             0.810      0.115      0.173           0.479            0.099
    # 40%             0.899      0.131      0.200           0.565            0.107
    # 45%             0.998      0.156      0.236           0.673            0.122
    # 50%             1.091      0.181      0.268           0.735            0.137
    # 55%             1.164      0.222      0.305           0.863            0.152
    # 60%             1.224      0.263      0.347           1.194            0.175
    # 65%             1.283      0.320      0.389           1.262            0.198
    # 70%             1.336      0.394      0.441           1.351            0.228
    # 75%             1.396      0.493      0.488           1.420            0.266
    # 80%             1.530      0.662      0.557           1.631            0.312
    # 85%             1.646      0.928      0.651           1.706            0.365
    # 90%             1.777      1.364      0.817           1.799            0.479
    # 95%             2.018      2.158      1.134           2.964            0.825
    # 100%            3.166     17.377     25.141           4.215           23.051
    
    #      br5207.dlpfc br5276.amy.neun br5400.amy.neun br5400.sacc br5701.amy
    # 0%          0.000           0.000           0.000       0.000      0.000
    # 5%          0.000           0.035           0.011       0.000      0.021
    # 10%         0.011           0.118           0.021       0.008      0.042
    # 15%         0.011           0.163           0.053       0.008      0.063
    # 20%         0.021           0.186           0.116       0.016      0.085
    # 25%         0.021           0.207           0.206       0.032      0.106
    # 30%         0.032           0.227           0.306       0.048      0.134
    # 35%         0.042           0.242           0.411       0.063      0.169
    # 40%         0.053           0.261           0.525       0.079      0.204
    # 45%         0.064           0.281           0.617       0.103      0.254
    # 50%         0.085           0.301           0.701       0.127      0.303
    # 55%         0.106           0.320           0.764       0.158      0.352
    # 60%         0.127           0.345           0.817       0.198      0.423
    # 65%         0.159           0.385           0.870       0.238      0.507
    # 70%         0.191           0.429           0.922       0.277      0.600
    # 75%         0.244           0.493           0.980       0.325      0.726
    # 80%         0.307           0.562           1.043       0.396      0.881
    # 85%         0.402           0.663           1.138       0.507      1.163
    # 90%         0.635           0.816           1.296       0.728      1.593
    # 95%         1.539           1.241           1.544       1.211      2.515
    # 100%       26.671           7.863           4.943       7.483     13.109

sapply(dbl.dens.focused, function(x) table(x >= 5))

# Percent that would be dropped at density score >= 5
round(sapply(names(dbl.dens.focused), function(x) {
  table(dbl.dens.focused[[x]] >= 5)["TRUE"] / ncol(pilot.data.2[[x]]) * 100
}), 3)
    #  br5276.sacc.neun.NA       br5400.nac.TRUE       br5276.nac.TRUE 
    #                   NA                 2.191                 1.028 
    #   br5701.nac.neun.NA br5701.sacc.neun.TRUE     br5207.dlpfc.TRUE 
    #                   NA                 1.367                 1.360 
    # br5276.amy.neun.TRUE    br5400.amy.neun.NA      br5400.sacc.TRUE 
    #                0.243                    NA                 0.177 
    #      br5701.amy.TRUE 
    #               2.043


    # --> Thresholding (this is arbitrary!) at a score >= 5 should be fair, but acknowledging
    #     there is no clear cut answer and some true doublets may remain in the dataset.
    #     -> see http://bioconductor.org/books/release/OSCA/doublet-detection.html#doublet-simulation
    
    #     Additionally: Will be good to just check downstream if higher scores are still associated
    #                   with any particular subcluster


# Add the doublet density scores to the colData
for(i in names(pilot.data.2)){
  pilot.data.2[[i]]$doubletScore <- dbl.dens.focused[[i]]
}
    # -> Will leave the thresholding at the region-specific level for flexibility
    #    since will save separate .rda for each of those

### Make/add some sample metadata ===
ref.sampleInfo <- data.frame(sampleID = names(pilot.data.2))
ref.sampleInfo$region <- ss(names(pilot.data.2),"\\.", 2)
ref.sampleInfo$donor <- ss(names(pilot.data.2),"\\.", 1)
ref.sampleInfo$sex <- ifelse(ref.sampleInfo$donor %in% c("br5400", "br5701"), "F", "M")

ref.sampleInfo$processBatch <- ifelse(ref.sampleInfo$sampleID %in% c("br5276.nac", "br5400.nac", "br5701.nac.neun",
                                                                     "br5276.sacc.neun", "br5701.sacc.neun"),
                                      "R6.10Feb2021", "R5.03Feb2021")


ref.sampleInfo$protocol <- "Frankenstein"
ref.sampleInfo$protocol[grep("neun", ref.sampleInfo$sampleID)] <- "Frank.NeuN.enriched"

ref.sampleInfo$sequencer <- "NovaSeq"

rownames(ref.sampleInfo) <- ref.sampleInfo$sampleID

## Add those to the colData:
for(i in names(pilot.data.2)){
  pilot.data.2[[i]]$sampleID <- i
  pilot.data.2[[i]]$region <- ref.sampleInfo[i, "region"]
  pilot.data.2[[i]]$donor <- ref.sampleInfo[i, "donor"]
  pilot.data.2[[i]]$sex <- ref.sampleInfo[i, "sex"]
  pilot.data.2[[i]]$processBatch <- ref.sampleInfo[i, "processBatch"]
  pilot.data.2[[i]]$protocol <- ref.sampleInfo[i, "protocol"]
  pilot.data.2[[i]]$sequencer <- ref.sampleInfo[i, "sequencer"]
}

## Save:
ref.sampleInfo.rev <- ref.sampleInfo
save(pilot.data.2, pilot.data.2.unfiltered, e.out.2, ref.sampleInfo.rev,
     file="rdas/revision/all-FACS-n10_2021rev_SCEs_processing-QC_MNTMar2021.rda")


        # === === === === === === === === === === ===
        # And end here -> proceed to 'step02' scripts
        # === === === === === === === === === === ===


sessionInfo()
### session info ====================================
# R version 4.0.4 RC (2021-02-08 r79975)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS:   /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/lib/libRblas.so
# LAPACK: /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/lib/libRlapack.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices datasets  utils    
# [8] methods   base     
# 
# other attached packages:
#   [1] gridExtra_2.3               Rtsne_0.15                 
# [3] jaffelab_0.99.30            rafalib_1.0.0              
# [5] DropletUtils_1.10.3         uwot_0.1.10                
# [7] Matrix_1.3-2                scran_1.18.5               
# [9] scater_1.18.6               ggplot2_3.3.3              
# [11] EnsDb.Hsapiens.v86_2.99.0   ensembldb_2.14.0           
# [13] AnnotationFilter_1.14.0     GenomicFeatures_1.42.3     
# [15] AnnotationDbi_1.52.0        batchelor_1.6.2            
# [17] scRNAseq_2.4.0              SingleCellExperiment_1.12.0
# [19] SummarizedExperiment_1.20.0 Biobase_2.50.0             
# [21] GenomicRanges_1.42.0        GenomeInfoDb_1.26.7        
# [23] IRanges_2.24.1              S4Vectors_0.28.1           
# [25] BiocGenerics_0.36.0         MatrixGenerics_1.2.1       
# [27] matrixStats_0.58.0         
# 
# loaded via a namespace (and not attached):
#   [1] AnnotationHub_2.22.0          BiocFileCache_1.14.0         
# [3] igraph_1.2.6                  lazyeval_0.2.2               
# [5] splines_4.0.4                 BiocParallel_1.24.1          
# [7] digest_0.6.27                 htmltools_0.5.1.1            
# [9] viridis_0.6.0                 fansi_0.4.2                  
# [11] magrittr_2.0.1                memoise_2.0.0                
# [13] limma_3.46.0                  Biostrings_2.58.0            
# [15] R.utils_2.10.1                askpass_1.1                  
# [17] prettyunits_1.1.1             colorspace_2.0-0             
# [19] blob_1.2.1                    rappdirs_0.3.3               
# [21] dplyr_1.0.5                   crayon_1.4.1                 
# [23] RCurl_1.98-1.3                glue_1.4.2                   
# [25] gtable_0.3.0                  zlibbioc_1.36.0              
# [27] XVector_0.30.0                DelayedArray_0.16.3          
# [29] BiocSingular_1.6.0            Rhdf5lib_1.12.1              
# [31] HDF5Array_1.18.1              scales_1.1.1                 
# [33] DBI_1.1.1                     edgeR_3.32.1                 
# [35] Rcpp_1.0.6                    viridisLite_0.4.0            
# [37] xtable_1.8-4                  progress_1.2.2               
# [39] dqrng_0.2.1                   bit_4.0.4                    
# [41] rsvd_1.0.3                    ResidualMatrix_1.0.0         
# [43] httr_1.4.2                    RColorBrewer_1.1-2           
# [45] ellipsis_0.3.1                pkgconfig_2.0.3              
# [47] XML_3.99-0.6                  R.methodsS3_1.8.1            
# [49] scuttle_1.0.4                 dbplyr_2.1.1                 
# [51] locfit_1.5-9.4                utf8_1.2.1                   
# [53] tidyselect_1.1.0              rlang_0.4.10                 
# [55] later_1.1.0.1                 munsell_0.5.0                
# [57] BiocVersion_3.12.0            tools_4.0.4                  
# [59] cachem_1.0.4                  generics_0.1.0               
# [61] RSQLite_2.2.6                 ExperimentHub_1.16.0         
# [63] stringr_1.4.0                 fastmap_1.1.0                
# [65] yaml_2.2.1                    bit64_4.0.5                  
# [67] purrr_0.3.4                   sparseMatrixStats_1.2.1      
# [69] mime_0.10                     R.oo_1.24.0                  
# [71] xml2_1.3.2                    biomaRt_2.46.3               
# [73] compiler_4.0.4                rstudioapi_0.13              
# [75] beeswarm_0.3.1                curl_4.3                     
# [77] interactiveDisplayBase_1.28.0 tibble_3.1.0                 
# [79] statmod_1.4.35                stringi_1.5.3                
# [81] lattice_0.20-41               bluster_1.0.0                
# [83] ProtGenerics_1.22.0           vctrs_0.3.6                  
# [85] pillar_1.6.0                  lifecycle_1.0.0              
# [87] rhdf5filters_1.2.0            BiocManager_1.30.12          
# [89] BiocNeighbors_1.8.2           bitops_1.0-6                 
# [91] irlba_2.3.3                   httpuv_1.5.5                 
# [93] rtracklayer_1.50.0            R6_2.5.0                     
# [95] promises_1.2.0.1              vipor_0.4.5                  
# [97] assertthat_0.2.1              rhdf5_2.34.0                 
# [99] openssl_1.4.3                 withr_2.4.1                  
# [101] GenomicAlignments_1.26.0      Rsamtools_2.6.0              
# [103] GenomeInfoDbData_1.2.4        hms_1.0.0                    
# [105] grid_4.0.4                    beachmat_2.6.4               
# [107] DelayedMatrixStats_1.12.3     googledrive_1.0.1            
# [109] segmented_1.3-3               shiny_1.6.0                  
# [111] ggbeeswarm_0.6.0

