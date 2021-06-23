################################################################################
### LIBD pilot 10x-Frankenstein (n=12) snRNA-seq samples
### STEP 01: Read in SCEs and perform nuclei calling and QC
### Initiated: MNT 29Jan2020
### Modified: MNT 03Mar2021
### Intention: To generate/have a streamlined, easy-to-follow pipeline
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
        
        ### Read in preprint 'samples.manifest.full' for streamlining
        samples.prepr <- read.table("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/FASTQ/samples.manifest.full",
                                       sep="\t", header=F)$V5
        # Drop Br5287-DLPFC (poor quality sample; dropped for preprint) and the test sucrose samples
        samples.prepr <- samples.prepr[-c(grep("Br5287_DLPFC", samples.prepr),
                                          grep("_suc", samples.prepr))]
        # Add '_NeuN' suffix to the two NeuN sorts
        samples.prepr[c(6,14)] <- paste0(samples.prepr[c(6,14)],"_NeuN")

        # Make list of paths
        paths.rawCounts <- c(paste0("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/",
                                    samples.prepr,"/outs/raw_feature_bc_matrix"))
        path.alt.5182nac <- "/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/Br5182_NAc_NeuN_reseq/outs/raw_feature_bc_matrix"
        # Make sure works
        sapply(paths.rawCounts, list.files) # good
        
        # Make names for individual SCEs
        names(paths.rawCounts) <- gsub("_",".", tolower(samples.prepr))
        
        
        ## Read in raw UMI x barcode matrix - **use pre-mRNA-aligned reads
        pilot.data <- lapply(paths.rawCounts, function(x){ read10xCounts(x, col.names=TRUE) })
        names(pilot.data) <- names(paths.rawCounts)
        
            # MNT 21Jun: 'pilot.data.alt' will be the stand-in for the re-sequenced Br5182-NAc sample,
            #            which we only realized in revision that had incorrectly been adapter-trimmed,
            #            scrambling all the R1 lengths
            br5182.nac.neun <- read10xCounts(path.alt.5182nac, col.names=TRUE)
        
        ### Gene annotation (from scater) ===
        # Pull in GTF information
        gtf = import("/dcl01/ajaffe/data/lab/singleCell/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf")
        gtf = gtf[gtf$type == "gene"]
        length(gtf)
        names(gtf) = gtf$gene_id
        table(names(gtf) == rowData(pilot.data[[1]])$ID)
        
        seqlevels(gtf)[1:25] = paste0("chr", seqlevels(gtf)[1:25])
        mcols(gtf) = mcols(gtf)[,c(5:9)]
        
        
        for(i in 1:length(pilot.data)){
          rowRanges(pilot.data[[i]]) <- gtf
          # Because some gene names are the same:
          rownames(pilot.data[[i]]) <- uniquifyFeatureNames(rowData(pilot.data[[i]])$gene_id, rowData(pilot.data[[i]])$gene_name)
          rowData(pilot.data[[i]])$Symbol.uniq <- rownames(pilot.data[[i]])
        }
        
            # And for re-seq'd 5182-NAc:
            rowRanges(br5182.nac.neun) <- gtf
            # Because some gene names are the same:
            rownames(br5182.nac.neun) <- uniquifyFeatureNames(rowData(br5182.nac.neun)$gene_id, rowData(br5182.nac.neun)$gene_name)
            rowData(br5182.nac.neun)$Symbol.uniq <- rownames(br5182.nac.neun)
        
        # store into List, just to reflect the other objects
        pilot.data.alt <- list(br5182.nac.neun)
        names(pilot.data.alt) <- "br5182.nac.neun"
        
        # Preliminary save
        save(pilot.data, pilot.data.alt,
             # Created/stored in this .rda downstream:
             pilot.data.unfiltered, e.out, ref.sampleInfo,
             file="rdas/revision/all-FACS-n14_preprint_SCEs_processing-QC_MNTMar2021.rda")
        
        
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
        
        e.out <- list()
        for(i in 1:length(pilot.data)){
          cat(paste0("Simulating empty drops for: ",names(pilot.data)[i],"... \n"))
          
          set.seed(109)
          e.out[[i]] <- emptyDrops(counts(pilot.data[[i]]), niters=20000)
          cat(paste0("\n\t...Simulations complete. \n\t", date(), "\n\n\n"))
          date()
        }
        
        names(e.out) <- names(pilot.data)

# This one done interactively:        
e.out.alt <- list(br5182.nac.neun = emptyDrops(counts(pilot.data.alt[["br5182.nac.neun"]]), niters=20000))
  
  save(pilot.data, pilot.data.alt, e.out, e.out.alt,
       # Created/stored in this .rda downstream:
       pilot.data.unfiltered, ref.sampleInfo,file="rdas/revision/all-FACS-n14_preprint_SCEs_processing-QC_MNTMar2021.rda")
  
        #### ** END JOB - pick up interactive assessment, below ** ====


### (Interactive:) Read in data with `emptyDrops` stats =====
load("rdas/revision/all-FACS-n14_preprint_SCEs_processing-QC_MNTMar2021.rda", verbose=T)
    # pilot.data, e.out
        
for(i in 1:length(e.out)){
  print(names(e.out)[[i]])
  print(table(Signif = e.out[[i]]$FDR <= 0.001, Limited = e.out[[i]]$Limited))
  cat("\n")
}
        ##[1] "br5161.amy"
        #       Limited
        # Signif  FALSE TRUE
        # FALSE  1019    0
        # TRUE    150 3537
        # 
        # [1] "br5161.dlpfc"
        #       Limited
        # Signif  FALSE TRUE
        #   FALSE  1229    0
        #   TRUE     75 4713
        # 
        # [1] "br5161.hpc"
        #       Limited
        # Signif  FALSE TRUE
        #   FALSE   735    0
        #   TRUE     91 4951
        # 
        # [1] "br5161.nac"
        #       Limited
        # Signif  FALSE TRUE
        #   FALSE    41    0
        #   TRUE     14 2280
        # 
        # [1] "br5161.sacc"
        #       Limited
        # Signif  FALSE TRUE
        #   FALSE  3762    0
        #   TRUE     93 3454
        # 
        # [1] "br5207.nac.neun"
        #       Limited
        # Signif  FALSE  TRUE
        #   FALSE 18112     0
        #   TRUE     10  4672
        # 
        # [1] "br5212.amy"
        #       Limited
        # Signif  FALSE  TRUE
        #   FALSE 16780     0
        #   TRUE     30  3631
        # 
        # [1] "br5212.dlpfc"
        #       Limited
        # Signif  FALSE TRUE
        #   FALSE   755    0
        #   TRUE     35 1971
        # 
        # [1] "br5212.hpc"
        #       Limited
        # Signif  FALSE TRUE
        #   FALSE   358    0
        #   TRUE     59 4479
        # 
        # [1] "br5212.nac"
        #       Limited
        # Signif  FALSE TRUE
        #   FALSE   456    0
        #   TRUE     37 1864
        # 
        # [1] "br5212.sacc"
        #       Limited
        # Signif  FALSE  TRUE
        #   FALSE 55255     0
        #   TRUE      0  4398
        # 
        # [1] "br5287.hpc"
        #       Limited
        # Signif  FALSE TRUE
        #   FALSE   325    0
        #   TRUE     35 2154
        # 
        # [1] "br5287.nac"
        #       Limited
        # Signif  FALSE TRUE
        #   FALSE   116    0
        #   TRUE     18  752
        # 
        # [1] "br5182.nac.neun"
        #       Limited
        # Signif  FALSE TRUE
        #   FALSE  3552    0
        #   TRUE     23 4579        - all are good and not lower-p-value-bound-limited


# Subset in for-loop:
for(i in 1:length(pilot.data)){
  pilot.data[[i]] <- pilot.data[[i]][ ,which(e.out[[i]]$FDR <= 0.001)]
}
# Check
sapply(pilot.data, dim)

# For re-sequenced sample:
    i <- "br5182.nac.neun"
    print(table(Signif = e.out.alt[[i]]$FDR <= 0.001, Limited = e.out.alt[[i]]$Limited))
        #       Limited
        # Signif  FALSE TRUE
        #   FALSE  6506    0
        #   TRUE     25 4560
    pilot.data.alt[[i]] <- pilot.data.alt[[i]][ ,which(e.out.alt[[i]]$FDR <= 0.001)]

# Save 
save(pilot.data, pilot.data.alt, e.out, e.out.alt,
     # Created/stored in this .rda downstream:
     pilot.data.unfiltered, ref.sampleInfo,
     file="rdas/revision/all-FACS-n14_preprint_SCEs_processing-QC_MNTMar2021.rda")



### Mito rate QC ==================
table(rownames(pilot.data[[1]])==rownames(pilot.data[[6]]))  # and checked various other pairs

location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(pilot.data[[1]])$gene_id, 
                   column="SEQNAME", keytype="GENEID")
    #Warning message:
    #Unable to map 312 of 33538 requested IDs. - ok bc none of these are MT genes (13 pt-coding; `table(location)`)

# ID those mito genes
stats <- list()
for(i in 1:length(pilot.data)){
  stats[[i]] <- perCellQCMetrics(pilot.data[[i]], subsets=list(Mito=which(location=="MT")))
}
names(stats) <- names(pilot.data)


### Trick: Add a pseudo-count==1 for a 'MT transcript' ===
  # Note: This was implemented because we realized samples with mito rate distributions that
  #       were 'clean' and tightly distributed about 0 would yield a 3x MAD = 0, thus over-penalizing
  #       nuclei even if they had a single MT transcript (throwing out upwards of 50% of the sample)

# First check computation of mito percent:
table(stats[[13]]$subsets_Mito_percent == (stats[[13]]$subsets_Mito_sum/stats[[13]]$sum)*100)
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
        #  br5161.amy.TRUE    br5161.dlpfc.TRUE      br5161.hpc.TRUE      br5161.nac.TRUE 
        #            0.107                0.120                0.123                0.104 
        # br5161.sacc.TRUE br5207.nac.neun.TRUE      br5212.amy.TRUE    br5212.dlpfc.TRUE 
        #            0.105                0.055                0.110                0.156 
        #  br5212.hpc.TRUE      br5212.nac.TRUE     br5212.sacc.TRUE      br5287.hpc.TRUE 
        #            0.124                0.067                0.118                0.146 
        #  br5287.nac.TRUE br5182.nac.neun.TRUE 
        #            0.116                0.073
    
    # Thresholds
    sapply(pseudo.high.mito, function(x){round(attributes(x)[["thresholds"]]["higher"], 4)})
        #  br5161.amy.higher    br5161.dlpfc.higher      br5161.hpc.higher      br5161.nac.higher 
        #             0.3214                 0.1186                 0.2361                 0.1500 
        # br5161.sacc.higher br5207.nac.neun.higher      br5212.amy.higher    br5212.dlpfc.higher 
        #             0.1688                 0.1187                 0.1110                 0.0607 
        #  br5212.hpc.higher      br5212.nac.higher     br5212.sacc.higher      br5287.hpc.higher 
        #             0.2072                 3.2164                 0.0704                 0.2527 
        #  br5287.nac.higher br5182.nac.neun.higher 
        #             0.1160                 0.3161

# Bind [true] stats to colData
for(i in 1:length(pilot.data)){
  colData(pilot.data[[i]]) <- cbind(colData(pilot.data[[i]]), stats[[i]],
                                    #high.mito[[i]]
                                    pseudo.high.mito[[i]]
                                    )
  colnames(colData(pilot.data[[i]]))[9] <- "high.mito"
}

# $sum vs. $total ?
for(i in 1:length(pilot.data)){
  print(table(pilot.data[[i]]$sum == pilot.data[[i]]$total))
}
    ## all TRUE so can remove this second column:
    for(i in 1:length(pilot.data)){
      pilot.data[[i]]$total <- NULL
    }

# Store original for comparison/plotting
pilot.data.unfiltered <- pilot.data

## Subset - remove those indexed as high.mito
for(i in 1:length(pilot.data)){
  pilot.data[[i]] <- pilot.data[[i]][ ,!pilot.data[[i]]$high.mito]
}
sapply(pilot.data, dim)


## Plot metrics === ===

mitoCutoffs <- unlist(lapply(pseudo.high.mito, function(x){attributes(x)$thresholds["higher"]}))
mean(mitoCutoffs)
    #[1] 0.3903217 (n=12);  [1] 0.364583 (n=14)
    ## With pseudo-MT count:
    # [1] 0.3902808
median(mitoCutoffs)
    #[1] 0.138046 (n=12);   [1] 0.138046 (n=14)
    ## With pseudo-MT count:
    # [1] 0.1594075
mitoCutoffs <- round(mitoCutoffs, 3)

#dir.create("pdfs/")
#pdf("pdfs/revision/all-FACS-n14_preprint_QCmetrics_high-mitoColored_MNT.pdf", height=4)
pdf("pdfs/revision/all-FACS-n14_preprint_QCmetrics_high-mitoColored_wPseudoMTcount_MNT.pdf", height=4)
for(i in 1:length(pilot.data.unfiltered)){
  grid.arrange(
    plotColData(pilot.data.unfiltered[[i]], y="sum", colour_by="high.mito") +
      scale_y_log10() + ggtitle(paste0("Total count: ", names(pilot.data.unfiltered)[[i]])),
    plotColData(pilot.data.unfiltered[[i]], y="detected", colour_by="high.mito") +
      scale_y_log10() + ggtitle("Detected features"),
    plotColData(pilot.data.unfiltered[[i]], y="subsets_Mito_percent",
                colour_by="high.mito") + ggtitle(paste0("Mito % (cutoff = ", mitoCutoffs[i],")")),
    ncol=3
  )
  # Mito rate vs n detected features
  print(
    plotColData(pilot.data.unfiltered[[i]], x="detected", y="subsets_Mito_percent",
                colour_by="high.mito", point_size=2.5, point_alpha=0.5) +
      ggtitle(paste0("Sample: ", names(pilot.data.unfiltered)[[i]],
                     ";   pre-QC nNuclei: ", ncol(pilot.data.unfiltered[[i]]),";    ",
                     "nNuclei kept: ", ncol(pilot.data[[i]])," (",
                     round(ncol(pilot.data[[i]]) / ncol(pilot.data.unfiltered[[i]]), 2), "%)"
      ))
  )
  # Detected features vs total count
  print(
    plotColData(pilot.data.unfiltered[[i]], x="sum", y="detected",
                colour_by="high.mito", point_size=2.5, point_alpha=0.5) +
      ggtitle(paste0("Sample: ", names(pilot.data.unfiltered)[[i]],
                     ";   pre-QC nNuclei: ", ncol(pilot.data.unfiltered[[i]]),";    ",
                     "nNuclei kept: ", ncol(pilot.data[[i]])," (",
                     round(ncol(pilot.data[[i]]) / ncol(pilot.data.unfiltered[[i]]), 2), "%)"
      ))
  )
}
dev.off()


## Save!
save(pilot.data, pilot.data.unfiltered, e.out,
     file="rdas/revision/all-FACS-n14_preprint_SCEs_processing-QC_MNTMar2021.rda")



### Doublet detection / removal ==============================
  # Use default params, because this is at the single-sample-level
  # (within-region normalization, PCA, etc. will be performed with corresponding samples)

#BiocManager::install("scDblFinder")
library(scDblFinder)

## To speed up, run on sample-level top-HVGs - just take top 1000 ===
pilot.data.normd <- lapply(pilot.data, logNormCounts)
geneVar.samples <- lapply(pilot.data.normd, modelGeneVar)
topHVGs <- lapply(geneVar.samples, function(x) {getTopHVGs(x, n=1000)})

# Generate doublet density scores
set.seed(109)
dbl.dens.focused <- lapply(names(pilot.data.normd), function(x) {
  computeDoubletDensity(pilot.data.normd[[x]], subset.row=topHVGs[[x]])})
names(dbl.dens.focused) <- names(pilot.data.normd)

sapply(dbl.dens.focused, function(x) round(quantile(x, probs=seq(0,1,by=0.05)),3))
    #     br5161.amy br5161.dlpfc br5161.hpc br5161.nac br5161.sacc br5207.nac.neun br5212.amy
    # 0%        0.000        0.000      0.000      0.004       0.000           0.000      0.000
    # 5%        0.013        0.017      0.018      0.029       0.013           0.027      0.007
    # 10%       0.033        0.034      0.035      0.055       0.019           0.053      0.020
    # 15%       0.053        0.051      0.053      0.086       0.032           0.080      0.026
    # 20%       0.072        0.059      0.071      0.111       0.044           0.106      0.046
    # 25%       0.099        0.076      0.088      0.144       0.057           0.150      0.059
    # 30%       0.125        0.101      0.115      0.177       0.083           0.195      0.078
    # 35%       0.152        0.118      0.150      0.222       0.102           0.239      0.098
    # 40%       0.184        0.160      0.177      0.263       0.127           0.292      0.124
    # 45%       0.231        0.202      0.221      0.292       0.165           0.354      0.156
    # 50%       0.277        0.253      0.265      0.345       0.209           0.416      0.196
    # 55%       0.329        0.312      0.318      0.407       0.267           0.487      0.235
    # 60%       0.389        0.396      0.380      0.497       0.330           0.575      0.293
    # 65%       0.468        0.481      0.460      0.600       0.413           0.673      0.363
    # 70%       0.573        0.607      0.592      0.699       0.521           0.779      0.460
    # 75%       0.698        0.742      0.778      0.851       0.678           0.903      0.606
    # 80%       0.896        0.927      0.973      1.037       0.876           1.062      0.763
    # 85%       1.219        1.197      1.211      1.261       1.156           1.274      1.023
    # 90%       1.456        1.568      1.539      1.575       1.479           1.611      1.349
    # 95%       2.108        2.467      2.290      2.205       2.012           2.434      1.851
    # 100%     12.148       15.781     15.606      7.028      20.326          24.444     17.286

    #      br5212.dlpfc br5212.hpc br5212.nac br5212.sacc br5287.hpc br5287.nac br5182.nac.neun
    # 0%          0.000      0.000      0.014       0.000      0.000      0.039           0.000
    # 5%          0.014      0.016      0.043       0.000      0.030      0.151           0.026
    # 10%         0.027      0.032      0.064       0.008      0.045      0.165           0.043
    # 15%         0.044      0.048      0.082       0.008      0.067      0.259           0.068
    # 20%         0.062      0.072      0.110       0.023      0.094      0.426           0.094
    # 25%         0.081      0.087      0.131       0.031      0.120      0.541           0.128
    # 30%         0.095      0.111      0.167       0.039      0.142      0.639           0.162
    # 35%         0.115      0.135      0.199       0.054      0.180      0.756           0.205
    # 40%         0.132      0.159      0.238       0.062      0.208      0.812           0.247
    # 45%         0.157      0.191      0.269       0.078      0.254      0.903           0.299
    # 50%         0.200      0.239      0.298       0.101      0.310      0.982           0.358
    # 55%         0.234      0.302      0.337       0.124      0.389      1.076           0.418
    # 60%         0.271      0.390      0.383       0.155      0.494      1.186           0.478
    # 65%         0.311      0.493      0.426       0.186      0.606      1.320           0.554
    # 70%         0.349      0.636      0.475       0.225      0.738      1.427           0.640
    # 75%         0.416      0.787      0.564       0.303      0.930      1.494           0.734
    # 80%         0.507      0.986      0.684       0.396      1.165      1.683           0.853
    # 85%         0.630      1.265      0.823       0.560      1.397      1.870           1.015
    # 90%         0.863      1.618      1.069       0.877      1.638      2.065           1.285
    # 95%         1.312      2.416      1.450       1.428      2.057      2.296           1.868
    # 100%        6.454     16.043      7.886       9.584      7.749      4.154          21.564

sapply(dbl.dens.focused, function(x) table(x >= 5))

# Percent that would be dropped at density score >= 5
round(sapply(names(dbl.dens.focused), function(x) {
  table(dbl.dens.focused[[x]] >= 5)["TRUE"] / ncol(pilot.data[[x]]) * 100
  }), 3)
    # br5161.amy.TRUE    br5161.dlpfc.TRUE      br5161.hpc.TRUE      br5161.nac.TRUE 
    #           1.154                1.447                1.448                0.487 
    #br5161.sacc.TRUE br5207.nac.neun.TRUE      br5212.amy.TRUE    br5212.dlpfc.TRUE 
    #           0.693                2.192                1.043                0.118 
    # br5212.hpc.TRUE      br5212.nac.TRUE     br5212.sacc.TRUE      br5287.hpc.TRUE 
    #           2.012                0.395                0.412                0.214 
    #   br5287.nac.NA br5182.nac.neun.TRUE 
    #              NA                2.251 

    # --> Thresholding (this is arbitrary!) at a score >= 5 should be fair, but acknowledging
    #     there is no clear cut answer and some true doublets may remain in the dataset.
    #     -> see http://bioconductor.org/books/release/OSCA/doublet-detection.html#doublet-simulation

    #     Additionally: Will be good to just check downstream if higher scores are still associated
    #                   with any particular subcluster


# Add the doublet density scores to the colData
for(i in names(pilot.data)){
  pilot.data[[i]]$doubletScore <- dbl.dens.focused[[i]]
}
    # -> Will leave the thresholding at the region-specific level for flexibility
    #    since will save separate .rda for each of those

### Make/add some sample metadata ===
ref.sampleInfo <- data.frame(sampleID = names(pilot.data))
ref.sampleInfo$region <- ss(names(pilot.data),"\\.", 2)
ref.sampleInfo$donor <- ss(names(pilot.data),"\\.", 1)
ref.sampleInfo$sex <- ifelse(ref.sampleInfo$donor %in% c("br5400", "br5701"), "F", "M")

ref.sampleInfo$processBatch <- ifelse(ref.sampleInfo$sampleID %in% c("br5212.dlpfc", "br5287.hpc", "br5161.nac", "br5212.nac"),
                                   "R2.23Jul2019", "R1.08May2019")
ref.sampleInfo$processBatch <- ifelse(ref.sampleInfo$sampleID %in% c("br5161.dlpfc", "br5212.hpc", "br5287.nac", "br5161.amy"),
                                   "R3.04Sep2019", ref.sampleInfo$processBatch)
ref.sampleInfo$processBatch <- ifelse(ref.sampleInfo$sampleID %in% c("br5212.amy", "br5161.sacc", "br5212.sacc", "br5182.nac.neun", "br5207.nac.neun"),
                                   "R4.25Sep2019", ref.sampleInfo$processBatch)

ref.sampleInfo$protocol <- ifelse(ref.sampleInfo$processBatch=="R2.23Jul2019", "pseudoSort", "Frankenstein")
ref.sampleInfo$protocol <- ifelse(ref.sampleInfo$sampleID %in% c("br5182.nac.neun", "br5207.nac.neun"),
                                  "Frank.NeuN.enriched", ref.sampleInfo$protocol)

ref.sampleInfo$sequencer <- "NextSeq"

rownames(ref.sampleInfo) <- ref.sampleInfo$sampleID

## Add those to the colData:
for(i in names(pilot.data)){
  pilot.data[[i]]$sampleID <- i
  pilot.data[[i]]$region <- ref.sampleInfo[i, "region"]
  pilot.data[[i]]$donor <- ref.sampleInfo[i, "donor"]
  pilot.data[[i]]$sex <- ref.sampleInfo[i, "sex"]
  pilot.data[[i]]$processBatch <- ref.sampleInfo[i, "processBatch"]
  pilot.data[[i]]$protocol <- ref.sampleInfo[i, "protocol"]
  pilot.data[[i]]$sequencer <- ref.sampleInfo[i, "sequencer"]
}

## Save:
save(pilot.data, pilot.data.unfiltered, e.out, ref.sampleInfo,
     file="rdas/revision/all-FACS-n14_preprint_SCEs_processing-QC_MNTMar2021.rda")

    # === === === === === === === === === === ===
    # And end here -> proceed to 'step02' scripts
    # === === === === === === === === === === ===


sessionInfo()
## session info =================================
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

