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
        # Make sure works
        sapply(paths.rawCounts, list.files) # good
        
        # Make names for individual SCEs
        names(paths.rawCounts) <- gsub("_",".", tolower(samples.prepr))
        
        
        ## Read in raw UMI x barcode matrix - **use pre-mRNA-aligned reads
        pilot.data <- lapply(paths.rawCounts, function(x){ read10xCounts(x, col.names=TRUE) })
        names(pilot.data) <- names(paths.rawCounts)
        
        
        
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
        
        
        # Preliminary save
        save(pilot.data, file="rdas/revision/all-FACS-n14_preprint_SCEs_processing-QC_MNTMar2021.rda")
        
        
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
        
        save(pilot.data, e.out, file="rdas/revision/all-FACS-n14_preprint_SCEs_processing-QC_MNTMar2021.rda")
        
        #### ** END JOB - pick up interactive assessment, below ** ====


### (Interactive:) Read in data with `emptyDrops` stats =====
load("rdas/revision/all-FACS-n14_preprint_SCEs_processing-QC_MNTMar2021.rda", verbose=T)
    # pilot.data, e.out
        
for(i in 1:length(e.out)){
  print(names(e.out)[[i]])
  print(table(Signif = e.out[[i]]$FDR <= 0.001, Limited = e.out[[i]]$Limited))
  cat("\n")
}
        ##[1] "amy.5161"
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


# Save 
save(pilot.data, e.out, file="rdas/revision/all-FACS-n14_preprint_SCEs_processing-QC_MNTMar2021.rda")



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

## No subset.row iteration:
set.seed(109)
dbl.dens <- lapply(pilot.data, computeDoubletDensity)
sapply(dbl.dens, function(x) round(quantile(x, probs=seq(0,1,by=0.05)),3))
    #     br5161.amy br5161.dlpfc br5161.hpc br5161.nac br5161.sacc br5207.nac.neun br5212.amy
    # 0%        0.000        0.000      0.000      0.000       0.000           0.000      0.000
    # 5%        0.000        0.000      0.000      0.016       0.000           0.009      0.000
    # 10%       0.000        0.008      0.009      0.025       0.000           0.018      0.007
    # 15%       0.007        0.017      0.009      0.033       0.006           0.027      0.007
    # 20%       0.007        0.017      0.018      0.045       0.013           0.035      0.013
    # 25%       0.013        0.025      0.027      0.058       0.019           0.044      0.020
    # 30%       0.020        0.042      0.035      0.074       0.025           0.062      0.026
    # 35%       0.026        0.051      0.044      0.086       0.038           0.080      0.033
    # 40%       0.033        0.067      0.053      0.103       0.051           0.097      0.046
    # 45%       0.040        0.084      0.071      0.120       0.070           0.115      0.065
    # 50%       0.053        0.110      0.088      0.140       0.102           0.142      0.085
    # 55%       0.079        0.135      0.115      0.164       0.140           0.177      0.104
    # 60%       0.112        0.169      0.141      0.197       0.197           0.212      0.130
    # 65%       0.145        0.211      0.186      0.226       0.279           0.257      0.163
    # 70%       0.204        0.278      0.239      0.291       0.382           0.319      0.196
    # 75%       0.288        0.362      0.318      0.349       0.495           0.398      0.248
    # 80%       0.435        0.489      0.451      0.452       0.641           0.504      0.319
    # 85%       0.626        0.632      0.699      0.604       0.819           0.637      0.463
    # 90%       1.015        0.902      1.070      0.898       1.073           0.938      0.763
    # 95%       1.660        1.669      1.671      1.434       1.596           1.673      1.323
    # 100%     13.295       13.657     18.551      7.793      19.723          25.010     13.570

    #      br5212.dlpfc br5212.hpc br5212.nac br5212.sacc br5287.hpc br5287.nac br5182.nac.neun
    # 0%          0.000      0.000      0.004       0.000      0.000      0.014           0.000
    # 5%          0.003      0.000      0.027       0.000      0.007      0.065           0.009
    # 10%         0.010      0.008      0.035       0.000      0.015      0.091           0.017
    # 15%         0.014      0.008      0.043       0.000      0.022      0.123           0.034
    # 20%         0.020      0.016      0.050       0.000      0.030      0.136           0.043
    # 25%         0.027      0.024      0.064       0.008      0.045      0.148           0.060
    # 30%         0.037      0.024      0.077       0.008      0.056      0.170           0.068
    # 35%         0.058      0.032      0.096       0.016      0.067      0.210           0.085
    # 40%         0.081      0.048      0.113       0.031      0.079      0.262           0.102
    # 45%         0.108      0.056      0.131       0.039      0.090      0.320           0.128
    # 50%         0.135      0.080      0.156       0.062      0.108      0.394           0.145
    # 55%         0.179      0.095      0.181       0.085      0.134      0.471           0.179
    # 60%         0.240      0.119      0.213       0.109      0.168      0.605           0.205
    # 65%         0.297      0.151      0.248       0.147      0.217      0.722           0.247
    # 70%         0.352      0.199      0.285       0.202      0.280      0.913           0.290
    # 75%         0.413      0.262      0.330       0.264      0.363      1.034           0.367
    # 80%         0.500      0.382      0.379       0.357      0.480      1.293           0.462
    # 85%         0.609      0.549      0.465       0.489      0.718      1.547           0.606
    # 90%         0.835      0.902      0.662       0.768      1.260      2.072           0.862
    # 95%         1.366      1.433      0.994       1.312      1.661      3.310           1.544
    # 100%        6.765     14.596      7.844      11.756      5.965      4.349          21.376

## What about on sample-level top-HVGs?  Just take top 1000 === ===
pilot.data.normd <- lapply(pilot.data, logNormCounts)
geneVar.samples <- lapply(pilot.data.normd, modelGeneVar)
topHVGs <- lapply(geneVar.samples, function(x) {getTopHVGs(x, n=1000)})
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

# How do they correlate?
round(sapply(names(dbl.dens), function(x){cor(dbl.dens[[x]], dbl.dens.focused[[x]])}), 3)
    #      br5161.amy    br5161.dlpfc      br5161.hpc      br5161.nac     br5161.sacc 
    #           0.892           0.834           0.918           0.722           0.794 
    # br5207.nac.neun      br5212.amy    br5212.dlpfc      br5212.hpc      br5212.nac 
    #           0.878           0.850           0.912           0.922           0.834 
    #     br5212.sacc      br5287.hpc      br5287.nac br5182.nac.neun 
    #           0.896           0.710           0.747           0.933 

sapply(dbl.dens, function(x) table(x >= 10))
sapply(dbl.dens.focused, function(x) table(x >= 10))
sapply(dbl.dens, function(x) table(x >= 5))
sapply(dbl.dens.focused, function(x) table(x >= 5))

    # Observation: seems that more NeuN-selected-sample nuclei are flagged with higher max scores,
    #              (and slightly more) which is unsurprising since they are more 'homogeneous'

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


# Add the doublet density scores to the colData
for(i in names(pilot.data)){
  pilot.data[[i]]$doubletScore <- dbl.dens.focused[[i]]
}

## Post-hoc: that there aren't any potential cell types associated with high scores === ===
new.hpc <- cbind(pilot.data[["br5161.hpc"]], pilot.data[["br5212.hpc"]], pilot.data[["br5287.hpc"]])


# Bring in preprint SCE for HPC:
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda",
     verbose=T)


# set a sampleID; make match style of that SCE:
new.hpc$sampleID <- c(rep("hpc.5161", ncol(pilot.data[["br5161.hpc"]])),
                    rep("hpc.5212", ncol(pilot.data[["br5212.hpc"]])),
                    rep("hpc.5287", ncol(pilot.data[["br5287.hpc"]])))

new.hpc$sampleBC <- paste0(new.hpc$sampleID, ".", new.hpc$Barcode)

# Do same for preprint SCE
sce.hpc$sampleBC <- paste0(sce.hpc$sample, ".", sce.hpc$Barcode)

length(intersect(sce.hpc$sampleBC, new.hpc$sampleBC))
    # 10213 out of 10268 for the 'new' iteration and 10444 for the preprint

intersectingBCs <- intersect(sce.hpc$sampleBC, new.hpc$sampleBC)
sce.hpc.sub <- sce.hpc[ ,sce.hpc$sampleBC %in% intersectingBCs]
new.hpc.sub <- new.hpc[ ,new.hpc$sampleBC %in% intersectingBCs]
table(colnames(sce.hpc.sub) == colnames(new.hpc.sub)) # all TRUE
sce.hpc.sub$doubletScore <- new.hpc.sub$doubletScore

table(sce.hpc$cellType.split)
table(sce.hpc.sub$cellType.split)
    ## Looks likes cell populations are mostly preserved--esp in the smaller subpops (Tcell, Inhib.1, etc.)

cellType.idx <- splitit(sce.hpc.sub$cellType.split)
sapply(cellType.idx, function(x){round(quantile(sce.hpc.sub[ ,x]$doubletScore), 3)})
    #     Ambig.lowNtrxts  Astro Excit.1 Excit.2 Excit.3 Excit.4 Excit.5 Inhib.1
    # 0%             0.032  0.000   0.026   0.018   0.016   0.414   0.477   1.167
    # 25%            0.043  0.062   0.052   0.071   0.087   0.479   0.813   1.441
    # 50%            0.057  0.120   0.097   0.124   0.151   0.497   0.920   1.490
    # 75%            0.692  0.231   0.608   0.760   0.207   0.678   0.981   1.539
    # 100%           1.361 16.043   4.080   3.210   7.739   2.087   6.499   2.962
    #      Inhib.2 Inhib.3 Inhib.4 Inhib.5  Micro Oligo    OPC Tcell
    # 0%     0.663   0.027   0.274   0.592  0.000 0.000  0.000 0.048
    # 25%    1.177   0.115   0.924   0.935  0.032 0.199  0.041 0.064
    # 50%    1.291   0.705   1.017   1.241  0.072 0.509  0.120 0.097
    # 75%    1.415   1.362   1.122   1.386  0.159 1.100  0.230 0.112
    # 100%   6.984   6.928   1.848   1.795 15.580 8.037 12.034 0.124
sapply(cellType.idx, function(x){table(sce.hpc.sub[ ,x]$doubletScore >= 5)["TRUE"]})
    # Ambig.lowNtrxts.NA         Astro.TRUE         Excit.1.NA         Excit.2.NA 
    #                 NA                 35                 NA                 NA 
    #       Excit.3.TRUE         Excit.4.NA       Excit.5.TRUE         Inhib.1.NA 
    #                 10                 NA                  1                 NA 
    #       Inhib.2.TRUE       Inhib.3.TRUE         Inhib.4.NA         Inhib.5.NA 
    #                  1                  1                 NA                 NA 
    #         Micro.TRUE         Oligo.TRUE           OPC.TRUE           Tcell.NA 
    #                 51                 22                 25                 NA 

    ### Conclusion: thresholding at 5 seems appropriate (e.g. = 10 doesn't flag any Oligos, here)
      #             and will be good to just check downstream if higher scores are still associated
      #             with any particular subcluster

    # === === === === === === === === === === ===
    # And end here -> proceed to 'step02' scripts
    # === === === === === === === === === === ===



