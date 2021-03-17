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



### Mito rate QC ===
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


## Lapply: MAD approach for mito rate thresholding
high.mito <- lapply(stats, function(x) isOutlier(x$subsets_Mito_percent, nmads=3, type="higher"))
high.mito.table <- lapply(high.mito, table)
# Percent dropped
sapply(high.mito.table, function(x) round(x[2]/sum(x), 3))  # about 10% across all (besides weird 5212dlpfc, ~50%)
    # br5161.amy.TRUE    br5161.dlpfc.TRUE      br5161.hpc.TRUE      br5161.nac.TRUE 
    #           0.107                0.082                0.104                0.099 
    # br5161.sacc.TRUE br5207.nac.neun.TRUE      br5212.amy.TRUE    br5212.dlpfc.TRUE 
    #           0.101                0.055                0.087                0.499 
    # br5212.hpc.TRUE      br5212.nac.TRUE     br5212.sacc.TRUE      br5287.hpc.TRUE 
    #           0.110                0.067                0.122                0.138 
    # br5287.nac.TRUE br5182.nac.neun.TRUE 
    #           0.082                0.073

# Variable thresholds
sapply(high.mito, function(x){round(attributes(x)[["thresholds"]]["higher"], 4)})
    #  br5161.amy.higher    br5161.dlpfc.higher      br5161.hpc.higher      br5161.nac.higher 
    #             0.2841                 0.0966                 0.2119                 0.1314 
    # br5161.sacc.higher br5207.nac.neun.higher      br5212.amy.higher    br5212.dlpfc.higher 
    #             0.1447                 0.1116                 0.1001                 0.0000 
    #  br5212.hpc.higher      br5212.nac.higher     br5212.sacc.higher      br5287.hpc.higher 
    #             0.1751                 3.1776                 0.0414                 0.2113 
    #  br5287.nac.higher br5182.nac.neun.higher 
    #             0.1096                 0.3087


## Trick: Add a pseudo-count==1 for a 'MT transcript' ===
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

# Bind stats to colData
for(i in 1:length(pilot.data)){
  colData(pilot.data[[i]]) <- cbind(colData(pilot.data[[i]]), stats[[i]],
                                    #high.mito[[i]]
                                    pseudo.high.mito[[i]]
                                    )
  colnames(colData(pilot.data[[i]]))[9] <- "high.mito"
}

# $sum vs. $total ??
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

#mitoCutoffs <- unlist(lapply(high.mito, function(x){attributes(x)$thresholds["higher"]}))
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
#save(pilot.data, pilot.data.unfiltered, e.out, file="rdas/revision/all-FACS-n14_preprint_SCEs_processing-QC_MNTMar2021.rda")


    
    # === === === === === === === === === === ===
    # And end here -> proceed to 'step02' scripts
    # === === === === === === === === === === ===



