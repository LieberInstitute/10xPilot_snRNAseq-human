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
        samples.revision <- read.table("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/FASTQ/Feb2021/samples.manifest",
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



### Mito rate QC ===
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


## Lapply: MAD approach for mito rate thresholding
high.mito <- lapply(stats, function(x) isOutlier(x$subsets_Mito_percent, nmads=3, type="higher"))
high.mito.table <- lapply(high.mito, table)
# Percept dropped
sapply(high.mito.table, function(x) round(x[2]/sum(x), 3))  # ~15-20% across all
  # br5276.sacc.neun.TRUE       br5400.nac.TRUE       br5276.nac.TRUE 
  #                 0.091                 0.126                 0.124 
  #  br5701.nac.neun.TRUE br5701.sacc.neun.TRUE     br5207.dlpfc.TRUE 
  #                 0.264                 0.153                 0.438 
  #  br5276.amy.neun.TRUE  br5400.amy.neun.TRUE      br5400.sacc.TRUE 
  #                 0.098                 0.257                 0.328 
  #       br5701.amy.TRUE 
  #                 0.159

# Variable thresholds
sapply(high.mito, function(x){round(attributes(x)[["thresholds"]]["higher"], 4)})
    # br5276.sacc.neun.higher       br5400.nac.higher       br5276.nac.higher 
    #                  5.1915                  0.0816                  1.2914 
    #  br5701.nac.neun.higher br5701.sacc.neun.higher     br5207.dlpfc.higher 
    #                  2.2890                  0.0863                  0.0000 
    #  br5276.amy.neun.higher  br5400.amy.neun.higher      br5400.sacc.higher 
    #                  5.2586                  0.1519                  0.0000 
    #       br5701.amy.higher 
    #                  0.1798


## Test: add a pseudo-count==1 for a 'MT transcript'
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

    
    
    
    
## Bind stats to colData
for(i in 1:length(pilot.data.2)){
  colData(pilot.data.2[[i]]) <- cbind(colData(pilot.data.2[[i]]), stats[[i]])
}
for(i in 1:length(pilot.data.2)){
  colData(pilot.data.2[[i]]) <- cbind(colData(pilot.data.2[[i]]), high.mito[[i]])
}
for(i in 1:length(pilot.data.2)){
  colnames(colData(pilot.data.2[[i]]))[13] <- "high.mito"
}

# $sum vs. $total ??
for(i in 1:length(pilot.data.2)){
  print(table(pilot.data.2[[i]]$sum == pilot.data.2[[i]]$total))
}
    ## all TRUE so can remove this second column:
for(i in 1:length(pilot.data.2)){
  pilot.data.2[[i]]$total <- NULL
}


# Store original for comparison/plotting
pilot.data.2.unfiltered <- pilot.data.2

## Subset - remove those indexed as high.mito
for(i in 1:length(pilot.data.2)){
  pilot.data.2[[i]] <- pilot.data.2[[i]][ ,!high.mito[[i]]]
}
sapply(pilot.data.2, dim)

## Plot metrics

mitoCutoffs <- unlist(lapply(high.mito, function(x){attributes(x)$thresholds["higher"]}))
mean(mitoCutoffs)
    # [1] 1.453033;;     0.3903217 for first batch (n=12)
median(mitoCutoffs)
    # [1] 0.1657229;;    0.138046 for first batch (n=12)
mitoCutoffs <- round(mitoCutoffs, 3)

#dir.create("pdfs/revision")
pdf("pdfs/revision/all-FACS-n10_2021rev_QCmetrics_high-mitoColored_MNT.pdf", height=5)
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
  # Total count vs mito rate
  print(
    plotColData(pilot.data.2.unfiltered[[i]], x="sum", y="subsets_Mito_percent",
              colour_by="high.mito", point_size=2.5, point_alpha=0.5) +
    ggtitle(paste0("Sample: ", names(pilot.data.2.unfiltered)[[i]],
                   ";    nNuclei (pre-QC): ", ncol(pilot.data.2.unfiltered[[i]]),";    ",
                   "Mito % (cutoff = ", mitoCutoffs[i],")"
                   ))
  )
}
dev.off()


## Save!
save(pilot.data.2, pilot.data.2.unfiltered, e.out.2, file="rdas/revision/all-FACS-n10_2021rev_SCEs_processing-QC_MNTMar2021.rda")


        # === === === === === === === === === === ===
        # And end here -> proceed to 'step02' scripts
        # === === === === === === === === === === ===




