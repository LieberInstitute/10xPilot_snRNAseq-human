### MNT 10x snRNA-seq workflow: step 02
###   **PAN-BRAIN analyses**
###     - (n=12) all regions from: Br5161 & Br5212 & Br5287
###     - Amyg, DLPFC, HPC, NAc, and sACC
### Initiated MNT 07Feb2020
### test.edit
#####################################################################

library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
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


load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_processing-QC_MNTJan2020.rda",
     verbose=T)
# pilot.data, pilot.data.unfiltered, e.out

### MNT comment: At this point, each sample (which is a SCE object in the list, 'pilot.data') has been
#              QC'd for cell/nucleus calling ('emptyDrops()' test) and mito rate thresholding


### Merging shared-region samples ============================================
  # Newest iterations for normalization: cbind, THEN take scaled LSFs computed on all nuclei
  # (i.e. no more MBN, because batch is so confounded with sample)

# Add $sample identity
for(i in 1:length(pilot.data)){
  pilot.data[[i]]$sample <- names(pilot.data)[i]
}

sce.all.n12 <- cbind(pilot.data[[1]], pilot.data[[2]], pilot.data[[3]], pilot.data[[4]],
                     pilot.data[[5]], pilot.data[[6]], pilot.data[[7]], pilot.data[[8]],
                     pilot.data[[9]], pilot.data[[10]], pilot.data[[11]], pilot.data[[12]])

# Remove $logcounts
assay(sce.all.n12, "logcounts") <- NULL
# Re-generate log-normalized counts
sce.all.n12 <- logNormCounts(sce.all.n12)

geneVar.all.n12 <- modelGeneVar(sce.all.n12)
chosen.hvgs.all.n12 <- geneVar.all.n12$bio > 0
sum(chosen.hvgs.all.n12)
    # [1] 8748    - kinda surprising...


### Dimensionality reduction ================================================================

# Run PCA, taking top 250 (instead of default 50 PCs)
set.seed(109)
sce.all.n12 <- runPCA(sce.all.n12, subset_row=chosen.hvgs.all.n12, ncomponents=250,
                  BSPARAM=BiocSingular::RandomParam())

# Save into a new data file, which will dedicate for pan-brain-analyses
save(sce.all.n12, chosen.hvgs.all.n12, file="rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda")



    # === === === === === === === ===
    ## 'getClusteredPCs()' evaluated in qsub mode (with 'R-batchJob_panBrain_optimalPCselxn_MNTFeb2020.R')
    #    --> saved into same .rda





### Picking up with optimally-defined PC space ===
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda",
     verbose=TRUE)
    # sce.all.n12, chosen.hvgs.all.n12, pc.choice.n12


# How many PCs is optimal?:
metadata(pc.choice.n12)$chosen
    ## 204

## Assign this chosen ( PCs) to 'PCA_opt'
reducedDim(sce.all.n12, "PCA_opt") <- reducedDim(sce.all.n12, "PCA")[ ,1:(metadata(pc.choice.n12)$chosen)]


## t-SNE
set.seed(109)
sce.all.n12 <- runTSNE(sce.all.n12, dimred="PCA_opt")


## UMAP
set.seed(109)
sce.all.n12 <- runUMAP(sce.all.n12, dimred="PCA_opt")



# Save for now
save(sce.all.n12, chosen.hvgs.all.n12, pc.choice.n12, ref.sampleInfo,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda")



### MNT 20Mar2020 === === ===
  # Add region-specific annotations

load("rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda", verbose=T)
    # sce.all.n12, chosen.hvgs.all.n12, pc.choice.n12, ref.sampleInfo, clusterRefTab.all.n12

## Amyg
load("rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    # sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo
    rm(chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy)

## DLPFC
# load("rdas/regionSpecific_DLPFC-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
#     # sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc, ref.sampleInfo
#     rm(chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc)
    #or
    load("rdas/regionSpecific_DLPFC-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
        # sce.dlpfc.st
    
## HPC
load("rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    # sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo
    rm(chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc)
# Change an annotation
sce.hpc$cellType <- factor(gsub(pattern="Ambig.glial", "Tcell", sce.hpc$cellType))

## NAc
#load("rdas/regionSpecific_NAc-n3_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
#    # sce.nac, chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac, ref.sampleInfo
#    rm(chosen.hvgs.nac, pc.choice.nac, clusterRefTab.nac)

    ## Actually use the annotations for these respective nuclei from the all-NAc analysis
    load("rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda", verbose=T)
        # sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all, ref.sampleInfo
        rm(chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all)
    
    # TEMPorarily set $cellType.final to $cellType and change casing of 'ambig.lowNtrxts'
    sce.nac.all$cellType <- sce.nac.all$cellType.final
    sce.nac.all$cellType <- ifelse(sce.nac.all$cellType=="ambig.lowNtrxts",
                                   "Ambig.lowNtrxts",
                                   as.character(sce.nac.all$cellType))
    # Then collapse D1 and D2 subclusters & otherwise-inhibitory
    sce.nac.all$cellType[grep("MSN.D1", sce.nac.all$cellType)] <- "MSN.D1"
    sce.nac.all$cellType[grep("MSN.D2", sce.nac.all$cellType)] <- "MSN.D2"
    sce.nac.all$cellType[grep("Inhib", sce.nac.all$cellType)] <- "Inhib"
    
    sce.nac.all$cellType <- factor(sce.nac.all$cellType)
    
## sACC
load("rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    # sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc, ref.sampleInfo
    rm(chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc)

    # Collapse as above
    sce.sacc$cellType <- as.character(sce.sacc$cellType)
    sce.sacc$cellType[grep("Inhib", sce.sacc$cellType)] <- "Inhib"
    sce.sacc$cellType[grep("Excit", sce.sacc$cellType)] <- "Excit"
    sce.sacc$cellType <- factor(sce.sacc$cellType)
    

## Assign to new colData column
regions <- list(sce.amy, sce.dlpfc, sce.hpc, sce.nac.all, sce.sacc)

BC.sample.ii <- paste0(colnames(sce.all.n12),".",sce.all.n12$sample)

#sce.all.n12$cellType.RS <- NA
    # or
    sce.all.n12$cellType.RS.sub <- NA

for(r in regions){
  matchTo <- paste0(colnames(r),".",r$sample)
  BC.sample.sub <- BC.sample.ii %in% matchTo
  sce.all.n12$cellType.RS[BC.sample.sub] <- as.character(r$cellType[match(BC.sample.ii[BC.sample.sub], matchTo)])
}

table(sce.all.n12$cellType.RS)
sce.all.n12$cellType.RS <- factor(sce.all.n12$cellType.RS)
table(sce.all.n12$cellType.RS)


## Added chunk 11May2020: add in NON-collapsed region-specific annotations ===
    # First, just load all objects but don't collapse, above

    # Temp - call everything $cellType.split (if not already)
    sce.nac.all$cellType.split <- sce.nac.all$cellType.final
    sce.sacc$cellType.split <- sce.sacc$cellType
    
    regions <- list(sce.amy, sce.dlpfc.st, sce.hpc, sce.nac.all, sce.sacc)
    names(regions) <- c("amy", "dlpfc", "hpc", "nac", "sacc")
    BC.sample.ii <- paste0(colnames(sce.all.n12),".",sce.all.n12$sample)
    
    # Add that annotation
    for(r in names(regions)){
      matchTo <- paste0(colnames(regions[[r]]),".",regions[[r]]$sample)
      BC.sample.sub <- BC.sample.ii %in% matchTo
      # Shorten
      regions[[r]]$cellType.split <- gsub("Excit", "Ex", regions[[r]]$cellType.split)
      regions[[r]]$cellType.split <- gsub("Inhib", "In", regions[[r]]$cellType.split)
      # Add region
      regions[[r]]$cellType.split <- paste0(regions[[r]]$cellType.split, "_", r)
      sce.all.n12$cellType.RS.sub[BC.sample.sub] <- as.character(regions[[r]]$cellType.split[match(BC.sample.ii[BC.sample.sub], matchTo)])
    }
    sce.all.n12$cellType.RS.sub <- factor(sce.all.n12$cellType.RS.sub)
    unique(sce.all.n12$cellType.RS.sub)    
        # 73 == 68 + 5 'Ambig.lowNtrxts' for each region.  good.



# Getting rid of sub-cell types in sACC sample and the pan-brain assignments to compare
table(gsub("MSN","Inhib",ss(as.character(sce.all.n12$cellType.RS),"\\.",1)) ==
        ss(as.character(sce.all.n12$cellType),"\\.",1))
    # FALSE  TRUE
    #   866 33204    - 33204/(866+ 33204) = 97.5% congruence

# Region specific
table(ss(as.character(sce.all.n12$cellType.RS),"\\.",1))
    # Ambig Astro Excit Inhib Micro   MSN Oligo   OPC Tcell
    #   445  3864  2927  2019  2956   642 18664  2527    26

table( ss(as.character(sce.all.n12$cellType),"\\.",1))
    #Ambig Astro Excit Inhib Micro Oligo   OPC
    #   32  3828  2848  3110  3077 18614  2561


## Pretty good - let's save
save(sce.all.n12, chosen.hvgs.all.n12, pc.choice.n12, ref.sampleInfo, clusterRefTab.all.n12,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda")


### Re-print some visualizations - MNT 10Apr2020 ====================
# With original coordinates and annotation
pdf("pdfs/panBrain-n12_reducedDims-with-collapsedClusters_Apr2020.pdf")
plotReducedDim(sce.all.n12, dimred="PCA", ncomponents=5, colour_by="cellType", point_alpha=0.5)
plotTSNE(sce.all.n12, colour_by="region", point_alpha=0.5, point_size=2.5) + ggtitle("t-SNE on opt PCs (d=204)")
plotTSNE(sce.all.n12, colour_by="processDate", point_alpha=0.5, point_size=2.5) + ggtitle("t-SNE on opt PCs (d=204)")
plotTSNE(sce.all.n12, colour_by="sample", point_alpha=0.5, point_size=2.5) + ggtitle("t-SNE on opt PCs (d=204)")
plotTSNE(sce.all.n12, colour_by="cellType", point_alpha=0.5, point_size=2.5) + ggtitle("t-SNE on opt PCs (d=204)")
# Region-specific annotation
plotTSNE(sce.all.n12, colour_by="cellType.RS", point_alpha=0.5, point_size=2.5) + ggtitle("t-SNE on opt PCs (d=204): region-specific annot")
plotTSNE(sce.all.n12, colour_by="sum", point_alpha=0.5, point_size=2.5) + ggtitle("t-SNE on opt PCs (d=204)")
plotUMAP(sce.all.n12, colour_by="cellType", point_alpha=0.5, point_size=2.5) + ggtitle("UMAP on opt PCs (d=204)")
dev.off()

# How many PCs?
head(attr(reducedDim(sce.all.n12, "PCA"), "percentVar"), n=75)
    # [1] 19.93610257  8.56551975  5.48499030  2.69104457  2.29910336  0.86823624  0.71376460
    # [8]  0.48377800  0.32700635  0.31751426  0.29901816  0.25911566  0.22686401  0.21400937
    # [15]  0.20303113  0.18768708  0.18136268  0.17008113  0.15380352  0.14366776  0.13666573
    # [22]  0.12527399  0.11727016  0.11141297  0.10514416  0.10317094  0.09743560  0.09296292
    # [29]  0.08796466  0.08415014  0.08363044  0.08132186  0.08054339  0.07687419  0.07353210
    # [36]  0.07313253  0.07035361  0.06841083  0.06810964  0.06737097  0.06464408  0.06338695
    # [43]  0.06261076  0.06170788  0.06062411  0.05978423  0.05923917  0.05845309  0.05738356
    # [50]  0.05623525  0.05588523  0.05546358  0.05420676  0.05329303  0.05254012  0.05199590
    # [57]  0.05148879  0.05073573  0.05023337  0.04904533  0.04896795  0.04883193  0.04804923
    # [64]  0.04767100  0.04691248  0.04680456  0.04633970  0.04582236  0.04560561  0.04480624
    # [71]  0.04469499  0.04421348  0.04352501  0.04339655  0.04307707

# 0.05% var or greater
reducedDim(sce.all.n12, "PCA_59") <- reducedDim(sce.all.n12, "PCA")[ ,c(1:59)]
# 0.1% var or greater
reducedDim(sce.all.n12, "PCA_26") <- reducedDim(sce.all.n12, "PCA")[ ,c(1:26)]

# First remove this reducedDim bc this has caused trouble previously
reducedDim(sce.all.n12, "TSNE") <- NULL

## 59 PCs tsNE === (this one looks better than 26 PCs actually)
set.seed(109)
sce.all.tsne.59pcs <- runTSNE(sce.all.n12, dimred="PCA_59")

save(sce.all.tsne.59pcs, file="rdas/ztemp_panBrain-n12_SCE-with-tSNEon59PCs_MNT.rda")
rm(sce.all.tsne.59pcs)


# MNT 16Apr: Deciding to remove the clusters won't focus on for plotting:
    #'Ambig.hiVCAN' & 'Excit.4' & those .RS-annot'd as 'Ambig.lowNtrxts'
sce.all.tsne.59pcs <- sce.all.tsne.59pcs[ ,sce.all.tsne.59pcs$cellType.RS != "Ambig.lowNtrxts"] # 445
sce.all.tsne.59pcs$cellType.RS <- droplevels(sce.all.tsne.59pcs$cellType.RS)

sce.all.tsne.59pcs <- sce.all.tsne.59pcs[ ,sce.all.tsne.59pcs$cellType != "Ambig.hiVCAN"] # 32 nuclei
sce.all.tsne.59pcs <- sce.all.tsne.59pcs[ ,sce.all.tsne.59pcs$cellType != "Excit.4"]  # 33 nuclei
sce.all.tsne.59pcs$cellType <- droplevels(sce.all.tsne.59pcs$cellType)

# Add broad cell type taken from pan-brain annotation
sce.all.tsne.59pcs$cellType.broad <- ss(as.character(sce.all.tsne.59pcs$cellType), "\\.", 1)

#pdf("pdfs/exploration/ztemp_panBrain-n12_TSNEon59PCs_MNT.pdf")
pdf("pdfs/pubFigures/panBrain-n12_tSNEon59PCs_3x3PCA_smallClustersDropped_MNTApr2020.pdf", width=9)
plotTSNE(sce.all.tsne.59pcs, colour_by="cellType", point_alpha=0.5, point_size=4.0,
         text_by="cellType.broad", text_size=8, theme_size=18) +
  ggtitle("t-SNE on top 59 PCs (pan-brain annot.)") + theme(plot.title = element_text(size=19))

plotTSNE(sce.all.tsne.59pcs, colour_by="cellType.broad", point_alpha=0.5, point_size=4.0,
         text_by="cellType.broad", text_size=7, theme_size=22) +
  ggtitle("t-SNE on top 59 PCs (broad pan-brain annot.)") + theme(plot.title = element_text(size=18))

plotTSNE(sce.all.tsne.59pcs, colour_by="cellType.RS", point_alpha=0.5, point_size=4.0,
         text_by="cellType.RS", text_size=5.5, theme_size=17) +
  ggtitle("t-SNE on top 59 PCs (region-specific annot.)") + theme(plot.title = element_text(size=18))

# Top 3 PCs
plotReducedDim(sce.all.tsne.59pcs, dimred="PCA", ncomponents=3, colour_by="cellType.broad",
               point_alpha=0.5, theme_size=15, add_legend=FALSE) + 
  ggtitle("Top three PCs (broad pan-brain annot.)") + theme(plot.title = element_text(size=18))
dev.off()


    ### MNT update 31Aug2020 === ===
      # Facet some different iterations of this 'best' tSNE: maybe by region
    
    pdf("pdfs/pubFigures/panBrain-n12_tSNEon59PCs_faceted_Aug2020.pdf", width=9)
    plotTSNE(sce.all.tsne.59pcs, colour_by="region", point_alpha=0.5, point_size=4.0, theme_size=22) +
      facet_wrap(~ sce.all.tsne.59pcs$region)
      ggtitle("t-SNE on top 59 PCs (broad pan-brain annot.)") + theme(plot.title = element_text(size=18))
    dev.off()
    
    ## More manually to have shadow of those for each region ======
    custom.cols <- c("DLPFC"=tableau20[1],
                  "sACC"=tableau20[3],
                  "HPC"=tableau20[5],
                  "AMY"=tableau20[7],
                  "NAc"=tableau20[9])
    
    sce.temp <- sce.all.tsne.59pcs
    
    ## DLPFC
    # Reorder to plot the region nuclei last
    sce.temp <- cbind(sce.temp[ ,sce.temp$region!="dlpfc"],
                      sce.temp[ ,sce.temp$region=="dlpfc"])
    # Grey out other regions
    sce.temp$reg.temp <- ifelse(sce.temp$region=="dlpfc", "DLPFC", NA)
    
    p.dlpfc <- plotTSNE(sce.temp, colour_by="reg.temp", point_alpha=0.9, point_size=3.5, theme_size=15,
             add_legend=FALSE) +
      scale_fill_manual(values=custom.cols["DLPFC"]) + ggtitle("DLPFC") +
      theme(plot.title = element_text(size=30),
            axis.title = element_text(size=0),
            axis.text = element_text(size=20))
    
    ## HPC
    # Reorder to plot the region nuclei last
    sce.temp <- cbind(sce.temp[ ,sce.temp$region!="hpc"],
                      sce.temp[ ,sce.temp$region=="hpc"])
    # Grey out other regions
    sce.temp$reg.temp <- ifelse(sce.temp$region=="hpc", "HPC", NA)
    
    p.hpc <- plotTSNE(sce.temp, colour_by="reg.temp", point_alpha=0.9, point_size=3.5, theme_size=15,
                        add_legend=FALSE) +
      scale_fill_manual(values=custom.cols["HPC"]) + ggtitle("HIPPO") +
      theme(plot.title = element_text(size=30),
            axis.title = element_text(size=0),
            axis.text = element_text(size=20))
    
    ## sACC
    # Reorder to plot the region nuclei last
    sce.temp <- cbind(sce.temp[ ,sce.temp$region!="sacc"],
                      sce.temp[ ,sce.temp$region=="sacc"])
    # Grey out other regions
    sce.temp$reg.temp <- ifelse(sce.temp$region=="sacc", "sACC", NA)
    
    p.sacc <- plotTSNE(sce.temp, colour_by="reg.temp", point_alpha=0.9, point_size=3.5, theme_size=15,
                      add_legend=FALSE) +
      scale_fill_manual(values=custom.cols["sACC"]) + ggtitle("sACC") +
      theme(plot.title = element_text(size=30),
            axis.title = element_text(size=0),
            axis.text = element_text(size=20))
    
    ## AMY
    # Reorder to plot the region nuclei last
    sce.temp <- cbind(sce.temp[ ,sce.temp$region!="amy"],
                      sce.temp[ ,sce.temp$region=="amy"])
    # Grey out other regions
    sce.temp$reg.temp <- ifelse(sce.temp$region=="amy", "AMY", NA)
    
    p.amy <- plotTSNE(sce.temp, colour_by="reg.temp", point_alpha=0.9, point_size=3.5, theme_size=15,
                       add_legend=FALSE) +
      scale_fill_manual(values=custom.cols["AMY"]) + ggtitle("AMY") +
      theme(plot.title = element_text(size=30),
            axis.title = element_text(size=0),
            axis.text = element_text(size=20))
    
    ## NAc
    # Reorder to plot the region nuclei last
    sce.temp <- cbind(sce.temp[ ,sce.temp$region!="nac"],
                      sce.temp[ ,sce.temp$region=="nac"])
    # Grey out other regions
    sce.temp$reg.temp <- ifelse(sce.temp$region=="nac", "NAc", NA)
    
    p.nac <- plotTSNE(sce.temp, colour_by="reg.temp", point_alpha=0.9, point_size=3.5, theme_size=15,
                      add_legend=FALSE) +
      scale_fill_manual(values=custom.cols["NAc"]) + ggtitle("NAc") +
      theme(plot.title = element_text(size=30),
            axis.title = element_text(size=0),
            axis.text = element_text(size=20))
    
            ## end region-colored t-SNEs ========
    
    
    ## All nuclei (pan-brain annotation, from above) ===
    p.full <- plotTSNE(sce.all.tsne.59pcs, colour_by="cellType", point_alpha=0.5, point_size=4.0,
             text_by="cellType.broad", text_size=8, theme_size=24) +
      ggtitle("t-SNE on top 59 PCs (pan-brain annot.)") + theme(plot.title = element_text(size=28))
    
    lay <- rbind(c(1,1,2),
                 c(1,1,3),
                 c(6,5,4))
    
    pdf("pdfs/pubFigures/panBrain-n12_tSNEon59PCs_faceted_v2_Aug2020.pdf", width=13.5, height=12.5)
    grid.arrange(grobs=list(p.full,
                         p.nac,
                         p.amy,
                         p.hpc,
                         p.dlpfc,
                         p.sacc),
                 layout_matrix=lay)
    dev.off()


        
## Print broad marker heatmap of pan-brain-defined clusters === === ===
load("rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda", verbose=T)
    #sce.all.n12, chosen.hvgs.all.n12, pc.choice.n12, ref.sampleInfo, clusterRefTab.all.n12

# As decided for marker detection, remove the clusters that won't focus on:
#     'Ambig.hiVCAN' & 'Excit.4' & those .RS-annot'd as 'Ambig.lowNtrxts'
sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType.RS != "Ambig.lowNtrxts"] # 445
sce.all.n12$cellType.RS <- droplevels(sce.all.n12$cellType.RS)

sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType != "Ambig.hiVCAN"] # 32 nuclei
sce.all.n12 <- sce.all.n12[ ,sce.all.n12$cellType != "Excit.4"]  # 33 nuclei
sce.all.n12$cellType <- droplevels(sce.all.n12$cellType)

cell.idx <- splitit(sce.all.n12$cellType)
dat <- as.matrix(assay(sce.all.n12, "logcounts"))


genes <- c('SNAP25','SLC17A6','SLC17A7','GAD1','GAD2','AQP4','GFAP','C3','CD74','MBP','PDGFRA','VCAN','CLDN5','FLT1','SKAP1','TRAC')
current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMeans(dat[genes, ii])))

current_dat <- current_dat[ ,c(3:14, 1,2, 15:17)]

pdf("pdfs/pubFigures/heatmap-geneExprs_panBrain-annot_mean-broadMarkers_MNT.pdf")
pheatmap(current_dat, cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(0.02, 4, length.out = 101),
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "OrRd"))(100),
         fontsize_row = 17.5, fontsize_col = 17.5,
         main="Broad cell type marker expression (pan-brain annot.)", fontsize=11)
dev.off()



