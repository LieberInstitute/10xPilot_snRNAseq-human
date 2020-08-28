### for MAGMA with LIBD 10x pilot analyses
  # MNT 14Aug2020 ========================


### Set up annotation === === === === === === === === === === === === 
library(rtracklayer)

## get annotation
map = read.delim("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/151675/151675_raw_feature_bc_matrix__features.tsv.gz",
                 as.is=TRUE, header=FALSE,col.names=c("EnsemblID", "Symbol", "Type"))
## get GTF, this seems like what they used
gtf = import("/dcl01/ajaffe/data/lab/singleCell/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf")
## of length 2565061
gtf = gtf[gtf$type == "gene"]
## of length 33538

table(seqnames(gtf) %in% c(1:21,"X","Y","MT"))
## 32672 TRUE - maybe just take these for the .loc file

# chr 1:21, X, Y, MT
annoTab <- as.data.frame(ranges(gtf[seqnames(gtf) %in% c(1:21,"X","Y","MT")]))

annoTab <- annoTab[ ,c("names", "start", "end")]
annoTab$seqlevel <- as.character(seqnames(gtf)[seqnames(gtf) %in% c(1:21,"X","Y","MT")])
annoTab$strand <- as.character(strand(gtf)[seqnames(gtf) %in% c(1:21,"X","Y","MT")])
annoTab$symbol <- mcols(gtf)$gene_name[as.logical(seqnames(gtf) %in% c(1:21,"X","Y","MT"))]
annoTab <- annoTab[ ,c("names", "seqlevel", "start", "end", "strand", "symbol")]

write.table(annoTab, file="./MAGMA/GRCh38_Ensembl-93_GENES_chr-x-y-mt.gene.loc", sep="\t",
            row.names=F, col.names=F, quote=F)

# Full version
annoTab.full <- as.data.frame(ranges(gtf))
annoTab.full <- annoTab.full[ ,c("names", "start", "end")]
annoTab.full$seqlevel <- as.character(seqnames(gtf))
annoTab.full$strand <- as.character(strand(gtf))
annoTab.full$symbol <- mcols(gtf)$gene_name
annoTab.full <- annoTab.full[ ,c("names", "seqlevel", "start", "end", "strand", "symbol")]

write.table(annoTab.full, file="./MAGMA/GRCh38_Ensembl-93_GENES_all-33538.gene.loc", sep="\t",
            row.names=F, col.names=F, quote=F)


### MNT update 21Aug2020: *** USE GRCh37 / hg19 GENE COORDINATES ***  ===
  # * ALSO only use those for expressed genes - otherwise penalizes stats for expressed genes

## Define all expressed genes (union)
load("../rdas/markers-stats_DLPFC_n2_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.dlpfc.t.1vAll, markers.dlpfc.t.design

load("../rdas/markers-stats_HPC-n3_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.hpc.t.1vAll, markers.hpc.t.design, markers.hpc.wilcox.block

load("../rdas/markers-stats_NAc-n5_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.nac.t.design, markers.nac.t.1vAll

load("../rdas/markers-stats_Amyg-n2_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.amy.t.1vAll, markers.amy.t.design, markers.amy.wilcox.block

load("../rdas/markers-stats_sACC-n2_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.sacc.t.design, markers.sacc.wilcox.block, markers.sacc.binom.block, markers.sacc.t.1vAll

    # Remove all pairwise stats
    rm(list=setdiff(ls(), ls(pattern=".1vAll")))


# Re-order rows in each list entry for each set of stats
expressedGenes.list <- list(amy=rownames(markers.amy.t.1vAll[["Astro"]]),
                            dlpfc=rownames(markers.dlpfc.t.1vAll[["Astro"]]),
                            hpc=rownames(markers.hpc.t.1vAll[["Astro"]]),
                            nac=rownames(markers.nac.t.1vAll[["Astro"]]),
                            sacc=rownames(markers.sacc.t.1vAll[["Astro"]])
)

expressedGenes <- unique(unlist(expressedGenes.list))
length(expressedGenes)  # 30252

# Convert to EnsemblID to subset
library(SingleCellExperiment)
load("../rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    # sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo

expressedGenes <- rowData(sce.amy)$ID[match(expressedGenes, rownames(sce.amy))]


## get GTF
gtf.full = import("/dcl01/ajaffe/data/lab/singleCell/refdata-cellranger-hg19-3.0.0/genes/genes.gtf")
    ## of length 2490452 (as opposed to 2565061, above)
gtf = gtf.full[gtf.full$type == "gene"]
    ## of length 32738 (as opposed to 33538, above)

    # table(seqnames(gtf) %in% c(1:21,"X","Y","MT"))
    # ## 32672 TRUE - maybe just take these for the .loc file
    # 
    # # chr 1:21, X, Y, MT
    # annoTab <- as.data.frame(ranges(gtf[seqnames(gtf) %in% c(1:21,"X","Y","MT")]))
    # 
    # annoTab <- annoTab[ ,c("names", "start", "end")]
    # annoTab$seqlevel <- as.character(seqnames(gtf)[seqnames(gtf) %in% c(1:21,"X","Y","MT")])
    # annoTab$strand <- as.character(strand(gtf)[seqnames(gtf) %in% c(1:21,"X","Y","MT")])
    # annoTab$symbol <- mcols(gtf)$gene_name[as.logical(seqnames(gtf) %in% c(1:21,"X","Y","MT"))]
    # annoTab <- annoTab[ ,c("names", "seqlevel", "start", "end", "strand", "symbol")]
    # 
    # write.table(annoTab, file="./MAGMA/GRCh38_Ensembl-93_GENES_chr-x-y-mt.gene.loc", sep="\t",
    #             row.names=F, col.names=F, quote=F)


table(expressedGenes %in% gtf$gene_id)
    # FALSE  TRUE
    #  1790 28462

notInHg19 <- expressedGenes[!(expressedGenes %in% gtf$gene_id)]
notInHg19.symbol <- rowData(sce.amy)$Symbol[match(notInHg19, rowData(sce.amy)$ID)]
    # Mostly 'AL' / 'AC' genes
    # * only thing that stood out on glance is seeing GABRQ would be dropped...
    #     - in hg19 this is 'ENSG00000147402'
    #     - in GRCh38 it's actually 'ENSG00000268089'

    ## Hmmm what does this look like?
    ens = read.delim("/dcl02/lieber/ajaffe/Nick_Clifton/magma/GRCh37_ensembl_GENES.gene.loc",
                     header=FALSE,as.is=TRUE)
    ens = ens[ens$V1 %in% rownames(stats),]
    write.table(ens, file = "MAGMA/GRCh37_ensembl_GENES_SpatialExprs.gene.loc",
                col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
    
    table(expressedGenes %in% ens$V6)
        # FALSE  TRUE
        # 10419 19833

    ###   *** --> go to 'lift.R' for lifting coordinates with R package `liftOver`, using
      #           this already-generated list of all 'expressedGenes'



## MNT 23Aug2020: Additional GWAS for alcohol/tobacco use === === ===
sumStats.ageInit <- read.table("./GWAS_Results/AgeofInitiation.txt", header=T)
class(sumStats.ageInit)
dim(sumStats.ageInit)
head(sumStats.ageInit)
unique(sumStats.ageInit$CHROM)  # only 1:22, interesting... (is this common?)
snploc.ageInit <- sumStats.ageInit[ ,c(3,1,2)]
colnames(snploc.ageInit) <- c("SNP", "CHR", "BP")
write.table(snploc.ageInit, file="./GWAS_Results/AgeofInitiation.snploc",
            sep="\t", col.names=T, row.names=F, quote=F)
rm(list=ls(pattern=".ageInit"))

sumStats.cigDay <- read.table("./GWAS_Results/CigarettesPerDay.txt", header=T)
class(sumStats.cigDay)
dim(sumStats.cigDay)
head(sumStats.cigDay)
unique(sumStats.cigDay$CHROM) # same
snploc.cigDay <- sumStats.cigDay[ ,c(3,1,2)]
colnames(snploc.cigDay) <- c("SNP", "CHR", "BP")
write.table(snploc.cigDay, file="./GWAS_Results/CigarettesPerDay.snploc",
            sep="\t", col.names=T, row.names=F, quote=F)
rm(list=ls(pattern=".cigDay"))

sumStats.driWk <- read.table("./GWAS_Results/DrinksPerWeek.txt", header=T)
class(sumStats.driWk)
dim(sumStats.driWk)
head(sumStats.driWk)
unique(sumStats.driWk$CHROM) # same
snploc.driWk <- sumStats.driWk[ ,c(3,1,2)]
colnames(snploc.driWk) <- c("SNP", "CHR", "BP")
write.table(snploc.driWk, file="./GWAS_Results/DrinksPerWeek.snploc",
            sep="\t", col.names=T, row.names=F, quote=F)
rm(list=ls(pattern=".driWk"))

sumStats.smokCess <- read.table("./GWAS_Results/SmokingCessation.txt", header=T)
class(sumStats.smokCess)
dim(sumStats.smokCess)
head(sumStats.smokCess)
unique(sumStats.smokCess$CHROM) # same
snploc.smokCess <- sumStats.smokCess[ ,c(3,1,2)]
colnames(snploc.smokCess) <- c("SNP", "CHR", "BP")
write.table(snploc.smokCess, file="./GWAS_Results/SmokingCessation.snploc",
            sep="\t", col.names=T, row.names=F, quote=F)
rm(list=ls(pattern=".smokCess"))

sumStats.smokInit <- read.table("./GWAS_Results/SmokingInitiation.txt", header=T)
class(sumStats.smokInit)
dim(sumStats.smokInit)
head(sumStats.smokInit)
unique(sumStats.smokInit$CHROM) # same
snploc.smokInit <- sumStats.smokInit[ ,c(3,1,2)]
colnames(snploc.smokInit) <- c("SNP", "CHR", "BP")
write.table(snploc.smokInit, file="./GWAS_Results/SmokingInitiation.snploc",
            sep="\t", col.names=T, row.names=F, quote=F)
rm(list=ls(pattern=".smokInit"))





### Gene marker sets (analagous to layer marker setup) ==============
  # We'll just use the 'enriched' stats--i.e. '1vAll'

## DLPFC ===
load("rdas/markers-stats_DLPFC_n2_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.dlpfc.t.1vAll, markers.dlpfc.t.design
    rm(markers.dlpfc.t.design)

sapply(markers.dlpfc.t.1vAll, function(x){table(x$log.FDR < log10(1e-12))})
    #      Oligo Astro Inhib.4 Excit.L4:5 Micro Inhib.6   OPC Excit.L2:3 Excit.ambig
    # FALSE 26195 25277   24745      23641 25565   24960 25912      25367       24707
    # TRUE   1916  2834    3366       4470  2546    3151  2199       2744        3404
    #       Excit.L3:4 Excit.L5:6 Excit.L6.broad Inhib.5 Excit.L5 Inhib.1 Inhib.2
    # FALSE      26329      25338          25069   25793    26300   27738   27663
    # TRUE        1782       2773           3042    2318     1811     373     448
    #       Inhib.3
    # FALSE   27281
    # TRUE      830     - this looks good

    ## Btw:
        layerMarkers <- read.table("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/MAGMA/laminar_sets.txt",
                                   sep="\t", header=T)
        table(layerMarkers$Set)
            #    All_Neuropil          Layer1 Layer1_Neuropil          Layer2 Layer2_Neuropil
            #             581            1876              27            1512              20
            #          Layer3 Layer3_Neuropil          Layer4 Layer4_Neuropil          Layer5
            #             270             344             610             296             794
            # Layer5_Neuropil          Layer6 Layer6_Neuropil              WM     WM_Neuropil
            #             253             432             246            5010             215



dlpfc.markerSet <- data.frame()
for(i in names(markers.dlpfc.t.1vAll)){
  dlpfc.i <- data.frame(Set=rep(i, sum(markers.dlpfc.t.1vAll[[i]]$log.FDR < log10(1e-12))),
             Gene=rownames(markers.dlpfc.t.1vAll[[i]])[markers.dlpfc.t.1vAll[[i]]$log.FDR < log10(1e-12)])
  dlpfc.markerSet <- rbind(dlpfc.markerSet, dlpfc.i)
}

head(dlpfc.markerSet)
table(dlpfc.markerSet$Set)  # looks good

# Need to convert to Ensembl ID (these symbols are uniquified, so have to refer to SCE info)
load("rdas/regionSpecific_DLPFC-n2_SCE_cellTypesSplit-fromST_Apr2020.rda")
    # sce.dlpfc.st, 

dlpfc.markerSet$Gene <- rowData(sce.dlpfc.st)$ID[match(dlpfc.markerSet$Gene, rownames(sce.dlpfc.st))]

# Write out
write.table(dlpfc.markerSet, file="./MAGMA/dlpfcMarkerSets_fdr1e-12.txt", sep="\t",
            row.names=F, col.names=T, quote=F)



## sACC ===
load("rdas/markers-stats_sACC-n2_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.sacc.t.design, markers.sacc.wilcox.block, markers.sacc.binom.block, markers.sacc.t.1vAll
    rm(markers.sacc.t.design, markers.sacc.wilcox.block, markers.sacc.binom.block)

sapply(markers.sacc.t.1vAll, function(x){table(x$log.FDR < log10(1e-12))})
    #       Astro Excit.1 Excit.2 Excit.3 Excit.4 Inhib.1 Inhib.2 Micro Oligo   OPC
    # FALSE 26058   21609   22359   25007   25887   22961   23800 26034 26789 26437
    # TRUE   2716    7165    6415    3767    2887    5813    4974  2740  1985  2337
        ## A lot of genes... but don't want to adjust too much the cutoff on a per-dataset-basis...

    
sacc.markerSet <- data.frame()
for(i in names(markers.sacc.t.1vAll)){
  sacc.i <- data.frame(Set=rep(i, sum(markers.sacc.t.1vAll[[i]]$log.FDR < log10(1e-12))),
                        Gene=rownames(markers.sacc.t.1vAll[[i]])[markers.sacc.t.1vAll[[i]]$log.FDR < log10(1e-12)])
  sacc.markerSet <- rbind(sacc.markerSet, sacc.i)
}

head(sacc.markerSet)
table(sacc.markerSet$Set)  # looks good


# Need to convert to Ensembl ID (these symbols are uniquified, so have to refer to SCE info)
load("rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    # sce.sacc, chosen.hvgs.sacc, pc.choice.sacc, clusterRefTab.sacc, ref.sampleInfo

sacc.markerSet$Gene <- rowData(sce.sacc)$ID[match(sacc.markerSet$Gene, rownames(sce.sacc))]
    # (btw:)
    length(unique(sacc.markerSet$Gene)) # [1] 16923 - though 40,799 entries

# Write out
write.table(sacc.markerSet, file="./MAGMA/saccMarkerSets_fdr1e-12.txt", sep="\t",
            row.names=F, col.names=T, quote=F)




## HPC ===
load("rdas/markers-stats_HPC-n3_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
  # markers.hpc.t.1vAll, markers.hpc.t.design, markers.hpc.wilcox.block

sapply(markers.hpc.t.1vAll, function(x){table(x$log.FDR < log10(1e-12))})
    # #      Astro Excit.1 Excit.2 Excit.3 Excit.4 Excit.5 Inhib.1 Inhib.2 Inhib.3
    # FALSE 24318   25940   25446   23167   26648   27400   27262   26239   25232
    # TRUE   4439    2817    3311    5590    2109    1357    1495    2518    3525
    #       Inhib.4 Inhib.5 Micro Oligo   OPC Tcell
    # FALSE   26884   26777 24765 25986 25602 27739
    # TRUE     1873    1980  3992  2771  3155  1018


hpc.markerSet <- data.frame()
for(i in names(markers.hpc.t.1vAll)){
  hpc.i <- data.frame(Set=rep(i, sum(markers.hpc.t.1vAll[[i]]$log.FDR < log10(1e-12))),
                       Gene=rownames(markers.hpc.t.1vAll[[i]])[markers.hpc.t.1vAll[[i]]$log.FDR < log10(1e-12)])
  hpc.markerSet <- rbind(hpc.markerSet, hpc.i)
}

head(hpc.markerSet)
table(hpc.markerSet$Set)  # looks good


# Need to convert to Ensembl ID (these symbols are uniquified, so have to refer to SCE info)
load("rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    # sce.hpc, chosen.hvgs.hpc, pc.choice.hpc, clusterRefTab.hpc, ref.sampleInfo

hpc.markerSet$Gene <- rowData(sce.hpc)$ID[match(hpc.markerSet$Gene, rownames(sce.hpc))]
    # (btw:)
    length(unique(hpc.markerSet$Gene)) # [1] 17552

# Write out
write.table(hpc.markerSet, file="./MAGMA/hpcMarkerSets_fdr1e-12.txt", sep="\t",
            row.names=F, col.names=T, quote=F)



## NAc ===
load("rdas/markers-stats_NAc-n5_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.nac.t.design, markers.nac.t.1vAll
    rm(markers.nac.t.design)

sapply(markers.nac.t.1vAll, function(x){table(x$log.FDR < log10(1e-12))})
    #       Astro Inhib.1 Inhib.2 Inhib.3 Inhib.4 Micro MSN.D1.1 MSN.D1.2 MSN.D1.3
    # FALSE 25362   28529   28310   25829   26081 26330    28352    27649    27571
    # TRUE   3874     707     926    3407    3155  2906      884     1587     1665
    #       MSN.D1.4 MSN.D2.1 MSN.D2.2 Oligo   OPC
    # FALSE    25006    27827    27042 27303 26849
    # TRUE      4230     1409     2194  1933  2387


nac.markerSet <- data.frame()
for(i in names(markers.nac.t.1vAll)){
  nac.i <- data.frame(Set=rep(i, sum(markers.nac.t.1vAll[[i]]$log.FDR < log10(1e-12))),
                      Gene=rownames(markers.nac.t.1vAll[[i]])[markers.nac.t.1vAll[[i]]$log.FDR < log10(1e-12)])
  nac.markerSet <- rbind(nac.markerSet, nac.i)
}

head(nac.markerSet)
table(nac.markerSet$Set)  # looks good


# Need to convert to Ensembl ID (these symbols are uniquified, so have to refer to SCE info)
load("rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda", verbose=T)
    # sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all, ref.sampleInfo

nac.markerSet$Gene <- rowData(sce.nac.all)$ID[match(nac.markerSet$Gene, rownames(sce.nac.all))]
    # (btw:)
    length(unique(nac.markerSet$Gene)) # [1] 15515

# Write out
write.table(nac.markerSet, file="./MAGMA/nacMarkerSets_fdr1e-12.txt", sep="\t",
            row.names=F, col.names=T, quote=F)





## AMY ===
load("rdas/markers-stats_Amyg-n2_findMarkers-SN-LEVEL_MNTMay2020.rda", verbose=T)
    # markers.amy.t.1vAll, markers.amy.t.design, markers.amy.wilcox.block
    rm(markers.amy.t.design, markers.amy.wilcox.block)

sapply(markers.amy.t.1vAll, function(x){table(x$log.FDR < log10(1e-12))})
    #       Astro Excit.1 Excit.2 Excit.3 Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5
    # FALSE 24426   22088   26685   26092   23941   25135   26957   27455   25612
    # TRUE   4038    6376    1779    2372    4523    3329    1507    1009    2852
    #       Micro Oligo   OPC
    # FALSE 25215 26127 25604
    # TRUE   3249  2337  2860   - also quite high numbers... just keep this in mind


amy.markerSet <- data.frame()
for(i in names(markers.amy.t.1vAll)){
  amy.i <- data.frame(Set=rep(i, sum(markers.amy.t.1vAll[[i]]$log.FDR < log10(1e-12))),
                      Gene=rownames(markers.amy.t.1vAll[[i]])[markers.amy.t.1vAll[[i]]$log.FDR < log10(1e-12)])
  amy.markerSet <- rbind(amy.markerSet, amy.i)
}

head(amy.markerSet)
table(amy.markerSet$Set)  # looks good


# Need to convert to Ensembl ID (these symbols are uniquified, so have to refer to SCE info)
load("rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    # sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo

amy.markerSet$Gene <- rowData(sce.amy)$ID[match(amy.markerSet$Gene, rownames(sce.amy))]
    # (btw:)
    length(unique(amy.markerSet$Gene)) # [1] 16533

# Write out
write.table(amy.markerSet, file="./MAGMA/amyMarkerSets_fdr1e-12.txt", sep="\t",
            row.names=F, col.names=T, quote=F)




