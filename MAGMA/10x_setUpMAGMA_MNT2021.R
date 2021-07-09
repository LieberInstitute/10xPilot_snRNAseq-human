### for MAGMA with LIBD 10x pilot analyses
  # MNT Jul2021 update ===================


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


## GWAS for coronary artery disease (2015, downloaded to /GWAS_Results/)
 #    added MNT 02Sep2020 for 'negative' control
sumStats.CAD <- read.table("./GWAS_Results/cad.add.160614.website.txt", header=T)
class(sumStats.CAD)
dim(sumStats.CAD)
head(sumStats.CAD)
unique(sumStats.CAD$CHROM)  # only 1:22 as well
snploc.CAD <- sumStats.CAD[ ,c(1:3)]
colnames(snploc.CAD) <- c("SNP", "CHR", "BP")
write.table(snploc.CAD, file="./GWAS_Results/CoronaryArteryDisease.snploc",
            sep="\t", col.names=T, row.names=F, quote=F)

# Add column for N_estim_mnt: using '60,801 CAD cases and 123,504 controls' across 48 studies
sumStats.CAD$N_estim_mnt <- round(sumStats.CAD$n_studies / 48 * (60801+123504), 0)
write.table(sumStats.CAD, file="./GWAS_Results/CoronaryArteryDisease_Nikpay2015_MNTaddn.txt",
            sep="\t", col.names=T, row.names=F, quote=F)

rm(list=ls(pattern=".CAD"))


## GWAS for ADHD (PGC x iPSYCH
#     added MNT Jun2021
sumStats.ADHD <- read.table("./GWAS_Results/daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta", header=T)
class(sumStats.ADHD)
dim(sumStats.ADHD)
head(sumStats.ADHD)
unique(sumStats.ADHD$CHR)  # 1:22 as well
snploc.ADHD <- sumStats.ADHD[ ,c(2,1,3)]
#colnames(snploc.ADHD) <- c("SNP", "CHR", "BP")
write.table(snploc.ADHD, file="./GWAS_Results/ADHD_PGC2018.snploc",
            sep="\t", col.names=T, row.names=F, quote=F)
rm(list=ls(pattern=".ADHD"))


## GWAS for AD (PGC-ALZ, IGAP, ADSP, UKB meta-meta-analysis: Jansen, et al.2019)
#     added MNT Jul2021
sumStats.AD <- read.table("./GWAS_Results/AD_sumstats_Jansenetal_2019sept.txt", header=T)
class(sumStats.AD)
dim(sumStats.AD)
head(sumStats.AD)
unique(sumStats.AD$CHR)  # 1:22 as well
snploc.AD <- sumStats.AD[ ,c(6,2,3)]
#colnames(snploc.AD) <- c("SNP", "CHR", "BP")
write.table(snploc.AD, file="./GWAS_Results/AD_PGC-IGAP-ADSP-UKB_2019.snploc",
            sep="\t", col.names=T, row.names=F, quote=F)
rm(list=ls(pattern=".AD"))



### Gene marker sets (analagous to layer marker setup) ==============
  # We'll just use the 'enriched' stats--i.e. '1vAll'
  # Update May 2021: Apply a filter for markers that median expression in respective subcluster != 0
  #                  (already computed and added to a different .rda file)
  # UPDATE FOR REVISION: the 'log.FDR' are on the natural log scale, NOT BASE 10
  #                      -> Will switch to log.FDR <= log(1e-6), since this is ~= log10(1e-12)
  #                         (log(1e-12) == -27.63, way too strict)
  
library(SingleCellExperiment)

## DLPFC ===
load("rdas/revision/markers-stats_DLPFC-n3_findMarkers-SN-LEVEL_LAHMay2021.rda", verbose=T)
    # markers.t.1vAll, markers.t.1vAll.db, markers.t.pw, markers.wilcox.block
    rm(markers.t.1vAll.db, markers.t.pw, markers.wilcox.block)

sapply(markers.t.1vAll, function(x){table(x$log.FDR < log(1e-6) & x$non0median==TRUE)})
    #       Oligo Astro Inhib.4 Excit.L4:5 Micro Inhib.6   OPC Excit.L2:3 Excit.ambig Excit.L3:4
    # FALSE 27306 27536   25959      25458 27527   26421 27198      26506       26040      26748
    # TRUE    805   575    2152       2653   584    1690   913       1605        2071       1363

    #       Excit.L5:6 Excit.L6.broad Inhib.5 Excit.L5 Inhib.1 Inhib.2 Inhib.3
    # FALSE      26325          26010   26964    26884   27917   27883   27509
    # TRUE        1786           2101    1147     1227     194     228     602

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
  dlpfc.i <- data.frame(Set=rep(i, sum(markers.dlpfc.t.1vAll[[i]]$log.FDR < log(1e-6) &
                                         markers.dlpfc.t.1vAll[[i]]$non0median==TRUE)),
             Gene=rownames(markers.dlpfc.t.1vAll[[i]])[markers.dlpfc.t.1vAll[[i]]$log.FDR < log(1e-6) &
                                                         markers.dlpfc.t.1vAll[[i]]$non0median==TRUE])
  dlpfc.markerSet <- rbind(dlpfc.markerSet, dlpfc.i)
}

head(dlpfc.markerSet)
table(dlpfc.markerSet$Set)  # looks good

dlpfc.markerSet$Gene <- rowData(sce.dlpfc)$gene_id[match(dlpfc.markerSet$Gene, rownames(sce.dlpfc))]

# Write out
write.table(dlpfc.markerSet, file="./MAGMA/dlpfcMarkerSets_fdr1e-6.txt", sep="\t",
            row.names=F, col.names=T, quote=F)



## sACC ===
load("rdas/revision/markers-stats_sACC-n5_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.sacc.t.pw, markers.sacc.wilcox.block, markers.sacc.t.1vAll, medianNon0.sacc
    rm(markers.sacc.t.pw, markers.sacc.wilcox.block, medianNon0.sacc)

# For EnsemblIDs
load("rdas/revision/regionSpecific_sACC-n5_cleaned-combined_SCE_MNT2021.rda", verbose=T)

## These stats now have both an '_enriched' & '_depleted' result - take the '_enriched'
names(markers.sacc.t.1vAll[[1]])
markers.sacc.enriched <- lapply(markers.sacc.t.1vAll, function(x){x[[2]]})

sapply(markers.sacc.enriched, function(x){table(x$log.FDR < log(1e-6) & x$non0median==TRUE)})
    #       Astro_A Astro_B Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Excit_G
    # FALSE   28657   29260   26458   26706   25799   27170   27236   27415   29329
    # TRUE      926     323    3125    2877    3784    2413    2347    2168     254
    
    #       Inhib_A Inhib_B Inhib_C Inhib_D Inhib_E Inhib_F Inhib_G Inhib_H Inhib_I
    # FALSE   27983   26930   28248   28231   28795   28525   28796   28454   29333
    # TRUE     1600    2653    1335    1352     788    1058     787    1129     250
    
    #       Inhib_J Inhib_K Micro Neu_FAT2.CDH15 Oligo_A Oligo_B   OPC
    # FALSE   29358   29389 28922          29361   28829   29470 28568
    # TRUE      225     194   661            222     754     113  1015

sacc.markerSet <- data.frame()
for(i in names(markers.sacc.enriched)){
  sacc.i <- data.frame(Set=rep(i, sum(markers.sacc.enriched[[i]]$log.FDR < log(1e-6) &
                                        markers.sacc.enriched[[i]]$non0median==TRUE)),
                        Gene=rownames(markers.sacc.enriched[[i]])[markers.sacc.enriched[[i]]$log.FDR < log(1e-6) &
                                                                   markers.sacc.enriched[[i]]$non0median==TRUE])
  sacc.markerSet <- rbind(sacc.markerSet, sacc.i)
}

head(sacc.markerSet)
table(sacc.markerSet$Set)  # looks good

sacc.markerSet$Gene <- rowData(sce.sacc)$gene_id[match(sacc.markerSet$Gene, rownames(sce.sacc))]
    # (btw:)
    length(unique(sacc.markerSet$Gene)) # [1] 7528

# Write out
write.table(sacc.markerSet, file="./MAGMA/saccMarkerSets_fdr1e-6.txt", sep="\t",
            row.names=F, col.names=T, quote=F)




## HPC ===
load("rdas/revision/markers-stats_HPC-n3_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.hpc.t.pw, markers.hpc.t.1vAll, medianNon0.hpc
    rm(markers.hpc.t.pw, medianNon0.hpc)

# For EnsemblIDs
load("rdas/revision/regionSpecific_HPC-n3_cleaned-combined_SCE_MNT2021.rda", verbose=T)

## These stats now have both an '_enriched' & '_depleted' result - take the '_enriched'
names(markers.hpc.t.1vAll[[1]])
markers.hpc.enriched <- lapply(markers.hpc.t.1vAll, function(x){x[[2]]})

sapply(markers.hpc.enriched, function(x){table(x$log.FDR < log(1e-6) & x$non0median==TRUE)})
    #       Astro_A Astro_B Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Excit_G
    # FALSE   27735   28127   26801   25912   28486   28186   28436   27461   28565
    # TRUE     1029     637    1963    2852     278     578     328    1303     199
    
    #       Excit_H Inhib_A Inhib_B Inhib_C Inhib_D Micro Mural Oligo   OPC OPC_COP
    # FALSE   28151   26181   27670   28568   27918 27975 28553 27853 27767   28570
    # TRUE      613    2583    1094     196     846   789   211   911   997     194
    
    #       Tcell
    # FALSE 28491
    # TRUE    273

# A priori: Remove rare / not-homeostatic cell pops
markers.hpc.enriched[["OPC_COP"]] <- NULL

hpc.markerSet <- data.frame()
for(i in names(markers.hpc.enriched)){
  hpc.i <- data.frame(Set=rep(i, sum(markers.hpc.enriched[[i]]$log.FDR < log(1e-6) &
                                       markers.hpc.enriched[[i]]$non0median==TRUE)),
                       Gene=rownames(markers.hpc.enriched[[i]])[markers.hpc.enriched[[i]]$log.FDR < log(1e-6) &
                                                                 markers.hpc.enriched[[i]]$non0median==TRUE])
  hpc.markerSet <- rbind(hpc.markerSet, hpc.i)
}

head(hpc.markerSet)
table(hpc.markerSet$Set)  # looks good

hpc.markerSet$Gene <- rowData(sce.hpc)$gene_id[match(hpc.markerSet$Gene, rownames(sce.hpc))]
    # (btw:)
    length(unique(hpc.markerSet$Gene)) # [1] 6517

# Write out
write.table(hpc.markerSet, file="./MAGMA/hpcMarkerSets_fdr1e-6.txt", sep="\t",
            row.names=F, col.names=T, quote=F)



## NAc ===
load("rdas/revision/markers-stats_NAc-n8_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.nac.t.pw, markers.nac.t.1vAll, medianNon0.nac
    rm(markers.nac.t.pw, medianNon0.nac)

# For EnsemblIDs
load("rdas/revision/regionSpecific_NAc-n8_cleaned-combined_MNT2021.rda", verbose=T)

## These stats now have both an '_enriched' & '_depleted' result - take the '_enriched'
markers.nac.enriched <- lapply(markers.nac.t.1vAll, function(x){x[[2]]})

sapply(markers.nac.enriched, function(x){table(x$log.FDR < log(1e-6) & x$non0median==TRUE)})
    #       Astro_A Astro_B Inhib_A Inhib_B Inhib_C Inhib_D Inhib_E Macro_infilt
    # FALSE   29115   28198   27938   29267   28462   28230   29161        29258
    # TRUE      565    1482    1742     413    1218    1450     519          422
    
    #       Micro Micro_resting MSN.D1_A MSN.D1_B MSN.D1_C MSN.D1_D MSN.D1_E MSN.D1_F
    # FALSE 28792         29615    26549    29203    28866    27408    28648    29317
    # TRUE    888            65     3131      477      814     2272     1032      363
    
    #       MSN.D2_A MSN.D2_B MSN.D2_C MSN.D2_D Oligo_A Oligo_B   OPC OPC_COP
    # FALSE    26742    28484    28741    29393   28377   29010 28461   29434
    # TRUE      2938     1196      939      287    1303     670  1219     246

# A priori: Remove rare / not-homeostatic cell pops
for(i in c("Macro_infilt", "Micro_resting", "OPC_COP")){
  markers.nac.enriched[[i]] <- NULL
}

nac.markerSet <- data.frame()
for(i in names(markers.nac.enriched)){
  nac.i <- data.frame(Set=rep(i, sum(markers.nac.enriched[[i]]$log.FDR <= log(1e-6) &
                                       markers.nac.enriched[[i]]$non0median==TRUE)),
                      Gene=rownames(markers.nac.enriched[[i]])[markers.nac.enriched[[i]]$log.FDR < log(1e-6) &
                                                                markers.nac.enriched[[i]]$non0median==TRUE])
  nac.markerSet <- rbind(nac.markerSet, nac.i)
}

head(nac.markerSet)
table(nac.markerSet$Set)  # looks good

nac.markerSet$Gene <- rowData(sce.nac)$gene_id[match(nac.markerSet$Gene, rownames(sce.nac))]
    # (btw:)
    length(unique(nac.markerSet$Gene)) # [1] 7224

# Write out
write.table(nac.markerSet, file="./MAGMA/nacMarkerSets_fdr1e-6.txt", sep="\t",
            row.names=F, col.names=T, quote=F)





## AMY ===
load("rdas/revision/markers-stats_Amyg-n5_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.amy.t.pw, markers.amy.wilcox.block, markers.amy.t.1vAll, medianNon0.amy
    rm(markers.amy.t.pw, markers.amy.wilcox.block, medianNon0.amy)

# For EnsemblIDs
load("rdas/revision/regionSpecific_Amyg-n5_cleaned-combined_SCE_MNT2021.rda", verbose=T)

## These stats now have both an '_enriched' & '_depleted' result - take the '_enriched'
names(markers.amy.t.1vAll[[1]]) # check
markers.amy.enriched <- lapply(markers.amy.t.1vAll, function(x){x[[2]]})

sapply(markers.amy.enriched, function(x){table(x$log.FDR < log(1e-6) & x$non0median==TRUE)})
    #       Astro_A Astro_B  Endo Excit_A Excit_B Excit_C Inhib_A Inhib_B Inhib_C
    # FALSE   28039   29314 29010   25936   28676   28182   26175   26980   27657
    # TRUE     1332      57   361    3435     695    1189    3196    2391    1714
    
    #       Inhib_D Inhib_E Inhib_F Inhib_G Inhib_H Micro Mural Oligo   OPC Tcell
    # FALSE   27091   28125   27808   29066   28864 28482 29000 28517 27994 29128
    # TRUE     2280    1246    1563     305     507   889   371   854  1377   243


amy.markerSet <- data.frame()
for(i in names(markers.amy.enriched)){
  amy.i <- data.frame(Set=rep(i, sum(markers.amy.enriched[[i]]$log.FDR < log(1e-6) &
                                       markers.amy.enriched[[i]]$non0median==TRUE)),
                      Gene=rownames(markers.amy.enriched[[i]])[markers.amy.enriched[[i]]$log.FDR < log(1e-6) &
                                                                markers.amy.enriched[[i]]$non0median==TRUE])
  amy.markerSet <- rbind(amy.markerSet, amy.i)
}

head(amy.markerSet)
table(amy.markerSet$Set)  # looks good

amy.markerSet$Gene <- rowData(sce.amy)$gene_id[match(amy.markerSet$Gene, rownames(sce.amy))]
    # (btw:)
    length(unique(amy.markerSet$Gene)) # [1] 7439

# Write out
write.table(amy.markerSet, file="./MAGMA/amyMarkerSets_fdr1e-6.txt", sep="\t",
            row.names=F, col.names=T, quote=F)




