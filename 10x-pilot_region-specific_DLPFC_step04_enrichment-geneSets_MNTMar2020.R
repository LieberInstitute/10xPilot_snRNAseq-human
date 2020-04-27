### MNT 10x snRNA-seq workflow: step 04 - enrichment testing
###   **Region-specific analyses**
###     - (2x) DLPFC samples from: Br5161 & Br5212
### MNT Feb-Mar2020
#####################################################################

library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(batchelor)
library(DropletUtils)
library(jaffelab)
library(limma)
library(lattice)
library(parallel)
library(RColorBrewer)
library(pheatmap)
library(fields)

### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

# ===

## From findMarkers() ====
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_DLPFC_n2_findMarkers_MNTMar2020.rda",
     verbose=T)
    # markers.dlpfc.t.design.log, markers.dlpfc.t.design.countsN
    rm(markers.dlpfc.t.design.log, markers.dlpfc.t.design.countsN)
  
sapply(markers.dlpfc.t.design.log, function(x){table(x$FDR<0.05)})
    ##       Astro Excit Inhib Micro Oligo   OPC
    # FALSE 28074 27912 28112 27219 27997 28118
    # TRUE     54   216    16   909   131    10


sapply(markers.dlpfc.t.design.countsN, function(x){table(x$FDR<0.05)})
    ##      Astro Excit Inhib Micro Oligo   OPC
    #FALSE 27470 27467 27980 26384 27428 27888
    #TRUE    658   661   148  1744   700   240


        # ## Actually what is this? ====
        # load("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/rda/tstats_Human_DLPFC_snRNAseq_Nguyen.Rdata",
        #      verbose=T)
        #     # tstats_Human_DLPFC_snRNAseq_Nguyen_topLayer
        # 
        # dim(tstats_Human_DLPFC_snRNAseq_Nguyen_topLayer)  # 692 x 31
        # ## not sure if applicable here ====



### The generation of these stats are adapted from:
  # `/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/layer_specificity.R`
  # (and performed in step03 script - "direct limma approach")

load("rdas/markers-stats_DLPFC_n2_manualContrasts_MNTMar2020.rda", verbose=T)
    # eb_contrasts.dlpfc.broad, eb_list.dlpfc.broad, sce.dlpfc.PB
        
        ## Tabulating significant hits for pairwise contrasts (won't look at really) ====
        pvals_contrasts <- eb_contrasts.dlpfc.broad$p.value
        
        data.frame(
          'FDRsig' = colSums(apply(pvals_contrasts, 2, p.adjust, 'fdr') < 0.05),
          'Pval10-6sig' = colSums(pvals_contrasts < 1e-6),
          'Pval10-8sig' = colSums(pvals_contrasts < 1e-8)
        )
            #             FDRsig Pval10.6sig Pval10.8sig
            #Astro-Excit   5636         247           7
            #Astro-Inhib   4573         183           4
            #Astro-Micro   6878         499          40
            #Astro-Oligo   4487         233           6
            #Astro-OPC     3211         131           6
            #Excit-Inhib    907          28           1
            #Excit-Micro   9029         586          43
            #Excit-Oligo   5984         306          11
            #Excit-OPC     4311         162           4
            #Inhib-Micro   8129         534          45
            #Inhib-Oligo   5061         261          11
            #Inhib-OPC     2971         113           3
            #Micro-Oligo   6063         464          37
            #Micro-OPC     6558         501          49
            #Oligo-OPC     3617         205           9 
        # ====


## Extract the p-values and compute fdrs
pvals0_broad <- sapply(eb_list.dlpfc.broad, function(x) {
  x$p.value[, 2, drop = FALSE]
})

fdrs0_broad = apply(pvals0_broad, 2, p.adjust, "fdr")

## Extract the tstats
t0_broad <- sapply(eb_list.dlpfc.broad, function(x) {
  x$t[, 2, drop = FALSE]
})


data.frame(
  'FDRsig' = colSums(apply(pvals0_broad, 2, p.adjust, 'fdr') < 0.05),
  'Pval10-6sig' = colSums(pvals0_broad < 1e-6),
  'Pval10-8sig' = colSums(pvals0_broad < 1e-8)
)

# For only (+) t-stats
data.frame(
  'FDRsig' = colSums(apply(pvals0_broad, 2, p.adjust, 'fdr') < 0.05 &
                       t0_broad > 0),
  'Pval10-6sig' = colSums(pvals0_broad < 1e-6 &
                            t0_broad > 0),
  'Pval10-8sig' = colSums(pvals0_broad < 1e-8 &
                            t0_broad > 0)
)

## Without t > 0 subset:
    #FDRsig Pval10.6sig Pval10.8sig
    #Astro    537          49          16
    #Excit    690          90          26
    #Inhib    123          28           5
    #Micro   3589         359         116
    #Oligo    770          71          13
    #OPC      113          24           4

## With t > 0
    #      FDRsig Pval10.6sig Pval10.8sig
    #Astro    461          47          16
    #Excit    635          89          26
    #Inhib    121          28           5
    #Micro   1909         264         101
    #Oligo    480          42          10
    #OPC       97          24           4




### Gene list enrichment test =============================================
  # Followwing approach taken in ST project in:
  # /dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/check_clinical_gene_sets.R
  # MNT re-created gene lists in [10x_pilot_FINAL]/ref_ST-check_clinical_gene_sets.R

## Load gene lists curated by AnJa
load("rdas/geneLists-fromSTpaper_forEnrichTests_MNT.rda", verbose=T)
    # geneLists.fromST
    ## 37 lists

# add and change rownames to Ensembl ID
rownames(t0_broad) <- rownames(sce.dlpfc.PB)
rownames(t0_broad) <- rowData(sce.dlpfc.PB)$ID[match(rownames(t0_broad), rownames(sce.dlpfc.PB))]

## filter for those present in stats
geneList_present = lapply(geneLists.fromST, function(x) {
  x = x[!is.na(x)]
  x[x %in% rownames(t0_broad)]
})

unname(data.frame(lengths(geneLists.fromST), lengths(geneList_present)))
    ## not too bad

## do enrichment
enrich_stat_list = eb_list.dlpfc.broad
for (i in seq(along = enrich_stat_list)) {
  #layer = t0_broad[, i] > 0 & fdrs0_broad[, i] < 0.1
  # or
  layer = t0_broad[, i] > 0 & fdrs0_broad[, i] < 0.05
  tabList = mclapply(geneList_present, function(g) {
    tt = table(Set = factor(names(layer) %in% g, c(FALSE, TRUE)),
               Layer = factor(layer, c(FALSE, TRUE)))
  }, mc.cores = 8)
  enrichList = lapply(tabList,fisher.test)
  
  o = data.frame(
    OR = sapply(enrichList, "[[", "estimate"),
    Pval = sapply(enrichList, "[[", "p.value"),
    NumSig = sapply(tabList, function(x) x[2,2])
  )
  rownames(o) = gsub(".odds ratio", "", rownames(o))
  enrich_stat_list[[i]] = o
}

enrichTab = do.call("cbind", enrich_stat_list)

#  name
enrichTab$Type = ss(rownames(enrichTab), "_", 1)
    enrichTab$Group = ss(rownames(enrichTab), "_", 2)
enrichTab$Type[enrichTab$Group == "Birnbaum"] = "Birnbaum"
enrichTab$Type[enrichTab$Type == "Gene"] = "ASD"
#enrichTab$Group = ss(rownames(enrichTab), "_", 2)
enrichTab$Set = ss(rownames(enrichTab), "_", 3)
enrichTab$ID = rownames(enrichTab)
enrichTab$SetSize = sapply(geneList_present, length)

### save a copy as a supp table
#enrichTabOut.fdr.10 = enrichTab[ ,c(22, 19:21, 23, 1:18)]
# or
enrichTabOut.fdr.05 = enrichTab[ ,c(22, 19:21, 23, 1:18)]

#dir.create("./tables/")
#write.csv(enrichTabOut.fdr.05, file = "tables/enrichTab_clinicalGeneLists_DLPFC-broadClusters-fdr05.csv", row.names=FALSE)
#write.csv(enrichTabOut.fdr.10, file = "tables/enrichTab_clinicalGeneLists_DLPFC-broadClusters-fdr10.csv", row.names=FALSE)


## look at enrichment
pMat = enrichTab[ , grep("Pval", colnames(enrichTab))]
orMat = enrichTab[ , grep("OR", colnames(enrichTab))]
colnames(pMat) = ss(colnames(pMat), "\\.")
colnames(orMat) = ss(colnames(orMat), "\\.")
pMat < 0.05 / nrow(pMat)
pMat < 0.001
round(-log10(pMat),1)

######################
## pull out results ##
######################
# 
# ## summary stats from genes
# enrichTab["Gene_SFARI_all",]
# enrichTab["Gene_Satterstrom_ASC102.2018",]
# enrichTab["Gene_Satterstrom_ASD53",]
# enrichTab["Gene_Satterstrom_DDID49",]
# 
# ## Satterstrom deep dive
# sat_102_l2= which(t0_contrasts[,"Layer2"] > 0 & fdrs0_contrasts[,"Layer2"] < 0.1 & 
#                     rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_ASC102.2018)
# rowData(sce_layer)$gene_name[sat_102_l2]
# sat_102_l5= which(t0_contrasts[,"Layer5"] > 0 & fdrs0_contrasts[,"Layer5"] < 0.1 & 
#                     rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_ASC102.2018)
# rowData(sce_layer)$gene_name[sat_102_l5]
# 
# sat_49_l2= which(t0_contrasts[,"Layer2"] > 0 & fdrs0_contrasts[,"Layer2"] < 0.1 & 
#                    rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_DDID49)
# cat(rowData(sce_layer)$gene_name[sat_49_l2], sep=", ")
# 
# sat_53_l5= which(t0_contrasts[,"Layer5"] > 0 & fdrs0_contrasts[,"Layer5"] < 0.1 & 
#                    rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_ASD53)
# cat(rowData(sce_layer)$gene_name[sat_53_l5], sep=", ")

## case control - asd
enrichTab["DE_PE_ASD.Up",]
enrichTab["DE_PE_ASD.Down",]

## case control - sczd
enrichTab[c("DE_PE_SCZ.Up","DE_BS2_SCZ.Up"),]
enrichTab[c("DE_PE_SCZ.Down","DE_BS2_SCZ.Down"),]

    enrichTabOut.fdr.05[c("DE_PE_SCZ.Up","DE_BS2_SCZ.Up", "TWAS_PE_SCZ.Up","TWAS_BS2_SCZ.Up"), ]
    enrichTabOut.fdr.05[c("DE_PE_SCZ.Down","DE_BS2_SCZ.Down", "TWAS_PE_SCZ.Down","TWAS_BS2_SCZ.Down"), ]

    # Additionally
    lengths(geneLists.fromST)[grep("SCZ", names(geneLists.fromST))]
    #          Gene_Birnbaum_SCZ.SNV      Gene_Birnbaum_SCZ.PGC.GWAS
    #                            212                             106
    #Gene_Birnbaum_SCZ.Meta.analysis           Gene_Birnbaum_SCZ.CNV
    #                             36                             113
    #                   DE_PE_SCZ.Up                  DE_PE_SCZ.Down
    #                           2452                            2373
    #                  DE_BS2_SCZ.Up                 DE_BS2_SCZ.Down
    #                            103                             142
    #                TWAS_BS2_SCZ.Up               TWAS_BS2_SCZ.Down
    #                            483                             505
    #                 TWAS_PE_SCZ.Up                TWAS_PE_SCZ.Down
    #                            561                             544
    #               TWAS_PE_SCZBD.Up              TWAS_PE_SCZBD.Down
    #                             76                              67
    


    
    
    
### Look at some new stats with 'cellType.split'-level clusters =================================
  # MNT 30Mar2020

load("rdas/markers-stats_DLPFC_n2_manualContrasts_neuronalSubs_MNTApr2020.rda", verbose=T)
#    # eb_list.dlpfc.neuronalSubs, sce.dlpfc.st.PB
load("rdas/markers-stats_DLPFC_n2_manualContrasts_neuronalSubs_noBroadTerm_MNTApr2020.rda", verbose=T)
    # eb_list.dlpfc.neuronalSubs.simple, sce.dlpfc.st.PB
    
# Previously
names(eb_list.dlpfc.broad)
    #[1] "Astro" "Excit" "Inhib" "Micro" "Oligo" "OPC"

names(eb_list.dlpfc.neuronalSubs.simple)
    # [1] "Excit.ambig"    "Excit.L2:3"     "Excit.L3:4"     "Excit.L4:5"
    # [5] "Excit.L5"       "Excit.L5:6"     "Excit.L6.broad" "Inhib.1"
    # [9] "Inhib.2"        "Inhib.3"        "Inhib.4"        "Inhib.5"
    # [13] "Inhib.6"

# Combine them
eb_list.dlpfc <- list(eb_list.dlpfc.broad[["Astro"]],
                      eb_list.dlpfc.broad[["Micro"]],
                      eb_list.dlpfc.broad[["Oligo"]],
                      eb_list.dlpfc.broad[["OPC"]])
names(eb_list.dlpfc) <- c("Astro","Micro","Oligo","OPC")
eb_list.dlpfc <- c(eb_list.dlpfc, eb_list.dlpfc.neuronalSubs.simple)
    # or, for stats including 'broad' [excit/inhib] term
eb_list.dlpfc <- c(eb_list.dlpfc, eb_list.dlpfc.neuronalSubs)




## Extract the p-values and compute fdrs
pvals0_full <- sapply(eb_list.dlpfc, function(x) {
  x$p.value[, 2, drop = FALSE]
})

fdrs0_full = apply(pvals0_full, 2, p.adjust, "fdr")

## Extract the tstats
t0_full <- sapply(eb_list.dlpfc, function(x) {
  x$t[, 2, drop = FALSE]
})

# How does this look?
data.frame(
  'FDRsig.05' = colSums(apply(pvals0_full, 2, p.adjust, 'fdr') < 0.05 &
                       t0_full > 0),
  'FDRsig.01' = colSums(apply(pvals0_full, 2, p.adjust, 'fdr') < 0.01 & 
                          t0_full > 0),
  'Pval10-6sig' = colSums(pvals0_full < 1e-6 &
                            t0_full > 0),
  'Pval10-8sig' = colSums(pvals0_full < 1e-8 &
                            t0_full > 0)
)
    #                FDRsig.05 FDRsig.01 Pval10.6sig Pval10.8sig
    # Astro                461       226          47          16
    # Micro               1909      1078         264         101
    # Oligo                480       210          42          10
    # OPC                   97        53          24           4
    # Excit.ambig          240       176          24          14
    # Excit.L2:3           140        60          33          13
    # Excit.L3:4           117        28          18           5
    # Excit.L4:5            96        44          29           6
    # Excit.L5             210        97          48          26
    # Excit.L5:6           117        49          29          13
    # Excit.L6.broad        71        17           8           5
    # Inhib.1              217       136          72          36
    # Inhib.2              602       289           9           0
    # Inhib.3              113        36          11           7
    # Inhib.4               91        41          19           6
    # Inhib.5              101        43          30          16
    # Inhib.6              182        94          58          30




## Load gene lists curated by AnJa
load("rdas/geneLists-fromSTpaper_forEnrichTests_MNT.rda", verbose=T)
    # geneLists.fromST
    ## 37 lists

# add and change rownames to Ensembl ID
rownames(t0_full) <- rownames(sce.dlpfc.st.PB)
rownames(t0_full) <- rowData(sce.dlpfc.st.PB)$ID[match(rownames(t0_full), rownames(sce.dlpfc.st.PB))]

## filter for those present in stats (this won't change b/tw broad & subtype-level stats)
geneList_present = lapply(geneLists.fromST, function(x) {
  x = x[!is.na(x)]
  x[x %in% rownames(t0_full)]
})

unname(data.frame(lengths(geneLists.fromST), lengths(geneList_present)))
    ## not too bad

## do enrichment
enrich_stat_list.full = eb_list.dlpfc
for (i in seq(along = enrich_stat_list.full)) {
  #layer = t0_full[, i] > 0 & fdrs0_full[, i] < 0.1
  # or
  layer = t0_full[, i] > 0 & fdrs0_full[, i] < 0.05
  tabList = mclapply(geneList_present, function(g) {
    tt = table(Set = factor(names(layer) %in% g, c(FALSE, TRUE)),
               Layer = factor(layer, c(FALSE, TRUE)))
  }, mc.cores = 8)
  enrichList = lapply(tabList,fisher.test)
  
  o = data.frame(
    OR = sapply(enrichList, "[[", "estimate"),
    Pval = sapply(enrichList, "[[", "p.value"),
    NumSig = sapply(tabList, function(x) x[2,2])
  )
  rownames(o) = gsub(".odds ratio", "", rownames(o))
  enrich_stat_list.full[[i]] = o
}

enrichTab.full = do.call("cbind", enrich_stat_list.full)

#  name
enrichTab.full$Type = ss(rownames(enrichTab.full), "_", 1)
enrichTab.full$Group = ss(rownames(enrichTab.full), "_", 2)
enrichTab.full$Type[enrichTab.full$Group == "Birnbaum"] = "Birnbaum"
enrichTab.full$Type[enrichTab.full$Type == "Gene"] = "ASD"
enrichTab.full$Set = ss(rownames(enrichTab.full), "_", 3)
enrichTab.full$ID = rownames(enrichTab.full)
enrichTab.full$SetSize = sapply(geneList_present, length)

### save a copy as a supp table
enrichTab.fullOut.fdr.05 = enrichTab.full[ ,c(55, 52:54, 56, 1:51)]
write.csv(enrichTab.fullOut.fdr.05, file = "tables/enrichTab_clinicalGeneLists_DLPFC-cellTypesSplit-fdr05.csv", row.names=FALSE)
# or
write.csv(enrichTab.fullOut.fdr.05, file = "tables/enrichTab_clinicalGeneLists_DLPFC-cellTypesSplit_noBroadTerm-fdr05.csv", row.names=FALSE)

## look at enrichment
pMat = enrichTab.fullOut.fdr.05[ , grep("Pval", colnames(enrichTab.fullOut.fdr.05))]
orMat = enrichTab.fullOut.fdr.05[ , grep("OR", colnames(enrichTab.fullOut.fdr.05))]
colnames(pMat) = ss(colnames(pMat), ".Pval")
colnames(orMat) = ss(colnames(orMat), ".OR")
pMat < 0.05 / nrow(pMat)
pMat < 0.001
round(-log10(pMat),1)


## SCZD gene sets
enrichTab.fullOut.fdr.05[c("DE_PE_SCZ.Up","DE_BS2_SCZ.Up", "TWAS_PE_SCZ.Up","TWAS_BS2_SCZ.Up"),
                         c(2:5,grep("Pval", colnames(enrichTab.fullOut.fdr.05)))]
        #                 Type Group    Set SetSize   Astro.Pval  Micro.Pval   Oligo.Pval
        # DE_PE_SCZ.Up      DE    PE SCZ.Up    2247 3.408401e-52 0.007603556 4.522922e-11
        # DE_BS2_SCZ.Up     DE   BS2 SCZ.Up      93 1.000000e+00 0.094115327 4.134514e-01
        # TWAS_PE_SCZ.Up  TWAS    PE SCZ.Up     456 8.511075e-01 0.132974645 4.110474e-02
        # TWAS_BS2_SCZ.Up TWAS   BS2 SCZ.Up     368 4.030077e-01 0.834424913 1.812590e-03
        #                   OPC.Pval Excit.ambig.Pval Excit.L2:3.Pval Excit.L3:4.Pval
        # DE_PE_SCZ.Up    0.02364615     6.356758e-06      0.87555871       0.1700458
        # DE_BS2_SCZ.Up   1.00000000     1.000000e+00      0.07859219       1.0000000
        # TWAS_PE_SCZ.Up  0.67241126     1.941352e-01      0.17900791       0.7152518
        # TWAS_BS2_SCZ.Up 0.36324032     8.113622e-02      0.26959825       0.4114099
        #                 Excit.L4:5.Pval Excit.L5.Pval Excit.L5:6.Pval
        # DE_PE_SCZ.Up       0.0005191715   0.009940076      0.02517485
        # DE_BS2_SCZ.Up      1.0000000000   1.000000000      1.00000000
        # TWAS_PE_SCZ.Up     1.0000000000   0.270739668      1.00000000
        # TWAS_BS2_SCZ.Up    0.6401488571   0.120999676      0.41140985
        #                 Excit.L6.broad.Pval Inhib.1.Pval Inhib.2.Pval Inhib.3.Pval
        # DE_PE_SCZ.Up             0.04402799    0.3156281   0.76112359    0.3835457
        # DE_BS2_SCZ.Up            1.00000000    0.5141599   1.00000000    1.0000000
        # TWAS_PE_SCZ.Up           1.00000000    0.5910278   0.25557934    1.0000000
        # TWAS_BS2_SCZ.Up          0.60810846    0.5377963   0.06893443    0.6615109
        #                 Inhib.4.Pval Inhib.5.Pval Inhib.6.Pval
        # DE_PE_SCZ.Up      0.01003077   0.02475195  0.001454071
        # DE_BS2_SCZ.Up     1.00000000   1.00000000  0.453960813
        # TWAS_PE_SCZ.Up    0.40911692   1.00000000  1.000000000
        # TWAS_BS2_SCZ.Up   0.63564520   0.64556513  1.000000000

enrichTab.fullOut.fdr.05[c("DE_PE_SCZ.Down","DE_BS2_SCZ.Down", "TWAS_PE_SCZ.Down","TWAS_BS2_SCZ.Down"),
                         c(2:5,grep("Pval", colnames(enrichTab.fullOut.fdr.05)))]
        #                  Type Group      Set SetSize  Astro.Pval   Micro.Pval
        # DE_PE_SCZ.Down      DE    PE SCZ.Down    2077 0.000214712 1.464783e-27
        # DE_BS2_SCZ.Down     DE   BS2 SCZ.Down     132 0.728452517 1.787082e-15
        # TWAS_PE_SCZ.Down  TWAS    PE SCZ.Down     463 0.263615712 1.000000e+00
        # TWAS_BS2_SCZ.Down TWAS   BS2 SCZ.Down     402 0.229405139 1.930231e-01
        #                     Oligo.Pval  OPC.Pval Excit.ambig.Pval Excit.L2:3.Pval
        # DE_PE_SCZ.Down    1.176215e-34 0.6963981      0.001742817     0.003043082
        # DE_BS2_SCZ.Down   2.604828e-02 1.0000000      0.632725573     1.000000000
        # TWAS_PE_SCZ.Down  1.000000e+00 0.6753821      0.195486991     0.732647496
        # TWAS_BS2_SCZ.Down 5.620994e-01 0.6516143      0.053818994     0.726697323
        #                   Excit.L3:4.Pval Excit.L4:5.Pval Excit.L5.Pval Excit.L5:6.Pval
        # DE_PE_SCZ.Down         0.04875122       0.3248871    0.00770607      0.04875122
        # DE_BS2_SCZ.Down        0.42410473       1.0000000    1.00000000      0.42410473
        # TWAS_PE_SCZ.Down       1.00000000       1.0000000    0.27166251      0.27068716
        # TWAS_BS2_SCZ.Down      0.68507946       0.6501990    0.08004020      0.41964416
        #                   Excit.L6.broad.Pval Inhib.1.Pval Inhib.2.Pval Inhib.3.Pval
        # DE_PE_SCZ.Down              0.4908060  0.005779045   0.13444415    0.5860449
        # DE_BS2_SCZ.Down             1.0000000  0.629934109   1.00000000    1.0000000
        # TWAS_PE_SCZ.Down            0.6339505  0.275086819   0.07137112    1.0000000
        # TWAS_BS2_SCZ.Down           0.6290907  0.381643105   0.07856522    0.4155913
        #                   Inhib.4.Pval Inhib.5.Pval Inhib.6.Pval
        # DE_PE_SCZ.Down      0.06705622    0.4467516    0.1529106
        # DE_BS2_SCZ.Down     1.00000000    1.0000000    1.0000000
        # TWAS_PE_SCZ.Down    1.00000000    0.4185697    0.7733553
        # TWAS_BS2_SCZ.Down   0.64363553    0.4089148    1.0000000






################   
## make plots ##
################


## ASD - skip for now
# ## make long
# enrichLong = reshape2::melt(enrichTab[,c(seq(1,19,by=3),22:26)],id.vars = 8:12)
# colnames(enrichLong)[6:7] = c("Layer", "OR")
# enrichLong_P = reshape2::melt(enrichTab[,c(seq(2,20,by=3),22:26)],id.vars = 8:12)
# identical(enrichLong$ID, enrichLong_P$ID)
# enrichLong$P = enrichLong_P$value
# enrichLong$Layer = ss(as.character(enrichLong$Layer), "\\.")
# enrichLong$ID = factor(enrichLong$ID, levels=rev(rownames(enrichTab)))
# enrichLong$Set = factor(enrichLong$Set, levels=unique(rev(enrichTab$Set)))
# enrichLong$FDR = p.adjust(enrichLong$P, "fdr")
# 
# ## what p-value controls FDR?
# enrichLongSort = enrichLong[order(enrichLong$P),]
# max(enrichLongSort$P[enrichLongSort$FDR < 0.05] )
# # 0.01009034
# 
# ## overall ##
# enrichLong$P_thresh = enrichLong$P
# enrichLong$P_thresh[enrichLong$P_thresh < 2.2e-16] = 2.2e-16
# 
# ### ASD focus
# enrichLong_ASD = enrichLong[enrichLong$ID %in% 
#                               c("Gene_SFARI_all", "Gene_Satterstrom_ASC102.2018",
#                                 "Gene_Satterstrom_ASD53", "Gene_Satterstrom_DDID49",
#                                 "DE_PE_ASD.Down", "DE_PE_ASD.Up",
#                                 "TWAS_PE_ASD.Up", "TWAS_PE_ASD.Down"),]
# enrichLong_ASD$ID2 =  as.character(droplevels(enrichLong_ASD$Set))
# enrichLong_ASD$ID2[enrichLong_ASD$ID2 == "all"] = "SFARI"
# enrichLong_ASD$ID2[enrichLong_ASD$ID2 == "ASC102.2018"] = "ASC102"
# enrichLong_ASD$ID2[enrichLong_ASD$ID == "DE_PE_ASD.Up"] = "DE.Up"
# enrichLong_ASD$ID2[enrichLong_ASD$ID == "DE_PE_ASD.Down"] = "DE.Down"
# enrichLong_ASD$ID2[enrichLong_ASD$ID == "TWAS_PE_ASD.Up"] = "TWAS.Up"
# enrichLong_ASD$ID2[enrichLong_ASD$ID == "TWAS_PE_ASD.Down"] = "TWAS.Down"
# enrichLong_ASD$ID2 = factor(enrichLong_ASD$ID2, unique(enrichLong_ASD$ID2))
# 
# enrichLong_ASD$LayerFac = factor(as.character(enrichLong_ASD$Layer), 
#                                  c("WM", paste0("Layer", 6:1)))
# enrichLong_ASD = enrichLong_ASD[order(enrichLong_ASD$ID2, enrichLong_ASD$LayerFac),]



### custom heatmap (SCZD)    - start here - *JUST USE MARGINAL p<0.05 here for plotting OR for now
midpoint = function(x) x[-length(x)] + diff(x)/2

customLayerEnrichment = function(enrichTab , groups, xlabs, 
                                 Pthresh = 12, ORcut = -log10(0.05), enrichOnly = FALSE,
                                 #layerHeights = c(0,40,55,75,85,110,120,135),
                                 layerHeights = seq(0,170, by=10),
                                 mypal = c("white", colorRampPalette(brewer.pal(9,"YlOrRd"))(50)), ...) {
  
  wide_p = -log10( enrichTab[groups,grep("Pval", colnames(enrichTab))])
  wide_p[wide_p > Pthresh] = Pthresh
  wide_p = t(round(wide_p[,
#                          c("WM.Pval", "Layer6.Pval", "Layer5.Pval", "Layer4.Pval", "Layer3.Pval","Layer2.Pval", "Layer1.Pval")],2))
                          ],2))
  
  wide_or = enrichTab[groups,grep("OR", colnames(enrichTab))]
  wide_or= round(t(wide_or[,
#                           c("WM.OR", "Layer6.OR", "Layer5.OR", "Layer4.OR", "Layer3.OR", "Layer2.OR", "Layer1.OR")]),1)
                           ]),1)

  if(enrichOnly) wide_p[wide_or < 1] = 0
  wide_or[wide_p < ORcut] = ""
      # or, if want to print -log10(p's)
      #wide_p_2plot <- wide_p
      #wide_p_2plot[wide_p < ORcut] = ""
  
  image.plot(x = seq(0,ncol(wide_p),by=1), y = layerHeights, z = as.matrix(t(wide_p)),
             col = mypal,xaxt="n", yaxt="n",xlab = "", ylab="", ...)
#  axis(2, c("WM", paste0("L", 6:1)), at = midpoint(layerHeights),las=1)
        axis(2, names(eb_list.dlpfc), at = midpoint(layerHeights),las=1)  # MNT add
      
  axis(1, rep("", ncol(wide_p)), at = seq(0.5,ncol(wide_p)-0.5))
  text(x = seq(0.5,ncol(wide_p)-0.5), y=-1*max(nchar(xlabs))/2, xlabs,
       xpd=TRUE, srt=45,cex=1.5,adj= 1)
  abline(h=layerHeights,v=0:ncol(wide_p))
  text(x = rep(seq(0.5,ncol(wide_p)-0.5),each = nrow(wide_p)), 
       y = rep(midpoint(layerHeights), ncol(wide_p)),
       as.character(wide_or),cex=1.5,font=2)
          #as.character(wide_p_2plot),cex=1.5,font=2)
}

## ASD - skip for now
# pdf("pdf/asd_geneSet_heatmap.pdf",w=6)
# par(mar=c(8,4.5,2.5,1), cex.axis=2,cex.lab=2)
# groups = unique(as.character(enrichLong_ASD$ID))[1:6]
# xlabs  = as.character(enrichLong_ASD$ID2[match(groups, enrichLong_ASD$ID)])
# customLayerEnrichment(enrichTab, groups,xlabs, enrichOnly=TRUE)
# abline(v=4,lwd=3)
# text(x = 3, y = 142, c("ASD"), xpd=TRUE,cex=2.5,font=2)
# 
# dev.off()


pdf("pdfs/exploration/enrichPlots_SCZD-geneSet_DLPFC-cellTypeSplit_heatmap_Apr2020.pdf",w=8)
par(mar=c(8,8,2.5,1), cex.axis=1.2, cex.lab=1.5)
groups =c("DE_PE_SCZ.Up", "DE_PE_SCZ.Down", 
          "DE_BS2_SCZ.Up", "DE_BS2_SCZ.Down", 
          "TWAS_BS2_SCZ.Up", "TWAS_BS2_SCZ.Down", "TWAS_PE_SCZ.Up",
          "TWAS_PE_SCZ.Down")
xlabs = ss(gsub("_SCZ", "", groups), "_", 2)
customLayerEnrichment(enrichTab.full, groups, xlabs, enrichOnly=TRUE)
abline(v=4,lwd=3)
text(x = c(2,6), y = 175, c("SCZD-DE", "SCZD-TWAS"), xpd=TRUE,cex=2,font=2)
dev.off()


pdf("pdfs/exploration/enrichPlots_birnbaum-geneSet_DLPFC-cellTypeSplit_heatmap_Apr2020.pdf",w=8)
par(mar=c(12,8,2.5,1), cex.axis=1, cex.lab=1.5)
groups =grep(enrichTab.full$ID, pattern = "Birnbaum", value=TRUE)
xlabs = ss(groups, "_", 3)
customLayerEnrichment(enrichTab.full, groups,xlabs, enrichOnly=TRUE,
                      breaks = seq(0,12,len = 52))
dev.off()






### 24Apr2020: using single-nucleus-level stats for enrichment ==========================================
  #   - as it's believed these are more 'real' markers

load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_DLPFC_n2_findMarkers-SN-LEVEL_MNTApr2020.rda",
     verbose=T)
    # markers.dlpfc.t.design

# What do these objects/DataFrames look like:
head(markers.dlpfc.t.design[["Astro"]])
head(markers.dlpfc.t.design[["Astro"]][ ,"stats.OPC"])

sapply(markers.dlpfc.t.design, function(x){table(x$FDR<0.05)})
    #       Astro Excit.ambig Excit.L2:3 Excit.L3:4 Excit.L4:5 Excit.L5 Excit.L5:6
    # FALSE 33308       33463      33536      33443      33513    33315      33490
    # TRUE    230          75          2         95         25      223         48
    #       Excit.L6.broad Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5 Inhib.6 Micro Oligo   OPC
    # FALSE          33501   33302   33347   33405   33513   33509   33512 33006 33390 33384
    # TRUE              37     236     191     133      25      29      26   532   148   154

## Load SCE for rowData
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_DLPFC-n2_SCE_cellTypesSplit-fromST_Apr2020.rda", verbose=T)
    # sce.dlpfc.st, clusterRefTab.dlpfc, chosen.hvgs.dlpfc, ref.sampleInfo
    rm(clusterRefTab.dlpfc, chosen.hvgs.dlpfc, ref.sampleInfo)

    # As for test, first drop "Ambig.lowNtrxts" (168 nuclei)
    sce.dlpfc.st <- sce.dlpfc.st[ ,sce.dlpfc.st$cellType.split != "Ambig.lowNtrxts"]
    sce.dlpfc.st$cellType.split <- droplevels(sce.dlpfc.st$cellType.split)
    
    # Remove 0 genes across all nuclei
    sce.dlpfc.st <- sce.dlpfc.st[!rowSums(assay(sce.dlpfc.st, "counts"))==0, ]  # keeps same 28111 genes



## Load gene lists curated by AnJa
load("rdas/geneLists-fromSTpaper_forEnrichTests_MNT.rda", verbose=T)
    # geneLists.fromST
    ## 37 lists

## filter for those present in stats (this won't change b/tw broad & subtype-level stats)
geneList_present = lapply(geneLists.fromST, function(x) {
  x = x[!is.na(x)]
  x[x %in% rowData(sce.dlpfc.st)$ID]
})

unname(data.frame(lengths(geneLists.fromST), lengths(geneList_present)))
## not too bad


## do enrichment ===
enrich_stat_list.full = markers.dlpfc.t.design
for (i in seq(along = enrich_stat_list.full)) {
  #cellType = t0_full[, i] > 0 & fdrs0_full[, i] < 0.1
  # or
  cellType = markers.dlpfc.t.design[[i]][, "FDR"] < 0.05 # a logical
      # don't need to subset for positive t's bc test alrdy does
  # Add names, then change to ensemblID
  names(cellType) <- rownames(markers.dlpfc.t.design[[i]])
  names(cellType) <- rowData(sce.dlpfc.st)$ID[match(names(cellType), rownames(sce.dlpfc.st))]
  tabList = mclapply(geneList_present, function(g) {
    tt = table(Set = factor(names(cellType) %in% g, c(FALSE, TRUE)),
               CellType = factor(cellType, c(FALSE, TRUE)))
  }, mc.cores = 8)
  enrichList = lapply(tabList,fisher.test)
  
  o = data.frame(
    OR = sapply(enrichList, "[[", "estimate"),
    Pval = sapply(enrichList, "[[", "p.value"),
    NumSig = sapply(tabList, function(x) x[2,2])
  )
  rownames(o) = gsub(".odds ratio", "", rownames(o))
  enrich_stat_list.full[[i]] = o
}

enrichTab.full = do.call("cbind", enrich_stat_list.full)



#  name
enrichTab.full$Type = ss(rownames(enrichTab.full), "_", 1)
enrichTab.full$Group = ss(rownames(enrichTab.full), "_", 2)
enrichTab.full$Type[enrichTab.full$Group == "Birnbaum"] = "Birnbaum"
enrichTab.full$Type[enrichTab.full$Type == "Gene"] = "ASD"
enrichTab.full$Set = ss(rownames(enrichTab.full), "_", 3)
enrichTab.full$ID = rownames(enrichTab.full)
enrichTab.full$SetSize = sapply(geneList_present, length)

### save a copy as a supp table
enrichTab.fullOut.fdr.05 = enrichTab.full[ ,c(55, 52:54, 56, 1:51)]
write.csv(enrichTab.fullOut.fdr.05, file = "tables/enrichTab_clinicalGeneLists_DLPFC-cellTypesSplit-SN-LEVEL-fdr05_Apr2020.csv", row.names=FALSE)

## look at enrichment
pMat = enrichTab.fullOut.fdr.05[ , grep("Pval", colnames(enrichTab.fullOut.fdr.05))]
orMat = enrichTab.fullOut.fdr.05[ , grep("OR", colnames(enrichTab.fullOut.fdr.05))]
colnames(pMat) = ss(colnames(pMat), ".Pval")
colnames(orMat) = ss(colnames(orMat), ".OR")
pMat < 0.05 / nrow(pMat)  #  0.001351351
pMat < 0.001
round(-log10(pMat),1)


## SCZD gene sets
enrichTab.fullOut.fdr.05[c("DE_PE_SCZ.Up","DE_BS2_SCZ.Up", "TWAS_PE_SCZ.Up","TWAS_BS2_SCZ.Up"),
                         c(2:5,grep("Pval", colnames(enrichTab.fullOut.fdr.05)))]
    #                 Type Group    Set SetSize   Astro.Pval Excit.ambig.Pval Excit.L2:3.Pval
    # DE_PE_SCZ.Up      DE    PE SCZ.Up    2247 8.222529e-59        1.0000000               1
    # DE_BS2_SCZ.Up     DE   BS2 SCZ.Up      93 4.278843e-02        1.0000000               1
    # TWAS_PE_SCZ.Up  TWAS    PE SCZ.Up     456 7.936539e-01        0.6392323               1
    # TWAS_BS2_SCZ.Up TWAS   BS2 SCZ.Up     368 2.389407e-01        0.6292215               1
    #                 Excit.L3:4.Pval Excit.L4:5.Pval Excit.L5.Pval Excit.L5:6.Pval
    # DE_PE_SCZ.Up          1.0000000       0.2591801     0.2213807       0.5884815
    # DE_BS2_SCZ.Up         1.0000000       1.0000000     0.1743730       1.0000000
    # TWAS_PE_SCZ.Up        1.0000000       1.0000000     0.7894808       0.5441878
    # TWAS_BS2_SCZ.Up       0.6401489       1.0000000     0.3786089       1.0000000
    #                 Excit.L6.broad.Pval Inhib.1.Pval Inhib.2.Pval Inhib.3.Pval Inhib.4.Pval
    # DE_PE_SCZ.Up                      1    0.1247846   0.03537435    0.1535045   0.47308656
    # DE_BS2_SCZ.Up                     1    0.5575512   0.14198264    1.0000000   0.08562455
    # TWAS_PE_SCZ.Up                    1    0.1232565   1.00000000    0.1787307   0.35710588
    # TWAS_BS2_SCZ.Up                   1    0.5667187   0.52785683    1.0000000   1.00000000
    #                 Inhib.5.Pval Inhib.6.Pval Micro.Pval   Oligo.Pval  OPC.Pval
    # DE_PE_SCZ.Up       0.1690452     0.723026  1.0000000 6.547380e-06 0.1813966
    # DE_BS2_SCZ.Up      1.0000000     1.000000  0.4292628 1.000000e+00 0.4045187
    # TWAS_PE_SCZ.Up     1.0000000     1.000000  0.7292709 3.903101e-02 0.5255598
    # TWAS_BS2_SCZ.Up    0.3266752     1.000000  0.8486287 1.417104e-01 0.7269583


enrichTab.fullOut.fdr.05[c("DE_PE_SCZ.Down","DE_BS2_SCZ.Down", "TWAS_PE_SCZ.Down","TWAS_BS2_SCZ.Down"),
                         c(2:5,grep("Pval", colnames(enrichTab.fullOut.fdr.05)))]
    #                   Type Group      Set SetSize   Astro.Pval Excit.ambig.Pval Excit.L2:3.Pval
    # DE_PE_SCZ.Down      DE    PE SCZ.Down    2077 0.0002128647        0.2771107               1
    # DE_BS2_SCZ.Down     DE   BS2 SCZ.Down     132 1.0000000000        1.0000000               1
    # TWAS_PE_SCZ.Down  TWAS    PE SCZ.Down     463 0.7963027289        0.6406838               1
    # TWAS_BS2_SCZ.Down TWAS   BS2 SCZ.Down     402 1.0000000000        0.6310541               1
    #                   Excit.L3:4.Pval Excit.L4:5.Pval Excit.L5.Pval Excit.L5:6.Pval
    # DE_PE_SCZ.Down          0.1195261       0.4268897     0.5266120      0.08631337
    # DE_BS2_SCZ.Down         1.0000000       0.1110539     1.0000000      0.20237169
    # TWAS_PE_SCZ.Down        0.4131792       1.0000000     0.5954856      1.00000000
    # TWAS_BS2_SCZ.Down       0.6501990       1.0000000     0.3890119      0.49940599
    #                   Excit.L6.broad.Pval Inhib.1.Pval Inhib.2.Pval Inhib.3.Pval Inhib.4.Pval
    # DE_PE_SCZ.Down             0.05826803   0.26788758    0.1347153            1   0.04521774
    # DE_BS2_SCZ.Down            1.00000000   1.00000000    1.0000000            1   1.00000000
    # TWAS_PE_SCZ.Down           1.00000000   0.44772716    0.3907894            1   0.36148717
    # TWAS_BS2_SCZ.Down          0.42172359   0.05421873    0.0683707            1   0.05661475
    #                   Inhib.5.Pval Inhib.6.Pval   Micro.Pval   Oligo.Pval     OPC.Pval
    # DE_PE_SCZ.Down      0.00117537    0.1485163 2.036640e-43 6.803246e-26 8.315935e-05
    # DE_BS2_SCZ.Down     1.00000000    1.0000000 1.094018e-28 3.560448e-02 1.666946e-01
    # TWAS_PE_SCZ.Down    1.00000000    0.3720139 2.305730e-01 1.875323e-01 5.262360e-01
    # TWAS_BS2_SCZ.Down   0.35100828    1.0000000 8.540537e-01 7.291697e-01 1.789745e-01



    # Marginal signal is interesting for 'TWAS_BS2_SCZ.Down'
    length(intersect(geneLists.fromST[["TWAS_PE_SCZ.Down"]], geneLists.fromST[["TWAS_BS2_SCZ.Down"]]))
        # 175 genes shared in these ~430 gene sets

# Plot into the custom heatmap
midpoint = function(x) x[-length(x)] + diff(x)/2

customLayerEnrichment = function(enrichTab , groups, xlabs, 
                                 Pthresh = 12, ORcut = -log10(0.05), enrichOnly = FALSE,
                                 #layerHeights = c(0,40,55,75,85,110,120,135),
                                 layerHeights = seq(0,153,by=9),
                                 mypal = c("white", colorRampPalette(brewer.pal(9,"YlOrRd"))(50)), ...) {
  
  wide_p = -log10(enrichTab[groups,grep("Pval", colnames(enrichTab))])
  wide_p[wide_p > Pthresh] = Pthresh
      colnames(wide_p) <- gsub(".Pval","",colnames(wide_p))
  wide_p = t(round(wide_p[,
                          # (From bottom to top on plot)
                          c(rev(names(markers.dlpfc.t.design)))
                          ],2))
  
  wide_or = enrichTab[groups,grep("OR", colnames(enrichTab))]
  colnames(wide_or) <- gsub(".OR","",colnames(wide_or))
  wide_or= round(t(wide_or[,
                           # (From bottom to top on plot)
                           c(rev(names(markers.dlpfc.t.design)))
                           ]),1)
  
  if(enrichOnly) wide_p[wide_or < 1] = 0
  wide_or[wide_p < ORcut] = ""
      # or, if want to print -log10(p's)
      wide_p_2plot <- wide_p
      wide_p_2plot[wide_p < ORcut] = ""
  
  image.plot(x = seq(0,ncol(wide_p),by=1), y = layerHeights, z = as.matrix(t(wide_p)),
             col = mypal,xaxt="n", yaxt="n",xlab = "", ylab="", ...)
  
  axis(2, rev(names(markers.dlpfc.t.design)), at = midpoint(layerHeights),las=1)  # MNT add
  axis(1, rep("", ncol(wide_p)), at = seq(0.5,ncol(wide_p)-0.5))
  
  text(x = seq(0.5,ncol(wide_p)-0.5), y=-1*max(nchar(xlabs))/2, xlabs,
       xpd=TRUE, srt=45,cex=1.5,adj= 1)
  abline(h=layerHeights,v=0:ncol(wide_p))
  text(x = rep(seq(0.5,ncol(wide_p)-0.5),each = nrow(wide_p)), 
       y = rep(midpoint(layerHeights), ncol(wide_p)),
       as.character(wide_or),cex=1.4,font=2)
       # or, for -log10(p)
       #as.character(wide_p_2plot),cex=1.5,font=2)
}


# Print this
pdf("pdfs/exploration/enrichPlots_SN-LEVEL-markers_SCZD-geneSet_DLPFC-cellTypeSplit_heatmap_Apr2020.pdf",w=8)
par(mar=c(8,8,2.5,1), cex.axis=1.2, cex.lab=1.5)
groups =c("DE_PE_SCZ.Up", "DE_PE_SCZ.Down", 
          "DE_BS2_SCZ.Up", "DE_BS2_SCZ.Down", 
          "TWAS_BS2_SCZ.Up", "TWAS_BS2_SCZ.Down", "TWAS_PE_SCZ.Up",
          "TWAS_PE_SCZ.Down")
xlabs = ss(gsub("_SCZ", "", groups), "_", 2)
customLayerEnrichment(enrichTab.full, groups, xlabs, enrichOnly=TRUE)
abline(v=4,lwd=3)
text(x = c(2,6), y = 160, c("SCZD-DE", "SCZD-TWAS"), xpd=TRUE,cex=2,font=2)
dev.off()


pdf("pdfs/exploration/enrichPlots_SN-LEVEL-markers_birnbaum-geneSet_DLPFC-cellTypeSplit_heatmap_Apr2020.pdf",w=8)
par(mar=c(12,8,2.5,1), cex.axis=1, cex.lab=1.5)
groups =grep(enrichTab.full$ID, pattern = "Birnbaum", value=TRUE)
xlabs = ss(groups, "_", 3)
customLayerEnrichment(enrichTab.full, groups,xlabs, enrichOnly=TRUE,
                      breaks = seq(0,12,len = 52))
dev.off()


