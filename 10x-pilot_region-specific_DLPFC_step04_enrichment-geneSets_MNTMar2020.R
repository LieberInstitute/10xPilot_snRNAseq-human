### MNT 10x snRNA-seq workflow: step 03 - marker detection
###   **Region-specific analyses**
###     - (2x) DLPFC samples from: Br5161 & Br5212
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
library(lattice)
library(RColorBrewer)
library(pheatmap)

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
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_DLPFC_n2_MNTFeb2020.rda",
     verbose=T)
    # markers.dlpfc.t.design.log, markers.dlpfc.t.design.countsN, markers.dlpfc.t.design.noAmbig.log, markers.dlpfc.t.design.noAmbig.countsN
    rm(markers.dlpfc.t.design.log, markers.dlpfc.t.design.countsN)
  
sapply(markers.dlpfc.t.design.noAmbig.log, function(x){table(x$FDR<0.05)})
    ##       Astro Excit Inhib Micro Oligo   OPC
    # FALSE 28074 27912 28112 27219 27997 28118
    # TRUE     54   216    16   909   131    10


sapply(markers.dlpfc.t.design.noAmbig.countsN, function(x){table(x$FDR<0.05)})
    ##      Astro Excit Inhib Micro Oligo   OPC
    #FALSE 27470 27467 27980 26384 27428 27888
    #TRUE    658   661   148  1744   700   240


        ## Actually what is this? ====
        load("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/rda/tstats_Human_DLPFC_snRNAseq_Nguyen.Rdata",
             verbose=T)
            # tstats_Human_DLPFC_snRNAseq_Nguyen_topLayer
        
        dim(tstats_Human_DLPFC_snRNAseq_Nguyen_topLayer)  # 692 x 31
        ## not sure if applicable here ====



### The generation of these stats are dapted from:
  # `/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/layer_specificity.R`
  # (and performed in step03 script - "direct limma approach")

load("rdas/markers-stats_DLPFC_n2_manualContrasts_MNTMar2020.rda", verbose=T)
    # eb_contrasts.dlpfc.broad, eb_list.dlpfc.broad, sce.dlpfc.PB
        
## Tabulating significant hits
pvals_contrasts <- eb_contrasts$p.value

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


## Extract the p-values and compute fdrs
pvals0_contrasts <- sapply(eb_list.dlpfc.broad, function(x) {
  x$p.value[, 2, drop = FALSE]
})

fdrs0_contrasts = apply(pvals0_contrasts, 2, p.adjust, "fdr")

## Extract the tstats
t0_contrasts <- sapply(eb_list.dlpfc.broad, function(x) {
  x$t[, 2, drop = FALSE]
})


data.frame(
  'FDRsig' = colSums(apply(pvals0_contrasts, 2, p.adjust, 'fdr') < 0.05),
  'Pval10-6sig' = colSums(pvals0_contrasts < 1e-6),
  'Pval10-8sig' = colSums(pvals0_contrasts < 1e-8)
)

# For only (+) t-stats
data.frame(
  'FDRsig' = colSums(apply(pvals0_contrasts, 2, p.adjust, 'fdr') < 0.05 &
                       t0_contrasts > 0),
  'Pval10-6sig' = colSums(pvals0_contrasts < 1e-6 &
                            t0_contrasts > 0),
  'Pval10-8sig' = colSums(pvals0_contrasts < 1e-8 &
                            t0_contrasts > 0)
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
rownames(t0_contrasts) <- rownames(eb_list.dlpfc.broad[["Astro"]])
rownames(t0_contrasts) <- rowData(sce.dlpfc.PB)$ID[match(rownames(t0_contrasts), rownames(sce.dlpfc.PB))]
table(rowData(sce.dlpfc.PB)$ID == rownames(t0_contrasts)) # all TRUE

## filter for those present in stats
geneList_present = lapply(geneLists.fromST, function(x) {
  x = x[!is.na(x)]
  x[x %in% rownames(t0_contrasts)]
})

unname(data.frame(lengths(geneLists.fromST), lengths(geneList_present)))
    ## not too bad

## do enrichment
enrich_stat_list = eb_list.dlpfc.broad
for (i in seq(along = enrich_stat_list)) {
  #layer = t0_contrasts[, i] > 0 & fdrs0_contrasts[, i] < 0.1
  # or
  layer = t0_contrasts[, i] > 0 & fdrs0_contrasts[, i] < 0.05
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
write.csv(enrichTabOut.fdr.05, file = "tables/enrichTab_clinicalGeneLists_DLPFC-broadClusters-fdr05.csv", row.names=FALSE)
write.csv(enrichTabOut.fdr.10, file = "tables/enrichTab_clinicalGeneLists_DLPFC-broadClusters-fdr10.csv", row.names=FALSE)


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

## summary stats from genes
enrichTab["Gene_SFARI_all",]
enrichTab["Gene_Satterstrom_ASC102.2018",]
enrichTab["Gene_Satterstrom_ASD53",]
enrichTab["Gene_Satterstrom_DDID49",]

## Satterstrom deep dive
sat_102_l2= which(t0_contrasts[,"Layer2"] > 0 & fdrs0_contrasts[,"Layer2"] < 0.1 & 
                    rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_ASC102.2018)
rowData(sce_layer)$gene_name[sat_102_l2]
sat_102_l5= which(t0_contrasts[,"Layer5"] > 0 & fdrs0_contrasts[,"Layer5"] < 0.1 & 
                    rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_ASC102.2018)
rowData(sce_layer)$gene_name[sat_102_l5]

sat_49_l2= which(t0_contrasts[,"Layer2"] > 0 & fdrs0_contrasts[,"Layer2"] < 0.1 & 
                   rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_DDID49)
cat(rowData(sce_layer)$gene_name[sat_49_l2], sep=", ")

sat_53_l5= which(t0_contrasts[,"Layer5"] > 0 & fdrs0_contrasts[,"Layer5"] < 0.1 & 
                   rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_ASD53)
cat(rowData(sce_layer)$gene_name[sat_53_l5], sep=", ")

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
    

################
## make plots ##
################

## make long
enrichLong = reshape2::melt(enrichTab[,c(seq(1,19,by=3),22:26)],id.vars = 8:12)
colnames(enrichLong)[6:7] = c("Layer", "OR")
enrichLong_P = reshape2::melt(enrichTab[,c(seq(2,20,by=3),22:26)],id.vars = 8:12)
identical(enrichLong$ID, enrichLong_P$ID)
enrichLong$P = enrichLong_P$value
enrichLong$Layer = ss(as.character(enrichLong$Layer), "\\.")
enrichLong$ID = factor(enrichLong$ID, levels=rev(rownames(enrichTab)))
enrichLong$Set = factor(enrichLong$Set, levels=unique(rev(enrichTab$Set)))
enrichLong$FDR = p.adjust(enrichLong$P, "fdr")

## what p-value controls FDR?
enrichLongSort = enrichLong[order(enrichLong$P),]
max(enrichLongSort$P[enrichLongSort$FDR < 0.05] )
# 0.01009034

## overall ##
enrichLong$P_thresh = enrichLong$P
enrichLong$P_thresh[enrichLong$P_thresh < 2.2e-16] = 2.2e-16

### ASD focus
enrichLong_ASD = enrichLong[enrichLong$ID %in% 
                              c("Gene_SFARI_all", "Gene_Satterstrom_ASC102.2018",
                                "Gene_Satterstrom_ASD53", "Gene_Satterstrom_DDID49",
                                "DE_PE_ASD.Down", "DE_PE_ASD.Up",
                                "TWAS_PE_ASD.Up", "TWAS_PE_ASD.Down"),]
enrichLong_ASD$ID2 =  as.character(droplevels(enrichLong_ASD$Set))
enrichLong_ASD$ID2[enrichLong_ASD$ID2 == "all"] = "SFARI"
enrichLong_ASD$ID2[enrichLong_ASD$ID2 == "ASC102.2018"] = "ASC102"
enrichLong_ASD$ID2[enrichLong_ASD$ID == "DE_PE_ASD.Up"] = "DE.Up"
enrichLong_ASD$ID2[enrichLong_ASD$ID == "DE_PE_ASD.Down"] = "DE.Down"
enrichLong_ASD$ID2[enrichLong_ASD$ID == "TWAS_PE_ASD.Up"] = "TWAS.Up"
enrichLong_ASD$ID2[enrichLong_ASD$ID == "TWAS_PE_ASD.Down"] = "TWAS.Down"
enrichLong_ASD$ID2 = factor(enrichLong_ASD$ID2, unique(enrichLong_ASD$ID2))

enrichLong_ASD$LayerFac = factor(as.character(enrichLong_ASD$Layer), 
                                 c("WM", paste0("Layer", 6:1)))
enrichLong_ASD = enrichLong_ASD[order(enrichLong_ASD$ID2, enrichLong_ASD$LayerFac),]

### custom heatmap

midpoint = function(x) x[-length(x)] + diff(x)/2

customLayerEnrichment = function(enrichTab , groups, xlabs, 
                                 Pthresh = 12, ORcut = 3, enrichOnly = FALSE,
                                 layerHeights = c(0,40,55,75,85,110,120,135),
                                 mypal = c("white", colorRampPalette(brewer.pal(9,"YlOrRd"))(50)), ...) {
  
  wide_p = -log10( enrichTab[groups,grep("Pval", colnames(enrichTab))])
  wide_p[wide_p > Pthresh] = Pthresh
  wide_p = t(round(wide_p[,
                          c("WM.Pval", "Layer6.Pval", "Layer5.Pval", "Layer4.Pval", "Layer3.Pval","Layer2.Pval", "Layer1.Pval")],2))
  
  wide_or = enrichTab[groups,grep("OR", colnames(enrichTab))]
  wide_or= round(t(wide_or[,
                           c("WM.OR", "Layer6.OR", "Layer5.OR", "Layer4.OR", "Layer3.OR", "Layer2.OR", "Layer1.OR")]),1)
  if(enrichOnly) wide_p[wide_or < 1] = 0
  wide_or[wide_p < ORcut] = ""
  
  image.plot(x = seq(0,ncol(wide_p),by=1), y = layerHeights, z = as.matrix(t(wide_p)),
             col = mypal,xaxt="n", yaxt="n",xlab = "", ylab="", ...)
  axis(2, c("WM", paste0("L", 6:1)), at = midpoint(layerHeights),las=1)
  axis(1, rep("", ncol(wide_p)), at = seq(0.5,ncol(wide_p)-0.5))
  text(x = seq(0.5,ncol(wide_p)-0.5), y=-1*max(nchar(xlabs))/2, xlabs,
       xpd=TRUE, srt=45,cex=2,adj= 1)
  abline(h=layerHeights,v=0:ncol(wide_p))
  text(x = rep(seq(0.5,ncol(wide_p)-0.5),each = nrow(wide_p)), 
       y = rep(midpoint(layerHeights), ncol(wide_p)),
       as.character(wide_or),cex=1.5,font=2)
}

pdf("pdf/asd_geneSet_heatmap.pdf",w=6)
par(mar=c(8,4.5,2.5,1), cex.axis=2,cex.lab=2)
groups = unique(as.character(enrichLong_ASD$ID))[1:6]
xlabs  = as.character(enrichLong_ASD$ID2[match(groups, enrichLong_ASD$ID)])
customLayerEnrichment(enrichTab, groups,xlabs, enrichOnly=TRUE)
abline(v=4,lwd=3)
text(x = 3, y = 142, c("ASD"), xpd=TRUE,cex=2.5,font=2)

dev.off()


pdf("pdf/sczd_geneSet_heatmap.pdf",w=8)
par(mar=c(8,4.5,2.5,1), cex.axis=2,cex.lab=2)

groups =c("DE_PE_SCZ.Up", "DE_PE_SCZ.Down", 
          "DE_BS2_SCZ.Up", "DE_BS2_SCZ.Down", 
          "TWAS_BS2_SCZ.Up", "TWAS_BS2_SCZ.Down", "TWAS_PE_SCZ.Up",
          "TWAS_PE_SCZ.Down")
xlabs = ss(gsub("_SCZ", "", groups), "_", 2)
customLayerEnrichment(enrichTab, groups,xlabs, enrichOnly=TRUE)
abline(v=4,lwd=3)
text(x = c(2,6), y = 142, c("SCZD-DE", "SCZD-TWAS"), xpd=TRUE,cex=2.5,font=2)
dev.off()


pdf("pdf/suppXX_birnbaum_geneSet_heatmap.pdf",w=8)
par(mar=c(12,5.5,2.5,1), cex.axis=2,cex.lab=2)
groups =grep(enrichTab$ID, pattern = "Birnbaum", value=TRUE)
xlabs = ss(groups, "_", 3)
customLayerEnrichment(enrichTab, groups,xlabs, enrichOnly=TRUE,
                      breaks = seq(0,12,len = 51))
dev.off()

















