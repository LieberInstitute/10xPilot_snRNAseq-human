### for MAGMA with LIBD 10x pilot analyses
#     - Comparing some gene-level stats vs.
#       those from Liu, et al.
#     - Color by cell class marker label
# MNT Jul2021 ============================

library(readr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(fields)
library(jaffelab)
library(SingleCellExperiment)
library(readxl)

source("../plotExpressionCustom.R")

### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

## ===

## Read in amy SCE as example for gene symbol
load("../rdas/revision/regionSpecific_Amyg-n5_cleaned-combined_SCE_MNT2021.rda", verbose=T)
    # sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo, annotationTab.amy, cell_colors.amy
    rm(chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, annotationTab.amy)

## Read in MAGMA gene-level stats
geneStats <- list()
for(i in c("AgeofInitiation","CigarettesPerDay","DrinksPerWeek","SmokingCessation","SmokingInitiation")){
  geneStats[[i]] <- read.table(paste0("./SNP_Data/Liu-etal_Addiction_",i,"_10xPilotGenes_snp-wise.genes.out"), header=T)
  rownames(geneStats[[i]]) <- geneStats[[i]]$GENE
  geneStats[[i]]$Symbol.uniq <- rowData(sce.amy)$Symbol.uniq[match(geneStats[[i]]$GENE,
                                                         rowData(sce.amy)$gene_id)]
}

colnames(geneStats[[1]])
    # [1] "GENE"        "CHR"         "START"       "STOP"        "NSNPS"      
    # [6] "NPARAM"      "N"           "ZSTAT"       "P"           "Symbol.uniq"


setdiff(rownames(geneStats[["CigarettesPerDay"]]), rownames(geneStats[["SmokingCess"]]))
    #[1] "ENSG00000186092" "ENSG00000238009" "ENSG00000224698" "ENSG00000278175"
    #[5] "ENSG00000275450" "ENSG00000279851"
    # which are, (and added others):
        #                    gene_id gene_version   gene_name    gene_source
        #                <character>  <character> <character>    <character>
        # OR4F5      ENSG00000186092            6       OR4F5 ensembl_havana
        # AL627309.1 ENSG00000238009            6  AL627309.1 ensembl_havana
        # AL390038.1 ENSG00000224698            1  AL390038.1         havana
        # GLIDR      ENSG00000278175            3       GLIDR         havana
        # AL845472.1 ENSG00000275450            1  AL845472.1         havana
        # CR382287.2 ENSG00000279851            2  CR382287.2         havana
        # AL512625.2 ENSG00000229422            2  AL512625.2         havana
        # CR382287.2 ENSG00000279851            2  CR382287.2         havana

# Clean up by taking intersecting and match rownames
keep <- intersect(rownames(geneStats[["AgeofInitiation"]]), rownames(geneStats[["CigarettesPerDay"]]))
keep <- intersect(rownames(geneStats[["DrinksPerWeek"]]), keep) # of length 29006
keep <- intersect(rownames(geneStats[["SmokingCessation"]]), keep) # of length 29003
keep <- intersect(rownames(geneStats[["SmokingInitiation"]]), keep) # of length 29001
    # 29000, final

for(i in names(geneStats)){
  geneStats[[i]] <- geneStats[[i]][keep, ]
}

# Shorten names for graphics
names(geneStats) <- c("AgeSmk", "CigDay", "DrnkWk", "SmkCes", "SmkInit")

## Read in PASCAL stats ===
liu.pascal.stats <- as.data.frame(read_xlsx(path="../MAGMA/GWAS_Results/Liu-etal_smoking-alcohol-GWAS_SuppTables.xlsx",
                                            sheet=20, skip=1, col_names=T))
colnames(liu.pascal.stats) <- gsub(" ",".", colnames(liu.pascal.stats))
# [1] "Phenotype"                           "Chromosome"                         
# [3] "Start"                               "End"                                
# [5] "Strand"                              "Gene.ID"                            
# [7] "Gene.Symbol"                         "Number.of.Variants.Included.in.Test"
# [9] "P-Value"                             "P-Value.Estimation.Outcome"   

dim(liu.pascal.stats)
    # [1] 1511   10

# Out of curiosity
table(liu.pascal.stats$Gene.Symbol %in% rowData(sce.amy)$gene_name)
    # FALSE  TRUE 
    #   289  1222

## Convert Entrez to EnsemblID ===
library(biomaRt)
# ensembl.db <- useMart("ensembl")
# listDatasets(ensembl.db)[grep("hsapiens", listDatasets(ensembl.db)$dataset), ]
#     #                  dataset              description    version
#     # 80 hsapiens_gene_ensembl Human genes (GRCh38.p13) GRCh38.p13
# ensembl.db <- useDataset("hsapiens_gene_ensembl",mart=ensembl.db)

# Simpler
ensembl.db <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

attributes.hs <- listAttributes(ensembl.db)
attributes.hs[grep("entrez",attributes.hs$name), ]
#                      name                                 description         page
# 61  entrezgene_trans_name               EntrezGene transcript name ID feature_page
# 81 entrezgene_description NCBI gene (formerly Entrezgene) description feature_page
# 82   entrezgene_accession   NCBI gene (formerly Entrezgene) accession feature_page
# 83          entrezgene_id          NCBI gene (formerly Entrezgene) ID feature_page

length(unique(liu.pascal.stats$Gene.ID))
    # 1280

genes.ref <- getBM(c("entrezgene_id", "ensembl_gene_id", "hgnc_symbol", "uniprot_gn_symbol"),
                   filters="entrezgene_id",
                   values=unique(liu.pascal.stats$Gene.ID),
                   mart=ensembl.db)
    # of dim 1367 x 4
length(unique(genes.ref$entrezgene_id))   # 1217
length(unique(genes.ref$ensembl_gene_id)) # 1283

table(unique(liu.pascal.stats$Gene.ID) %in% genes.ref$entrezgene_id)  # 1217

# Lost from conversion:
    liu.pascal.stats$Gene.Symbol[!(liu.pascal.stats$Gene.ID %in% genes.ref$entrezgene_id)]
        # [1] "HTT-AS1"        "NOP14-AS1"      "LOC646938"      "MIA-RAB4B"      "RAB4B-EGLN2"   
        # [6] "ARIH2OS"        "LAMB2P1"        "C16orf97"       "TSPY26P"        "NOP14-AS1"     
        # [11] "MIA-RAB4B"      "RAB4B-EGLN2"    "LAMB2P1"        "RTEL1-TNFRSF6B" "ARIH2OS"       
        # [16] "LOC646938"      "KDM4A-AS1"      "C10orf32-AS3MT" "LINC00678"      "BRWD1-IT2"     
        # [21] "TRAF3IP2-AS1"   "TFAMP1"         "C22orf26"       "LOC100131289"   "TOB2P1"        
        # [26] "ZSCAN12P1"      "DLX6-AS1"       "MYH16"          "LINC00899"      "LOC388906"     
        # [31] "PCCA-AS1"       "FLJ30838"       "NARR"           "LOC285819"      "DIAPH3-AS2"    
        # [36] "LOC100270746"   "CELF2-AS2"      "VN1R10P"        "SNORD47"        "SNORD80"       
        # [41] "OTX2-AS1"       "SNORD79"        "GAS5-AS1"       "SNORD81"        "SNORD78"       
        # [46] "SNORD44"        "SNORD76"        "SNORD77"        "SNORD74"        "LAMB2P1"       
        # [51] "SNORD75"        "TRHDE-AS1"      "LOC286186"      "USP3-AS1"       "RTEL1-TNFRSF6B"
        # [56] "LOC646903"      "HNRNPUL2-BSCL2" "TSPY26P"        "MAPT-AS1"       "LOC401127"     
        # [61] "LOC440354"      "BRE-AS1"        "LOC613038"      "LOC606724"      "LOC644172"     
        # [66] "LRRC37A4P"      "FTH1P3"         "MYH16"          "SNORD47"        "SNORD80"       
        # [71] "SNORD79"        "SNORD78"        "SNORD44"        "GAS5-AS1"       "SNORD81"       
        # [76] "SNORD77"        "SNORD76"        "SNORD74"        "SNORD75"        "SNORD91B"      
        # [81] "SNORD91A"       "MIR181A2HG"     "SNORD10"        "ZRANB2-AS1"

# Subset those so can match to ensembl ID in MAGMA stats
liu.pascal.stats <- liu.pascal.stats[liu.pascal.stats$Gene.ID %in% genes.ref$entrezgene_id, ]
    # 1511 rows to 1427
liu.pascal.stats$ensemblID <- genes.ref$ensembl_gene_id[match(liu.pascal.stats$Gene.ID,
                                                              genes.ref$entrezgene_id)]

# Now how many in MAGMA gene-level results?
table(unique(liu.pascal.stats$ensemblID) %in% geneStats[[1]]$GENE)
    # FALSE  TRUE 
    #   191  1025 
    
    # For reference maybe:
    notKept.more <- unique(liu.pascal.stats$Gene.Symbol[!(liu.pascal.stats$ensemblID %in% geneStats[[1]]$GENE)])
    # a lot of miR's:
    length(grep("MIR", notKept.more)) # 82
    
# Subset that too
liu.pascal.stats <- liu.pascal.stats[liu.pascal.stats$ensemblID %in% geneStats[[1]]$GENE, ]

liu.pascal.list <- splitit(liu.pascal.stats$Phenotype)
    # Age of Initiation of Regular Smoking                   Cigarettes per Day 
    #                                   18                                  121 
    #                      Drinks per Week                    Smoking Cessation 
    #                                  266                                  120 
    #                   Smoking Initiation 
    #                                  675


names(liu.pascal.list) <- c("AgeSmk", "CigDay", "DrnkWk", "SmkCes", "SmkInit")
liu.pascal.list <- lapply(liu.pascal.list, function(x){liu.pascal.stats[x, ]})

# Pull in MAGMA z-scores for those that are matching and probit-transform PASCAL p-values
for(i in names(liu.pascal.list)){
  # Take that 'P-Value.Estimation.Outcome' == "DAVIES_SUCCESS" 
  liu.pascal.list[[i]] <- liu.pascal.list[[i]][liu.pascal.list[[i]][ ,"P-Value.Estimation.Outcome"] == "DAVIES_SUCCESS", ]
  
  # First make PASCAL z's with probit
  liu.pascal.list[[i]]$z.pascal <- qnorm(liu.pascal.list[[i]][ ,"P-Value"], mean=0, sd=1, lower.tail=F)
}

sapply(liu.pascal.list, dim)
    # AgeSmk CigDay DrnkWk SmkCes SmkInit
    # [1,]     18     96    238    112     601
    # [2,]     12     12     12     12      12    good

## Save this for now ===
Readme <- "These stats have been subsetted for intersecting genes and can now be cross-compared."
save(geneStats, liu.pascal.list, Readme,
     file="zForRef_Liu-etal-PASCAL-gene-results_subsetted_MNT.rda")

# Now bring in corresponding MAGMA gene-level stats for 
for(i in names(geneStats)){
  liu.pascal.list[[i]]$z.magma.gene <- geneStats[[i]]$ZSTAT[match(liu.pascal.list[[i]]$ensemblID,
                                                                  geneStats[[i]]$GENE)]
}



plot(liu.pascal.list[["AgeSmk"]]$z.pascal, liu.pascal.list[["AgeSmk"]]$z.magma.gene, pch=16)
    # pearson's cor: [1] 0.2902735
plot(liu.pascal.list[["CigDay"]]$z.pascal, liu.pascal.list[["CigDay"]]$z.magma.gene, pch=16)
    # pearson's cor: [1] 0.6945395
plot(liu.pascal.list[["DrnkWk"]]$z.pascal, liu.pascal.list[["DrnkWk"]]$z.magma.gene, pch=16)
    # pearson's cor: [1] 0.6740618
plot(liu.pascal.list[["SmkCes"]]$z.pascal, liu.pascal.list[["SmkCes"]]$z.magma.gene, pch=16)
    # pearson's cor: [1] 0.4160359
plot(liu.pascal.list[["SmkInit"]]$z.pascal, liu.pascal.list[["SmkInit"]]$z.magma.gene, pch=16)
    # pearson's cor: [1] 0.4755282




## Print these ===
pdf("graphics/MAGMA-vs-PASCAL-gene-levelResults-from-Liu-etal2019_MNT2021.pdf")
for(i in names(geneStats)){
  p <- ggplot(data=liu.pascal.list[[i]], aes(x=z.pascal, y=z.magma.gene)) + geom_point(size=3.5, alpha=0.8) +
    geom_smooth(method = "lm", se=FALSE) +
    labs(title = paste0(i,": MAGMA vs PASCAL-significant loci (from Liu, et al. 2019)"),
         subtitle = paste0("Pearson's r = ", round(cor(liu.pascal.list[[i]]$z.pascal, liu.pascal.list[[i]]$z.magma.gene, method="pearson"),3)),
         x = "PASCAL gene-level z-transformed (significant) p-values",
         y = "Corresponding MAGMA z-scores") +
    theme(plot.title = element_text(size=12, face="bold"),
          axis.title=element_text(size=10,face="bold"))
  print(p)
}
dev.off()


### Any markers in the PASCAL x MAGMA gene sets? =============
  #   Want to use the strict, pairwise marker tests, so that can say that any given [marker]
  #   gene is specific to X cell class

### NAc ===
load("../rdas/revision/markers-stats_NAc-n8_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.nac.t.pw, markers.nac.t.1vAll, medianNon0.nac

sapply(markers.nac.t.pw, function(x){table(x$FDR<0.05 & x$non0median==TRUE)["TRUE"]})
    #  Astro_A.TRUE     Astro_B.TRUE     Inhib_A.TRUE     Inhib_B.TRUE     Inhib_C.TRUE 
    #            97              250               78               40               69 
    #  Inhib_D.TRUE     Inhib_E.TRUE  Macrophage.TRUE       Micro.TRUE Micro_resting.NA 
    #            40               70              154              227               NA 
    # MSN.D1_A.TRUE    MSN.D1_B.TRUE    MSN.D1_C.TRUE    MSN.D1_D.TRUE    MSN.D1_E.TRUE 
    #            12               16               15               21               38 
    # MSN.D1_F.TRUE    MSN.D2_A.TRUE    MSN.D2_B.TRUE    MSN.D2_C.TRUE    MSN.D2_D.TRUE 
    #            33               16                6               99               41 
    #  Oligo_A.TRUE     Oligo_B.TRUE         OPC.TRUE     OPC_COP.TRUE 
    #           213               59               67               23 
sum(sapply(markers.nac.t.pw, function(x){table(x$FDR<0.05 & x$non0median==TRUE)["TRUE"]}), na.rm=T)
    # [1] 1684

markers.nac.specific <- data.frame()
for(i in names(markers.nac.t.pw)){
  markers.temp <- as.data.frame(markers.nac.t.pw[[i]])[ ,c("p.value", "FDR", "non0median")]
  markers.temp$gene <- rownames(markers.temp)
  markers.temp$cellType <- i
  markers.nac.specific <- rbind(markers.nac.specific,
                                markers.temp[markers.temp$FDR<0.05 & markers.temp$non0median==TRUE,
                                             c("p.value", "FDR", "cellType","gene")])
}

table(markers.nac.specific$cellType)
length(unique(markers.nac.specific$gene)) # [1] 1684 - good.

# Add ensemblID so can match, and color
markers.nac.specific$geneID <- rowData(sce.amy)$gene_id[match(markers.nac.specific$gene,
                                                              rowData(sce.amy)$Symbol.uniq)]

# Make same plots, coloring genes by their cell-type-specific color
pdf("graphics/MAGMA-vs-PASCAL-gene-levelResults-from-Liu-etal2019_NAcMarkers_MNT2021.pdf")
for(i in names(geneStats)){
  
  pascal.temp <- liu.pascal.list[[i]]
  pascal.temp$cellType <- markers.nac.specific$cellType[match(pascal.temp$ensemblID, markers.nac.specific$geneID)]
  pascal.temp$cellType <- ifelse(is.na(pascal.temp$cellType), "nonspecific", pascal.temp$cellType)
  
  print(
    ggplot(data=pascal.temp, aes(x=z.pascal, y=z.magma.gene)) +
      geom_point(aes(color=cellType), size=ifelse(pascal.temp$cellType=="nonspecific",2.5,5.0),
                 alpha=ifelse(pascal.temp$cellType=="nonspecific",0.1,0.8)) +
      scale_color_manual(values=c(cell_colors.nac, nonspecific="#333333")) +
      geom_smooth(method = "lm", se=FALSE) +
      labs(title = paste0(i," (Liu, et al. 2019): NAc cell class markers"),
           subtitle = paste0("Pearson's r = ", round(cor(liu.pascal.list[[i]]$z.pascal, liu.pascal.list[[i]]$z.magma.gene, method="pearson"),3)),
           x = "PASCAL gene-level z-scores",
           y = "Corresponding MAGMA z-scores") +
      theme(plot.title = element_text(size=16, face="bold"),
            axis.title=element_text(size=15,face="bold"),
            axis.text = element_text(size=12))
  )
  
  # With gene labels
  print(
    ggplot(data=pascal.temp, aes(x=z.pascal, y=z.magma.gene)) +
      geom_point(aes(color=cellType), size=ifelse(pascal.temp$cellType=="nonspecific",2.5,5.0),
                 alpha=ifelse(pascal.temp$cellType=="nonspecific",0.1,0.8)) +
      geom_text(label=ifelse(!pascal.temp$cellType=="nonspecific",pascal.temp$Gene.Symbol,NA), hjust=0, check_overlap=F,
                fontface="italic", size=3) + #, position=position_jitter(width=0.05,height=0.05)) +
      scale_color_manual(values=c(cell_colors.nac, nonspecific="#333333")) +
      geom_smooth(method = "lm", se=FALSE) +
      labs(title = paste0(i," (Liu, et al. 2019): NAc cell class markers"),
           subtitle = paste0("Pearson's r = ", round(cor(liu.pascal.list[[i]]$z.pascal, liu.pascal.list[[i]]$z.magma.gene, method="pearson"),3)),
           x = "PASCAL gene-level z-scores",
           y = "Corresponding MAGMA z-scores") +
      theme(plot.title = element_text(size=16, face="bold"),
            axis.title=element_text(size=15,face="bold"),
            axis.text = element_text(size=12))
  )
}
dev.off()



### MSN.D1_C is the strongest associating with 'SmkInit'
  # There're only four pairwise-specific D1_C markers printed - might be too strict
markers.nac.enriched <- lapply(markers.nac.t.1vAll, function(x){x[[2]]})
markerList.nac <- lapply(markers.nac.enriched, function(x){rownames(x)[x$log.FDR < log(1e-6) & x$non0median==T]})

D1_C.markers.enriched <- markerList.nac[["MSN.D1_C"]]

table(D1_C.markers.enriched %in% liu.pascal.list[["SmkInit"]]$Gene.Symbol &
        D1_C.markers.enriched %in% geneStats[["SmkInit"]]$Symbol.uniq[geneStats[["SmkInit"]]$ZSTAT>3])
    #FALSE  TRUE 
    #  777    37

printThese <- D1_C.markers.enriched[D1_C.markers.enriched %in% liu.pascal.list[["SmkInit"]]$Gene.Symbol &
                                      D1_C.markers.enriched %in% geneStats[["SmkInit"]]$Symbol.uniq[geneStats[["SmkInit"]]$ZSTAT>3]]
    # [1] "KCNJ6"      "RPL23A"     "KLHL29"     "DAB1"       "TRPC4"      "NYAP2"     
    # [7] "GRM8"       "RPL6"       "ZFHX3"      "CHD3"       "OPCML"      "PPP1R1B"   
    # [13] "FAT3"       "SMAD3"      "CACNA1D"    "AMBRA1"     "NBEAL1"     "CNTN4"     
    # [19] "CAMKV"      "CELF2"      "AFF3"       "DCC"        "RBFOX1"     "GRID2"     
    # [25] "RUNX1T1"    "NOL4"       "CDH8"       "DGKZ"       "CAMTA1"     "NECAB1"    
    # [31] "XKR6"       "SSBP4"      "ST6GALNAC3" "BCL11B"     "INPP4B"     "CTNNA2"    
    # [37] "WDR12" 

plotExpressionCustom(sce.nac, anno_name="cellType", features_name="MSN.D1_C-containing 'SmkInit'", ncol=3,
                     features=head(printThese,9), scales="free_y") +
  scale_color_manual(values=cell_colors.nac)


# Print an iteration for Fig 4C, highlighting the Bonf-significant classes' PW markers
i <- "SmkInit"
pascal.temp <- liu.pascal.list[[i]]
pascal.temp$cellType <- markers.nac.specific$cellType[match(pascal.temp$ensemblID, markers.nac.specific$geneID)]
pascal.temp$cellType <- ifelse(is.na(pascal.temp$cellType), "nonspecific", pascal.temp$cellType)

pdf("graphics/MAGMA-vs-PASCAL-gene-levelResults-from-Liu-etal2019_NAc-SmkInit-highlight_MNT2021.pdf")
set.seed(901)
print(
  ggplot(data=pascal.temp, aes(x=z.pascal, y=z.magma.gene)) +
    geom_point(aes(color=cellType), size=ifelse(pascal.temp$cellType=="nonspecific",2.5,5.0),
               alpha=ifelse(pascal.temp$cellType=="nonspecific",0.1,0.8)) +
    geom_text(label=ifelse(pascal.temp$cellType %in% c("Inhib_A","MSN.D1_C","MSN.D2_C","OPC"),
                           pascal.temp$Gene.Symbol,NA), hjust=0, check_overlap=F,
              fontface="italic", size=3, position=position_jitter(width=0.03,height=0.07)) +
    scale_color_manual(values=c(cell_colors.nac, nonspecific="#333333")) +
    geom_smooth(method = "lm", se=FALSE) +
    labs(title = paste0(i," (Liu, et al. 2019): gene-level scores"),
         subtitle = "Pairwise-significant markers for Bonferroni-significant NAc cell classes",
         x = "PASCAL gene-level z-scores",
         y = "Corresponding MAGMA z-scores",
         colour="Cell type") +
    theme(plot.title = element_text(size=16, face="bold"),
          axis.title=element_text(size=15,face="bold"),
          axis.text = element_text(size=12))
)
dev.off()



### MSN.D1_E for 'CigDay'
D1_E.markers.enriched <- markerList.nac[["MSN.D1_E"]]

table(D1_E.markers.enriched %in% liu.pascal.list[["CigDay"]]$Gene.Symbol &
        D1_E.markers.enriched %in% geneStats[["CigDay"]]$Symbol.uniq[geneStats[["CigDay"]]$ZSTAT>3])
    # FALSE  TRUE 
    #  1023    9

printThese <- D1_E.markers.enriched[D1_E.markers.enriched %in% liu.pascal.list[["CigDay"]]$Gene.Symbol &
                                      D1_E.markers.enriched %in% geneStats[["CigDay"]]$Symbol.uniq[geneStats[["CigDay"]]$ZSTAT>3]]

plotExpressionCustom(sce.nac, anno_name="cellType", features_name="MSN.D1_E-containing 'CigDay'", ncol=3,
                     features=head(printThese,9), scales="free_y") +
  scale_color_manual(values=cell_colors.nac)

    ## Not including genes in the PASCAL subsets looked better for MSN.D1_C
     #    -> Might look better here too


### AMY cell classes markers' enrichment ===================================
load("../rdas/revision/markers-stats_Amyg-n5_findMarkers-SN-LEVEL_MNT2021.rda", verbose=T)
    # markers.amy.t.pw, markers.amy.wilcox.block, markers.amy.t.1vAll, medianNon0.amy
    rm(markers.amy.wilcox.block, medianNon0.amy)

sapply(markers.amy.t.pw, function(x){table(x$FDR<0.05 & x$non0median==TRUE)["TRUE"]})
# Astro_A.TRUE Astro_B.TRUE    Endo.TRUE Excit_A.TRUE Excit_B.TRUE Excit_C.TRUE 
#          326            3          180           59          173          142 
# Inhib_A.TRUE Inhib_B.TRUE Inhib_C.TRUE Inhib_D.TRUE Inhib_E.TRUE Inhib_F.TRUE 
#           26           88           14          131          147           43 
# Inhib_G.TRUE Inhib_H.TRUE   Micro.TRUE   Mural.TRUE   Oligo.TRUE     OPC.TRUE 
#            4           79          284          172          306          172 
#   Tcell.TRUE 
#          112
sum(sapply(markers.amy.t.pw, function(x){table(x$FDR<0.05 & x$non0median==TRUE)["TRUE"]}), na.rm=T)
# [1] 2461

markers.amy.specific <- data.frame()
for(i in names(markers.amy.t.pw)){
  markers.temp <- as.data.frame(markers.amy.t.pw[[i]])[ ,c("p.value", "FDR", "non0median")]
  markers.temp$gene <- rownames(markers.temp)
  markers.temp$cellType <- i
  markers.amy.specific <- rbind(markers.amy.specific,
                                markers.temp[markers.temp$FDR<0.05 & markers.temp$non0median==TRUE,
                                             c("p.value", "FDR", "cellType","gene")])
}

table(markers.amy.specific$cellType)
length(unique(markers.amy.specific$gene)) # [1] 2461 - good.

# Add ensemblID so can match, and color
markers.amy.specific$geneID <- rowData(sce.amy)$gene_id[match(markers.amy.specific$gene,
                                                              rowData(sce.amy)$Symbol.uniq)]

# Make same plots, coloring genes by their cell-type-specific color
pdf("graphics/MAGMA-vs-PASCAL-gene-levelResults-from-Liu-etal2019_AMY-Markers_MNT2021.pdf")
for(i in names(geneStats)){
  
  pascal.temp <- liu.pascal.list[[i]]
  pascal.temp$cellType <- markers.amy.specific$cellType[match(pascal.temp$ensemblID, markers.amy.specific$geneID)]
  pascal.temp$cellType <- ifelse(is.na(pascal.temp$cellType), "nonspecific", pascal.temp$cellType)
  
  print(
    ggplot(data=pascal.temp, aes(x=z.pascal, y=z.magma.gene)) +
      geom_point(aes(color=cellType), size=ifelse(pascal.temp$cellType=="nonspecific",2.5,5.0),
                 alpha=ifelse(pascal.temp$cellType=="nonspecific",0.1,0.8)) +
      scale_color_manual(values=c(cell_colors.amy, nonspecific="#333333")) +
      geom_smooth(method = "lm", se=FALSE) +
      labs(title = paste0(i," (Liu, et al. 2019): AMY cell class markers"),
           subtitle = paste0("Pearson's r = ", round(cor(liu.pascal.list[[i]]$z.pascal, liu.pascal.list[[i]]$z.magma.gene, method="pearson"),3)),
           x = "PASCAL gene-level z-scores",
           y = "Corresponding MAGMA z-scores") +
      theme(plot.title = element_text(size=16, face="bold"),
            axis.title=element_text(size=15,face="bold"),
            axis.text = element_text(size=12))
  )
  
  # With gene labels
  print(
    ggplot(data=pascal.temp, aes(x=z.pascal, y=z.magma.gene)) +
      geom_point(aes(color=cellType), size=ifelse(pascal.temp$cellType=="nonspecific",2.5,5.0),
                 alpha=ifelse(pascal.temp$cellType=="nonspecific",0.1,0.8)) +
      geom_text(label=ifelse(!pascal.temp$cellType=="nonspecific",pascal.temp$Gene.Symbol,NA), hjust=0, check_overlap=F,
                fontface="italic", size=3) + #, position=position_jitter(width=0.05,height=0.05)) +
      scale_color_manual(values=c(cell_colors.amy, nonspecific="#333333")) +
      geom_smooth(method = "lm", se=FALSE) +
      labs(title = paste0(i," (Liu, et al. 2019): AMY cell class markers"),
           subtitle = paste0("Pearson's r = ", round(cor(liu.pascal.list[[i]]$z.pascal, liu.pascal.list[[i]]$z.magma.gene, method="pearson"),3)),
           x = "PASCAL gene-level z-scores",
           y = "Corresponding MAGMA z-scores") +
      theme(plot.title = element_text(size=16, face="bold"),
            axis.title=element_text(size=15,face="bold"),
            axis.text = element_text(size=12))
  )
}
dev.off()


## 'CigDay'
sce.amy <- sce.amy[ ,-grep("drop.", sce.amy$cellType)]
sce.amy$cellType <- droplevels(sce.amy$cellType)

plotExpressionCustom(sce.amy, anno_name="cellType", features_name="'CigDay'-containing AMY cell class", ncol=3,
                     features=c("CLU","GRK4","SEZ6","BTBD3","RASGRF1","ADD1","DRD2","SOX6"), scales="free_y") +
  scale_color_manual(values=cell_colors.amy)


## 'DrnkWk' - highlight Inhib_E & Inhib_G '_enriched' markers
markers.amy.enriched <- lapply(markers.amy.t.1vAll, function(x){x[[2]]})
markerList.amy <- lapply(markers.amy.enriched, function(x){rownames(x)[x$log.FDR < log(1e-6) & x$non0median==T]})

Inhib_E.enriched <- markerList.amy[["Inhib_E"]]
Inhib_G.enriched <- markerList.amy[["Inhib_G"]]

# Which are more enriched in one or the other?
intersect(setdiff(Inhib_G.enriched, Inhib_E.enriched), liu.pascal.list[["DrnkWk"]]$Gene.Symbol)
    # [1] "RUNX1T1"
intersect(setdiff(Inhib_E.enriched, Inhib_G.enriched), liu.pascal.list[["DrnkWk"]]$Gene.Symbol)
    #  [1] "DRD2"      "TENM2"     "SLC24A2"   "RAB11FIP4" "RPS6KA5"   "CELF1"    
    # [7] "POU2F2"    "RABGAP1L"  "PHC2"      "RAI1"      "AUTS2"     "FTO"      
    # [13] "INPP4B"    "NR6A1"

intersect(intersect(Inhib_G.enriched, Inhib_E.enriched), liu.pascal.list[["DrnkWk"]]$Gene.Symbol)
    # [1] "MYT1L"   "GRIK5"   "ST8SIA3" "SULT1A1"

    # The above are those that are _in_ the Liu et al. PASCAL stats - what about original input?
    length(setdiff(Inhib_G.enriched, Inhib_E.enriched)) # 38 specific to Inhib_G
    length(setdiff(Inhib_E.enriched, Inhib_G.enriched)) # 979 specific to Inhib_E
    length(intersect(Inhib_G.enriched, Inhib_E.enriched)) # 267 shared enriched markers
        # So most of 'Inhib_G' were captured in the input for 'Inhib_E'
    head(intersect(Inhib_G.enriched, Inhib_E.enriched),n=20)
        #  [1] "AC008574.1"  "CRYM"        "STARD4-AS1"  "C3orf67-AS1" "AC093523.1" 
        #  [6] "KCNH4"       "YPEL1"       "GPR149"      "MME"         "ST6GALNAC5" 
        # [11] "GABRA5"      "LAMB1"       "LDB2"        "CNTNAP3B"    "LINC01250"  
        # [16] "GNAL"        "ITPR1"       "CDH9"        "GRM1"        "SLC47A1"

i <- "DrnkWk"
pascal.temp <- liu.pascal.list[[i]]

    inhib.g <- intersect(setdiff(Inhib_G.enriched, Inhib_E.enriched), liu.pascal.list[["DrnkWk"]]$Gene.Symbol)
    inhib.e <- intersect(setdiff(Inhib_E.enriched, Inhib_G.enriched), liu.pascal.list[["DrnkWk"]]$Gene.Symbol)
    both <- intersect(intersect(Inhib_G.enriched, Inhib_E.enriched), liu.pascal.list[["DrnkWk"]]$Gene.Symbol)

    pascal.temp$cellType <- ifelse(pascal.temp$Gene.Symbol %in% inhib.g, "Inhib_G", NA)
    pascal.temp$cellType <- ifelse(pascal.temp$Gene.Symbol %in% inhib.e, "Inhib_E", pascal.temp$cellType)
    pascal.temp$cellType <- ifelse(pascal.temp$Gene.Symbol %in% both, "shared", pascal.temp$cellType)
    
pascal.temp$cellType <- ifelse(is.na(pascal.temp$cellType), "neither", pascal.temp$cellType)

#pdf("graphics/MAGMA-vs-PASCAL-gene-levelResults-from-Liu-etal2019_AMY-DrnkWk-highlight_MNT2021.pdf")
# set.seed(901)
# print(
#   ggplot(data=pascal.temp, aes(x=z.pascal, y=z.magma.gene)) +
#     geom_point(aes(color=cellType), size=ifelse(pascal.temp$cellType=="nonspecific",2.5,5.0),
#                alpha=ifelse(pascal.temp$cellType=="neither",0.1,0.8)) +
#     geom_text(label=ifelse(pascal.temp$cellType %in% c("Inhib_G","Inhib_E","shared"),
#                            pascal.temp$Gene.Symbol,NA), hjust=0, check_overlap=F,
#               fontface="italic", size=3, position=position_jitter(width=0.03,height=0.07)) +
#     scale_color_manual(values=c(cell_colors.amy, shared="#669966", neither="#333333")) +
#     geom_smooth(method = "lm", se=FALSE) +
#     labs(title = paste0(i," (Liu, et al. 2019): gene-level scores"),
#          subtitle = "Enriched markers for FDR-significant AMY classes Inhib_E/_G",
#          x = "PASCAL gene-level z-scores",
#          y = "Corresponding MAGMA z-scores",
#          colour="Enriched in:") +
#     theme(plot.title = element_text(size=16, face="bold"),
#           axis.title=element_text(size=15,face="bold"),
#           axis.text = element_text(size=12))
# )
#dev.off()
    # Nvm don't like this


## For 4D ===
# Print an iteration for Fig 4C, highlighting the Bonf-significant classes' PW markers
i <- "DrnkWk"
pascal.temp <- liu.pascal.list[[i]]
pascal.temp$cellType <- markers.amy.specific$cellType[match(pascal.temp$ensemblID, markers.amy.specific$geneID)]
pascal.temp$cellType <- ifelse(is.na(pascal.temp$cellType), "nonspecific", pascal.temp$cellType)

pdf("graphics/MAGMA-vs-PASCAL-gene-levelResults-from-Liu-etal2019_AMY-DrnkWk-highlight_MNT2021.pdf")
set.seed(720)
print(
  ggplot(data=pascal.temp, aes(x=z.pascal, y=z.magma.gene)) +
    geom_point(aes(color=cellType), size=ifelse(pascal.temp$cellType=="nonspecific",2.5,5.0),
               alpha=ifelse(pascal.temp$cellType=="nonspecific",0.1,0.8)) +
    geom_text(label=ifelse(!pascal.temp$cellType=="nonspecific",pascal.temp$Gene.Symbol,NA), hjust=0, check_overlap=F,
              fontface="italic", size=2.7, position=position_jitter(width=0.07,height=0.07)) +
    scale_color_manual(values=c(cell_colors.amy, nonspecific="#333333")) +
    geom_smooth(method = "lm", se=FALSE) + xlim(4.4,7.1) +
    labs(title = paste0(i," (Liu, et al. 2019): gene-level scores"),
         subtitle = "Pairwise-significant markers in AMY cell classes",
         x = "PASCAL gene-level z-scores",
         y = "Corresponding MAGMA z-scores",
         colour="Cell type") +
    theme(plot.title = element_text(size=16, face="bold"),
          axis.title=element_text(size=15,face="bold"),
          axis.text = element_text(size=12))
)
dev.off()



### Session info for 20Jul2021 ====================================
sessionInfo()




