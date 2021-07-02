### MNT 10x snRNA-seq workflow: (step 04?:)
### Miscellaneous lookings-into / finalizing graphics for manuscript
###   - Brief neuron-specific clustering for DLPFC (Maynard-Collado-Torres et al.)
###   - 10x pilot snRNA-seq paper (Tran-Maynard et al.)
##################################################################################

library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(batchelor)
library(DropletUtils)
library(jaffelab)


### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")



#### Does sex impact cluster proportions? ================
  ## --> Logistic regression: Sex ~ cell type proportion [for each cell type]

### Velmeshev, et al (PFC & ACC) ===
load("/dcl01/ajaffe/data/lab/singleCell/velmeshev2019/analysis_MNT/SCE_asd-velmeshev-etal_MNT.rda",
     verbose=T)
    # *Since these regions were analyzed together, just keep them together

table(sce.asd$individual, sce.asd$sex)
    # 31 individuals, 7 female
table(sce.asd$individual, sce.asd$BAregion)
    # Basically individuals could have been sampled multiple times - BA24[ACC] and/or a PFC BAregion

apply(table(sce.asd$cluster, sce.asd$sex), 2, function(x){round(prop.table(x),3)})
    #                      F     M
    # AST-FB           0.042 0.044
    # AST-PP           0.098 0.067
    # Endothelial      0.019 0.027
    # IN-PV            0.028 0.041
    # IN-SST           0.035 0.046
    # IN-SV2C          0.018 0.018
    # IN-VIP           0.056 0.058
    # L2/3             0.133 0.126
    # L4               0.069 0.065
    # L5/6             0.032 0.035
    # L5/6-CC          0.043 0.044
    # Microglia        0.034 0.031
    # Neu-mat          0.037 0.040
    # Neu-NRGN-I       0.043 0.032
    # Neu-NRGN-II      0.082 0.079
    # Oligodendrocytes 0.130 0.151
    # OPC              0.101 0.096    - broadly it looks pretty much the same

# Make cell proportion table per individual
propDat <- apply(table(sce.asd$cluster, sce.asd$individual), 2, function(x){signif(prop.table(x),6)})

# Matching sex, binarized
ind.idx <- colnames(propDat)
sex.idx <- sce.asd$sex[match(ind.idx, sce.asd$individual)]
sex.idx <- as.numeric(sex.idx=="F")

regrStats.sex <- list()
for(i in levels(sce.asd$cluster)){
  regrStats.sex[[i]] <- glm(sex.idx ~ propDat[i, ], family="binomial")
}

# P-values for regression on cell type prop.
sapply(regrStats.sex, function(x){coef(summary(x))["propDat[i, ]", "Pr(>|z|)"]})
    #    AST-FB           AST-PP      Endothelial            IN-PV 
    # 0.8598569        0.6190065        0.2571417        0.1365944 
    #    IN-SST          IN-SV2C           IN-VIP             L2/3 
    # 0.2652023        0.7407465        0.8813374        0.5832419 
    #        L4             L5/6          L5/6-CC        Microglia 
    # 0.7443358        0.7411243        0.9355193        0.9933819 
    #   Neu-mat       Neu-NRGN-I      Neu-NRGN-II Oligodendrocytes 
    # 0.6790237        0.9039778        0.7741251        0.7656764 
    #       OPC 
    # 0.7381165 


# (Btw marker statistics computed on:
    # mod <- with(colData(sce.asd.pfc), model.matrix(~ diagnosis + age + sex + Capbatch + RIN))
    # mod <- mod[ ,-1])
    #   --> don't think I want to include these additional predictors tho





### And for Mathys, et al ============
load("rdas/referenceDatasets/SCE_mathys-PFC-BA10_MNT.rda", verbose=T)
    # sce.mathys

colnames(colData(sce.mathys))
    # [1] "projid"          "pre.cluster"     "broad.cell.type" "Subcluster"     
    # [5] "individualID"    "Dx"              "age_death"       "msex"           
    # [9] "race"            "sizeFactor"  

table(sce.mathys$Subcluster, sce.mathys$msex)
    # Pretty even distribution, by eye

# Do the same thing as above
propDat <- apply(table(sce.mathys$Subcluster, sce.mathys$individualID), 2, function(x){signif(prop.table(x),6)})

# Matching sex, binarized
ind.idx <- colnames(propDat)
sex.idx <- sce.mathys$msex[match(ind.idx, sce.mathys$individualID)]
sex.idx <- as.numeric(sex.idx=="0")

regrStats.sex <- list()
for(i in sort(unique(sce.mathys$Subcluster))){
  regrStats.sex[[i]] <- glm(sex.idx ~ propDat[i, ], family="binomial")
}

# P-values for regression on cell type prop.
sapply(regrStats.sex, function(x){signif(coef(summary(x)),6)["propDat[i, ]", "Pr(>|z|)"]})
    #      Ast0      Ast1      Ast2      Ast3      End1      End2       Ex0       Ex1 
    # 0.7961580 0.8174390 0.7589670 0.9583130 0.5789900 0.4910990 0.6330960 0.0545296 
    #      Ex11      Ex12      Ex14       Ex2       Ex3       Ex4       Ex5       Ex6 
    # 0.4111320 0.6637870 0.2281180 0.7249140 0.1487160 0.0754658 0.1397110 0.4928440 
    #       Ex7       Ex8       Ex9       In0       In1      In10      In11       In2 
    # 0.3383150 0.1877990 0.8266440 0.3886890 0.0700323 0.0788781 0.4623090 0.1488410 
    #       In3       In4       In5       In6       In7       In8       In9      Mic0 
    # 0.2427360 0.2654090 0.1671240 0.7676490 0.2821590 0.6259190 0.9835870 0.7097980 
    #      Mic1      Mic2      Mic3      Oli0      Oli1      Oli3      Oli4      Oli5 
    # 0.9702780 0.5932430 0.2823240 0.3986670 0.1447250 0.1929830 0.5829920 0.9261180 
    #      Opc0      Opc1      Opc2       Per 
    # 0.8286980 0.1100380 0.5692740 0.7238730 
ps.logisticReg <- sapply(regrStats.sex, function(x){signif(coef(summary(x)),6)["propDat[i, ]", "Pr(>|z|)"]})


## By simpler `lm(prop ~ sex)` ===
regrStats.sex.lm <- list()
for(i in sort(unique(sce.mathys$Subcluster))){
  regrStats.sex.lm[[i]] <- lm(propDat[i, ] ~ sex.idx)
}
# P-values for regression on cell type prop.
sapply(regrStats.sex.lm, function(x){signif(coef(summary(x)),6)["sex.idx", "Pr(>|t|)"]})
ps.lm <- sapply(regrStats.sex.lm, function(x){signif(coef(summary(x)),6)["sex.idx", "Pr(>|t|)"]})
    # All > 0.05 except 'Ex1': 0.0484782

cor(ps.lm, ps.logisticReg)
    # [1] 0.999268

## What if add 'age_death'? ===
ageDeath.idx <- sce.mathys$age_death[match(ind.idx, sce.mathys$individualID)]

regrStats.sex.lm <- list()
for(i in sort(unique(sce.mathys$Subcluster))){
  regrStats.sex.lm[[i]] <- lm(propDat[i, ] ~ sex.idx + ageDeath.idx)
}
# P-values for sex term, with 'age_death' included
sapply(regrStats.sex.lm, function(x){signif(coef(summary(x)),6)["sex.idx", "Pr(>|t|)"]})
    #      Ast0      Ast1      Ast2      Ast3      End1      End2       Ex0       Ex1 
    # 0.8945570 0.7448440 0.7069090 0.8059920 0.4090890 0.5901650 0.5657990 0.0831644 
    #      Ex11      Ex12      Ex14       Ex2       Ex3       Ex4       Ex5       Ex6 
    # 0.5430500 0.5407110 0.3459960 0.5924880 0.2069800 0.0685688 0.2072260 0.5739530 
    #       Ex7       Ex8       Ex9       In0       In1      In10      In11       In2 
    # 0.5152840 0.1394130 0.7542460 0.3782150 0.0766577 0.1220210 0.6955990 0.1371930 
    #       In3       In4       In5       In6       In7       In8       In9      Mic0 
    # 0.1957020 0.3230890 0.0927227 0.9837080 0.3858550 0.6341780 0.9307060 0.7890800 
    #      Mic1      Mic2      Mic3      Oli0      Oli1      Oli3      Oli4      Oli5 
    # 0.9976360 0.5730220 0.2560560 0.4411790 0.1849180 0.1724730 0.6886770 0.8794820 
    #      Opc0      Opc1      Opc2       Per 
    # 0.9913400 0.1193400 0.5304300 0.6870260 





## Median gene capture per cell type? ===

# NAc
load("rdas/regionSpecific_NAc-ALL-n5_cleaned-combined_SCE_MNTMar2020.rda", verbose=T)
    # sce.nac.all, chosen.hvgs.nac.all, pc.choice.nac.all, clusterRefTab.nac.all, ref.sampleInfo

sce.nac.all <- sce.nac.all[ ,-which(sce.nac.all$cellType.final=="ambig.lowNtrxts")]
sce.nac.all$cellType.final <- droplevels(sce.nac.all$cellType.final)

cell.idx <- splitit(sce.nac.all$cellType.final)
sapply(cell.idx, function(x){quantile(sce.nac.all$detected[x])})
    #      Astro Inhib.1 Inhib.2 Inhib.3 Inhib.4 Micro MSN.D1.1 MSN.D1.2 MSN.D1.3
    # 0%     473    2782 2917.00 1980.00  1047.0   836  1021.00     1508   1632.0
    # 25%   2729    5205 4066.75 4872.50  4280.5  1706  4513.25     4385   4919.0
    # 50%   3147    5974 4723.50 5480.00  4896.0  2190  5005.50     4933   5363.0
    # 75%   3587    6823 6219.00 6226.25  5445.0  2720  5577.00     5342   5783.5
    # 100%  6526    9332 8298.00 9561.00  9021.0  3792  7846.00     7743   8378.0
    #      MSN.D1.4 MSN.D2.1 MSN.D2.2  Oligo    OPC
    # 0%       1372   800.00   1240.0  303.0  910.0
    # 25%      5480  5104.75   5258.0 1689.5 2884.0
    # 50%      5972  5603.50   5724.0 2216.0 3300.0
    # 75%      6468  6188.25   6186.5 2767.5 3708.5
    # 100%     9874  9108.00   9917.0 6129.0 6325.0

apply(table(sce.nac.all$cellType.final, sce.nac.all$donor), 2, function(x){round(prop.table(x),2)})
    #          Br5161 Br5182 Br5207 Br5212 Br5287
    # Astro      0.07   0.00   0.00   0.22   0.02
    # Inhib.1    0.00   0.00   0.00   0.00   0.00
    # Inhib.2    0.00   0.01   0.00   0.00   0.00
    # Inhib.3    0.00   0.02   0.04   0.00   0.01
    # Inhib.4    0.00   0.02   0.01   0.00   0.01
    # Micro      0.04   0.00   0.00   0.04   0.05
    # MSN.D1.1   0.00   0.03   0.00   0.00   0.00
    # MSN.D1.2   0.00   0.07   0.00   0.00   0.00
    # MSN.D1.3   0.01   0.09   0.07   0.00   0.01
    # MSN.D1.4   0.09   0.35   0.41   0.10   0.11
    # MSN.D2.1   0.00   0.03   0.03   0.00   0.00
    # MSN.D2.2   0.02   0.38   0.42   0.07   0.01
    # Oligo      0.71   0.00   0.00   0.49   0.73
    # OPC        0.05   0.00   0.00   0.06   0.05


# AMY
load("rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose=T)
    # sce.amy, chosen.hvgs.amy, pc.choice.amy, clusterRefTab.amy, ref.sampleInfo

sce.amy <- sce.amy[ ,-which(sce.amy$cellType.split=="Ambig.lowNtrxts")]
sce.amy$cellType.split <- droplevels(sce.amy$cellType.split)

cell.idx <- splitit(sce.amy$cellType.split)
sapply(cell.idx, function(x){quantile(sce.amy$detected[x])})
    #        Astro  Excit.1 Excit.2 Excit.3 Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5
    # 0%    460.00   634.00 1990.00    2722  1197.0    2566  3025.0 3233.00     629
    # 25%  2926.25  7426.00 3257.75    4506  6623.5    6458  6809.0 5208.00    6480
    # 50%  3567.00  8389.00 3672.50    4876  7420.0    7277  7378.0 6098.00    7342
    # 75%  4383.50  9096.25 3992.00    6054  8177.5    7994  7946.5 7019.25    7853
    # 100% 6609.00 11625.00 5072.00    8118 10176.0    9158  8734.0 7936.00    9368
    #      Micro Oligo    OPC
    # 0%     217   600 1223.0
    # 25%   1800  2044 3382.5
    # 50%   2370  2555 3921.0
    # 75%   2749  3138 4422.5
    # 100%  4271  5607 6378.0

# Broadly
quantile(sce.amy$detected)
    #     0%      25%      50%      75%     100% 
    # 217.00  2225.00  2934.50  3911.75 11625.00 




## AMY vs. mouse MeA ===
load("/dcl01/ajaffe/data/lab/singleCell/ucla_mouse-MeA/2019Cell/SCE_mouse-MeA_downstream-processing_MNT.rda", verbose=T)
    # sce.amy.mm, chosen.hvgs.amy.mm

table(sce.amy.mm$subCluster)

## Dist of genes expressed?
sce.amy.mm$detected <- apply(assay(sce.amy.mm,"logcounts"), 2, function(x){table(x > 0)["TRUE"]})

subClust.idx <- splitit(sce.amy.mm$subCluster)
sapply(subClust.idx, function(x){quantile(sce.amy.mm$detected[x])[c("25%", "75%")]})
    #       AS     EN   MG      MU    N.1 N.2    N.3    N.4     N.5    N.6    N.7
    # 25%  273  353.5  312  412.25 241.00 249 285.00 264.00  505.25  265.5  834.5
    # 75% 1043 1106.5 1123 1210.75 813.25 705 888.75 765.75 1450.00 1148.0 2178.5
    #
    #         N.8    N.9 N.10    N.11 N.12 N.13  N.14 N.15 N.16   OL     OPC OPC.OL
    # 25%  365.75 320.50  226  379.00  514  228 245.5  230  253  696  639.00  454.0
    # 75% 1245.50 961.75  479 1802.25 2001  374 452.0  354  474 1864 2181.75 2708.5

# Broadly
quantile(sce.amy.mm$detected)
    #  0%  25%  50%  75% 100% 
    # 177  299  612 1299 8547 


    # For some reason, these were actually somewhat ordered by n Genes detected...
        head(sce.amy.mm$detected, n=40)
        tail(sce.amy.mm$detected, n=40)
        
        apply(assay(sce.amy.mm,"logcounts")[ ,1:20], 2, function(x){table(x > 0)["TRUE"]})
        apply(assay(sce.amy.mm,"logcounts")[ ,43316:43345], 2, function(x){table(x > 0)["TRUE"]})

    # oh they're not.  They're ordered by $sample.  So they're quite batch effect-y.....
        table(sce.amy.mm$sample)
        sce.amy.mm$sample[1:40]
        sce.amy.mm$sample[3258:3269]
        3261+ 3540 +2188 +4913+ 4149+ 2960  # n males
            #[1] 21011
        sce.amy.mm$sample[21000:21040]
            # [1] "M6" "M6" "M6" "M6" "M6" "M6" "M6" "M6" "M6" "M6" "M6" "M6" "F1" "F1" "F1"
            # [16] "F1" "F1" "F1" "F1" "F1" "F1" "F1" "F1" "F1" "F1" "F1" "F1" "F1" "F1" "F1"
            # [31] "F1" "F1" "F1" "F1" "F1" "F1" "F1" "F1" "F1" "F1" "F1"
    
    # This is strange:
        sample.idx <- splitit(sce.amy.mm$sample)
        sapply(sample.idx, function(x){quantile(sce.amy.mm$detected[x])})
        
    # AHhhhhh:
    sce.amy.mm$detected[21000:21040]
        
        #   ** so these data are ordered by $sample, and then ~ in descending order of detected # genes
        #      (or maybe capture).  Curious.
    
    
    

## Concordance b/tw region-specific and pan-brain broad annotation ===
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/all-FACS-homogenates_n12_PAN-BRAIN-Analyses_MNTFeb2020.rda",
     verbose=T)
    # sce.all.n12, chosen.hvgs.all.n12, pc.choice.n12, ref.sampleInfo, clusterRefTab.all.n12

# Getting rid of sub-cell types in sACC sample and the pan-brain assignments to compare
table(gsub("MSN","Inhib",ss(as.character(sce.all.n12$cellType.RS),"\\.",1)) ==
        ss(as.character(sce.all.n12$cellType),"\\.",1))
    # FALSE  TRUE
    #   866 33204    - 33204/(866+ 33204) = 97.5% concordant

# Region specific
table(ss(as.character(sce.all.n12$cellType.RS),"\\.",1))
    # Ambig Astro Excit Inhib Micro   MSN Oligo   OPC Tcell
    #   445  3864  2927  2019  2956   642 18664  2527    26

    table( ss(as.character(sce.all.n12$cellType),"\\.",1))
    #Ambig Astro Excit Inhib Micro Oligo   OPC
    #   32  3828  2848  3110  3077 18614  2561

# What's the broad assignment of the .RS 'Ambig's??
table(ss(as.character(sce.all.n12$cellType.RS),"\\.",1))
ambig.idx <- which(ss(as.character(sce.all.n12$cellType.RS),"\\.",1)=="Ambig")
table(sce.all.n12$cellType[ambig.idx])
    # Ambig.hiVCAN      Astro.1      Astro.2      Excit.1      Excit.2      Excit.3 
    #            0            0            2            0            0            0 
    #      Excit.4      Excit.5      Excit.6      Excit.7      Excit.8      Inhib.1 
    #            0            0            0            0            0          368     - that's interesting lol
    #      Inhib.2      Inhib.3      Inhib.4      Inhib.5        Micro        Oligo       (out of 950)
    #            2            0            0            0           66            1 
    #          OPC 
    #            6

panClust.idx <- splitit(sce.all.n12$cellType)
sapply(panClust.idx, function(x){quantile(sce.all.n12$sum[x])})
    # Inhib.1 median/IQR is lower than all the other neuronal subclusters...
    #     start to wonder if a lot of this cluster is driven by this technical artifact

# Let's remove those
sce.all.sub <- sce.all.n12[ ,-ambig.idx]
table(gsub("MSN","Inhib",ss(as.character(sce.all.sub$cellType.RS),"\\.",1)) ==
        ss(as.character(sce.all.sub$cellType),"\\.",1))
    # FALSE  TRUE 
    #   421 33204 
33204/(33204+421) # 98.7%   (and if you remove that weird 32 'Ambig.hiVCAN' -> 98.8% concordant)



