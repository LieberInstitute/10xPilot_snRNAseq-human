### for MAGMA with LIBD 10x pilot analyses
  #     - plotting Results_rev/heatmaps
  # UPDATE: re-running with new v1.08
  # MNT 14Jul2021: plotting all stats from GSA tests with 102 revision
  #                cell class markers ================================

library(readr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(fields)
library(jaffelab)

# ===

regions <- c("dlpfc","sacc","hpc","nac","amy")
magmaStats <- list()

for(i in regions){
  magmaStats[[i]][["SCZ.PGC2"]] <- read.table(paste0("./Results_rev/",i,"_clozuk_pgc2.gsa.out"), header=T)
  #magmaStats[[i]][["PGC3.SCZ"]] <- read.table(paste0("./Results_rev/",i,"_pgc3_scz.gsa.out"), header=T)
  magmaStats[[i]][["ASD.PGC"]] <- read.table(paste0("./Results_rev/",i,"_PGC_ASD.gsa.out"), header=T)
  magmaStats[[i]][["BIP.PGC"]] <- read.table(paste0("./Results_rev/",i,"_PGC_BIP.gsa.out"), header=T)
  magmaStats[[i]][["MDD.PGC"]] <- read.table(paste0("./Results_rev/",i,"_MDD29_23andMe.gsa.out"), header=T)
  magmaStats[[i]][["PTSD.PGC2"]] <- read.table(paste0("./Results_rev/",i,"_PTSD.gsa.out"), header=T)
  magmaStats[[i]][["ADHD.PGC.iPsych"]] <- read.table(paste0("./Results_rev/",i,"_ADHD.gsa.out"), header=T)
  magmaStats[[i]][["AD.metaPh3"]] <- read.table(paste0("./Results_rev/",i,"_AD.gsa.out"), header=T)
  
  # Addiction GWAS
  magmaStats[[i]][["addxn.AgeSmk"]] <- read.table(paste0("./Results_rev/",i,"_AgeSmk.gsa.out"), header=T)
  magmaStats[[i]][["addxn.CigDay"]] <- read.table(paste0("./Results_rev/",i,"_CigDay.gsa.out"), header=T)
  magmaStats[[i]][["addxn.DrnkWk"]] <- read.table(paste0("./Results_rev/",i,"_DrnkWk.gsa.out"), header=T)
  magmaStats[[i]][["addxn.SmkInit"]] <- read.table(paste0("./Results_rev/",i,"_SmkInit.gsa.out"), header=T)
  magmaStats[[i]][["addxn.SmkCes"]] <- read.table(paste0("./Results_rev/",i,"_SmkCes.gsa.out"), header=T)
  # Added for control but not including in multiple test corrxn, and just for supplement:
  #   - CARDIoGRAM+C4D coronary artery dx meta-GWAS
  #magmaStats[[i]][["CoronArtDx"]] <- read.table(paste0("./Results_rev/",i,"_CAD.gsa.out"), header=T)
}

## merge to assess significance thresholds (did this before adding CAD meta-GWAS) ===
magmaStats_list = lapply(magmaStats, function(m) {
	z = sapply(m, "[[", "P")
	rownames(z) = m[[1]]$VARIABLE
	z
})
magmaStats_wide = as.data.frame(do.call("rbind", magmaStats_list))
magmaStats_wide$Region = rep(names(magmaStats_list), 
	times = sapply(magmaStats_list, nrow))
magmaStats_wide$CellType = rownames(magmaStats_wide)
## reshape to long
magmaStats_long = reshape2::melt(magmaStats_wide)
colnames(magmaStats_long)[3:4] = c("GWAS", "P")
dim(magmaStats_long)
    # [1] 1284    4     - 1248 bc 107 reported cell classes x 12 GWAS'

table(p.adjust(magmaStats_long$P, "fdr") < 0.05)
    #FALSE  TRUE 
    #  892   392
betacut.fdr <- max(magmaStats_long$P[p.adjust(magmaStats_long$P, "fdr") < 0.05])
    # [1] 0.015109
table(p.adjust(magmaStats_long$P, "bonf") < 0.05)
    # FALSE  TRUE 
    #  1168   116
betacut.bonf <- max(magmaStats_long$P[p.adjust(magmaStats_long$P, "bonf") < 0.05])
    #[1] 3.7819e-05

magmaStats_long$P.adj.fdr <- p.adjust(magmaStats_long$P, "fdr")
          
# (for interactive exploration:)
magmaStats_long[which(magmaStats_long$P.adj.fdr < 0.05 &
                        magmaStats_long$GWAS=="AD.metaPh3"), ]


    
    # For betas (to include in supplementary table) ===
    magmaStats_list.beta = lapply(magmaStats, function(m) {
      z = sapply(m, "[[", "BETA")
      rownames(z) = m[[1]]$VARIABLE
      z
    })
    magmaStats_wide.beta = as.data.frame(do.call("rbind", magmaStats_list.beta))
    magmaStats_wide.beta$Region = rep(names(magmaStats_list.beta), 
                                      times = sapply(magmaStats_list.beta, nrow))
    magmaStats_wide.beta$CellType = rownames(magmaStats_wide.beta)
    ## reshape to long
    magmaStats_long.beta = reshape2::melt(magmaStats_wide.beta)
    colnames(magmaStats_long.beta)[3:4] = c("GWAS", "Beta")

    #Check
    table(magmaStats_long$CellType == magmaStats_long.beta$CellType) # (& $GWAS) - all TRUE

# cbind Beta
magmaStats_long$Beta <- magmaStats_long.beta$Beta
# Reorder
magmaStats_long <- magmaStats_long[ ,c("Region", "CellType", "GWAS", "Beta", "P", "P.adj.fdr")]

## Supplementary table: All GWAS included in main section (12 phenotypes, no PGC3-SCZD) ===
write.csv(magmaStats_long, file = "../tables/revision/suppTable_magma-v1.08_byCluster_12xGWAS_MNT2021.csv",
          row.names=F)

# CAD results, separately (without computing 'P.adj.fdr') because not part of main section
#   -> first include all GWAS in 'magmaStats' & remake 'magmaStats_long'/'magmaStats_long.beta'
magmaStats.sup2 <- magmaStats_long[magmaStats_long$GWAS %in% c("CoronArtDx"), ]
write.csv(magmaStats.sup2, file = "../tables/revision/suppTable_magma-v1.08_byCluster_CAD-GWAS_MNT2021.csv",
          row.names=F)



## SKIP BELOW (v2 heatmap) Adapted from Spatial_dlpfc paper - custom heatmap =====

midpoint = function(x) x[-length(x)] + diff(x)/2

customMAGMAplot = function(region, Pthresh, BETAcut, ...) {
      ## Set up -log10(p's)
      wide_p = sapply(magmaStats[[region]], function(x){cbind(-log10(x$P))})
      rownames(wide_p) <- magmaStats[[region]][[1]]$VARIABLE
      wide_p[wide_p > Pthresh] = Pthresh
      wide_p <- round(wide_p[rev(sort(rownames(wide_p))), ], 3)

      
      ## Set up betas
      wide_beta <- sapply(magmaStats[[region]], function(x){cbind(x$BETA)})
      rownames(wide_beta) <- magmaStats[[region]][[1]]$VARIABLE
      wide_beta <- round(wide_beta[rev(sort(rownames(wide_beta))), ], 2)
      
      # Use FDR/Bonf cutoffs for printing betas, set above
      wide_beta[wide_p < -log10(BETAcut)] = ""
      
      
      ## Plot
      clusterHeights <- seq(0,160,length.out=nrow(wide_p)+1)
      mypal = c("white", colorRampPalette(brewer.pal(9,"YlOrRd"))(60))[1:50]
      xlabs <- colnames(wide_p)
      
      # Heatmap of p's
      image.plot(x = seq(0,ncol(wide_p),by=1), y = clusterHeights, z = as.matrix(t(wide_p)),
                 col = mypal,xaxt="n", yaxt="n",xlab = "", ylab="", ...)
      axis(2, rownames(wide_p), at=midpoint(clusterHeights), las=1)
      axis(1, rep("", ncol(wide_p)), at = seq(0.5,ncol(wide_p)-0.5))
      text(x = seq(0.5,ncol(wide_p)-0.5), y=-1*max(nchar(xlabs))/2, xlabs,
           xpd=TRUE, srt=45,cex=1.25,adj= 1)
      abline(h=clusterHeights,v=0:ncol(wide_p))
      
      # Print top decile of betas
      text(x = rep(seq(0.5,ncol(wide_p)-0.5),each = nrow(wide_p)), 
           y = rep(midpoint(clusterHeights), ncol(wide_p)),
           as.character(wide_beta),cex=1.0,font=2)
    }

## Plot
#dir.create("./graphics/")

# DLPFC
pdf("graphics/heatmap_dlpfc-magma-v.08_PGC-and-addiction_GWAS_MNT09Sep2020.pdf", w=8)
par(mar=c(8,7.5,6,1), cex.axis=1.0, cex.lab=0.5)
customMAGMAplot(region="dlpfc", Pthresh=12, BETAcut=betacut.fdr)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.0, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: DLPFC subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests at FDR < 0.05)", xpd=TRUE, cex=1, font=1)

customMAGMAplot(region="dlpfc", Pthresh=12, BETAcut=betacut.bonf)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.0, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: DLPFC subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests with Bonf. < 0.05)", xpd=TRUE, cex=1, font=1)
dev.off()


# sACC
pdf("graphics/heatmap_sacc-magma-v.08_PGC-and-addiction_GWAS_MNT09Sep2020.pdf", w=8)
par(mar=c(8,7.5,6,1), cex.axis=1.3, cex.lab=0.5)
customMAGMAplot(region="sacc", Pthresh=12, BETAcut=betacut.fdr)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.0, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: sACC subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests at FDR < 0.05)", xpd=TRUE, cex=1, font=1)

customMAGMAplot(region="sacc", Pthresh=12, BETAcut=betacut.bonf)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.0, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: sACC subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests with Bonf. < 0.05)", xpd=TRUE, cex=1, font=1)
dev.off()


# HPC
pdf("graphics/heatmap_hpc-magma-v.08_PGC-and-addiction_GWAS_MNT09Sep2020.pdf", w=8)
par(mar=c(8,7.5,5,1), cex.axis=1.1, cex.lab=0.5)
customMAGMAplot(region="hpc", Pthresh=12, BETAcut=betacut.fdr)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.0, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: HIPPO subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests at FDR < 0.05)", xpd=TRUE, cex=1, font=1)

customMAGMAplot(region="hpc", Pthresh=12, BETAcut=betacut.bonf)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.0, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: HIPPO subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests with Bonf. < 0.05)", xpd=TRUE, cex=1, font=1)
dev.off()


# NAc
pdf("graphics/heatmap_nac-magma-v.08_PGC-and-addiction_GWAS_MNT09Sep2020.pdf", w=8)
par(mar=c(8,7.5,6,1), cex.axis=1.0, cex.lab=0.5)
customMAGMAplot(region="nac", Pthresh=12, BETAcut=betacut.fdr)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.0, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: NAc subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests at FDR < 0.05)", xpd=TRUE, cex=1, font=1)

customMAGMAplot(region="nac", Pthresh=12, BETAcut=betacut.bonf)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.0, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: NAc subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests with Bonf. < 0.05)", xpd=TRUE, cex=1, font=1)
dev.off()


# Amy
pdf("graphics/heatmap_amy-magma-v.08_PGC-and-addiction_GWAS_MNT09Sep2020.pdf", w=8)
par(mar=c(8,7.5,5,1), cex.axis=1.2, cex.lab=0.5)
customMAGMAplot(region="amy", Pthresh=12, BETAcut=betacut.fdr)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.0, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: AMY subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests at FDR < 0.05)", xpd=TRUE, cex=1, font=1)

customMAGMAplot(region="amy", Pthresh=12, BETAcut=betacut.bonf)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.0, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: AMY subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests with Bonf. < 0.05)", xpd=TRUE, cex=1, font=1)
dev.off()




## Version 2: printing all FDR-signif, bolding where Bonf.-signif ==============

midpoint = function(x) x[-length(x)] + diff(x)/2

customMAGMAplot.b = function(region, Pthresh, fdrThresh, bonfThresh, ...) {
  ## Set up -log10(p's)
  wide_p = sapply(magmaStats[[region]], function(x){cbind(-log10(x$P))})
  rownames(wide_p) <- magmaStats[[region]][[1]]$VARIABLE
  wide_p[wide_p > Pthresh] = Pthresh
  wide_p <- round(wide_p[rev(sort(rownames(wide_p))), ], 3)
  
  
  ## Set up betas
  wide_beta <- sapply(magmaStats[[region]], function(x){cbind(x$BETA)})
  rownames(wide_beta) <- magmaStats[[region]][[1]]$VARIABLE
  wide_beta <- round(wide_beta[rev(sort(rownames(wide_beta))), ], 2)
  
  # Use FDR cutoff for printing betas
  wide_beta[wide_p < -log10(fdrThresh)] = ""
  # and Bonf. cutoff for bolding
  customFont <- ifelse(wide_p < -log10(bonfThresh), 1, 2)
  customCex <- ifelse(wide_p < -log10(bonfThresh), 0.9, 1.0)
  
  
  ## Plot
  clusterHeights <- seq(0,160,length.out=nrow(wide_p)+1)
  mypal = c("white", colorRampPalette(brewer.pal(9,"YlOrRd"))(60))[1:50]
  xlabs <- colnames(wide_p)
  
  # Heatmap of p's
  image.plot(x = seq(0,ncol(wide_p),by=1), y = clusterHeights, z = as.matrix(t(wide_p)),
             col = mypal,xaxt="n", yaxt="n",xlab = "", ylab="", ...)
  axis(2, rownames(wide_p), at=midpoint(clusterHeights), las=1)
  axis(1, rep("", ncol(wide_p)), at = seq(0.5,ncol(wide_p)-0.5))
  text(x = seq(0.5,ncol(wide_p)-0.5), y=-1*max(nchar(xlabs))/2, xlabs,
       xpd=TRUE, srt=45, cex=1.2, adj= 1)
  abline(h=clusterHeights,v=0:ncol(wide_p))
  
  # Print top decile of betas
  text(x = rep(seq(0.5,ncol(wide_p)-0.5),each = nrow(wide_p)), 
       y = rep(midpoint(clusterHeights), ncol(wide_p)),
       as.character(wide_beta),
       # If Bonf, a little bigger
       cex=customCex,
       # If Bonf, 2 (bold)
       font=customFont)
}

## Now just need one plot per region :)
# DLPFC
pdf("graphics/heatmap-rev_dlpfc-magma-v.08_PGC-and-addiction_GWAS_MNT2021.pdf", w=8)
par(mar=c(8,7.5,6,1), cex.axis=1.0, cex.lab=0.5)
customMAGMAplot.b(region="dlpfc", Pthresh=12, fdrThresh=betacut.fdr, bonfThresh=betacut.bonf)
abline(v=7,lwd=3)
text(x = c(3.5,9.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=0.9, font=2)
text(x = 6, y=185, "MAGMA gene set analyses: DLPFC cell classes", xpd=TRUE, cex=1.5, font=2)
text(x = 6, y=175, "(Betas controlling all tests at FDR < 0.05)", xpd=TRUE, cex=1, font=1)
grid::grid.text(label="-log10(p-value)", x=0.93, y=0.825, gp=gpar(fontsize=9))
dev.off()

# sACC
pdf("graphics/heatmap-rev_sacc-magma-v.08_PGC-and-addiction_GWAS_MNT2021.pdf", w=8, h=8)
par(mar=c(8,7.5,6,1), cex.axis=0.9, cex.lab=0.5)
customMAGMAplot.b(region="sacc", Pthresh=12, fdrThresh=betacut.fdr, bonfThresh=betacut.bonf)
abline(v=7,lwd=3)
text(x = c(3.5,9.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=0.9, font=2)
text(x = 6, y=185, "MAGMA gene set analyses: sACC cell classes", xpd=TRUE, cex=1.5, font=2)
text(x = 6, y=175, "(Betas controlling all tests at FDR < 0.05)", xpd=TRUE, cex=1, font=1)
grid::grid.text(label="-log10(p-value)", x=0.93, y=0.835, gp=gpar(fontsize=9))
dev.off()

# HIPPO
pdf("graphics/heatmap-rev_hpc-magma-v.08_PGC-and-addiction_GWAS_MNT2021.pdf", w=8)
par(mar=c(8,7.5,6,1), cex.axis=1.0, cex.lab=0.5)
customMAGMAplot.b(region="hpc", Pthresh=12, fdrThresh=betacut.fdr, bonfThresh=betacut.bonf)
abline(v=7,lwd=3)
text(x = c(3.5,9.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=0.9, font=2)
text(x = 6, y=185, "MAGMA gene set analyses: HPC cell classes", xpd=TRUE, cex=1.5, font=2)
text(x = 6, y=175, "(Betas controlling all tests at FDR < 0.05)", xpd=TRUE, cex=1, font=1)
grid::grid.text(label="-log10(p-value)", x=0.93, y=0.825, gp=gpar(fontsize=9))
dev.off()

# NAc
pdf("graphics/heatmap-rev_nac-magma-v.08_PGC-and-addiction_GWAS_MNT2021.pdf", w=8)
par(mar=c(8,7.5,6,1), cex.axis=1.0, cex.lab=0.5)
customMAGMAplot.b(region="nac", Pthresh=12, fdrThresh=betacut.fdr, bonfThresh=betacut.bonf)
abline(v=7,lwd=3)
text(x = c(3.5,9.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=0.9, font=2)
text(x = 6, y=185, "MAGMA gene set analyses: NAc cell classes", xpd=TRUE, cex=1.5, font=2)
text(x = 6, y=175, "(Betas controlling all tests at FDR < 0.05)", xpd=TRUE, cex=1, font=1)
grid::grid.text(label="-log10(p-value)", x=0.93, y=0.825, gp=gpar(fontsize=9))
dev.off()

# AMY
pdf("graphics/heatmap-rev_amy-magma-v.08_PGC-and-addiction_GWAS_MNT2021.pdf", w=8)
par(mar=c(8,7.5,6,1), cex.axis=1.0, cex.lab=0.5)
customMAGMAplot.b(region="amy", Pthresh=12, fdrThresh=betacut.fdr, bonfThresh=betacut.bonf)
abline(v=7,lwd=3)
text(x = c(3.5,9.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=0.9, font=2)
text(x = 6, y=185, "MAGMA gene set analyses: AMY cell classes", xpd=TRUE, cex=1.5, font=2)
text(x = 6, y=175, "(Betas controlling all tests at FDR < 0.05)", xpd=TRUE, cex=1, font=1)
grid::grid.text(label="-log10(p-value)", x=0.93, y=0.825, gp=gpar(fontsize=9))
dev.off()



### Revision exploration =======
  #   ?: Cortical vs subcortical for substance use phenotypes?

# (First re-create 'magmaStats_long', through line 95)
magmaStats_cort.liu <- magmaStats_long[magmaStats_long$Region %in% c("dlpfc","sacc"), ]
magmaStats_cort.liu <- magmaStats_cort.liu[grep("addxn.", magmaStats_cort.liu$GWAS), ]

magmaStats_sub.liu <- magmaStats_long[magmaStats_long$Region %in% c("hpc","nac","amy"), ]
magmaStats_sub.liu <- magmaStats_sub.liu[grep("addxn.", magmaStats_sub.liu$GWAS), ]

quantile(magmaStats_cort.liu$Beta[magmaStats_cort.liu$P.adj.fdr < 0.05])
    #        0%        25%        50%        75%       100% 
    #0.04411300 0.05888325 0.07044550 0.08067500 0.20759000
quantile(magmaStats_sub.liu$Beta[magmaStats_sub.liu$P.adj.fdr < 0.05])
    #      0%      25%      50%      75%     100% 
    #0.040552 0.061550 0.078397 0.106460 0.305230   -mean is greater too -> t-test?

t.test(magmaStats_sub.liu$Beta[magmaStats_sub.liu$P.adj.fdr < 0.05],
       magmaStats_cort.liu$Beta[magmaStats_cort.liu$P.adj.fdr < 0.05],
       alternative="greater")
    # t = 2.5416, df = 138.5, p-value = 0.006068
    # alternative hypothesis: true difference in means is greater than 0
    # 95 percent confidence interval:
    #   0.005195748         Inf
    # sample estimates:
    #   mean of x  mean of y 
    # 0.09095591 0.07604576


### Session info for 20Aug2021 ================================================
sessionInfo()
# R version 4.0.4 RC (2021-02-08 r79975)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS:   /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/lib/libRblas.so
# LAPACK: /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/lib/libRlapack.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices datasets  utils     methods   base     
# 
# other attached packages:
#   [1] jaffelab_0.99.30   rafalib_1.0.0      fields_11.6        spam_2.6-0        
# [5] dotCall64_1.0-1    ggplot2_3.3.3      RColorBrewer_1.1-2 stringr_1.4.0     
# [9] readr_1.4.0       
# 
# loaded via a namespace (and not attached):
#   [1] SummarizedExperiment_1.20.0 tidyselect_1.1.1            reshape2_1.4.4             
# [4] purrr_0.3.4                 gargle_1.2.0                splines_4.0.4              
# [7] lattice_0.20-41             colorspace_2.0-0            vctrs_0.3.8                
# [10] generics_0.1.0              stats4_4.0.4                utf8_1.2.1                 
# [13] rlang_0.4.11                pillar_1.6.0                glue_1.4.2                 
# [16] withr_2.4.2                 DBI_1.1.1                   BiocGenerics_0.36.1        
# [19] segmented_1.3-4             plyr_1.8.6                  matrixStats_0.58.0         
# [22] GenomeInfoDbData_1.2.4      lifecycle_1.0.0             zlibbioc_1.36.0            
# [25] MatrixGenerics_1.2.1        munsell_0.5.0               gtable_0.3.0               
# [28] Biobase_2.50.0              IRanges_2.24.1              GenomeInfoDb_1.26.7        
# [31] parallel_4.0.4              fansi_0.4.2                 Rcpp_1.0.6                 
# [34] scales_1.1.1                limma_3.46.0                DelayedArray_0.16.3        
# [37] S4Vectors_0.28.1            XVector_0.30.0              fs_1.5.0                   
# [40] hms_1.0.0                   stringi_1.5.3               dplyr_1.0.5                
# [43] GenomicRanges_1.42.0        tools_4.0.4                 bitops_1.0-7               
# [46] magrittr_2.0.1              maps_3.3.0                  RCurl_1.98-1.3             
# [49] tibble_3.1.1                crayon_1.4.1                pkgconfig_2.0.3            
# [52] ellipsis_0.3.2              Matrix_1.3-4                googledrive_2.0.0          
# [55] assertthat_0.2.1            R6_2.5.0                    compiler_4.0.4
