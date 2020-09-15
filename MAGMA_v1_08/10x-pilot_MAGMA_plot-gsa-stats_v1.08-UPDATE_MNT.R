### for MAGMA with LIBD 10x pilot analyses
  #     - plotting results/heatmaps
  # UPDATE: re-running with new v1.08
  # nvm - don't include this 13Sep2020: PGC3-SCZD GWAS now published. lol. 
  # MNT 15Sep2020 ======================

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
  magmaStats[[i]][["PGC2.SCZ"]] <- read.table(paste0("./Results/",i,"_clozuk_pgc2.gsa.out"), header=T)
  #magmaStats[[i]][["PGC3.SCZ"]] <- read.table(paste0("./Results/",i,"_pgc3_scz.gsa.out"), header=T)
  magmaStats[[i]][["PGC.ASD"]] <- read.table(paste0("./Results/",i,"_PGC_ASD.gsa.out"), header=T)
  magmaStats[[i]][["PGC.BIP"]] <- read.table(paste0("./Results/",i,"_PGC_BIP.gsa.out"), header=T)
  magmaStats[[i]][["PGC.MDD"]] <- read.table(paste0("./Results/",i,"_MDD29_23andMe.gsa.out"), header=T)
  # PTSD - added 09Sep
  magmaStats[[i]][["PGCFr2.PTSD"]] <- read.table(paste0("./Results/",i,"_PTSD.gsa.out"), header=T)
  
  # Addiction GWAS
  magmaStats[[i]][["addxn.AgeSmk"]] <- read.table(paste0("./Results/",i,"_AgeSmk.gsa.out"), header=T)
  magmaStats[[i]][["addxn.CigDay"]] <- read.table(paste0("./Results/",i,"_CigDay.gsa.out"), header=T)
  magmaStats[[i]][["addxn.DrnkWk"]] <- read.table(paste0("./Results/",i,"_DrnkWk.gsa.out"), header=T)
  magmaStats[[i]][["addxn.SmkInit"]] <- read.table(paste0("./Results/",i,"_SmkInit.gsa.out"), header=T)
  magmaStats[[i]][["addxn.SmkCes"]] <- read.table(paste0("./Results/",i,"_SmkCes.gsa.out"), header=T)
  # Added for control but not including in multiple test corrxn, and just for supplement:
  #   - CARDIoGRAM+C4D coronary artery dx meta-GWAS
  #magmaStats[[i]][["CoronArtDx"]] <- read.table(paste0("./Results/",i,"_CAD.gsa.out"), header=T)
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


table(p.adjust(magmaStats_long$P, "fdr") < 0.05)
    # FALSE  TRUE
    #   444   236 (PGC2-SCZ-only + PTSD, no CAD; v1.08)
betacut.fdr <- max(magmaStats_long$P[p.adjust(magmaStats_long$P, "fdr") < 0.05])
    # [1] 0.017309 (PGC2-SCZ-only + PTSD, no CAD; v1.08)
table(p.adjust(magmaStats_long$P, "bonf") < 0.05)
    # FALSE  TRUE
    #   620    60 (PGC2-SCZ-only + PTSD, no CAD; v1.08)
betacut.bonf <- max(magmaStats_long$P[p.adjust(magmaStats_long$P, "bonf") < 0.05])
    #[1] 6.7732e-05 (PGC2-SCZ-only + PTSD, no CAD; v1.08)


    ## Btw - for Discussion:
              ## Thresholds are (PGC3-SCZ-only + PTSD, no CAD; v1.08)
              # betacut.bonf
              #   [1] 6.7732e-05
              # betacut.fdr
              #   [1] 0.017309
    
    magmaStats_long$P.adj.fdr <- p.adjust(magmaStats_long$P, "fdr")
              
    magmaStats_long[which(magmaStats_long$P.adj.fdr < 0.05 & magmaStats_long$GWAS=="PGCFr2.PTSD"), ]
        #     Region       CellType        GWAS         P   P.adj.fdr
        # 279  dlpfc Excit.L6.broad PGCFr2.PTSD 0.0048336 0.020935338
        # 289  dlpfc          Oligo PGCFr2.PTSD 0.0035806 0.016340993
        # 291   sacc        Excit.1 PGCFr2.PTSD 0.0017134 0.009488194
        # 303    hpc      Excit.3.1 PGCFr2.PTSD 0.0154070 0.046153128
        # 310    hpc      Inhib.5.1 PGCFr2.PTSD 0.0015284 0.008660933
        # 313    hpc        Oligo.2 PGCFr2.PTSD 0.0093868 0.033244917
        # 316    nac      Inhib.1.3 PGCFr2.PTSD 0.0084156 0.031442901
        # 317    nac      Inhib.2.3 PGCFr2.PTSD 0.0169410 0.049230256
        # 319    nac      Inhib.4.2 PGCFr2.PTSD 0.0037616 0.017052587
        # 320    nac       MSN.D1.1 PGCFr2.PTSD 0.0100760 0.034957551
        # 324    nac       MSN.D2.1 PGCFr2.PTSD 0.0024500 0.012250000
        # 328    nac        Oligo.3 PGCFr2.PTSD 0.0146640 0.045027027
        # 330    amy      Excit.1.2 PGCFr2.PTSD 0.0115380 0.038871089
        # 337    amy      Inhib.5.2 PGCFr2.PTSD 0.0141300 0.044075229 *TLL1/NPFFR2 (HPA stress axis-involved)
        # 340    amy        Oligo.4 PGCFr2.PTSD 0.0143040 0.044414247      - none meet Bonf. for PTSD


    
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

## Supplementary table: PGC3.SCZ in lieu of PGC2 & PTSD ===
write.csv(magmaStats_long, file = "../tables/suppTable_magma-v1.08_byCluster_GWASinMainSection.csv",
          row.names=F)

# PTSD results, separately because not part of main section
#   -> first include all GWAS in 'magmaStats' & remake 'magmaStats_long'/'magmaStats_long.beta'
magmaStats.sup2 <- magmaStats_long[magmaStats_long$GWAS %in% c("PGC2.SCZ","PGCFr2.PTSD"), ]
write.csv(magmaStats.sup2, file = "../tables/suppTable_magma-v1.08_byCluster_GWAS-PGC2-CAD_separately.csv",
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
       xpd=TRUE, srt=45,cex=1.25,adj= 1)
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
pdf("graphics/heatmap-v2_dlpfc-magma-v.08_PGC-and-addiction_GWAS_MNT15Sep2020.pdf", w=8)
par(mar=c(8,7.5,6,1), cex.axis=1.0, cex.lab=0.5)
customMAGMAplot.b(region="dlpfc", Pthresh=12, fdrThresh=betacut.fdr, bonfThresh=betacut.bonf)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.0, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: DLPFC subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests at FDR < 0.05)", xpd=TRUE, cex=1, font=1)
dev.off()

# sACC
pdf("graphics/heatmap-v2_sacc-magma-v.08_PGC-and-addiction_GWAS_MNT15Sep2020.pdf", w=8)
par(mar=c(8,7.5,6,1), cex.axis=1.3, cex.lab=0.5)
customMAGMAplot.b(region="sacc", Pthresh=12, fdrThresh=betacut.fdr, bonfThresh=betacut.bonf)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.0, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: sACC subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests at FDR < 0.05)", xpd=TRUE, cex=1, font=1)
dev.off()

# HIPPO
pdf("graphics/heatmap-v2_hpc-magma-v.08_PGC-and-addiction_GWAS_MNT15Sep2020.pdf", w=8)
par(mar=c(8,7.5,5,1), cex.axis=1.1, cex.lab=0.5)
customMAGMAplot.b(region="hpc", Pthresh=12, fdrThresh=betacut.fdr, bonfThresh=betacut.bonf)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.0, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: HIPPO subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests at FDR < 0.05)", xpd=TRUE, cex=1, font=1)
dev.off()

# NAc
pdf("graphics/heatmap-v2_nac-magma-v.08_PGC-and-addiction_GWAS_MNT15Sep2020.pdf", w=8)
par(mar=c(8,7.5,6,1), cex.axis=1.0, cex.lab=0.5)
customMAGMAplot.b(region="nac", Pthresh=12, fdrThresh=betacut.fdr, bonfThresh=betacut.bonf)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.0, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: NAc subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests at FDR < 0.05)", xpd=TRUE, cex=1, font=1)
dev.off()

# AMY
pdf("graphics/heatmap-v2_amy-magma-v.08_PGC-and-addiction_GWAS_MNT15Sep2020.pdf", w=8)
par(mar=c(8,7.5,5,1), cex.axis=1.2, cex.lab=0.5)
customMAGMAplot.b(region="amy", Pthresh=12, fdrThresh=betacut.fdr, bonfThresh=betacut.bonf)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.0, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: AMY subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests at FDR < 0.05)", xpd=TRUE, cex=1, font=1)
dev.off()


