### for MAGMA with LIBD 10x pilot analyses
  #     - plotting results/heatmaps
  # MNT 24Aug2020 ======================

library(readr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(fields)
library(jaffelab)



# 
# x <- paste0("./Results/",i,"_clozuk_pgc2.gsa.out")
# 
# ## From AnJa
#     parse_magma = function(x) {
#       y = read_delim(x, delim = " ", skip=4)
#       colnames(y) = str_trim(colnames(y))
#       y = apply(y, 2, str_trim)
#       y = as.data.frame(y)
#       for(i in 3:6) y[,i] = as.numeric(y[,i])
#       return(y)
#     }



regions <- c("dlpfc","sacc","hpc","nac","amy")
magmaStats <- list()

for(i in regions){
  magmaStats[[i]][["PGC2.SCZD"]] <- read.table(paste0("./Results/",i,"_clozuk_pgc2.gsa.out"), header=T)
  magmaStats[[i]][["PGC3.SCZD"]] <- read.table(paste0("./Results/",i,"_pgc3_scz.gsa.out"), header=T)
  magmaStats[[i]][["PGC.ASD"]] <- read.table(paste0("./Results/",i,"_PGC_ASD.gsa.out"), header=T)
  magmaStats[[i]][["PGC.BPD"]] <- read.table(paste0("./Results/",i,"_PGC_BIP.gsa.out"), header=T)
  magmaStats[[i]][["PGC.MDD"]] <- read.table(paste0("./Results/",i,"_MDD29_23andMe.gsa.out"), header=T)
  # Addiction GWAS
  magmaStats[[i]][["addxn.AgeSmk"]] <- read.table(paste0("./Results/",i,"_AgeSmk.gsa.out"), header=T)
  magmaStats[[i]][["addxn.CigDay"]] <- read.table(paste0("./Results/",i,"_CigDay.gsa.out"), header=T)
  magmaStats[[i]][["addxn.DrnkWk"]] <- read.table(paste0("./Results/",i,"_DrnkWk.gsa.out"), header=T)
  magmaStats[[i]][["addxn.SmkInit"]] <- read.table(paste0("./Results/",i,"_SmkInit.gsa.out"), header=T)
  magmaStats[[i]][["addxn.SmkCes"]] <- read.table(paste0("./Results/",i,"_SmkCes.gsa.out"), header=T)
}

## merge to assess significance thresholds
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
  # 384   296
betacut.fdr <- max(magmaStats_long$P[p.adjust(magmaStats_long$P, "fdr") < 0.05])
# [1] 0.021686

table(p.adjust(magmaStats_long$P, "bonf") < 0.05)
# FALSE  TRUE
  # 562   118
betacut.bonf <- max(magmaStats_long$P[p.adjust(magmaStats_long$P, "bonf") < 0.05])
# [1] 6.7218e-05




## Adapted from Spatial_dlpfc paper - custom heatmap ===

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
      
      # # Make "" if for any p's < 80% decile
      # BETAcut <- quantile(wide_p, probs=seq(0.1,1,by=0.1))["80%"]
      # Update MNT 27Aug: Use FDR/Bonf cutoffs, set above
      wide_beta[wide_p < -log10(BETAcut)] = ""
      
      
      ## Plot
      clusterHeights <- seq(0,160,length.out=nrow(wide_p)+1)
      mypal = c("white", colorRampPalette(brewer.pal(9,"YlOrRd"))(50))
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
pdf("graphics/heatmap_dlpfc-magma_PGC-and-addiction_GWAS_MNT27Aug2020.pdf", w=8)
par(mar=c(8,7.5,6,1), cex.axis=1.0, cex.lab=0.5)
customMAGMAplot(region="dlpfc", Pthresh=12, BETAcut=betacut.fdr)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.05, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: DLPFC subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests at FDR < 0.05)", xpd=TRUE, cex=1, font=1)

customMAGMAplot(region="dlpfc", Pthresh=12, BETAcut=betacut.bonf)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.05, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: DLPFC subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests with Bonf. < 0.05)", xpd=TRUE, cex=1, font=1)
dev.off()
    ## BETAcut (80% percentile -log10(p-value) = 3.8836)


# sACC
pdf("graphics/heatmap_sacc-magma_PGC-and-addiction_GWAS_MNT27Aug2020.pdf", w=8)
par(mar=c(8,7.5,6,1), cex.axis=1.3, cex.lab=0.5)
customMAGMAplot(region="sacc", Pthresh=12, BETAcut=betacut.fdr)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.05, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: sACC subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests at FDR < 0.05)", xpd=TRUE, cex=1, font=1)

customMAGMAplot(region="sacc", Pthresh=12, BETAcut=betacut.bonf)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.05, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: sACC subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests with Bonf. < 0.05)", xpd=TRUE, cex=1, font=1)
dev.off()
    ## BETAcut (80% percentile -log10(p-value) = 3.7592)


# HPC
pdf("graphics/heatmap_hpc-magma_PGC-and-addiction_GWAS_MNT27Aug2020.pdf", w=8)
par(mar=c(8,7.5,5,1), cex.axis=1.1, cex.lab=0.5)
customMAGMAplot(region="hpc", Pthresh=12, BETAcut=betacut.fdr)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.05, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: HIPPO subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests at FDR < 0.05)", xpd=TRUE, cex=1, font=1)

customMAGMAplot(region="hpc", Pthresh=12, BETAcut=betacut.bonf)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.05, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: HIPPO subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests with Bonf. < 0.05)", xpd=TRUE, cex=1, font=1)
dev.off()
    ## BETAcut (80% percentile -log10(p-value) = 4.3684)


# NAc
pdf("graphics/heatmap_nac-magma_PGC-and-addiction_GWAS_MNT27Aug2020.pdf", w=8)
par(mar=c(8,7.5,6,1), cex.axis=1.0, cex.lab=0.5)
customMAGMAplot(region="nac", Pthresh=12, BETAcut=betacut.fdr)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.05, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: NAc subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests at FDR < 0.05)", xpd=TRUE, cex=1, font=1)

customMAGMAplot(region="nac", Pthresh=12, BETAcut=betacut.bonf)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.05, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: NAc subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests with Bonf. < 0.05)", xpd=TRUE, cex=1, font=1)
dev.off()
    ## BETAcut (80% percentile -log10(p-value) = 2.8812)


# Amy
pdf("graphics/heatmap_amy-magma_PGC-and-addiction_GWAS_MNT27Aug2020.pdf", w=8)
par(mar=c(8,7.5,5,1), cex.axis=1.2, cex.lab=0.5)
customMAGMAplot(region="amy", Pthresh=12, BETAcut=betacut.fdr)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.05, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: AMY subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests at FDR < 0.05)", xpd=TRUE, cex=1, font=1)

customMAGMAplot(region="amy", Pthresh=12, BETAcut=betacut.bonf)
abline(v=5,lwd=3)
text(x = c(2.5,7.5), y=165, c("Psychiatric disorder GWAS", "Alcohol/substance use GWAS"), xpd=TRUE, cex=1.05, font=2)
text(x = 5, y=185, "MAGMA gene set analyses: AMY subclusters", xpd=TRUE, cex=1.5, font=2)
text(x = 5, y=175, "(Betas controlling all tests with Bonf. < 0.05)", xpd=TRUE, cex=1, font=1)
dev.off()
    ## BETAcut (80% percentile -log10(p-value) = 3.8006)




