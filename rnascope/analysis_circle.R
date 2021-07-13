####
library(readxl)
library(jaffelab)
library(readr)
library(pracma)
library(RColorBrewer)

dir.create("circle_pdfs")

## read in reference data
ref = read_excel("raw_data/Circle_expected_expression_revisionMNT.xlsx")
ref = as.data.frame(ref[,1:5])
ref_mat = as.matrix(ref[,2:5])
rownames(ref_mat) = ref$Population

## read in long data
dat =read.csv("raw_data/long_data_Circle.csv",as.is=TRUE)
colnames(dat)[1]= "Image"
dat$Image =gsub(" unmix", "Umix", dat$Image)

## add proportions of rois covered by dots
dat$M_PP_Opal520_Lp20 = dat$MP_Opal520_Lp20 / dat$RVolume
dat$M_PP_Opal570Lp1_0 = dat$MP_Opal570Lp1_0 / dat$RVolume
dat$M_PP_Opal620_LP10 = dat$MP_Opal620_LP10 / dat$RVolume
dat$M_PP_Opal690Lp30 = dat$MP_Opal690Lp30 / dat$RVolume

############
## make pheno
pheno = data.frame(Image = dat$Image, stringsAsFactors=FALSE)
pheno$BrNum = ss(pheno$Image, "_", 1)
pheno$Region = ss(pheno$Image, "_", 2)
pheno$Combo_Code = ss(pheno$Image, "_", 3)
pheno$Section = ss(pheno$Image, "_", 4)
pheno$Image_Number = ss(ss(pheno$Image, "690_", 2), "_Linear")
pheno$Image_Number[grep("_", pheno$Image_Number)] = ss(
	pheno$Image_Number[grep("_", pheno$Image_Number)], "_", 2)
pheno$Image_Number = as.numeric(pheno$Image_Number)

## by channel
pheno$gene_520 = gsub("520", "", ss(pheno$Image, "_", 5))
pheno$gene_570 = gsub("570", "", ss(pheno$Image, "_", 6))
pheno$gene_620 = gsub("620", "", ss(pheno$Image, "_", 7))
pheno$gene_690 = gsub("690", "", ss(pheno$Image, "_", 8))

## make shorter combo variable
pheno$Combo_Genes = paste(pheno$gene_520, pheno$gene_570,
				pheno$gene_620, pheno$gene_690, sep = "_")

######################################
## determine drd1 expression cutoff ##
######################################

table(dat$MD_Opal620_LP10 > 3)

pdf("circle_pdfs/DRD1_dotCount_lipoMasked_allROIs_hist.pdf",h=5,w=10)
par(mar= c(5,6,2,2), cex.axis=1.7, cex.lab=2, cex.main = 2)
hist(log2(dat$MD_Opal620_LP10+1),xaxt = "n",
	xlab = "# Dots (Lipofuscin-Masked)",
	main = "DRD1 (Channel 620)")
axis(1,at=0:8,label = 2^(0:8)-1) 
abline(v=log2(3-1),col="red",lwd=3)
dev.off()

## filter to >3 dots
drd1_index = which(dat$MD_Opal620_LP10 > 3)
dat_d1 = dat[drd1_index,]
pheno_d1 = pheno[drd1_index,]

######################
## other channels w/in drd1

### 520 
table(dat_d1$MD_Opal520_Lp20 > 3)

pdf("circle_pdfs/RXFP2_dotCount_lipoMasked_D1ROIs_hist.pdf",h=5,w=10)
par(mar= c(5,6,2,2), cex.axis=1.7, cex.lab=2, cex.main = 2)
hist(log2(dat_d1$MD_Opal520_Lp20+1),xaxt = "n",breaks=20,
	xlab = "# Dots (Lipofuscin-Masked)",
	main = "RXFP2 (Channel 520)")
axis(1,at=0:8,label = 2^(0:8)-1) 
abline(v=log2(3-1),col="red",lwd=3)
dev.off()



### 570 

table(dat_d1$MD_Opal570Lp1_0 > 3)

pdf("circle_pdfs/CRHR2_dotCount_lipoMasked_D1ROIs_hist.pdf",h=5,w=10)
par(mar= c(5,6,2,2), cex.axis=1.7, cex.lab=2, cex.main = 2)
hist(log2(dat_d1$MD_Opal570Lp1_0+1),breaks=20, xaxt = "n",
	xlab = "# Dots (Lipofuscin-Masked)",
	main = "CRHR2 (Channel 570)")
axis(1,at=0:8,label = 2^(0:8)-1) 
abline(v=log2(3-1),col="red",lwd=3)
dev.off()


### 690 

table(dat_d1$MD_Opal690Lp30 > 6)

pdf("circle_pdfs/TAC1_dotCount_lipoMasked_D1ROIs_hist.pdf",h=5,w=10)
par(mar= c(5,6,2,2), cex.axis=1.7, cex.lab=2, cex.main = 2)
hist(log2(dat_d1$MD_Opal690Lp30+1),breaks=20, xaxt = "n",
	xlab = "# Dots (Lipofuscin-Masked)",
	main = "TAC1 (Channel 690)")
axis(1,at=0:8,label = 2^(0:8)-1) 
abline(v=log2(6-1),col="red",lwd=3)
dev.off()


dat_d1$RXFP1 = as.numeric(dat_d1$MD_Opal520_Lp20 > 3)
dat_d1$CRHR2 = as.numeric(dat_d1$MD_Opal570Lp1_0 > 3)
dat_d1$TAC1 = as.numeric(dat_d1$MD_Opal690Lp30 > 6)

###############################
## compare to reference #######
roi_mat = as.matrix(dat_d1[,c("RXFP1", "CRHR2", "TAC1")])
colMeans(roi_mat)

ref_mat_binary = ref_mat
ref_mat_binary[ref_mat_binary > 0] = 1
dist_mat = t(distmat(ref_mat_binary[,-3], roi_mat))

## label
d1_match = cbind(roi_mat, round(dist_mat,2), 
	data.frame(cell_type = apply(dist_mat,1, which.min),
		dist =apply(dist_mat,1, min)))
d1_match$cell_name = rownames(ref_mat)[d1_match$cell_type]
d1_match$cell_name = factor(d1_match$cell_name, rownames(ref_mat))

####
## metrics

table(d1_match$dist)
  # 0   1
# 204  47

mean(d1_match$dist==0)
# [1] 0.812749

table(d1_match$cell_name, d1_match$dist)

         # 0   1
  # D1.1   6   5
  # D1.2   2  38
  # D1.3  64   4
  # D1.4 132   0

prop.table(table(d1_match$cell_name, d1_match$dist), 2)

                 # 0           1
  # D1.1 0.029411765 0.106382979
  # D1.2 0.009803922 0.808510638
  # D1.3 0.313725490 0.085106383
  # D1.4 0.647058824 0.000000000

###############
## plotting ###

## circle  was about distinguishing the 4 different D1 subtypes.

pdf("circle_pdfs/supp_figure_circle_rnascope.pdf")
palette(brewer.pal(6,"Dark2"))
par(mar=c(5,6,2,2), cex.axis=1.5,cex.lab=1.5,cex.main=2)
plot(log2(MD_Opal520_Lp20+1) ~ log2(MD_Opal690Lp30+1),
	data = dat_d1,pch = 21, bg = CRHR2+1,main = "Circle",
	ylab = "RXFP1 (log2 [masked dot count + 1])",
	xlab = "TAC1 (log2 [masked dot count + 1])",
	cex=2-d1_match$dist)
abline(h=log2(4.5), v = log2(7.5),lwd=2,lty=2) 
legend("topleft", c("CRHR2-", "CRHR2+"), 
	col = 1:2, pch = 15,cex=1.5,bty="n")

text(x=c(1,1,4,4), y = c(5,1,6.7,6.2), c("D1.1","D1.2","D1.3","D1.4"),
	col = c(1,1,2,1), cex=2)
dev.off()

## cross-tab
with(d1_match,table(RXFP1 = RXFP1, TAC1=TAC1,CRHR2 = CRHR2))

## numbers for paper
nrow(pheno)
length(table(pheno$Image))
length(table(pheno$BrNum))
length(unique(paste0(pheno$Section,":", pheno$BrNum)))

dim(dat_d1)

### MNT / ABS update May2021 =============
# objective: make more-interpretable graphics with RNAscope quantification data
#(Abby can add/save this to the end of Andrew's script - analysis_square.R:)
# Some ideas ===
### 1) boxplot log2-transformed (/maybe add RVolume normalization)
## For PTHLH (Opal690)
pdf("circle_pdfs/TAC1_circle_rnascope.pdf")
boxplot(log2(dat_d1$MD_Opal690Lp30 + 1) ~ d1_match$cell_name)
dev.off()

# or (with $RVolume normalization)
pdf("circle_pdfs/TAC1_norm_circle_rnascope.pdf")
boxplot(log2(dat_d1$MD_Opal690Lp30/dat_d1$RVolume + 1) ~ d1_match$cell_name)
dev.off()

## For KIT
pdf("circle_pdfs/DRD1_circle_rnascope.pdf")
boxplot(log2(dat_d1$MD_Opal620_LP10 + 1) ~ d1_match$cell_name)
dev.off()

# or (with $RVolume normalization)
pdf("circle_pdfs/DRD1_norm_circle_rnascope.pdf")
boxplot(log2(dat_d1$MD_Opal620_LP10/dat_d1$RVolume + 1) ~ d1_match$cell_name)
dev.off()

## For PVALB
pdf("circle_pdfs/CRHR2circle_rnascope.pdf")
boxplot(log2(dat_d1$MD_Opal570Lp1_0 + 1) ~ d1_match$cell_name)
dev.off()

# or (with $RVolume normalization)
pdf("circle_pdfs/CRHR2_norm_circle_rnascope.pdf")
boxplot(log2(dat_d1$MD_Opal570Lp1_0/dat_d1$RVolume + 1) ~ d1_match$cell_name)
dev.off()

## For GAD1
pdf("circle_pdfs/RXFP1_circle_rnascope.pdf")
boxplot(log2(dat_d1$MD_Opal520_Lp20  + 1) ~ d1_match$cell_name)
dev.off()

# or (with $RVolume normalization)
pdf("circle_pdfs/RXFP1_norm_circle_rnascope.pdf")
boxplot(log2(dat_d1$MD_Opal520_Lp20/dat_d1$RVolume + 1) ~ d1_match$cell_name)
dev.off()


### MNT add 12Jul2021 ===
  # Above boxplots with predictions to the new cell classes. Prediction groups combine
  #   cell classes where not able to be differentiated, given only four probes
quantile(dat_d1$RVolume)
    #    0%      25%      50%      75%     100% 
    # 480.0  11062.5  15648.0  20667.5 128348.0 
    
    # -> Normalization scheme: by '10k ROI pixels' keeps relatively on the same scale
    #    (and is more interpretable than 0.003, e.g.)

pdf("circle_pdfs/quantified-dotsPerROI_all-circle-probes_MNT2021.pdf")
# TAC1
boxplot(log2(dat_d1$MD_Opal690Lp30 + 1) ~ d1_match$cell_name,
        ylab="log2( TAC1 dots per ROI )",
        main="TAC1, quantified by predicted cell class")
boxplot(log2(dat_d1$MD_Opal690Lp30/dat_d1$RVolume*10000 + 1) ~ d1_match$cell_name,
        ylab="log2( TAC1 dots per 10k ROI pixels )",
        main="TAC1, quantified by predicted cell class \n (ROI volume-normalized)")
# DRD1
boxplot(log2(dat_d1$MD_Opal620_LP10 + 1) ~ d1_match$cell_name,
        ylab="log2( DRD1 dots per ROI )",
        main="DRD1, quantified by predicted cell class")
boxplot(log2(dat_d1$MD_Opal620_LP10/dat_d1$RVolume*10000 + 1) ~ d1_match$cell_name,
        ylab="log2( DRD1 dots per 10k ROI pixels )",
        main="DRD1, quantified by predicted cell class \n (ROI volume-normalized)")
# CRHR2
boxplot(log2(dat_d1$MD_Opal570Lp1_0 + 1) ~ d1_match$cell_name,
        ylab="log2( CRHR2 dots per ROI )",
        main="CRHR2, quantified by predicted cell class")
boxplot(log2(dat_d1$MD_Opal570Lp1_0/dat_d1$RVolume*10000 + 1) ~ d1_match$cell_name,
        ylab="log2( CRHR2 dots per 10k ROI pixels )",
        main="CRHR2, quantified by predicted cell class \n (ROI volume-normalized)")
# RXFP1
boxplot(log2(dat_d1$MD_Opal520_Lp20  + 1) ~ d1_match$cell_name,
        ylab="log2( RXFP1 dots per ROI )",
        main="RXFP1, quantified by predicted cell class")
boxplot(log2(dat_d1$MD_Opal520_Lp20/dat_d1$RVolume*10000  + 1) ~ d1_match$cell_name,
        ylab="log2( RXFP1 dots per 10k ROI pixels )",
        main="RXFP1, quantified by predicted cell class \n (ROI volume-normalized)")
dev.off()


## Forcing to use like SCE data, for aesthetics ===
library(scater)
library(SingleCellExperiment)
source('../plotExpressionCustom.R')
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")
blgngy <- tableau20[c(1:2, 19:20, 5:6, 17:18, 15:16)]

MDcounts <- t(dat_d1[ ,colnames(dat_d1)[grep('^MD', colnames(dat_d1))]])
rownames(MDcounts) <- c("rnascope_RXFP1", "rnascope_CRHR2",
                        "rnascope_DRD1","rnascope_TAC1")
# check
table(colnames(MDcounts) == rownames(d1_match))
    # all 251 TRUE
colnames(d1_match)[colnames(d1_match)=="cell_name"] <- "class_predict"

# Create SCE
sce.circle <- SingleCellExperiment(list(counts=MDcounts), colData=d1_match)
# a couple normalizations:
assay(sce.circle, "logcounts") <- log2(assay(sce.circle, "counts") + 1)
assay(sce.circle, "logcounts.RVnorm") <- t(apply(assay(sce.circle, "counts"),
                                               1,function(x){
  log2( x/dat_d1$RVolume*10000 + 1 )
}))


pdf("circle_pdfs/quantified-dotsPerROI_all-circle-probes_betterVlnPlots_MNT2021.pdf", width=2)
# Log2-transform, only
plotExpressionCustom(sce.circle, exprs_values="logcounts", scales="free_y",
                     features=c("rnascope_DRD1", "rnascope_TAC1",
                                "rnascope_RXFP1", "rnascope_CRHR2"),
                     anno_name="class_predict", point_alpha=0.8, point_size=1.75, ncol=1,
                     features_name = "RNAscope quantified by predicted cell class",
                     xlab="Distance-predicted cell class [group]") +
  scale_color_manual(values = blgngy)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        axis.title.x = element_text(size = 7),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(angle = 90, size = 12),
        plot.title = element_text(size = 6),
        panel.grid.major=element_line(colour="grey95", size=0.8),
        panel.grid.minor=element_line(colour="grey95", size=0.4))

# + RVolume normalization
plotExpressionCustom(sce.circle, exprs_values="logcounts.RVnorm", scales="free_y",
                     features=c("rnascope_DRD1", "rnascope_TAC1",
                                "rnascope_RXFP1", "rnascope_CRHR2"),
                     anno_name="class_predict", point_alpha=0.8, point_size=1.75, ncol=1,
                     features_name = "RNAscope quantified by predicted cell class (ROI volume-normalized)",
                     xlab="Distance-predicted cell class [group]") +
  scale_color_manual(values = blgngy)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        axis.title.x = element_text(size = 7),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(angle = 90, size = 12),
        plot.title = element_text(size = 6),
        panel.grid.major=element_line(colour="grey95", size=0.8),
        panel.grid.minor=element_line(colour="grey95", size=0.4))
dev.off()

