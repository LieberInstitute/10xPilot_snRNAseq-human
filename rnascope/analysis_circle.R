####
library(readxl)
library(jaffelab)
library(readr)
library(pracma)
library(RColorBrewer)

dir.create("circle_pdfs")

## read in reference data
ref = read_excel("raw_data/Circle_expected_expression.xlsx")
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