####
library(readxl)
library(jaffelab)
library(readr)
library(pracma)
library(RColorBrewer)
options(width=100)

dir.create("triangle_pdfs")

## read in reference data
ref = read_excel("raw_data/Triangle_expected_expression.xlsx")
ref = as.data.frame(ref[,1:5])
ref_mat = as.matrix(ref[,2:5])
rownames(ref_mat) = ref$Population

head(ref_mat)

image(ref_mat)

## read in long data
dat =read.csv("raw_data/long_data_Triangle.csv",as.is=TRUE)
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

table(DRD1=dat$MD_Opal570Lp1_0 > 3, 
		DRD2=dat$MD_Opal620_LP10 > 3)
		

pdf("triangle_pdfs/DRD1_dotCount_lipoMasked_allROIs_hist.pdf",h=5,w=10)
par(mar= c(5,6,2,2), cex.axis=1.7, cex.lab=2, cex.main = 2)
hist(log2(dat$MD_Opal570Lp1_0+1),xaxt = "n",
	xlab = "# Dots (Lipofuscin-Masked)",
	main = "DRD1 (Channel 570)")
axis(1,at=0:8,label = 2^(0:8)-1) 
abline(v=log2(3-1),col="red",lwd=3)
dev.off()

pdf("triangle_pdfs/DRD2_dotCount_lipoMasked_allROIs_hist.pdf",h=5,w=10)
par(mar= c(5,6,2,2), cex.axis=1.7, cex.lab=2, cex.main = 2)
hist(log2(dat$MD_Opal620_LP10+1),xaxt = "n",
	xlab = "# Dots (Lipofuscin-Masked)",
	main = "DRD1 (Channel 620)")
axis(1,at=0:8,label = 2^(0:8)-1) 
abline(v=log2(3-1),col="red",lwd=3)
dev.off()

## filter to >3 dots in either d1 or d2
drd_index = which(dat$MD_Opal620_LP10 > 3 | dat$MD_Opal570Lp1_0 > 3)
dat_d = dat[drd_index,]
pheno_d = pheno[drd_index,]

######################
## other channels w/in drd1|2

### 520 
table(dat_d$MD_Opal520_Lp20 > 6)

pdf("triangle_pdfs/TAC1_dotCount_lipoMasked_D1or2ROIs_hist.pdf",h=5,w=10)
par(mar= c(5,6,2,2), cex.axis=1.7, cex.lab=2, cex.main = 2)
hist(log2(dat_d$MD_Opal520_Lp20+1),xaxt = "n",breaks=20,
	xlab = "# Dots (Lipofuscin-Masked)",
	main = "TAC1 (Channel 520)")
axis(1,at=0:8,label = 2^(0:8)-1) 
abline(v=log2(6-1),col="red",lwd=3)
dev.off()

### 690 

table(dat_d$MD_Opal690Lp30 > 6)

pdf("triangle_pdfs/PENK_dotCount_lipoMasked_D1or2ROIs_hist.pdf",h=5,w=10)
par(mar= c(5,6,2,2), cex.axis=1.7, cex.lab=2, cex.main = 2)
hist(log2(dat_d$MD_Opal690Lp30+1),breaks=20, xaxt = "n",
	xlab = "# Dots (Lipofuscin-Masked)",
	main = "PENK (Channel 690)")
axis(1,at=0:8,label = 2^(0:8)-1) 
abline(v=log2(6-1),col="red",lwd=3)
dev.off()


dat_d$DRD1 = as.numeric(dat_d$MD_Opal570Lp1_0  > 3)
dat_d$DRD2 = as.numeric(dat_d$MD_Opal620_LP10  > 3)
dat_d$TAC1 = as.numeric(dat_d$MD_Opal520_Lp20 > 6)
dat_d$PENK = as.numeric(dat_d$MD_Opal690Lp30 > 6)

###############################
## compare to reference #######
roi_mat = as.matrix(dat_d[,c("TAC1", "DRD1", "DRD2", "PENK")])
colMeans(roi_mat)

ref_mat_binary = ref_mat
ref_mat_binary[ref_mat_binary > 0] = 1
dist_mat = t(distmat(ref_mat_binary, roi_mat))

## label
d_match = cbind(roi_mat, round(dist_mat,2), 
	data.frame(cell_type = apply(dist_mat,1, which.min),
		dist =apply(dist_mat,1, min)))
d_match$cell_name = rownames(ref_mat)[d_match$cell_type]
d_match$cell_name = factor(d_match$cell_name, rownames(ref_mat))
table(d_match$dist)
  # 0   1
# 133 138

mean(d_match$dist==0)
# [1] 0.4907749

table(d_match$cell_name, d_match$dist)
        # 0  1
  # D1.1 17 38
  # D1.2  0  0
  # D1.3 37 82
  # D1.4 26  0
  # D2.1  7 18
  # D2.2 46  0


prop.table(table(d_match$cell_name, d_match$dist),2)

                # 0          1
  # D1.1 0.12781955 0.27536232
  # D1.2 0.00000000 0.00000000
  # D1.3 0.27819549 0.59420290
  # D1.4 0.19548872 0.00000000
  # D2.1 0.05263158 0.13043478
  # D2.2 0.34586466 0.00000000

########################
## make figure for paper

pdf("triangle_pdfs/supp_figure_triangle_rnascope.pdf")
par(mar=c(5,6,2,2), cex.axis=1.5,cex.lab=1.5,cex.main=2)
plot(log2(MD_Opal520_Lp20+1) ~ log2(MD_Opal690Lp30+1),
	data = dat_d,pch = 21+2*DRD1 , bg = DRD2 + 1,
	xlab = "PENK (log2 [masked dot count + 1])",
	ylab = "TAC1 (log2 [masked dot count + 1])",
	cex=2-d_match$dist, main = "Triangle")
abline(h=log2(7.5), v = log2(7.5),lwd=2,lty=2) 
legend("topleft", c("DRD2-", "DRD2+"), col = 1:2, pch = 15)
legend("bottomright", c("DRD1-", "DRD1+"), pt.bg = "black", pch = c(21,23))
text(x = c(1,4,4,1,6), y = c(0.5,6.4,6.1,5.5,0.5), 
	c("D1.1/2","D1.3","D1.4","D2.1","D2.2"), 
	col = c(1,2,1,2,2),cex=2)

dev.off()


## cross-tab
with(d_match[d_match$DRD2 == 1,],
	table(PENK = PENK,  TAC1 = TAC1))

with(d_match[d_match$DRD2 == 1 & d_match$DRD1== 0,],
	table(PENK = PENK,  TAC1 = TAC1))
	
		
## numbers for paper
nrow(pheno)
length(table(pheno$Image))
length(table(pheno$BrNum))
length(unique(paste0(pheno$Section,":", pheno$BrNum)))

dim(dat_d)