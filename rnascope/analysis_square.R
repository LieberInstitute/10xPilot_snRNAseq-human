####
library(readxl)
library(jaffelab)
library(readr)
library(pracma)
library(RColorBrewer)

dir.create("square_pdfs")

## read in reference data
ref = read_excel("raw_data/Square_expected_expression.xlsx")
ref = as.data.frame(ref[,1:5])
ref_mat = as.matrix(ref[,2:5])
rownames(ref_mat) = ref$Population

## read in long data
dat =read.csv("raw_data/long_data_Square.csv",as.is=TRUE)
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

pdf("square_pdfs/DRD1_dotCount_lipoMasked_allROIs_hist.pdf",h=5,w=10)
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
table(dat_d1$MD_Opal520_Lp20 > 6)

pdf("square_pdfs/RXFP2_dotCount_lipoMasked_D1ROIs_hist.pdf",h=5,w=10)
par(mar= c(5,6,2,2), cex.axis=1.7, cex.lab=2, cex.main = 2)
hist(log2(dat_d1$MD_Opal520_Lp20+1),xaxt = "n",
	xlab = "# Dots (Lipofuscin-Masked)",
	main = "RXFP2 (Channel 520)")
axis(1,at=0:8,label = 2^(0:8)-1) 
abline(v=log2(6-1),col="red",lwd=3)
dev.off()

### 570 
table(dat_d1$MD_Opal570Lp1_0 > 6)

pdf("square_pdfs/GABRQ_dotCount_lipoMasked_D1ROIs_hist.pdf",h=5,w=10)
par(mar= c(5,6,2,2), cex.axis=1.7, cex.lab=2, cex.main = 2)
hist(log2(dat_d1$MD_Opal570Lp1_0+1),breaks=20, xaxt = "n",
	xlab = "# Dots (Lipofuscin-Masked)",
	main = "GABRQ (Channel 570)")
axis(1,at=0:8,label = 2^(0:8)-1) 
abline(v=log2(6-1),col="red",lwd=3)
dev.off()

### 690 
table(dat_d1$MD_Opal690Lp30 > 6)

pdf("square_pdfs/TAC1_dotCount_lipoMasked_D1ROIs_hist.pdf",h=5,w=10)
par(mar= c(5,6,2,2), cex.axis=1.7, cex.lab=2, cex.main = 2)
hist(log2(dat_d1$MD_Opal690Lp30+1),breaks=20, xaxt = "n",
	xlab = "# Dots (Lipofuscin-Masked)",
	main = "TAC1 (Channel 690)")
axis(1,at=0:8,label = 2^(0:8)-1) 
abline(v=log2(6-1),col="red",lwd=3)
dev.off()

## recode
dat_d1$RXFP2 = as.numeric(dat_d1$MD_Opal520_Lp20 > 6)
dat_d1$GABRQ = as.numeric(dat_d1$MD_Opal570Lp1_0 > 6)
dat_d1$TAC1 = as.numeric(dat_d1$MD_Opal690Lp30 > 6)

###############################
## compare to reference #######
roi_mat = as.matrix(dat_d1[,c("RXFP2", "GABRQ", "TAC1")])
colMeans(roi_mat)

ref_mat_binary = ref_mat
ref_mat_binary[ref_mat_binary > 0] = 1
dist_mat = t(distmat(ref_mat_binary[,-3], roi_mat))

## label
d1_match = cbind(roi_mat, round(dist_mat,2), 
	data.frame(cell_type = apply(dist_mat,1, which.min),
		dist =apply(dist_mat,1, min)))
d1_match$cell_name = rownames(ref_mat)[d1_match$cell_type]

table(d1_match$dist)
  # 0   1
# 166 175
mean(d1_match$dist==0)
# [1] 0.4868035

table(d1_match$cell_name, d1_match$dist)
         # 0   1
  # D1.1   9  19
  # D1.2   6   0
  # D1.3  45 156
  # D1.4 106   0
  

prop.table(table(d1_match$cell_name, d1_match$dist), 2)

                # 0          1
  # D1.1 0.05421687 0.10857143
  # D1.2 0.03614458 0.00000000
  # D1.3 0.27108434 0.89142857
  # D1.4 0.63855422 0.00000000

#######################
## numbers for paper ##
#######################

## Square was about distinguishing the 4 different D1 subtypes.

pdf("square_pdfs/supp_figure_square_rnascope.pdf")
palette(brewer.pal(6,"Dark2"))
par(mar=c(5,6,2,2), cex.axis=1.5,cex.lab=1.5,cex.main=2)
plot(log2(MD_Opal520_Lp20+1) ~ log2(MD_Opal690Lp30+1),
	data = dat_d1,pch = 21, bg = GABRQ+1,main = "Square",
	ylab = "RXFP2 (log2 [masked dot count + 1])",
	xlab = "TAC1 (log2 [masked dot count + 1])",
	cex=2-d1_match$dist)
abline(h=log2(7.5), v = log2(7.5),lwd=2,lty=2) 
legend("topleft", c("GABRQ-", "GABRQ+"), col = 1:2, pch = 15,cex=1.5)
text(x=c(1,1,6,6), y = c(5,1,0.5,6.4), c("D1.2","D1.1","D1.3","D1.4"),
	col = c(1,2,2,1), cex=2)
dev.off()

## cross-tab
with(d1_match,table(RXFP2 = RXFP2,  GABRQ = GABRQ, TAC1=TAC1))
with(d1_match,table(RXFP2 = RXFP2,  TAC1=TAC1, GABRQ = GABRQ))


######################
## BELOW NOT USED ####
######################

## kmeans #####

dat_mask = dat_d1[,grep("M", colnames(dat_d1))]
# dat_mask = dat_mask[,-grep("MI", colnames(dat_mask))]

colIndex = splitit(ss(colnames(dat_mask), "Opal",2))
colIndex = colIndex[-3] # drop drd1

## run kmeans
kList= lapply(colIndex, function(ii) {
	d = as.matrix(dat_mask[,ii])
	d[is.nan(d)] = 0
	kmeans(d, centers = 3)
})

lapply(lapply(kList, "[[", "cluster"),table)
lapply(kList, "[[", "centers")

## relevel based on 
mat_qual_k3 = sapply(kList, function(x) {
	m = x$centers
	new = factor(x$cluster, 
		levels = order(m[,grep("^MD", colnames(m))]))
	as.numeric(new) - 1
})

## quick view
boxplot(log2(dat_d1$MD_Opal520_Lp20+1) ~ mat_qual_k3[,1])
boxplot(log2(dat_d1$MD_Opal570Lp1_0+1) ~ mat_qual_k3[,2])
boxplot(log2(dat_d1$MD_Opal690Lp30+1) ~ mat_qual_k3[,3])

apply(mat_qual_k3, 2, table)


################
## make boxplot for paper ##
plot_df = data.frame(dots_masked = 
	c(dat_d1$MD_Opal520_Lp20, dat_d1$MD_Opal570Lp1_0, dat_d1$MD_Opal690Lp30),
		gene = rep(c("RXFP2", "GABRQ", "TAC1"), each = nrow(dat_d1)),
		cluster = as.numeric(mat_qual_k3),
		exprs_group = as.numeric(roi_mat),
		dist = rep(d1_match$dist, times=3))
## make label
plot_df$plot_label = paste0(plot_df$gene, ":", plot_df$cluster)
plot_df$plot_label = factor(plot_df$plot_label ,
	c("RXFP2:0", "RXFP2:1", "RXFP2:2",
		"GABRQ:0", "GABRQ:1", "GABRQ:2", 
		"TAC1:0", "TAC1:1", "TAC1:2"))

par(mar= c(10,6,2,2), cex.axis=1.7, cex.lab=2, cex.main = 2)
boxplot(log2(dots_masked+1) ~ plot_label, data=plot_df,
	las = 3, outline = FALSE, yaxt= "n",xlab="",varwidth=TRUE)
axis(2,at=0:8,label = 2^(0:8)-1) 
	
	
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
pdf("square_pdfs/TAC1_square_rnascope.pdf")
boxplot(log2(dat_inhib$MD_Opal690Lp30 + 1) ~ d1_match$cell_name)
dev.off()

# or (with $RVolume normalization)
pdf("square_pdfs/TAC1_norm_square_rnascope.pdf")
boxplot(log2(dat_d1$MD_Opal690Lp30/dat_d1$RVolume + 1) ~ d1_match$cell_name)
dev.off()

## For KIT
pdf("square_pdfs/DRD1_square_rnascope.pdf")
boxplot(log2(dat_d1$MD_Opal620_LP10 + 1) ~ d1_match$cell_name)
dev.off()

# or (with $RVolume normalization)
pdf("square_pdfs/DRD1_norm_square_rnascope.pdf")
boxplot(log2(dat_d1$MD_Opal620_LP10/dat_d1$RVolume + 1) ~ d1_match$cell_name)
dev.off()

## For PVALB
pdf("square_pdfs/GABRQ_square_rnascope.pdf")
boxplot(log2(dat_d1$MD_Opal570Lp1_0 + 1) ~ d1_match$cell_name)
dev.off()

# or (with $RVolume normalization)
pdf("square_pdfs/GABRQ_norm_square_rnascope.pdf")
boxplot(log2(dat_d1$MD_Opal570Lp1_0/dat_d1$RVolume + 1) ~ d1_match$cell_name)
dev.off()

## For GAD1
pdf("square_pdfs/RXFP2_square_rnascope.pdf")
boxplot(log2(dat_d1$MD_Opal520_Lp20  + 1) ~ d1_match$cell_name)
dev.off()

# or (with $RVolume normalization)
pdf("square_pdfs/RXFP2_norm_square_rnascope.pdf")
boxplot(log2(dat_d1$MD_Opal520_Lp20/dat_d1$RVolume + 1) ~ d1_match$cell_name)
dev.off()

### 2) Hierarchically cluster the 212 (GAD1+) ROIs on normalized (i.e. `/RVolume`)
#    n transcript dots? (so 212 x 3 probes as input matrix)
#    (MD_* : total transcript dots in that ROI after Lipofuscin masking)
# - hopefully would see four (for four inhib. subpops) main branches??