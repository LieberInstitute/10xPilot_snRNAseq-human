####
library(readxl)
library(jaffelab)
library(readr)
library(pracma)
library(RColorBrewer)
options(width=100)

dir.create("triangle_pdfs")

## read in reference data
ref = read_excel("raw_data/Triangle_expected_expression_revisionMNT.xlsx")
ref = as.data.frame(ref[,1:5])
ref_mat = as.matrix(ref[,2:5])
rownames(ref_mat) = ref$Population
ref_mat

# Swap 'D1_A/D' & 'D1_E' row order
ref_mat <- ref_mat[c(1,3,2,4:5), ]
# Call 'D1_E' expression of DRD2 == 1 because the median is non-0
ref_mat["D1_E", "Channel_1_620_DRD2"] <- 1
ref_mat

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
    #  0   1 
    # 86 185

mean(d_match$dist==0)
   # [1] 0.3173432

table(d_match$cell_name, d_match$dist)
    #           0  1
    # D1_B/C/F 17 30
    # D1_A/D    8 26
    # D1_E      8 37
    # D2_C      7 48
    # D2_A/B/D 46 44


prop.table(table(d_match$cell_name, d_match$dist),2)

    #                   0          1
    # D1_B/C/F 0.19767442 0.16216216
    # D1_A/D   0.09302326 0.14054054
    # D1_E     0.09302326 0.20000000
    # D2_C     0.08139535 0.25945946
    # D2_A/B/D 0.53488372 0.23783784

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

## Forcing to use like SCE data, for aesthetics ===
library(scater)
library(SingleCellExperiment)
source('../plotExpressionCustom.R')
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")
blgngy <- tableau20[c(1:2, 19:20, 5:6, 17:18, 15:16)]

MDcounts <- t(dat_d[ ,colnames(dat_d)[grep('^MD', colnames(dat_d))]])
rownames(MDcounts)
    #[1] "MD_Opal520_Lp20" "MD_Opal570Lp1_0" "MD_Opal620_LP10" "MD_Opal690Lp30"
# From ref_mat:
rownames(MDcounts) <- c("rnascope_TAC1", "rnascope_DRD1",
                        "rnascope_DRD2","rnascope_PENK1")

# check
table(colnames(MDcounts) == rownames(d_match))
    # all 271 TRUE
colnames(d_match)[colnames(d_match)=="cell_name"] <- "class_predict"

# Create SCE
sce.triangle <- SingleCellExperiment(list(counts=MDcounts), colData=d_match)
# a couple normalizations:
assay(sce.triangle, "logcounts") <- log2(assay(sce.triangle, "counts") + 1)
assay(sce.triangle, "logcounts.RVnorm") <- t(apply(assay(sce.triangle, "counts"),
                                                 1,function(x){
                                                   log2( x/dat_d$RVolume*10000 + 1 )
                                                 }))


pdf("triangle_pdfs/quantified-dotsPerROI_all-triangle-probes_betterVlnPlots_MNT2021.pdf", width=2)
# Log2-transform, only
plotExpressionCustom(sce.triangle, exprs_values="logcounts", scales="free_y",
                     features=c("rnascope_DRD1", "rnascope_DRD2",
                                "rnascope_TAC1","rnascope_PENK1"),
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
plotExpressionCustom(sce.triangle, exprs_values="logcounts.RVnorm", scales="free_y",
                     features=c("rnascope_DRD1", "rnascope_DRD2",
                                "rnascope_TAC1","rnascope_PENK1"),
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

