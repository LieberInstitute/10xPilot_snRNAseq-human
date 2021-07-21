####
library(readxl)
library(jaffelab)
library(readr)
library(pracma)

dir.create("star_pdfs")
options(width=100)

## read in reference data
ref = read_excel("raw_data/Star_expected_expression_revisionMNT.xlsx")
ref = as.data.frame(ref[,1:5])
ref_mat = as.matrix(ref[,2:5])
rownames(ref_mat) = ref$Population

## read in long data
dat =read.csv("raw_data/long_data_Star.csv",as.is=TRUE)
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

table(DRD1=dat$MD_Opal520_Lp20 > 3, 
		DRD2=dat$MD_Opal620_LP10 > 3)
		
pdf("star_pdfs/DRD1_dotCount_lipoMasked_allROIs_hist.pdf",h=5,w=10)
par(mar= c(5,6,2,2), cex.axis=1.7, cex.lab=2, cex.main = 2)
hist(log2(dat$MD_Opal520_Lp20+1),xaxt = "n",
	xlab = "# Dots (Lipofuscin-Masked)",
	main = "DRD1 (Channel 520)")
axis(1,at=0:8,label = 2^(0:8)-1) 
abline(v=log2(3-1),col="red",lwd=3)
dev.off()

pdf("star_pdfs/DRD2_dotCount_lipoMasked_allROIs_hist.pdf",h=5,w=10)
par(mar= c(5,6,2,2), cex.axis=1.7, cex.lab=2, cex.main = 2)
hist(log2(dat$MD_Opal620_LP10+1),xaxt = "n",
	xlab = "# Dots (Lipofuscin-Masked)",
	main = "DRD1 (Channel 620)")
axis(1,at=0:8,label = 2^(0:8)-1) 
abline(v=log2(3-1),col="red",lwd=3)
dev.off()

## filter to >3 dots in either d1 or d2
drd_index = which(dat$MD_Opal520_Lp20 > 3 | dat$MD_Opal620_LP10 > 3)
dat_d = dat[drd_index,]
pheno_d = pheno[drd_index,]

######################
## other channels w/in drd1|2

### 570 
table(dat_d$MD_Opal570Lp1_0 > 6)

pdf("star_pdfs/HTR7_dotCount_lipoMasked_D1or2ROIs_hist.pdf",h=5,w=10)
par(mar= c(5,6,2,2), cex.axis=1.7, cex.lab=2, cex.main = 2)
hist(log2(dat_d$MD_Opal570Lp1_0+1),xaxt = "n",breaks=20,
	xlab = "# Dots (Lipofuscin-Masked)",
	main = "HTR7 (Channel 520)")
axis(1,at=0:8,label = 2^(0:8)-1) 
abline(v=log2(6-1),col="red",lwd=3)
dev.off()

### 690 
table(dat_d$MD_Opal690Lp30 > 6)

pdf("star_pdfs/CRHR21_dotCount_lipoMasked_D1or2ROIs_hist.pdf",h=5,w=10)
par(mar= c(5,6,2,2), cex.axis=1.7, cex.lab=2, cex.main = 2)
hist(log2(dat_d$MD_Opal690Lp30+1),breaks=20, xaxt = "n",
	xlab = "# Dots (Lipofuscin-Masked)",
	main = "CRHR21 (Channel 690)")
axis(1,at=0:8,label = 2^(0:8)-1) 
abline(v=log2(6-1),col="red",lwd=3)
dev.off()


dat_d$DRD1 = as.numeric(dat_d$MD_Opal520_Lp20  > 3)
dat_d$DRD2 = as.numeric(dat_d$MD_Opal620_LP10  > 3)
dat_d$HTR7 = as.numeric(dat_d$MD_Opal570Lp1_0 > 6)
dat_d$CRHR21 = as.numeric(dat_d$MD_Opal690Lp30 > 6)

###############################
## compare to reference #######
roi_mat = as.matrix(dat_d[,c( "DRD1","HTR7", "DRD2", "CRHR21")])
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
d_match$dist = signif(d_match$dist, 2)

table(d_match$dist)
    #   0   1 1.4
    # 157 324   1



mean(d_match$dist==0)
    # [1] 0.3257261

table(d_match$cell_name, d_match$dist)
    #                0   1 1.4
    # D1_A/B/C/D/F  61  15   1
    # D1_E          70 144   0
    # D2_C          12 165   0
    # D2_A/B/D      14   0   0




prop.table(table(d_match$cell_name, d_match$dist),2)
    #                       0          1        1.4
    # D1_A/B/C/D/F 0.38853503 0.04629630 1.00000000
    # D1_E         0.44585987 0.44444444 0.00000000
    # D2_C         0.07643312 0.50925926 0.00000000
    # D2_A/B/D     0.08917197 0.00000000 0.00000000

##############
## plots #####

pdf("star_pdfs/supp_figure_star_rnascope.pdf")
palette(brewer.pal(6,"Dark2"))
par(mar=c(5,6,2,2), cex.axis=1.5,cex.lab=1.5,cex.main=2)

plot(log2(MD_Opal570Lp1_0+1) ~ log2(MD_Opal690Lp30+1),
	data = dat_d,pch = 21+2*DRD1 , bg = DRD2 + 1,
	xlab = "CRHR21 (log2 [masked dot count + 1])",
	ylab = "HTR7 (log2 [masked dot count + 1])",
	cex=2-d_match$dist, main = "Star")
abline(h=log2(7.5), v = log2(7.5),lwd=2,lty=2) 
legend("topleft", c("DRD2-", "DRD2+"), col = 1:2, pch = 15,cex=1.5)
legend("top", c("DRD1-", "DRD1+"), pt.bg = "black", pch = c(21,23),cex=1.5)

text(x=c(1,7,7,1), y = c(0.3,0.5, 6.5 , 0.7), 
	c("D1.1/2/4", "D1.3", "D2.1", "D2.2"),
	col = c(1,2,1,2), cex=2)
dev.off()
	
	
		
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
rownames(MDcounts) <- c("rnascope_DRD1", "rnascope_HTR7",
                        "rnascope_DRD2","rnascope_CRHR2")

# check
table(colnames(MDcounts) == rownames(d_match))
    # all 482 TRUE
colnames(d_match)[colnames(d_match)=="cell_name"] <- "class_predict"

# Create SCE
sce.star <- SingleCellExperiment(list(counts=MDcounts), colData=d_match)
# a couple normalizations:
assay(sce.star, "logcounts") <- log2(assay(sce.star, "counts") + 1)
assay(sce.star, "logcounts.RVnorm") <- t(apply(assay(sce.star, "counts"),
                                                 1,function(x){
                                                   log2( x/dat_d$RVolume*10000 + 1 )
                                                 }))


pdf("star_pdfs/quantified-dotsPerROI_all-star-probes_betterVlnPlots_MNT2021.pdf", width=2)
# Log2-transform, only
plotExpressionCustom(sce.star, exprs_values="logcounts", scales="free_y",
                     features=c("rnascope_DRD1", "rnascope_DRD2",
                                "rnascope_CRHR2", "rnascope_HTR7"),
                     anno_name="class_predict", point_alpha=0.8, point_size=1.5, ncol=1,
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
plotExpressionCustom(sce.star, exprs_values="logcounts.RVnorm", scales="free_y",
                     features=c("rnascope_DRD1", "rnascope_DRD2",
                                "rnascope_CRHR2", "rnascope_HTR7"),
                     anno_name="class_predict", point_alpha=0.8, point_size=1.5, ncol=1,
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
	