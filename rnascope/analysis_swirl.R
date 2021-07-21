####
library(readxl)
library(jaffelab)
library(readr)
library(pracma)
library(RColorBrewer)

options(width=100)

dir.create("swirl_pdfs")

## read in reference data
ref = read_excel("raw_data/Swirl_expected_expression_revisionMNT.xlsx")
ref = as.data.frame(ref[,1:5])
ref_mat = as.matrix(ref[,2:5])
rownames(ref_mat) = ref$Population

## read in long data
dat =read.csv("raw_data/long_data_Swirl.csv",as.is=TRUE)
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


######################
## all channels w/o filtering


### 520 
table(dat$MD_Opal520_Lp20 > 6)

pdf("swirl_pdfs/GAD1_dotCount_lipoMasked_allROIs_hist.pdf",h=5,w=10)
par(mar= c(5,6,2,2), cex.axis=1.7, cex.lab=2, cex.main = 2)
hist(log2(dat$MD_Opal520_Lp20+1),xaxt = "n",breaks=20,
	xlab = "# Dots (Lipofuscin-Masked)",
	main = "GAD1 (Channel 520)")
axis(1,at=0:8,label = 2^(0:8)-1) 
abline(v=log2(6-1),col="red",lwd=3)
dev.off()


### filter to inhib
inhib_index = which(dat$MD_Opal520_Lp20 > 6)
dat_inhib = dat[inhib_index,]
pheno_inhib = pheno[inhib_index,]

######################
## other channels w/in inhib
#####################

### 570 

table(dat_inhib$MD_Opal570Lp1_0 > 6)

pdf("swirl_pdfs/PVALB_dotCount_lipoMasked_inhibROIs_hist.pdf",h=5,w=10)
par(mar= c(5,6,2,2), cex.axis=1.7, cex.lab=2, cex.main = 2)
hist(log2(dat_inhib$MD_Opal570Lp1_0+1),breaks=20, xaxt = "n",
	xlab = "# Dots (Lipofuscin-Masked)",
	main = "PVALB (Channel 570)")
axis(1,at=0:8,label = 2^(0:8)-1) 
abline(v=log2(6-1),col="red",lwd=3)
dev.off()


### 620 
table(dat_inhib$MD_Opal620_LP10 > 6)

pdf("swirl_pdfs/KIT_dotCount_lipoMasked_inhibROIs_hist.pdf",h=5,w=10)
par(mar= c(5,6,2,2), cex.axis=1.7, cex.lab=2, cex.main = 2)
hist(log2(dat_inhib$MD_Opal620_LP10+1),xaxt = "n",
	xlab = "# Dots (Lipofuscin-Masked)",breaks=20,
	main = "KIT (Channel 620)")
axis(1,at=0:8,label = 2^(0:8)-1) 
abline(v=log2(6-1),col="red",lwd=3)
dev.off()

### 690 
table(dat_inhib$MD_Opal690Lp30 > 6)
pdf("swirl_pdfs/PTHLH_dotCount_lipoMasked_inhibROIs_hist.pdf",h=5,w=10)
par(mar= c(5,6,2,2), cex.axis=1.7, cex.lab=2, cex.main = 2)
hist(log2(dat_inhib$MD_Opal690Lp30+1),breaks=20, xaxt = "n",
	xlab = "# Dots (Lipofuscin-Masked)",
	main = "PTHLH (Channel 690)")
axis(1,at=0:8,label = 2^(0:8)-1) 
abline(v=log2(6-1),col="red",lwd=3)
dev.off()

## get counts
dat_inhib$PVALB = as.numeric(dat_inhib$MD_Opal570Lp1_0 > 6)
dat_inhib$KIT = as.numeric(dat_inhib$MD_Opal620_LP10 > 6)
dat_inhib$PTHLH = as.numeric(dat_inhib$MD_Opal690Lp30 > 6)

###############################
## compare to reference #######
roi_mat = as.matrix(dat_inhib[,c("PVALB", "KIT", "PTHLH")])
colMeans(roi_mat)

ref_mat_binary = as.matrix(ref_mat[,-1])
ref_mat_binary[ref_mat_binary > 0] = 1
dist_mat = t(distmat(ref_mat_binary, roi_mat))

## label
dat_match = cbind(roi_mat, round(dist_mat,2), 
	data.frame(cell_type = apply(dist_mat,1, which.min),
		dist =apply(dist_mat,1, min)))
dat_match$cell_name = rownames(ref_mat)[dat_match$cell_type]
dat_match$cell_name = factor(dat_match$cell_name, rownames(ref_mat))
dat_match$dist = signif(dat_match$dist, 2) 
####
## metrics

table(dat_match$dist)
    #  0   1 
    # 76 136

mean(dat_match$dist==0)
    # [1] 0.3584906

table(dat_match$cell_name, dat_match$dist)
      #            0  1
      # Inhib_A    5 92
      # Inhib_B/E 21 32
      # Inhib_C/D 50 12

round(prop.table(table(dat_match$cell_name, dat_match$dist), 2),3)
    #               0     1
    # Inhib_A   0.066 0.676
    # Inhib_B/E 0.276 0.235
    # Inhib_C/D 0.658 0.088
  
## figure for paper
pdf("swirl_pdfs/supp_figure_swirl_rnascope.pdf")
palette(brewer.pal(6,"Dark2"))
par(mar=c(5,6,2,2), cex.axis=1.5,cex.lab=1.5,cex.main=2)

plot(log2(MD_Opal570Lp1_0+1) ~ log2(MD_Opal690Lp30+1),
	data = dat_inhib,pch = 21, bg = KIT + 1,
	xlab = "PTHLH (log2 [masked dot count + 1])",
	ylab = "PVALB (log2 [masked dot count + 1])",
	cex=2-dat_match$dist, main = "Swirl")
abline(h=log2(7.5), v = log2(7.5),lwd=2,lty=2) 
legend("topleft", c("KIT-", "KIT+"), col = 1:2, pch = 15,cex=1.5)

text(x=c(1,1,6,5), y = c(0.3,0.7, 0.5, 6.5), 
	c("Inhib1", "Inhib4", "Inhib2", "Inhib3"),
	col = c(1,2,1,2), cex=2)
dev.off()


## numbers for paper
nrow(pheno)
length(table(pheno$Image))
length(table(pheno$BrNum))
length(unique(paste0(pheno$Section,":", pheno$BrNum)))

dim(dat_inhib)

### MNT / ABS update May2021 =============
# objective: make more-interpretable graphics with RNAscope quantification data
#(Abby can add/save this to the end of Andrew's script - analysis_swirl.R:)
	# Some ideas ===
	### 1) boxplot log2-transformed (/maybe add RVolume normalization)
	## For PTHLH (Opal690)
  pdf("swirl_pdfs/PTHLH_swirl_rnascope.pdf")
	boxplot(log2(dat_inhib$MD_Opal690Lp30 + 1) ~ dat_match$cell_name)
	dev.off()
	
	# or (with $RVolume normalization)
	pdf("swirl_pdfs/PTHLH_norm_swirl_rnascope.pdf")
	boxplot(log2(dat_inhib$MD_Opal690Lp30/dat_inhib$RVolume + 1) ~ dat_match$cell_name)
	dev.off()
	
	## For KIT
	pdf("swirl_pdfs/KIT_swirl_rnascope.pdf")
	boxplot(log2(dat_inhib$MD_Opal620_LP10 + 1) ~ dat_match$cell_name)
	dev.off()
	
	# or (with $RVolume normalization)
	pdf("swirl_pdfs/KIT_norm_swirl_rnascope.pdf")
	boxplot(log2(dat_inhib$MD_Opal620_LP10/dat_inhib$RVolume + 1) ~ dat_match$cell_name)
	dev.off()
	
	## For PVALB
	pdf("swirl_pdfs/PVALB_swirl_rnascope.pdf")
	boxplot(log2(dat_inhib$MD_Opal570Lp1_0 + 1) ~ dat_match$cell_name)
	dev.off()
	
	# or (with $RVolume normalization)
	pdf("swirl_pdfs/PVALB_norm_swirl_rnascope.pdf")
	boxplot(log2(dat_inhib$MD_Opal570Lp1_0/dat_inhib$RVolume + 1) ~ dat_match$cell_name)
	dev.off()
	
	## For GAD1
	pdf("swirl_pdfs/GAD1_swirl_rnascope.pdf")
	boxplot(log2(dat_inhib$MD_Opal520_Lp20  + 1) ~ dat_match$cell_name)
	dev.off()
	
	# or (with $RVolume normalization)
	pdf("swirl_pdfs/GAD1_norm_swirl_rnascope.pdf")
	boxplot(log2(dat_inhib$MD_Opal520_Lp20/dat_inhib$RVolume + 1) ~ dat_match$cell_name)
	dev.off()
	
	### 2) Hierarchically cluster the 212 (GAD1+) ROIs on normalized (i.e. `/RVolume`)
		#    n transcript dots? (so 212 x 3 probes as input matrix)
		#    (MD_* : total transcript dots in that ROI after Lipofuscin masking)
		# - hopefully would see four (for four inhib. subpops) main branches??
	
	## Forcing to use like SCE data, for aesthetics ===
	library(scater)
	library(SingleCellExperiment)
	source('../plotExpressionCustom.R')
	tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
	              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
	              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
	              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")
	blgngy <- tableau20[c(1:2, 19:20, 5:6, 17:18, 15:16)]
	
	MDcounts <- t(dat_inhib[ ,colnames(dat_inhib)[grep('^MD', colnames(dat_inhib))]])
	rownames(MDcounts)
	    #[1] "MD_Opal520_Lp20" "MD_Opal570Lp1_0" "MD_Opal620_LP10" "MD_Opal690Lp30" 
	    # (check these with 'ref_mat' - good)
	rownames(MDcounts) <- c("rnascope_GAD1", "rnascope_PVALB",
	                        "rnascope_KIT","rnascope_PTHLH")
	# check
	table(colnames(MDcounts) == rownames(dat_match))

	colnames(dat_match)[colnames(dat_match)=="cell_name"] <- "class_predict"
	
	# Create SCE
	sce.swirl <- SingleCellExperiment(list(counts=MDcounts), colData=dat_match)
	# a couple normalizations:
	assay(sce.swirl, "logcounts") <- log2(assay(sce.swirl, "counts") + 1)
	assay(sce.swirl, "logcounts.RVnorm") <- t(apply(assay(sce.swirl, "counts"),
	                                                 1,function(x){
	                                                   log2( x/dat_inhib$RVolume*10000 + 1 )
	                                                 }))
	
	
	pdf("swirl_pdfs/quantified-dotsPerROI_all-swirl-probes_betterVlnPlots_MNT2021.pdf", width=2)
	# Log2-transform, only
	plotExpressionCustom(sce.swirl, exprs_values="logcounts", scales="free_y",
	                     features=c("rnascope_GAD1", "rnascope_KIT",
	                                "rnascope_PTHLH","rnascope_PVALB"),
	                     anno_name="class_predict", point_alpha=0.7, point_size=1.5, ncol=1,
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
	plotExpressionCustom(sce.swirl, exprs_values="logcounts.RVnorm", scales="free_y",
	                     features=c("rnascope_GAD1", "rnascope_KIT",
	                                "rnascope_PTHLH","rnascope_PVALB"),
	                     anno_name="class_predict", point_alpha=0.7, point_size=1.5, ncol=1,
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
	
	
	