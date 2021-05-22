####
library(readxl)
library(jaffelab)
library(readr)
library(pracma)

## read in reference data
ref = read_csv("raw_data/square_MSN_exp_rnascope.csv")
ref = as.data.frame(ref[,1:5])
ref_mat = as.matrix(ref[,2:5])
rownames(ref_mat) = ref$Population
ref_mat = ref_mat[,c(1,2,4,3)]

## read in long data
dat =read.csv("raw_data/long_data.csv",as.is=TRUE)
colnames(dat)[1]= "Image"
dat$Image =gsub(" unmix", "Umix", dat$Image)

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

#######				
## make a few matrices
####

## proportions of rois covered by dots
dat$M_PP_Opal520_Lp20 = dat$MP_Opal520_Lp20 / dat$RVolume
dat$M_PP_Opal570Lp1_0 = dat$MP_Opal570Lp1_0 / dat$RVolume
dat$M_PP_Opal620_LP10 = dat$MP_Opal620_LP10 / dat$RVolume
dat$M_PP_Opal690Lp30 = dat$MP_Opal690Lp30 / dat$RVolume

## check for ROIs with high DRD1 expression
## channel 620
hist(dat$MI_P_Opal620_LP10/dat$RVolume)
hist(dat$MI_MP_Opal620_LP10)
table(dat$MD_Opal620_LP10 > 0)
table(dat$MD_Opal620_LP10 > 10)

## filter to any dots
dat_d1 = dat[which(dat$MD_Opal620_LP10 > 6),]
pheno_d1 = pheno[which(dat$MD_Opal620_LP10 > 6),]

# ## cluster ROI that express d1
# 
# pp_mat = as.matrix(dat_d1[,grep("^M_PP", colnames(dat_d1))])
# k = kmeans(pp_mat, 4)
# round(k$centers*100,2)
# table(k$cluster)

################
#### kmeans
set.seed(2131)
dat_mask = dat_d1[,grep("M", colnames(dat_d1))]
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

## relevel
mat_qual_k3 = sapply(kList, function(x) {
	m = x$centers
	new = factor(x$cluster, levels = order(m[,4]))
	as.numeric(new) - 1
})

## compare to reference
dist_mat = t(distmat(ref_mat[,-3], mat_qual_k3))
tt = table(round(apply(dist_mat,1, min),2))
tt

## and binary
ref_mat_binary = ref_mat
ref_mat_binary[ref_mat_binary > 0] = 1
mat_qual_k3_binary = mat_qual_k3
mat_qual_k3_binary[mat_qual_k3_binary > 0] = 1
dist_mat_binary = t(distmat(ref_mat_binary[,-3], 
	mat_qual_k3_binary))

## make tables
table(round(apply(dist_mat_binary,1, min),2))

## cross tab above
table(cat = round(apply(dist_mat,1, min),2),
	bin = round(apply(dist_mat_binary,1, min),2))


d1_match = data.frame(cell_type = apply(dist_mat_binary,1, which.min),
	dist =apply(dist_mat_binary,1, min))
d1_match$cell_name = rownames(ref_mat)[d1_match$cell_type]
table(d1_match$cell_name)
prop.table(table(d1_match$cell_name))
table(d1_match$cell_name[d1_match$dist == 0])
prop.table(table(d1_match$cell_name[d1_match$dist == 0]))

# ### conditional on DRD1
# d1_index_lo = mat_qual_k3[,"620_LP10"] >0
# d1_index_hi = mat_qual_k3[,"620_LP10"] >1
# table(round(apply(dist_mat[d1_index_lo,],1, min),2))
# table(round(apply(dist_mat_binary[d1_index_lo,],1, min),2))
# table(round(apply(dist_mat[d1_index_hi,],1, min),2))
