### for MAGMA with LIBD 10x pilot analyses
# UPDATE: re-running with new v1.08
# MNT 09:10Sep2020 ======================

sumStats.ptsd <- read.table("/dcl02/lieber/ptsd/GWAS/ncbi-23818/files/release/submission/sub20190709/MVP.TE.REX.txt.gz",
                            header=T)
class(sumStats.ptsd)
dim(sumStats.ptsd)
head(sumStats.ptsd)
unique(sumStats.ptsd$CHR) # same
snploc.ptsd <- sumStats.ptsd[ ,c(1:3)]
colnames(snploc.ptsd) <- c("SNP", "CHR", "BP")
# Just write into where the others are stored (in the initial MAGMA/ dir, with v1.07)
write.table(snploc.ptsd, file="../MAGMA/GWAS_Results/PTSD_Nievergelt2019.snploc",
            sep="\t", col.names=T, row.names=F, quote=F)
rm(list=ls(pattern=".ptsd"))

