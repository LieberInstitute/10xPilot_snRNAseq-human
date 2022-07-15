#!/bin/bash
#$ -cwd
#$ -N magma-gsa_step1-annot_MNT
#$ -o ./logs/magma-gsa_step1-annot_MNT23Aug2020.o
#$ -e ./logs/magma-gsa_step1-annot_MNT23Aug2020.e
#$ -l bluejay,mem_free=16G,h_vmem=20G

echo "**** Job starts ****"
date

model="snp-wise"
ANNO=/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/GRCh38-ensembl93_to_hg19-lifted_30k-expressing-GENES.gene.loc
MAGMA=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/MAGMA/magma
BFILE=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur
PREVPATH=/dcl02/lieber/ajaffe/Nick_Clifton/magma
setcol=1
genecol=2


mkdir -p SNP_Data
mkdir -p Results


### Step 1 - Annotation (SNP : gene mapping)

# pgc clozuk2 schizophrenia
$MAGMA --annotate window=35,10 --snp-loc $PREVPATH/SNP_Data/clozuk_pgc2_pardinas2018.snploc --gene-loc $ANNO --out SNP_Data/clozuk_pgc2_pardinas2018_ensembl_10xPilotGenes

# PGC3 schizophrenia
$MAGMA --annotate window=35,10 --snp-loc $PREVPATH/SNP_Data/daner_PGC_SCZ_w3_90_0418b_INF06.snploc --gene-loc $ANNO --out SNP_Data/pgc3_scz_ensembl_10xPilotGenes

# bipolar disorder
$MAGMA --annotate window=35,10 --snp-loc $PREVPATH/SNP_Data/daner_PGC_BIP32b_mds7a_0416a_INF06.snploc --gene-loc $ANNO --out SNP_Data/daner_PGC_BIP32b_mds7a_0416a_ensembl_10xPilotGenes

# depression
$MAGMA --annotate window=35,10 --snp-loc $PREVPATH/SNP_Data/MDD29_23andMe_Meta_Analysed1_CHR_POS.snploc --gene-loc $ANNO --out SNP_Data/MDD29_23andMe_Meta_Analysed1_ensembl_10xPilotGenes

# autism
$MAGMA --annotate window=35,10 --snp-loc $PREVPATH/SNP_Data/iPSYCH-PGC_ASD_Nov2017.snploc --gene-loc $ANNO --out SNP_Data/iPSYCH-PGC_ASD_Nov2017_ensembl_10xPilotGenes


## GWAS for alcohol/tobacco use
$MAGMA --annotate window=35,10 --snp-loc ./GWAS_Results/AgeofInitiation.snploc --gene-loc $ANNO --out SNP_Data/Liu-etal_Addiction_AgeofInitiation_10xPilotGenes

$MAGMA --annotate window=35,10 --snp-loc ./GWAS_Results/CigarettesPerDay.snploc --gene-loc $ANNO --out SNP_Data/Liu-etal_Addiction_CigarettesPerDay_10xPilotGenes

$MAGMA --annotate window=35,10 --snp-loc ./GWAS_Results/DrinksPerWeek.snploc --gene-loc $ANNO --out SNP_Data/Liu-etal_Addiction_DrinksPerWeek_10xPilotGenes

$MAGMA --annotate window=35,10 --snp-loc ./GWAS_Results/SmokingCessation.snploc --gene-loc $ANNO --out SNP_Data/Liu-etal_Addiction_SmokingCessation_10xPilotGenes

$MAGMA --annotate window=35,10 --snp-loc ./GWAS_Results/SmokingInitiation.snploc --gene-loc $ANNO --out SNP_Data/Liu-etal_Addiction_SmokingInitiation_10xPilotGenes


echo "**** Job ends ****"
date
