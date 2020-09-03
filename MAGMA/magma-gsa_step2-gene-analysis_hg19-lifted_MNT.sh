#!/bin/bash
#$ -cwd
#$ -N magma-gsa_step2-gene_MNT
#$ -o ./logs/magma-gsa_step2-gene_MNT23Aug2020.o
#$ -e ./logs/magma-gsa_step2-gene_MNT23Aug2020.e
#$ -l bluejay,mem_free=32G,h_vmem=40G


echo "**** Job starts ****"
date

model="snp-wise"
ANNO=/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/GRCh38-ensembl93_to_hg19-lifted_30k-expressing-GENES.gene.loc
MAGMA=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/MAGMA/magma
BFILE=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur
PREVPATH=/dcl02/lieber/ajaffe/Nick_Clifton/magma
setcol=1
genecol=2

mkdir -p SNP_Data
mkdir -p Results


### Step 2 - Gene analysis (from SNP-wise summary stats)

# pgc clozuk2 schizophrenia
$MAGMA --bfile $BFILE --gene-annot SNP_Data/clozuk_pgc2_pardinas2018_ensembl_10xPilotGenes.genes.annot --pval $PREVPATH/SNP_Data/clozuk_pgc2_pardinas2018.meta.sumstats.txt use=SNP,P N=35802 --gene-model ${model} --out SNP_Data/clozuk_pgc2_ensembl_10xPilotGenes_${model}

# PGC3 schizophrenia
$MAGMA --bfile $BFILE --gene-annot SNP_Data/pgc3_scz_ensembl_10xPilotGenes.genes.annot --pval $PREVPATH/SNP_Data/daner_PGC_SCZ_w3_90_0418bN use=SNP,P ncol=Nsum --gene-model ${model} --out SNP_Data/pgc3_scz_ensembl_10xPilotGenes_${model}

# bipolar disorder
$MAGMA --bfile $BFILE --gene-annot SNP_Data/daner_PGC_BIP32b_mds7a_0416a_ensembl_10xPilotGenes.genes.annot --pval $PREVPATH/SNP_Data/daner_PGC_BIP32b_mds7a_0416aN use=SNP,P ncol=Nsum --gene-model ${model} --out SNP_Data/daner_PGC_BIP32b_mds7a_0416a_ensembl_10xPilotGenes_${model}

# depression
$MAGMA --bfile $BFILE --gene-annot SNP_Data/MDD29_23andMe_Meta_Analysed1_ensembl_10xPilotGenes.genes.annot --pval $PREVPATH/SNP_Data/MDD29_23andMe_Meta_Analysed1_CHR_POS.metal use=MarkerName,P.value N=480359 --gene-model ${model} --out SNP_Data/MDD29_23andMe_ensembl_10xPilotGenes_${model}

# autism
$MAGMA --bfile $BFILE --gene-annot SNP_Data/iPSYCH-PGC_ASD_Nov2017_ensembl_10xPilotGenes.genes.annot --pval $PREVPATH/SNP_Data/iPSYCH-PGC_ASD_Nov2017 use=SNP,P N=46350 --gene-model ${model} --out SNP_Data/PGC_ASD_ensembl_10xPilotGenes_${model}


## GWAS for alcohol/tobacco use
$MAGMA --bfile $BFILE --gene-annot SNP_Data/Liu-etal_Addiction_AgeofInitiation_10xPilotGenes.genes.annot --pval GWAS_Results/AgeofInitiation.txt use=RSID,PVALUE ncol=N --gene-model ${model} --out SNP_Data/Liu-etal_Addiction_AgeofInitiation_10xPilotGenes_${model}

$MAGMA --bfile $BFILE --gene-annot SNP_Data/Liu-etal_Addiction_CigarettesPerDay_10xPilotGenes.genes.annot --pval GWAS_Results/CigarettesPerDay.txt use=RSID,PVALUE ncol=N --gene-model ${model} --out SNP_Data/Liu-etal_Addiction_CigarettesPerDay_10xPilotGenes_${model}

$MAGMA --bfile $BFILE --gene-annot SNP_Data/Liu-etal_Addiction_DrinksPerWeek_10xPilotGenes.genes.annot --pval GWAS_Results/DrinksPerWeek.txt use=RSID,PVALUE ncol=N --gene-model ${model} --out SNP_Data/Liu-etal_Addiction_DrinksPerWeek_10xPilotGenes_${model}

$MAGMA --bfile $BFILE --gene-annot SNP_Data/Liu-etal_Addiction_SmokingCessation_10xPilotGenes.genes.annot --pval GWAS_Results/SmokingCessation.txt use=RSID,PVALUE ncol=N --gene-model ${model} --out SNP_Data/Liu-etal_Addiction_SmokingCessation_10xPilotGenes_${model}

$MAGMA --bfile $BFILE --gene-annot SNP_Data/Liu-etal_Addiction_SmokingInitiation_10xPilotGenes.genes.annot --pval GWAS_Results/SmokingInitiation.txt use=RSID,PVALUE ncol=N --gene-model ${model} --out SNP_Data/Liu-etal_Addiction_SmokingInitiation_10xPilotGenes_${model}



echo "**** Job ends ****"
date
