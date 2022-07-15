#!/bin/bash
#$ -cwd
#$ -N magma-gsa_step3-geneSet_MNT
#$ -o ./logs/magma-gsa_step3-geneSet_MNT24Aug2020.o
#$ -e ./logs/magma-gsa_step3-geneSet_MNT24Aug2020.e
#$ -l bluejay,mem_free=32G,h_vmem=40G

echo "**** Job starts ****"
date

model="snp-wise"

ANNO=/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/GRCh38-ensembl93_to_hg19-lifted_30k-expressing-GENES.gene.loc
MAGMA=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/MAGMA/magma
BFILE=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur
PREVPATH=/dcl02/lieber/ajaffe/Nick_Clifton/magma

setcol=1
genecol=2

gs_dlpfc=/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/dlpfcMarkerSets_fdr1e-12.txt
gs_sacc=/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/saccMarkerSets_fdr1e-12.txt
gs_hpc=/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/hpcMarkerSets_fdr1e-12.txt
gs_nac=/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/nacMarkerSets_fdr1e-12.txt
gs_amy=/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/amyMarkerSets_fdr1e-12.txt


mkdir -p Results

### Step 3 - Gene set analyses (using gene-level output)

## DLPFC ==================
# pgc clozuk2 schizophrenia
$MAGMA --gene-results SNP_Data/clozuk_pgc2_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results/dlpfc_clozuk_pgc2

# PGC3 schizophrenia
$MAGMA --gene-results SNP_Data/pgc3_scz_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results/dlpfc_pgc3_scz

# bipolar disorder
$MAGMA --gene-results SNP_Data/daner_PGC_BIP32b_mds7a_0416a_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results/dlpfc_PGC_BIP

# depression
$MAGMA --gene-results SNP_Data/MDD29_23andMe_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results/dlpfc_MDD29_23andMe

# autism
$MAGMA --gene-results SNP_Data/PGC_ASD_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results/dlpfc_PGC_ASD

## Addiction GWAS set
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_AgeofInitiation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results/dlpfc_AgeSmk
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_CigarettesPerDay_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results/dlpfc_CigDay
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_DrinksPerWeek_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results/dlpfc_DrnkWk
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_SmokingCessation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results/dlpfc_SmkCes
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_SmokingInitiation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results/dlpfc_SmkInit


echo "** DLPFC MAGMA gene sets analyses complete **"



## sACC ==================
# pgc clozuk2 schizophrenia
$MAGMA --gene-results SNP_Data/clozuk_pgc2_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results/sacc_clozuk_pgc2

# PGC3 schizophrenia
$MAGMA --gene-results SNP_Data/pgc3_scz_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results/sacc_pgc3_scz

# bipolar disorder
$MAGMA --gene-results SNP_Data/daner_PGC_BIP32b_mds7a_0416a_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results/sacc_PGC_BIP

# depression
$MAGMA --gene-results SNP_Data/MDD29_23andMe_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results/sacc_MDD29_23andMe

# autism
$MAGMA --gene-results SNP_Data/PGC_ASD_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results/sacc_PGC_ASD

## Addiction GWAS set
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_AgeofInitiation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results/sacc_AgeSmk
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_CigarettesPerDay_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results/sacc_CigDay
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_DrinksPerWeek_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results/sacc_DrnkWk
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_SmokingCessation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results/sacc_SmkCes
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_SmokingInitiation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results/sacc_SmkInit

echo "** sACC MAGMA gene sets analyses complete **"



## HPC ==================
# pgc clozuk2 schizophrenia
$MAGMA --gene-results SNP_Data/clozuk_pgc2_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results/hpc_clozuk_pgc2

# PGC3 schizophrenia
$MAGMA --gene-results SNP_Data/pgc3_scz_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results/hpc_pgc3_scz

# bipolar disorder
$MAGMA --gene-results SNP_Data/daner_PGC_BIP32b_mds7a_0416a_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results/hpc_PGC_BIP

# depression
$MAGMA --gene-results SNP_Data/MDD29_23andMe_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results/hpc_MDD29_23andMe

# autism
$MAGMA --gene-results SNP_Data/PGC_ASD_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results/hpc_PGC_ASD

## Addiction GWAS set
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_AgeofInitiation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results/hpc_AgeSmk
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_CigarettesPerDay_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results/hpc_CigDay
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_DrinksPerWeek_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results/hpc_DrnkWk
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_SmokingCessation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results/hpc_SmkCes
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_SmokingInitiation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results/hpc_SmkInit

echo "** HPC MAGMA gene sets analyses complete **"



## NAc ==================
# pgc clozuk2 schizophrenia
$MAGMA --gene-results SNP_Data/clozuk_pgc2_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results/nac_clozuk_pgc2

# PGC3 schizophrenia
$MAGMA --gene-results SNP_Data/pgc3_scz_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results/nac_pgc3_scz

# bipolar disorder
$MAGMA --gene-results SNP_Data/daner_PGC_BIP32b_mds7a_0416a_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results/nac_PGC_BIP

# depression
$MAGMA --gene-results SNP_Data/MDD29_23andMe_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results/nac_MDD29_23andMe

# autism
$MAGMA --gene-results SNP_Data/PGC_ASD_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results/nac_PGC_ASD

## Addiction GWAS set
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_AgeofInitiation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results/nac_AgeSmk
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_CigarettesPerDay_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results/nac_CigDay
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_DrinksPerWeek_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results/nac_DrnkWk
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_SmokingCessation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results/nac_SmkCes
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_SmokingInitiation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results/nac_SmkInit

echo "** NAc MAGMA gene sets analyses complete **"



## AMY ==================
# pgc clozuk2 schizophrenia
$MAGMA --gene-results SNP_Data/clozuk_pgc2_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results/amy_clozuk_pgc2

# PGC3 schizophrenia
$MAGMA --gene-results SNP_Data/pgc3_scz_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results/amy_pgc3_scz

# bipolar disorder
$MAGMA --gene-results SNP_Data/daner_PGC_BIP32b_mds7a_0416a_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results/amy_PGC_BIP

# depression
$MAGMA --gene-results SNP_Data/MDD29_23andMe_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results/amy_MDD29_23andMe

# autism
$MAGMA --gene-results SNP_Data/PGC_ASD_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results/amy_PGC_ASD

## Addiction GWAS set
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_AgeofInitiation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results/amy_AgeSmk
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_CigarettesPerDay_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results/amy_CigDay
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_DrinksPerWeek_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results/amy_DrnkWk
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_SmokingCessation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results/amy_SmkCes
$MAGMA --gene-results SNP_Data/Liu-etal_Addiction_SmokingInitiation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results/amy_SmkInit

echo "** AMY MAGMA gene sets analyses complete **"


echo "**** Job ends ****"
date
