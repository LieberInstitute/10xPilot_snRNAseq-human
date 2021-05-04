#!/bin/bash
#$ -cwd
#$ -N magma-gsa_step3-GSA_v2
#$ -o ./logs/magma-gsa_step3-GSA_v2_MNT03May2021.o
#$ -e ./logs/magma-gsa_step3-GSA_v2_MNT03May2021.e
#$ -l bluejay,mem_free=32G,h_vmem=40G

echo "**** Job starts ****"
date

model="snp-wise"

ANNO=/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/GRCh38-ensembl93_to_hg19-lifted_30k-expressing-GENES.gene.loc
MAGMA=/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA_v1_08/magma
BFILE=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur
PREVPATH=/dcl02/lieber/ajaffe/Nick_Clifton/magma

setcol=1
genecol=2

gs_dlpfc=/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/dlpfcMarkerSets_fdr1e-12.txt
gs_sacc=/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/saccMarkerSets_fdr1e-12.txt
gs_hpc=/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/hpcMarkerSets_fdr1e-12.txt
gs_nac=/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/nacMarkerSets_fdr1e-12.txt
gs_amy=/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/amyMarkerSets_fdr1e-12.txt


echo "This update (rolled out Sep2020) uses a new test statistic and approach that better handles increased SNP density and LD, for Step 2 (gene analysis). Step 3 approach remains the same."
echo "Update MNT 03May2021 - running an iteration (into 'Results_v2') with the preprint stats, applying a filter for non-0-median expression, in each respective subcluster."

### Step 3 - Gene set analyses (using gene-level output)

## DLPFC ==================
# pgc clozuk2 schizophrenia
$MAGMA --gene-results SNP_Data_v2/clozuk_pgc2_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results_v2/dlpfc_clozuk_pgc2

# PGC3 schizophrenia
$MAGMA --gene-results SNP_Data_v2/pgc3_scz_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results_v2/dlpfc_pgc3_scz

# bipolar disorder
$MAGMA --gene-results SNP_Data_v2/daner_PGC_BIP32b_mds7a_0416a_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results_v2/dlpfc_PGC_BIP

# depression
$MAGMA --gene-results SNP_Data_v2/MDD29_23andMe_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results_v2/dlpfc_MDD29_23andMe

# autism
$MAGMA --gene-results SNP_Data_v2/PGC_ASD_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results_v2/dlpfc_PGC_ASD

## Addiction GWAS set
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_AgeofInitiation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results_v2/dlpfc_AgeSmk
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_CigarettesPerDay_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results_v2/dlpfc_CigDay
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_DrinksPerWeek_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results_v2/dlpfc_DrnkWk
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_SmokingCessation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results_v2/dlpfc_SmkCes
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_SmokingInitiation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results_v2/dlpfc_SmkInit


echo "** DLPFC MAGMA gene sets analyses complete **"



## sACC ==================
# pgc clozuk2 schizophrenia
$MAGMA --gene-results SNP_Data_v2/clozuk_pgc2_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results_v2/sacc_clozuk_pgc2

# PGC3 schizophrenia
$MAGMA --gene-results SNP_Data_v2/pgc3_scz_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results_v2/sacc_pgc3_scz

# bipolar disorder
$MAGMA --gene-results SNP_Data_v2/daner_PGC_BIP32b_mds7a_0416a_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results_v2/sacc_PGC_BIP

# depression
$MAGMA --gene-results SNP_Data_v2/MDD29_23andMe_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results_v2/sacc_MDD29_23andMe

# autism
$MAGMA --gene-results SNP_Data_v2/PGC_ASD_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results_v2/sacc_PGC_ASD

## Addiction GWAS set
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_AgeofInitiation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results_v2/sacc_AgeSmk
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_CigarettesPerDay_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results_v2/sacc_CigDay
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_DrinksPerWeek_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results_v2/sacc_DrnkWk
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_SmokingCessation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results_v2/sacc_SmkCes
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_SmokingInitiation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results_v2/sacc_SmkInit

echo "** sACC MAGMA gene sets analyses complete **"



## HPC ==================
# pgc clozuk2 schizophrenia
$MAGMA --gene-results SNP_Data_v2/clozuk_pgc2_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results_v2/hpc_clozuk_pgc2

# PGC3 schizophrenia
$MAGMA --gene-results SNP_Data_v2/pgc3_scz_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results_v2/hpc_pgc3_scz

# bipolar disorder
$MAGMA --gene-results SNP_Data_v2/daner_PGC_BIP32b_mds7a_0416a_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results_v2/hpc_PGC_BIP

# depression
$MAGMA --gene-results SNP_Data_v2/MDD29_23andMe_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results_v2/hpc_MDD29_23andMe

# autism
$MAGMA --gene-results SNP_Data_v2/PGC_ASD_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results_v2/hpc_PGC_ASD

## Addiction GWAS set
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_AgeofInitiation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results_v2/hpc_AgeSmk
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_CigarettesPerDay_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results_v2/hpc_CigDay
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_DrinksPerWeek_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results_v2/hpc_DrnkWk
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_SmokingCessation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results_v2/hpc_SmkCes
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_SmokingInitiation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results_v2/hpc_SmkInit

echo "** HPC MAGMA gene sets analyses complete **"



## NAc ==================
# pgc clozuk2 schizophrenia
$MAGMA --gene-results SNP_Data_v2/clozuk_pgc2_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results_v2/nac_clozuk_pgc2

# PGC3 schizophrenia
$MAGMA --gene-results SNP_Data_v2/pgc3_scz_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results_v2/nac_pgc3_scz

# bipolar disorder
$MAGMA --gene-results SNP_Data_v2/daner_PGC_BIP32b_mds7a_0416a_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results_v2/nac_PGC_BIP

# depression
$MAGMA --gene-results SNP_Data_v2/MDD29_23andMe_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results_v2/nac_MDD29_23andMe

# autism
$MAGMA --gene-results SNP_Data_v2/PGC_ASD_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results_v2/nac_PGC_ASD

## Addiction GWAS set
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_AgeofInitiation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results_v2/nac_AgeSmk
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_CigarettesPerDay_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results_v2/nac_CigDay
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_DrinksPerWeek_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results_v2/nac_DrnkWk
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_SmokingCessation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results_v2/nac_SmkCes
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_SmokingInitiation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results_v2/nac_SmkInit

echo "** NAc MAGMA gene sets analyses complete **"



## AMY ==================
# pgc clozuk2 schizophrenia
$MAGMA --gene-results SNP_Data_v2/clozuk_pgc2_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results_v2/amy_clozuk_pgc2

# PGC3 schizophrenia
$MAGMA --gene-results SNP_Data_v2/pgc3_scz_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results_v2/amy_pgc3_scz

# bipolar disorder
$MAGMA --gene-results SNP_Data_v2/daner_PGC_BIP32b_mds7a_0416a_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results_v2/amy_PGC_BIP

# depression
$MAGMA --gene-results SNP_Data_v2/MDD29_23andMe_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results_v2/amy_MDD29_23andMe

# autism
$MAGMA --gene-results SNP_Data_v2/PGC_ASD_ensembl_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results_v2/amy_PGC_ASD

## Addiction GWAS set
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_AgeofInitiation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results_v2/amy_AgeSmk
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_CigarettesPerDay_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results_v2/amy_CigDay
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_DrinksPerWeek_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results_v2/amy_DrnkWk
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_SmokingCessation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results_v2/amy_SmkCes
$MAGMA --gene-results SNP_Data_v2/Liu-etal_Addiction_SmokingInitiation_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results_v2/amy_SmkInit

echo "** AMY MAGMA gene sets analyses complete **"


## 'Negative control': Meta-GWAS for CAD (CARDIoGRAM+C4D) === ===
$MAGMA --gene-results SNP_Data_v2/CoronaryArteryDisease_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results_v2/dlpfc_CAD
$MAGMA --gene-results SNP_Data_v2/CoronaryArteryDisease_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results_v2/sacc_CAD
$MAGMA --gene-results SNP_Data_v2/CoronaryArteryDisease_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results_v2/hpc_CAD
$MAGMA --gene-results SNP_Data_v2/CoronaryArteryDisease_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results_v2/nac_CAD
$MAGMA --gene-results SNP_Data_v2/CoronaryArteryDisease_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results_v2/amy_CAD


echo "**** Job ends ****"
date
