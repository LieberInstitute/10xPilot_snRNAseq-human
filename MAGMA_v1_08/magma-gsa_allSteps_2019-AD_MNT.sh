#!/bin/bash
#$ -cwd
#$ -N magma-gsa_step3-gsa_AD
#$ -o ./logs/magma-gsa_step3-gsa_AD_MNT18Jul2021.o
#$ -e ./logs/magma-gsa_step3-gsa_AD_MNT18Jul2021.e
#$ -l bluejay,mem_free=16G,h_vmem=20G

echo "**** Job starts ****"
date

model="snp-wise"
ANNO=/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/GRCh38-ensembl93_to_hg19-lifted_30k-expressing-GENES.gene.loc
MAGMA=/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA_v1_08/magma
BFILE=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur
PREVPATH=/dcl02/lieber/ajaffe/Nick_Clifton/magma
setcol=1
genecol=2

gs_dlpfc=/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/dlpfcMarkerSets_fdr1e-6.txt
gs_sacc=/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/saccMarkerSets_fdr1e-6.txt
gs_hpc=/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/hpcMarkerSets_fdr1e-6.txt
gs_nac=/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/nacMarkerSets_fdr1e-6.txt
gs_amy=/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/amyMarkerSets_fdr1e-6.txt

SUMMSTATS=/dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/GWAS_Results/AD_sumstats_Jansenetal_2019sept.txt

## Step 1 - Annotation (SNP : gene mapping)
# $MAGMA --annotate window=35,10 --snp-loc ../MAGMA/GWAS_Results/AD_PGC-IGAP-ADSP-UKB_2019.snploc --gene-loc $ANNO --out SNP_Data/AD_Jansen2019_10xPilotGenes


## Step 2 - Gene analysis (from SNP-wise summary stats)
# $MAGMA --bfile $BFILE --gene-annot SNP_Data/AD_Jansen2019_10xPilotGenes.genes.annot --pval $SUMMSTATS use=SNP,P ncol=Neff --gene-model ${model} --out SNP_Data/AD_Jansen2019_10xPilotGenes_${model}


## Step 3 - Gene set analyses (using gene-level output)  -  for five brain regions

$MAGMA --gene-results SNP_Data/AD_Jansen2019_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results_rev/dlpfc_AD
# $MAGMA --gene-results SNP_Data/AD_Jansen2019_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results_rev/sacc_AD
$MAGMA --gene-results SNP_Data/AD_Jansen2019_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results_rev/hpc_AD
$MAGMA --gene-results SNP_Data/AD_Jansen2019_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results_rev/nac_AD
# $MAGMA --gene-results SNP_Data/AD_Jansen2019_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results_rev/amy_AD


echo "**** Job ends ****"
date
