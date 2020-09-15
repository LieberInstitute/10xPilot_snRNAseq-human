#!/bin/bash
#$ -cwd
#$ -N magma-gsa_all_PTSD_updated_MNT
#$ -o ./logs/magma-gsa_all_PTSD_v1.08_MNT10Sep2020.o
#$ -e ./logs/magma-gsa_all_PTSD_v1.08_MNT10Sep2020.e
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


## Step 1 - Annotation (SNP : gene mapping)
$MAGMA --annotate window=35,10 --snp-loc ../MAGMA/GWAS_Results/PTSD_Nievergelt2019.snploc --gene-loc $ANNO --out SNP_Data/PTSD_Nievergelt2019_10xPilotGenes

## Step 2 - Gene analysis (from SNP-wise summary stats)
$MAGMA --bfile $BFILE --gene-annot SNP_Data/PTSD_Nievergelt2019_10xPilotGenes.genes.annot --pval ../MAGMA/GWAS_Results/MVP.TE.REX.txt use=SNP,P ncol=N --gene-model ${model} --out SNP_Data/PTSD_Nievergelt2019_10xPilotGenes_${model}

## Step 3 - Gene set analyses (using gene-level output)  -  for five brain regions
$MAGMA --gene-results SNP_Data/PTSD_Nievergelt2019_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results/dlpfc_PTSD
$MAGMA --gene-results SNP_Data/PTSD_Nievergelt2019_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results/sacc_PTSD
$MAGMA --gene-results SNP_Data/PTSD_Nievergelt2019_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results/hpc_PTSD
$MAGMA --gene-results SNP_Data/PTSD_Nievergelt2019_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results/nac_PTSD
$MAGMA --gene-results SNP_Data/PTSD_Nievergelt2019_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results/amy_PTSD


echo "**** Job ends ****"
date
