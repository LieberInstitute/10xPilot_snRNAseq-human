#!/bin/bash
#$ -cwd
#$ -N magma-gsa_all_CAD_MNT
#$ -o ./logs/magma-gsa_all_CAD_MNT02Sep2020.o
#$ -e ./logs/magma-gsa_all_CAD_MNT02Sep2020.e
#$ -l bluejay,mem_free=24G,h_vmem=32G

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


## Step 1 - Annotation (SNP : gene mapping)
# $MAGMA --annotate window=35,10 --snp-loc ./GWAS_Results/CoronaryArteryDisease.snploc --gene-loc $ANNO --out SNP_Data/CoronaryArteryDisease_2015_10xPilotGenes

## Step 2 - Gene analysis (from SNP-wise summary stats)
$MAGMA --bfile $BFILE --gene-annot SNP_Data/CoronaryArteryDisease_2015_10xPilotGenes.genes.annot --pval GWAS_Results/CoronaryArteryDisease_Nikpay2015_MNTaddn.txt use=markername,p_dgc ncol=N_estim_mnt --gene-model ${model} --out SNP_Data/CoronaryArteryDisease_10xPilotGenes_${model}

## Step 3 - Gene set analyses (using gene-level output)  -  for five brain regions
$MAGMA --gene-results SNP_Data/CoronaryArteryDisease_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results/dlpfc_CAD
$MAGMA --gene-results SNP_Data/CoronaryArteryDisease_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results/sacc_CAD
$MAGMA --gene-results SNP_Data/CoronaryArteryDisease_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results/hpc_CAD
$MAGMA --gene-results SNP_Data/CoronaryArteryDisease_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results/nac_CAD
$MAGMA --gene-results SNP_Data/CoronaryArteryDisease_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results/amy_CAD


echo "**** Job ends ****"
date
