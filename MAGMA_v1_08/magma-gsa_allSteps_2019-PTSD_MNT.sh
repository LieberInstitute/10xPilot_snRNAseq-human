#!/bin/bash
#$ -cwd
#$ -N magma-gsa_step3_PTSD-rev_MNT
#$ -o ./logs/magma-gsa_step3-PTSD_rev_MNT18Jul2021.o
#$ -e ./logs/magma-gsa_step3-PTSD_rev_MNT18Jul2021.e
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


## Step 1 - Annotation (SNP : gene mapping)
# $MAGMA --annotate window=35,10 --snp-loc ../MAGMA/GWAS_Results/PTSD_Nievergelt2019.snploc --gene-loc $ANNO --out SNP_Data/PTSD_Nievergelt2019_10xPilotGenes

## Step 2 - Gene analysis (from SNP-wise summary stats)
# $MAGMA --bfile $BFILE --gene-annot SNP_Data/PTSD_Nievergelt2019_10xPilotGenes.genes.annot --pval ../MAGMA/GWAS_Results/MVP.TE.REX.txt use=SNP,P ncol=N --gene-model ${model} --out SNP_Data/PTSD_Nievergelt2019_10xPilotGenes_${model}

## Step 3 - Gene set analyses (using gene-level output)  -  for five brain regions
echo "Update MNT 05May2021 - running an iteration (into 'Results_v2') with the preprint stats, applying a filter for non-0-median expression, in each respective subcluster."
echo "Update MNT 14Jul2021 - running an iteration (into 'Results_rev') following the above approach, finally with final revision cell classes."


$MAGMA --gene-results SNP_Data/PTSD_Nievergelt2019_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results_rev/dlpfc_PTSD
# $MAGMA --gene-results SNP_Data/PTSD_Nievergelt2019_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results_rev/sacc_PTSD
$MAGMA --gene-results SNP_Data/PTSD_Nievergelt2019_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results_rev/hpc_PTSD
$MAGMA --gene-results SNP_Data/PTSD_Nievergelt2019_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results_rev/nac_PTSD
# $MAGMA --gene-results SNP_Data/PTSD_Nievergelt2019_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results_rev/amy_PTSD


echo "**** Job ends ****"
date
