#!/bin/bash
#$ -cwd
#$ -N magma-gsa_step1-2_ADHD
#$ -o ./logs/magma-gsa_step1-2-ADHD_MNT22Jun2021.o
#$ -e ./logs/magma-gsa_step1-2-ADHD_MNT22Jun2021.e
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

SUMMSTATS=/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/MAGMA/GWAS_Results/daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta

## Step 1 - Annotation (SNP : gene mapping)
$MAGMA --annotate window=35,10 --snp-loc ../MAGMA/GWAS_Results/ADHD_PGC2018.snploc --gene-loc $ANNO --out SNP_Data/ADHD_PGC2018_10xPilotGenes

## Step 2 - Gene analysis (from SNP-wise summary stats)
$MAGMA --bfile $BFILE --gene-annot SNP_Data/ADHD_PGC2018_10xPilotGenes.genes.annot --pval $SUMMSTATS use=SNP,P ncol=Neff --gene-model ${model} --out SNP_Data/ADHD_PGC2018_10xPilotGenes_${model}


## Step 3 - Gene set analyses (using gene-level output)  -  for five brain regions

# $MAGMA --gene-results SNP_Data/ADHD_PGC2018_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_dlpfc gene-col=${genecol} set-col=${setcol} --out Results_v2/dlpfc_ADHD
# $MAGMA --gene-results SNP_Data/ADHD_PGC2018_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_sacc gene-col=${genecol} set-col=${setcol} --out Results_v2/sacc_ADHD
# $MAGMA --gene-results SNP_Data/ADHD_PGC2018_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_hpc gene-col=${genecol} set-col=${setcol} --out Results_v2/hpc_ADHD
# $MAGMA --gene-results SNP_Data/ADHD_PGC2018_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_nac gene-col=${genecol} set-col=${setcol} --out Results_v2/nac_ADHD
# $MAGMA --gene-results SNP_Data/ADHD_PGC2018_10xPilotGenes_snp-wise.genes.raw --set-annot $gs_amy gene-col=${genecol} set-col=${setcol} --out Results_v2/amy_ADHD


echo "**** Job ends ****"
date
