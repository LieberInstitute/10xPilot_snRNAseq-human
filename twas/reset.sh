#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=4G,h_vmem=4G
#$ -o logs/reset.txt
#$ -e logs/reset.txt
#$ -m e

echo "**** Job starts ****"
date

rm -fr /dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/twas/NAc_gene
rm -fr /dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/twas/logs/build_bims_NAc_genes_*.txt

echo "**** Job ends ****"
date
