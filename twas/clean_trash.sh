#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=1G,h_vmem=1G
#$ -o logs/clean_trash.txt
#$ -e logs/clean_trash.txt
#$ -m e

echo "**** Job starts ****"
date

rm -fr /dcs04/lieber/marmaypag/Tran_LIBD001/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/twas/trash/*

echo "**** Job ends ****"
date
