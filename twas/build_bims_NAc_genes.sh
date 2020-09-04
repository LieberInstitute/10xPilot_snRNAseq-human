#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -pe local 10
#$ -N "build_bims_NAc_genes"
#$ -m e
#$ -j y
#$ -o logs/build_bims_NAc_genes_$JOB_ID.txt

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load dependencies
module load plink/1.90b6.6
module load fusion_twas/github
module load conda_R/4.0

## List current modules
module list

## Compute weights for the given region/feature pair
Rscript build_bims.R -c 10 -t FALSE

echo "**** Job ends ****"
date
