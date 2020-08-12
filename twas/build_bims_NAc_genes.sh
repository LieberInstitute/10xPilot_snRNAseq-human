#!/bin/bash
#$ -cwd
#$ -l mem_free=7G,h_vmem=7G,h_fsize=100G
#$ -pe local 10
#$ -N "build_bims_NAc_genes"
#$ -o ./logs/build_bims_NAc_genes.txt
#$ -e ./logs/build_bims_NAc_genes.txt
#$ -m e

mkdir -p logs

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
Rscript build_bims.R -c 10 -p FALSE

echo "**** Job ends ****"
date
