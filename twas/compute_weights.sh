#!/bin/bash
#$ -cwd
#$ -l mem_free=4G,h_vmem=4G,h_fsize=100G
#$ -pe local 10
#$ -N compute_weights_full_NAc_genes
#$ -o ./logs/compute_weights_full_NAc_genes.txt
#$ -e ./logs/compute_weights_full_NAc_genes.txt
#$ -m a

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${TASK_ID}"

## Load dependencies
module load plink/1.90b6.6
module load fusion_twas/github
module load conda_R/3.6

## List current modules
module list

## Compute weights for the given region/feature pair
Rscript compute_weights.R -c 10 

echo "**** Job ends ****"
date
