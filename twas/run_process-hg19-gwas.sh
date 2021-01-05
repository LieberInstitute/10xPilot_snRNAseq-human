#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -pe local 6
#$ -N process-hg19-gwas
#$ -j y
#$ -o logs/process-hg19-gwas_$JOB_ID.txt

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${TASK_ID}"

module load conda_R/4.0.x

## List current modules
module list

Rscript process-hg19-gwas.R

echo "**** Job ends ****"
date
