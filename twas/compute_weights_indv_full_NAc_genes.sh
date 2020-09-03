#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=4G,h_vmem=4G,h_fsize=100G
#$ -N compute_weights_full_NAc_genes
#$ -o logs/NAc_genes/compute_weights_indv_full_NAc_genes.txt
#$ -e logs/NAc_genes/compute_weights__indv_full_NAc_genes.txt
#$ -t 1-${NUM}
#$ -tc 40
#$ -m a

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${SGE_TASK_ID}"

## Load dependencies
module load plink/1.90b6.6
module load fusion_twas/github
module load conda_R/4.0

## List current modules
module list

# relative path for FILELIST
FILELIST=$(echo "twas/NAc_gene/input_ids.txt")

## File id and feature name
FEATURENUM=\$(awk 'BEGIN {FS="\t"} {print \$1}' ${FILELIST} | awk "NR==\${SGE_TASK_ID}")
FEATUREID=\$(awk 'BEGIN {FS="\t"} {print \$2}' ${FILELIST} | awk "NR==\${SGE_TASK_ID}")

## Change directories
cd NAc/genes

## Define files
FILTBIM="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/twas/filter_data/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38_filtered_NAc_Nicotine"
TMPFILES="tmp_files/genes_\${FEATURENUM}"
OUTFILES="out_files/genes_\${FEATURENUM}"

Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.compute_weights.R \
    --bfile \${FILTBIM} \
    --tmp \${TMPFILES} \
    --out \${OUTFILES} \
    --PATH_gemma /dcl01/lieber/ajaffe/lab/twas/software/gemma-0.98.1-linux-static \
    --models top1,blup,lasso,enet --hsq_p 1.0001 --verbose 1 --save_hsq

## If you want to keep the temp files
# --noclean TRUE
#
# echo "**** Job ends ****"
# date
# EOF
#         #call="qsub .compute_weights_full_NAc_genes.sh"
#         call="Rscript -e \"sgejobs::array_submit_num(job_bash = '.compute_weights_full_NAc_genes_indv.sh', array_num = ${NUM}, submit = TRUE); options(width = 120); sessioninfo::session_info()\""
#         echo $call
#         #$call
#         Rscript -e "sgejobs::array_submit_num(job_bash = '.compute_weights_full_NAc_genes_indv.sh', array_num = ${NUM}, submit = TRUE); options(width = 120); sessioninfo::session_info()"
#     done
# done
