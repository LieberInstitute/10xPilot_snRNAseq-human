#!/bin/bash

## Usage:
# sh compute_weights_indv.sh

mkdir -p logs
module load conda_R/3.6


for region in DLPFC HIPPO DentateGyrus
do
    
    for feature in gene exon jxn tx
    do

        SHORT="compute_weights_full_${region}_${feature}"
        
        FILELIST="/dcl01/lieber/ajaffe/lab/twas/bsp2/${region}/${feature}/input_ids.txt"
        NUM=$(cat $FILELIST | awk '{print $NF}' | uniq | wc -l)
		
		echo ${NUM}
        
        ## Create logs dir
        mkdir -p logs/${region}/${feature}
        
        ## Create output directories
        mkdir -p ${region}/${feature}/tmp_files
        mkdir -p ${region}/${feature}/out_files
        
        ## For GEMMA
        ln -s /dcl01/lieber/ajaffe/lab/twas/bsp2/${region}/${feature}/ ${region}/${feature}/output

        # Construct shell file
        echo "Creating script for chromosome ${region} at the ${feature} level"
        cat > .${SHORT}_indv.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=4G,h_vmem=4G,h_fsize=100G
#$ -N ${SHORT}
#$ -o ./logs/${region}/${feature}/${SHORT}.\$TASK_ID.txt
#$ -e ./logs/${region}/${feature}/${SHORT}.\$TASK_ID.txt
# -t 300001-375000
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
module load conda_R/3.6

## List current modules
module list

## File id and feature name
FEATURENUM=\$(awk 'BEGIN {FS="\t"} {print \$1}' ${FILELIST} | awk "NR==\${SGE_TASK_ID}")
FEATUREID=\$(awk 'BEGIN {FS="\t"} {print \$2}' ${FILELIST} | awk "NR==\${SGE_TASK_ID}")

## Change directories
cd ${region}/${feature}

## Define files
FILTBIM="bim_files/${region}_${feature}_\${FEATURENUM}/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_LDfiltered_${region}_${feature}_\${FEATURENUM}"
TMPFILES="tmp_files/${feature}_\${FEATURENUM}"
OUTFILES="out_files/${feature}_\${FEATURENUM}"

Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.compute_weights.R \
    --bfile \${FILTBIM} \
    --tmp \${TMPFILES} \
    --out \${OUTFILES} \
    --PATH_gemma /dcl01/lieber/ajaffe/lab/twas/software/gemma-0.98.1-linux-static \
    --models top1,blup,lasso,enet --hsq_p 1.0001 --verbose 1 --save_hsq

## If you want to keep the temp files
# --noclean TRUE

echo "**** Job ends ****"
date
EOF
        #call="qsub .${SHORT}.sh"
        call="Rscript -e \"sgejobs::array_submit_num(job_bash = '.${SHORT}_indv.sh', array_num = ${NUM}, submit = TRUE); options(width = 120); sessioninfo::session_info()\""
        echo $call
        #$call
        Rscript -e "sgejobs::array_submit_num(job_bash = '.${SHORT}_indv.sh', array_num = ${NUM}, submit = TRUE); options(width = 120); sessioninfo::session_info()"
    done
done
