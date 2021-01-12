#!/bin/bash

## Usage:
# sh apply_weights.sh

umask u=rwx,g=rwx,o=

mkdir -p logs

for region in NAc
do

    for feature in gene
    # for feature in gene exon jxn tx
    do

        # set of summary stats
        for summstats in si sc dpw cpd aoi
        do

            SHORT="apply_weights_full_${region}_${feature}_${summstats}"

            # Construct shell file
            echo "Creating script for ${region} at the ${feature} level for ${summstats}"
            cat > .${SHORT}.sh <<EOF

#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=30G,h_vmem=30G,h_fsize=100G
#$ -N ${SHORT}
#$ -o ./logs/${SHORT}.txt
#$ -e ./logs/${SHORT}.txt
#$ -m e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${TASK_ID}"

## Load dependencies
module load fusion_twas/github
module load conda_R/4.0.x

## List current modules
module list

## Choose the correct GWAS summary statistics file
if [ "${summstats}" == "si" ]
then
    summstatsfile="clean_gwas/SmokingInitiation_Clean.txt"
elif [ "${summstats}" == "sc" ]
then
    summstatsfile="clean_gwas/SmokingCessation_Clean.txt"
elif [ "${summstats}" == "dpw" ]
then
    summstatsfile="clean_gwas/DrinksPerWeek_Clean.txt"
elif [ "${summstats}" == "cpd" ]
then
    summstatsfile="clean_gwas/CigarettesPerDay_Clean.txt"
elif [ "${summstats}" == "aoi" ]
then
    summstatsfile="clean_gwas/AgeofInitiation_Clean.txt"
else
    echo "Unexpected ${summstats} input"
fi

## Apply weights for the given region/feature pair and the given GWAS summary statistics
mkdir -p "${region}_${feature}/${summstats}"


for chr in {1..22}
do
    echo "*************************"
    echo ""
    echo "processing chromosome \${chr}"
    date
    echo ""

## Create summarized analysis
Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.assoc_test.R \
    --sumstats \${summstatsfile} \
    --weights ${region}_${feature}/${region}_${feature}.pos \
    --weights_dir ${region}_${feature}/ \
    --ref_ld_chr /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38/1000G.EUR. \
    --chr \${chr} \
    --out ${region}_${feature}/${summstats}/${summstats}.\${chr}.dat

    echo ""
    echo "making plots for chromosome \${chr}"
    date
    echo ""

## companion post-processing step (plots only for genes)
if [ "$feature" == "gene" ]
then
    Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.post_process.R \
        --sumstats \${summstatsfile} \
        --input ${region}_${feature}/${summstats}/${summstats}.\${chr}.dat \
        --out ${region}_${feature}/${summstats}/${summstats}.\${chr}.analysis \
        --ref_ld_chr /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38/1000G.EUR. \
        --chr \${chr} \
        --plot --locus_win 100000 --verbose 2 --plot_individual --plot_eqtl --plot_corr \
        --glist_path "/jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/glist-hg38"
fi

done

echo "**** Job ends ****"
date
EOF
            call="qsub .${SHORT}.sh"
            echo $call
            $call
        done
    done
done
