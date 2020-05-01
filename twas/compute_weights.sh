#!/bin/bash

## Usage:
# sh compute_weights.sh

mkdir -p logs

CORES=10

for region in DLPFC HIPPO #DentateGyrus
do
    
    # for feature in gene exon jxn tx
    for feature in exon jxn
    do
        echo "$feature"
        if [ "${feature}" == "gene" ] 
        then
            COREMEM=4
        elif [ "${feature}" == "exon" ] 
        then
            COREMEM=10
        elif [ "$feature" == 'jxn' ]
        then
            COREMEM=10
        else
            COREMEM=6
        fi
        
        SHORT="compute_weights_full_${region}_${feature}"

        # Construct shell file
        echo "Creating script for chromosome ${region} at the ${feature} level for ${CORES} cores with ${COREMEM}G of mem per core"
        cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=${COREMEM}G,h_vmem=${COREMEM}G,h_fsize=100G
#$ -pe local ${CORES}
#$ -N ${SHORT}
#$ -o ./logs/${SHORT}.txt
#$ -e ./logs/${SHORT}.txt
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
Rscript compute_weights.R -r ${region} -f ${feature} -c ${CORES} -p FALSE

echo "**** Job ends ****"
date
EOF
        call="qsub .${SHORT}.sh"
        echo $call
        $call
    done
done