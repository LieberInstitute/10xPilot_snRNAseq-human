#!/bin/bash
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=100G
#$ -pe local 1
#$ -N munge_gwas
#$ -j y
#$ -o logs/munge-gwas_$JOB_ID.txt

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${TASK_ID}"

module load conda/3-4.6.14
module load git

## List current modules
module list

if [ -d "ldsc/" ]
then
    echo "LDSC exists. Reformatting GWAS"
    cd ldsc/
    mkdir -p ../NAc_GWAS/munge_gwas/
    ./munge_sumstats.py --sumstats ../NAc_GWAS/AgeofInitiation.txt --out ../NAc_GWAS/munge_gwas/AgeofInitiation_hg19_munge.txt --a1 REF --a2 ALT --p PVALUE --snp RSID --nstudy Number_of_Studies --chunksize 512
    ./munge_sumstats.py --sumstats ../NAc_GWAS/CigarettesPerDay.txt --out ../NAc_GWAS/munge_gwas/CigarettesPerDay_hg19_munge.txt --a1 REF --a2 ALT --p PVALUE --snp RSID --nstudy Number_of_Studies --chunksize 512
    ./munge_sumstats.py --sumstats ../NAc_GWAS/SmokingCessation.txt --out ../NAc_GWAS/munge_gwas/SmokingCessation_hg19_munge.txt --a1 REF --a2 ALT --p PVALUE --snp RSID --nstudy Number_of_Studies --chunksize 512
    ./munge_sumstats.py --sumstats ../NAc_GWAS/SmokingInitiation.txt --out ../NAc_GWAS/munge_gwas/SmokingInitiation_hg19_munge.txt --a1 REF --a2 ALT --p PVALUE --snp RSID --nstudy Number_of_Studies --chunksize 512
    ./munge_sumstats.py --sumstats ../NAc_GWAS/DrinksPerWeek.txt --out ../NAc_GWAS/munge_gwas/DrinksPerWeek_hg19_munge.txt --a1 REF --a2 ALT --p PVALUE --snp RSID --nstudy Number_of_Studies --chunksize 512
else
    echo "No LDSC/, Gitcloning it now"
    git clone https://github.com/bulik/ldsc.git
    cd ldsc
    conda env create --file environment.yml
    source activate ldsc
    echo "Reformatting GWAS"
    mkdir -p ../NAc_GWAS/munge_gwas/
    ./munge_sumstats.py --sumstats ../NAc_GWAS/AgeofInitiation.txt --out ../NAc_GWAS/munge_gwas/AgeofInitiation_hg19_munge.txt --a1 REF --a2 ALT --p PVALUE --snp RSID --nstudy Number_of_Studies --chunksize 512
    ./munge_sumstats.py --sumstats ../NAc_GWAS/CigarettesPerDay.txt --out ../NAc_GWAS/munge_gwas/CigarettesPerDay_hg19_munge.txt --a1 REF --a2 ALT --p PVALUE --snp RSID --nstudy Number_of_Studies --chunksize 512
    ./munge_sumstats.py --sumstats ../NAc_GWAS/SmokingCessation.txt --out ../NAc_GWAS/munge_gwas/SmokingCessation_hg19_munge.txt --a1 REF --a2 ALT --p PVALUE --snp RSID --nstudy Number_of_Studies --chunksize 512
    ./munge_sumstats.py --sumstats ../NAc_GWAS/SmokingInitiation.txt --out ../NAc_GWAS/munge_gwas/SmokingInitiation_hg19_munge.txt --a1 REF --a2 ALT --p PVALUE --snp RSID --nstudy Number_of_Studies --chunksize 512
    ./munge_sumstats.py --sumstats ../NAc_GWAS/DrinksPerWeek.txt --out ../NAc_GWAS/munge_gwas/DrinksPerWeek_hg19_munge.txt --a1 REF --a2 ALT --p PVALUE --snp RSID --nstudy Number_of_Studies --chunksize 512
fi

python

echo "**** Job ends ****"
date
