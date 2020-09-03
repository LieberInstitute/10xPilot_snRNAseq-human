#!/bin/bash

## Clean up older files
# rm -fr trash

## Could be fancier and use the date or something
mkdir -p trash
mv logs/build_bims_NAc_genes.txt trash/
mv NAc_genes trash/

## Create logs dir if needed
mkdir -p logs
## Submit new job

# Arg1: test
# Arg2: degradation
qsub build_bims_NAc_genes.sh TRUE TRUE
