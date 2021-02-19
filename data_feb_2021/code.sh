#!/bin/bash

## Increase the h_fsize just in case
#  qrsh -l h_fsize=200G
module load bs/1.3.0

## Follow the instructions from 
## https://support.illumina.com/bulletins/2017/02/options-for-downloading-run-folders-from-basespace-sequence-hub.html
## and
## https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview
bs auth ## Involves opening a link in your browser

## Unlike the instructions, I couldn't find the run ID numbers. yet
## bs download project --help
## gave me the info I needed to run this command

bs download project -n KMay021021 -o FASTQ

## Gives output like this:

# 5701SAcc_S38_L002_R1_001.fastq.gz 635.00 MiB / 1.10 GiB [==================================>---------------------------]  56.16% 03m07s
# Br5207DLPFC_S3_L002_R1_001.fastq.gz 605.57 MiB / 1.05 GiB [=================================>--------------------------]  56.38% 03m06s
# Br5400SACC_S8_L002_R1_001.fastq.gz 610.00 MiB / 1.03 GiB [===================================>-------------------------]  57.59% 02m57s
