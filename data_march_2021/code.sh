#!/bin/bash

## Increase the h_fsize just in case
#  qrsh -l h_fsize=200G
module load bs/1.3.0

## Follow the instructions from 
## https://support.illumina.com/bulletins/2017/02/options-for-downloading-run-folders-from-basespace-sequence-hub.html
## and
## https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview
# bs auth ## Involves opening a link in your browser
## Not needed if ~/.basepace/default.cfg exists already

## From Linda:
# I transferred 3 runs to Leo.  One is the reprocessed 5701NAc, one is the new 3 samples and one, just in case, is the re-analysis of the original KMay021021.  If you don't want to review those, you can feel free to discard, but I just wanted to make sure there were no issues since discovering the sample sheet was wrong.  Let me know if anyone has questions!
# 
# Thanks
# Linda
# 
# Linda D. Orzolek, M.S.
# Director
# Johns Hopkins University
# Transcriptomics & Deep Sequencing Core
# 733 North Broadway, MRB 363
# Baltimore, MD 21205
# 410-502-6658

bs download project -n KMay021021b -o KMay021021b
bs download project -n KMay022421 -o KMay022421
bs download project -n KMay021021_01 -o KMay021021_01
