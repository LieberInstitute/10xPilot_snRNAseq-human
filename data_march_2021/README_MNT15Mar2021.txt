## MNT comment 15Mar2021
# All of these sample data have been moved to their respective 'donor_tissue' dir in
MOVETO=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-10-15_Tran2021_published/Feb2021
# (where the n14+ pilot expansion data FASTQs are)

## 3 samples (these were Frankenstein-prepped and sorted/run through 10x Chromium in-house
#mv KMay022421/Br2743_DLPFC_mid*/* $MOVETO/Br2743_DLPFC-mid/
#mv KMay022421/Br1204_Hb*/* $MOVETO/Br1204_Hb/
#mv KMay022421/Br5558_Hb*/* $MOVETO/Br5558_Hb/

# 1 sample: 5701_Nac_neun, which was re-library-prepped & sequnced (low QC and "sharp peak" in Bioanalyzer trace for first run)
mv KMay021021b/5701NAC*/* $MOVETO/Br5701_NAc_neun_v2/

# done MNT 22:27 15Mar2021
