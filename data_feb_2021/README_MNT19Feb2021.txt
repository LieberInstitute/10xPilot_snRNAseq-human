## MNT comment 19Feb2021
# All of these sample data have been moved to their respective 'donor_tissue' dir in
MOVETO=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-10-15_Tran2021_published/Feb2021
# (where the n14+ pilot data FASTQs are)

## 10 samples
#mv 5276SAcc*/* $MOVETO/Br5276_sACC_neun/
#mv 5400NAC*/* $MOVETO/Br5400_NAc/
#mv 5476NAC*/* $MOVETO/Br5276_NAc/
#mv 5701Acc*/* $MOVETO/Br5701_NAc_neun/
	# Looks like the file names in these folders were prefixed '5701SAcc'...
	#	(which is the same prefix as those in the next FASTQ dirs)
	# -> in that dir they were moved to, run `rename 5701SAcc 5701NAc 5701SAcc*`
#mv 5701SAcc*/* $MOVETO/Br5701_sACC_neun/
#mv Br5207DLPFC*/* $MOVETO/Br5207_DLPFC/
#mv Br5276AMY*/* $MOVETO/Br5276_Amy_neun/
#mv Br5400AMY*/* $MOVETO/Br5400_Amy_neun/
#mv Br5400SACC*/* $MOVETO/Br5400_sACC/
#mv Br5701AMY*/* $MOVETO/Br5701_Amy/

# done MNT ~21:00 19Feb2021
