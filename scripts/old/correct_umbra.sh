#! /bin/bash

# I used this script, because I falsely merged samples of D. umbra, due to bad sample naming!!!

cd /home/uibk/c7701178/scratch/DAPHNIA/daphnia_snakemake_pbs/bbmap/rapid

for i in FP2_* F3J_*; do
	#echo $i
	id=$(echo $i | cut -f1,2 -d '_')
	#echo $id
	#cp $i $id.merged.bb.RAPID.bam
done

## remove falsely merged bam
#rm FP2.merged.bb.RAPID.bam
#rm F3J.merged.bb.RAPID.bam
