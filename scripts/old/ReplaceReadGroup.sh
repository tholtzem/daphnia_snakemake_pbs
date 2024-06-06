#! /bin/bash

BASEDIR=/home/uibk/c7701178/scratch/DAPHNIA/daphnia_snakemake_pbs

UMBRALIST=$BASEDIR/list/read_groups_umbra.list

BAMS=$BASEDIR/bbmap/rapid

for i in `cat $UMBRALIST`; do
	ID=$( echo $i | sed 's/.bb.RAPID.bam//')
	#echo $ID
	SM=$( echo $ID | cut -f1,2 -d'_')
	#echo $SM
	LB=$( echo $ID | cut -f3 -d'_')
	#echo $LB
	PL=ILLUMINA
	#echo $PL
	#ls $BAMS/$SM*
	picard AddOrReplaceReadGroups --INPUT $BAMS/$i --OUTPUT $BAMS/${ID}.bb.RAPID.RG.bam --RGLB $LB --RGPL $PL --RGSM $SM --RGID $ID --RGPU $LB
done
