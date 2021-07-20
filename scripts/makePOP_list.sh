#! /bin/bash


SAMPLELIST=/home/uibk/c7701178/scratch/DAPHNIA/daphnia_snakemake_pbs/list/prefixes_20210709.list # Path to a list of prefixes of the raw fastq files. It should be a subset of the the 1st column of the fastq table.

SAMPLETABLE=/home/uibk/c7701178/scratch/DAPHNIA/daphnia_snakemake_pbs/list/meta_data_20210709.tsv # Path to a fastq table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID. 


for SAMPLEFILE in `cat $SAMPLELIST | tail -n+2`; do
	# For each prefix, extract the associated species (column 2) from the table
	SPECIES=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
	#echo $SAMPLEFILE
	#echo $SPECIES
	#echo mkdir deDup/$SPECIES
	#echo rsync -avP deDup/$SAMPLEFILE* deDup/$SPECIES/
	ls deDup/$SAMPLEFILE*bam >> list/$SPECIES.list
done

