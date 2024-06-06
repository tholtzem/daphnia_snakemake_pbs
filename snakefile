include: "rules/common.smk"



# -----------------------------------------------

rule all:
	input:
		#expand('ref/bb_indexRef.done'),
		#expand('bbmap/HiC/minid76/{sample}.bam', sample=fastq_prefix),
		#expand('/home/uibk/c7701178/local/bbmap/HiC/minid76/{sample}.bam', sample=['Pfaeff2608_FDSW202515190-1r_HHG5FDSXY_L2']),
		#expand('/home/uibk/c7701178/local/bbmap/HiC/minid76/{sample}.rsync.done', sample=fastq_prefix)
		#expand('bbmap/HiC/minid95/{sample}.{ext}', sample=fastq_prefix, ext=['bam']),
		#expand('bbmap/HiC/minid95/{sample}.{ext}', sample=fastq_prefix, ext=['bam.bai']),
		#expand('bbmap/HiC/minid95/merged/{sample}.{ext}', sample=sample_ID, ext=['merged.bb.HiC.bam']),
		#expand('bbmap/HiC/minid95/merged/{sample}.{ext}', sample=sample_ID, ext=['merged.bb.HiC.bam.bai']),
		#expand('deDup/HiC/minid95/{sample}.{ext}', sample=sample_ID, ext=['bb.HiC.dedup.bam', 'bb.HiC.dedup.metrics.txt']),
                #expand('deDup/HiC/minid95/{sample}.{ext}', sample=sample_ID, ext=['bb.HiC.overlapclipped.bam']),
                #expand('deDup/HiC/minid95/{sample}.{ext}', sample=sample_ID, ext=['bb.HiC.overlapclipped.bam.bai']),
                #expand('ref/Dgaleata_M5_PBasm.FINAL.{ext}', sample=sample_ID, ext=['fasta.fai', 'dict']),
                #expand('deDup/HiC/minid95/{sample}.{ext}', sample=sample_ID, ext=['added2ClippedList.done']),
                #'list/HiC_indels.list',
                #expand('realigned/HiC/minid95/{sample}.{ext}', sample=sample_ID, ext=['bb.HiC.realigned.bam']),
		expand('/scratch/c7701178/mach2/DAPHNIA/daphnia_snakemake_pbs/realigned/HiC/minid95/{sample}.{ext}', sample=sample_ID, ext=['added2RealignedList.list']),
		#expand('bcf/daphnia_init_{chrom}.{ext}', chrom=chromosom_names[0:10], ext=['bcf']),
                #expand('depth/HiC/minid95/{sample}.{ext}', sample=sample_ID, ext=['bb.HiC.realigned.bam.depth.gz']),
                #expand('depth/HiC/minid95/{sample}.{ext}', sample=sample_ID, ext=['bb.HiC.realigned.bam.coverage.hist']),
		#expand('depth/HiC/minid95/stats/references_{ext}', ext=['depth_statistics.txt']),
                #expand('depth/HiC/minid95/stats/references_{ext}', ext=['depth_hist_dfAll.pdf', 'depth_hist_df10.pdf', 'depth_df1_boxplot.pdf', 'depth_df10_boxplot.pdf', 'depthFilter.list', 'depth_NonZero_hist_perSample.pdf', 'dedupBAM_depth1.list', 'dedupBAM_depth10.list', 'dedupBAM_depth_under10.list']),

		#expand('qmap/HiC/minid95/{sample}.{ext}', sample=sample_ID, ext=['bb.HiC.realigned'])
		#expand('qmap/HiC/minid95/{sample}.{ext}', sample=sample=sample_ID, ext=['bb.HiC.realigned']),
		#expand('bedtools/HiC/minid95/{sample}.{ext}', sample=sample_ID, ext=['bb.HiC.realigned.genomecov.bed']),
		#expand('bedtools/HiC/minid95/plots/{sample}.{ext}', sample=sample_ID, ext=['bb.HiC.pdf']),
		#expand('depth/stats/genome_stats.done'),
		##expand('plot_summary.done'),
		##expand('list/makePOPlist.done', species=species), # execute only once!
		##expand('bcftools/daphnia_init_dgal39_cucullata.vcf.gz.tbi')
		#### --- works --- ##
		
		#expand('deDup/{sample}.{ext}', sample=sample_ID, ext=['bb.RAPID.dedup.bam.bai']),
		#expand('list/BCF_{species}.list',  species=['lacustris', 'umbra']),
		#expand('bcf/concat/daphnia_{species}.{ext}', species=['cucullata', 'galeata', 'longispina', 'longispinaFIN', 'mendotae', 'dentifera', 'curvirostris', 'lacustris', 'umbra', 'tanakai', 'zschokkei'], ext=['bcf', 'bcf.csi']),
		#expand('bcf/GRE14/GRE1408_{chrom}.bcf', chrom=chromosom_names),
		#expand('list/BCF_GRE1408.list'),
		#expand('bcf/concat/GRE1408.{ext}', ext=['bcf', 'bcf.csi']),
		#expand('bcf/concat/daphnia_galeata2.{ext}', ext=['bcf']),
		#expand('bcf/merged/daphnia_10pops.{ext}' , ext=['bcf']),
		#expand('vcf/daphnia_10pops.{ext}', ext=['vcf.gz']),
		#expand('vcf/daphnia_{species}.{ext}', species=['cucullata', 'galeata2', 'longispina', 'longispinaFIN', 'mendotae', 'dentifera', 'curvirostris', 'lacustris', 'umbra', 'tanakai', 'zschokkei'], ext=['vcf.gz'])
		#expand('bcf/concat/daphnia_{species}.{ext}', species=['lacustris', 'umbra'], ext=['bcf.csi']),
		#expand('bcf/merged/daphnia_9pops.{ext}', ext=['bcf']),
		#expand('vcf/daphnia_{species}.{ext}', species=['cucullata', 'galeata', 'longispina', 'longispinaFIN', 'mendotae', 'dentifera', 'curvirostris', 'lacustris', 'umbra', 'tanakai', 'zschokkei', '9pops'], ext=['vcf.gz']),
		#expand('vcf/stats/daphnia_{species}_nbrSites.txt', species=['cucullata', 'galeata', 'longispina', 'longispinaFIN', 'mendotae', 'dentifera', 'curvirostris', 'lacustris', 'umbra', 'tanakai', 'zschokkei'], ext=['vcf.gz']), 
                #expand('vcf/subset/daphnia_{species}_100Ksubset.{ext}', species=['cucullata', 'galeata', 'longispina', 'longispinaFIN', 'mendotae', 'dentifera', 'curvirostris', 'lacustris', 'umbra', 'tanakai', 'zschokkei'], ext=['vcf.gz']),
                #expand('vcf/stats/daphnia_{species}_100Ksubset.{ext}', species=['cucullata', 'galeata', 'longispina', 'longispinaFIN', 'mendotae', 'dentifera', 'curvirostris', 'lacustris', 'umbra', 'tanakai', 'zschokkei'], ext=['frq.done', 'idepth.done', 'ldepth.mean.done', 'lqual.done', 'imiss.done', 'lmiss.done', 'MQ_DP.txt'])
    #-------------------------------------------------filtering-----------------------------------------------------------------#
    #expand('vcf/concat/tags/daphnia_{species}_alltags.{ext}', species=['galeata', 'longispina', 'cucullata'], ext=['vcf.gz.tbi']),
    #expand('vcf/filtered/daphnia_{species}_Allsites_hardfiltered_1.{ext}', species=['galeata', 'longispina', 'cucullata'], ext=['vcf.gz']),
    #expand('vcf/filtered/daphnia_{species}_Allsites_hardfiltered_2.{ext}', species=['galeata', 'longispina', 'cucullata'], ext=['vcf.gz']),
    #expand('vcf/filtered/daphnia_{species}_{ext}', species=['galeata', 'longispina', 'cucullata'], ext=['invariant.vcf.gz', 'variant.vcf.gz', 'invariant.vcf.gz.tbi', 'variant.vcf.gz.tbi']),
    #expand('vcf/filtered/daphnia_{species}_{ext}', species=['galeata', 'longispina', 'cucullata'], ext=['final.vcf.gz']) 

# -----------------------------------------------


#include: "rules/hts.smk"
#include: "rules/markdups.smk"
#include: "rules/gencov.smk"
#include: "rules/samtools_depth.smk"
#include: "rules/read_depth.smk"
#include: "rules/bcftools_GRE.smk"
include: "rules/bcftools.smk"
#include: "rules/vcf_stats.smk"
#include: "rules/filterVCF.smk"
#include: "rules/tabix.smk"
#include: "rules/vars.smk"
#include: "rules/qualimap.smk"
