include: "rules/common.smk"



# -----------------------------------------------

rule all:
	input:
		#expand('bbmap/rapid/{sample}.{ext}', sample=sample_names, ext=['bam']),
		##expand("mytask.done"),
		#expand('deDup/{sample}.{ext}', sample=sample_ID, ext=['bb.RAPID.dedup.bam', 'bb.RAPID.dedup.metrics.txt']),
		##expand('markDup/{sample}.{ext}', sample=sample_ID, ext=['bb.RAPID.mrkdup.bam', 'bb.RAPID.mrkdup.metrics.txt']),
		##expand('bedtools/{sample}.{ext}', sample=sample_ID, ext=['bb.RAPID.dedup.genomecov.bed', 'bb.RAPID.dedup.genomecov.bedgraph']), 
		##expand('bedtools/plots/{sample}.{ext}', sample=sample_ID, ext=['bb.RAPID.pdf']),
		##expand('depth/{sample}.{ext}', sample=sample_ID, ext=['bb.RAPID.dedup.bam.depth.gz']),
		##expand('list/depth_dedup.list'),
		##expand('genome_stats.done'),
		##expand('plot_summary.done'),
		##expand('list/makePOPlist.done', species=species), # execute only once!
		##expand('bcftools/daphnia_init_dgal39_cucullata.vcf.gz.tbi')
		#### --- works --- ##
		##expand('deDup/{sample}.{ext}', sample=sample_ID, ext=['bb.RAPID.dedup.bam.bai']),
		##expand('vars/daphnia_init_{species}.{ext}', species=['galeata', 'longispina', 'cucullata', 'mendotae', 'dentifera'], ext=['vcf.gz']),
		expand('bcf/daphnia_init_{chrom}_{species}.{ext}', chrom=chromosom_names[0:100], species=['cucullata', 'galeata', 'longispina'], ext=['bcf']),
		expand('bcf/all/daphnia_init_{chrom}.{ext}', chrom=chromosom_names[0:2], ext=['mpileup.bcf']),
		expand('bcf/all/daphnia_init_{chrom}.{ext}', chrom=chromosom_names[0:2], ext=['call.bcf'])
    #expand('list/VCF_{species}.list', species=species[2]),
    #expand('vcf/concat/daphnia_{species}.{ext}', species=species[2], ext=['vcf.gz'])
    #expand('vcf/concat/tags/daphnia_{species}_alltags.{ext}', species=species, ext=['vcf.gz']),
    #expand('vcf/concat/tags/daphnia_{species}_alltags_rndsubset.{ext}', species=species, ext=['vcf.gz']),
    #expand('vcf/concat/stats/daphnia_{species}_alltags_rndsubset.{ext}', species=species, ext=['frq.done', 'idepth.done', 'ldepth.mean.done', 'lqual.done', 'imiss.done', 'lmiss.done']),
    #-------------------------------------------------filtering-----------------------------------------------------------------#
    #expand('vcf/concat/tags/daphnia_{species}_alltags.{ext}', species=['galeata', 'longispina', 'cucullata'], ext=['vcf.gz.tbi']),
    #expand('vcf/filtered/daphnia_{species}_Allsites_hardfiltered_1.{ext}', species=['galeata', 'longispina', 'cucullata'], ext=['vcf.gz']),
    #expand('vcf/filtered/daphnia_{species}_Allsites_hardfiltered_2.{ext}', species=['galeata', 'longispina', 'cucullata'], ext=['vcf.gz']),
    #expand('vcf/filtered/daphnia_{species}_{ext}', species=['galeata', 'longispina', 'cucullata'], ext=['invariant.vcf.gz', 'variant.vcf.gz', 'invariant.vcf.gz.tbi', 'variant.vcf.gz.tbi']),
    #expand('vcf/filtered/daphnia_{species}_{ext}', species=['galeata', 'longispina', 'cucullata'], ext=['final.vcf.gz']) 

# -----------------------------------------------


#include: "rules/hts.smk"
#include: "rules/s2b_mergeBam.smk"
#include: "rules/markdups.smk"
#include: "rules/gencov.smk"
#include: "rules/samtools_depth.smk"
#include: "rules/read_depth.smk"
include: "rules/bcftools.smk"
#include: "rules/callvars.smk"
#include: "rules/vcf_stats.smk"
#include: "rules/filterVCF.smk"
#include: "rules/tabix.smk"
#include: "rules/vars.smk"
