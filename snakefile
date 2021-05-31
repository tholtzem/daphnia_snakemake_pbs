include: "rules/common.smk"



# -----------------------------------------------

rule all:
  input:
    expand('bbmap/rapid/{sample}.{ext}', sample=sample_names, ext=['bam']),
    expand("mytask.done"),
    expand('deDup/{sample}.{ext}', sample=sample_ID, ext=['bb.RAPID.dedup.bam', 'bb.RAPID.dedup.metrics.txt']),
    expand('markDup/{sample}.{ext}', sample=sample_ID, ext=['bb.RAPID.mrkdup.bam', 'bb.RAPID.mrkdup.metrics.txt']),
    expand('bedtools/{sample}.{ext}', sample=sample_ID, ext=['bb.RAPID.dedup.genomecov.bed', 'bb.RAPID.dedup.genomecov.bedgraph']), 
    expand('bedtools/plots/{sample}.{ext}', sample=sample_ID, ext=['bb.RAPID.pdf']),
    #expand('depth/{sample}.{ext}', sample=sample_ID, ext=['bb.RAPID.dedup.bam.depth.gz']),
    #expand('list/depth_dedup.list'),
    #expand('genome_stats.done'),
    #expand('plot_summary.done'),
    expand('list/dedupBAM.list'),
    expand('bcftools/daphnia_init_{chrom}.{ext}', chrom=chromosom_names, ext=['vcf.gz'])


# -----------------------------------------------


#include: "rules/hts.smk"
include: "rules/s2b_mergeBam.smk"
include: "rules/markdups.smk"
include: "rules/gencov.smk"
include: "rules/samtools_depth.smk"
include: "rules/read_depth.smk"
include: "rules/bcftools.smk"
#include: "rules/vars.smk"
