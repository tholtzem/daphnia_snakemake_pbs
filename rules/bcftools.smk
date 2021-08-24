rule ls_dedupBAM:
  input:
    #'depth/{sample}.depth.gz'
  output:
    'list/dedupBAM.list'
  log: 'log/dedupBAM.list.log'
  message: """--- Creating a sample list of deduplicated bam files ---"""
  shell:
    """
    ls deDup/*dedup.bam > {output} 2> {log}
    """

rule makePop_lists:
  input:#
  output:
    touch('list/makePOPlist.done')
  log: 'log/makePOPlist.log'
  message: """--- Creating population lists of deduplicated bam files ---"""
  shell:
    """
    scripts/makePOP_list.sh
    """

rule index_bam:
  input:
    'deDup/{sample}.bb.RAPID.dedup.bam'
  output:
    'deDup/{sample}.bb.RAPID.dedup.bam.bai'
  log: 'log/{sample}.index_bam.log'
  message: """ --- Indexing bam files --- """
  shell:
    """
    samtools index -b {input} 2> {log}
    """


rule bcf_mpileup_call:
  input:
    ref = config['ref'],
    touched = 'list/makePOPlist.done'
  output:
    'bcftools/daphnia_init_{chrom}_{species}.vcf.gz'
  log:
    'log/daphnia_init_{chrom}_{species}.log'
    #'log/daphnia_init_dgal1_{species}.log'
  threads: 48
  message: """--- Generate AllSites vcf using BCFtools (mpileup/call) ---"""
  shell:
    """
    bcftools mpileup -f {input.ref} -b list/{wildcards.species}.list -r {wildcards.chrom} | bcftools call -m -Oz -f GQ -o {output} --threads {threads} 2> {log}
    """

