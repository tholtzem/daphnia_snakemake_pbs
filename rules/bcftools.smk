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

rule bcf_mpileup_call:
    input:
      ref = config['ref'],
      #bamlist = 'list/dedupBAM.list',
      #region_file = 'list/dgal_rapid_regionFile.tsv'
    output:
      protected('bcftools/daphnia_init_{chrom}.vcf.gz')
    log: 'log/bcftools_daphnia_init_{chrom}.log'
    threads: 24
    message: """--- Generate AllSites vcf using BCFtools (mpileup/call) ---"""
    shell:
      """
      bcftools mpileup -f {input.ref} -b list/dedupBAM.list -r {wildcards.chrom} | bcftools call -m -Oz -f GQ -o {output} 2> {log}
      """
