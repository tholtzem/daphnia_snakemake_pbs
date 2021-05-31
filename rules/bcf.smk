
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

rule bcf_mpileup:
    input:
      ref = config['ref'],
      #bamlist = 'list/dedupBAM.list',
      region_file = 'list/dgal_rapid_regionFile.tsv'
    output:
      'bcftools/daphnia_init.vcf.gz'
      #touch('bcftools/bcftools_done')
    log: 'log/bcftools_done.log'
    threads: 24
    message: """--- Generate AllSites vcf using BCFtools (mpileup/call) ---"""
    shell:
      """
      bcftools mpileup -f {input.ref} -b list/dedupBAM.list -R {input.region_file} | bcftools call -m -Oz -f GQ -o {output} 2> {log}
      """
