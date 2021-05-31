rule samtools_depth:
  input:
    'deDup/{sample}.dedup.bam'
  output:
    'depth/{sample}.dedup.bam.depth.gz'
  log: 'log/{sample}.dedup.bam.depth.log'
  threads: 12
  message:
    """ Count per position depth per sample using samtools depth """
  shell:
    """
    samtools depth -aa {input} | cut -f3 | gzip > {output} 2> {log}
    """
