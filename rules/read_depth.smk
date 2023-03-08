rule ls_depth:
  input:
    #'depth/{sample}.depth.gz'
  output:
    #touch('bbmap/all/{sample}.ls.done'),
    'list/depth_dedup.list'
    #'list/depth_outgroup.list'
  log: 'log/depth_dedup.list.log'
  message: """--- Creating a sample list of samtools-depth files for Rscript ---"""
  shell:
    """
    ls depth/*dedup.bam.depth.gz | cut -f2 -d '/' > {output} 2> {log}
    """

rule read_depth:
  #input:
    #test.list
  output:
    #pdf = "some.pdf"
    touch('depth/stats/genome_stats.done')
  log: 'log/genome_stats.done.log'
  threads: 12
  message:
    """ Running Rscript to plot the genome-wide distribution of coverage """
  shell:
    """
    Rscript scripts/read_depth_R.R 2> {log}
    """

rule plot_summary:
  #input:
    #test.list
  output:
    #pdf = "some.pdf"
    touch('plot_summary.done')
  log: 'log/plot_summary.done.log'
  threads: 12
  message:
    """ Running Rscript to plot the genome-wide distribution of coverage """
  shell:
    """
    Rscript scripts/plot_summary.R 2> {log}
    """

