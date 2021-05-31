
rule genome_coverage_bed:
  input:
    'deDup/{ID}.bb.RAPID.dedup.bam'
  output:
    'bedtools/{ID}.bb.RAPID.dedup.genomecov.bed'
  threads: 12
  message:
    """ Computes BED summaries using bedtools """
  shell:
    """
    bedtools genomecov -ibam {input} > {output}
    """

rule plot_gencov:
  input:
    bed= "bedtools/{ID}.bb.RAPID.dedup.genomecov.bed"
  output:
    pdf = "bedtools/plots/{ID}.bb.RAPID.pdf"
  threads: 4
  message:
    """ Running Rscript to plot the genome-wide distribution of coverage """
  shell:
    """
    Rscript scripts/plot_gene_covs.R {input.bed} {output.pdf}
    """
    

