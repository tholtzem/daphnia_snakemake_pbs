rule list_speciesVcf:
  input:
    #'vcf/daphnia_init_{chrom}_{species}.vcf.gz'
  output:
    'vcf/concat/daphnia_{species}.vcf.list'
  shell:
    """
      ls vcf/daphnia_init_*_{wildcards.species}.vcf.gz
    """
