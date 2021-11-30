rule tabix:
  input:
    "bcftools/daphnia_init_dgal39_cucullata.vcf.gz"
  output:
    "bcftools/daphnia_init_dgal39_cucullata.vcf.gz.tbi"
  log:
    "log/tabix/daphnia_init_dgal39_cucullata.log"
  shell:
    """
    tabix -p vcf {input}
    """
