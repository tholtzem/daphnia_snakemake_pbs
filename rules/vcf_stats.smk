rule allel_frequency:
  input: 
    #'vcf/concat/daphnia_{species}_alltags.vcf.gz'
    'vcf/concat/tags/daphnia_{species}_alltags_rndsubset.vcf.gz'
  output:
    touch('vcf/concat/stats/daphnia_{species}_alltags_rndsubset.frq.done')
  log: 'log/daphnia_{species}_alltags.frq.log'
  threads: 24
  message: """--- Calculate allele frequency ---"""
  shell:
    """
    vcftools --gzvcf {input} --freq2 --out vcf/concat/stats/daphnia_{wildcards.species}_alltags_rndsubset --max-alleles 2 2> {log}
    """

rule mean_depth_individual:
  input:
    #'vcf/concat/daphnia_{species}_alltags.vcf.gz'
    'vcf/concat/tags/daphnia_{species}_alltags_rndsubset.vcf.gz'
  output:
    touch('vcf/concat/stats/daphnia_{species}_alltags_rndsubset.idepth.done')
  log: 'log/daphnia_{species}_alltags.idepth.log'
  threads: 24
  message: """--- Calculate mean depth per individual ---"""
  shell:
    """
    vcftools --gzvcf {input} --depth --out vcf/concat/stats/daphnia_{wildcards.species}_alltags_rndsubset 2> {log}
    """

rule mean_depth_site:
  input:
    #'vcf/concat/daphnia_{species}_alltags.vcf.gz'
    'vcf/concat/tags/daphnia_{species}_alltags_rndsubset.vcf.gz'
  output:
    touch('vcf/concat/stats/daphnia_{species}_alltags_rndsubset.ldepth.mean.done')
  log: 'log/daphnia_{species}_alltags.ldepth.mean.log'
  threads: 24
  message: """--- Calculate mean depth per site ---"""
  shell:
    """
    vcftools --gzvcf {input} --site-mean-depth --out vcf/concat/stats/daphnia_{wildcards.species}_alltags_rndsubset 2> {log}

    """

rule site_quality:
  input:
    #'vcf/concat/daphnia_{species}_alltags.vcf.gz'
    'vcf/concat/tags/daphnia_{species}_alltags_rndsubset.vcf.gz'
  output:
    touch('vcf/concat/stats/daphnia_{species}_alltags_rndsubset.lqual.done')
  log: 'log/daphnia_{species}_alltags.lqual.log'
  threads: 24
  message: """--- Calculate site quality  ---"""
  shell:
    """
    vcftools --gzvcf {input} --site-quality --out vcf/concat/stats/daphnia_{wildcards.species}_alltags_rndsubset 2> {log}
    """

rule missingdata_individual:
  input:
    #'vcf/concat/daphnia_{species}_alltags.vcf.gz'
    'vcf/concat/tags/daphnia_{species}_alltags_rndsubset.vcf.gz'
  output:
    touch('vcf/concat/stats/daphnia_{species}_alltags_rndsubset.imiss.done')
  log: 'log/daphnia_{species}_alltags.imiss.log'
  threads: 24
  message: """--- Calculate proportion of missing data per individual ---"""
  shell:
    """
    vcftools --gzvcf {input} --missing-indv --out vcf/concat/stats/daphnia_{wildcards.species}_alltags_rndsubset 2> {log}
    """

rule missingdata_site:
  input:
    #'vcf/concat/daphnia_{species}_alltags.vcf.gz'
    'vcf/concat/tags/daphnia_{species}_alltags_rndsubset.vcf.gz'
  output:
    touch('vcf/concat/stats/daphnia_{species}_alltags_rndsubset.lmiss.done')
  log: 'log/daphnia_{species}_alltags.lmiss.log'
  threads: 24
  message: """--- Calculate proportion of missing data per site ---"""
  shell:
    """
    vcftools --gzvcf {input} --missing-site --out vcf/concat/stats/daphnia_{wildcards.species}_alltags_rndsubset 2> {log}
    """
