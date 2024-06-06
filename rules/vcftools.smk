rule count_variants:
  input:
    'vcf/daphnia_{species}.vcf.gz'
  output:
    'vcf/stats/daphnia_{species}_nbrSites.txt'
  log: 'log/daphnia_{species}_nbrSites.log'
  threads: 12
  message: """--- Count variants ---"""
  shell:
    """
    bcftools view -H {input} | wc -l > {output} 2> {log}
    """


rule vcf_randomsample:
  input:
    'vcf/daphnia_{species}.vcf.gz'
  output:
    'vcf/subset/daphnia_{species}_100Ksubset.vcf.gz'
  log: 'log/daphnia_{species}_100Ksubset.log'
  threads: 12
  message: """--- Count variants and randomly subsample vcf ---"""
  shell:
    """
    bcftools view {input} | vcfrandomsample -r 0.0007 | bgzip -c > {output} 2> {log}
    """


rule allel_frequency:
  input: 
    'vcf/subset/daphnia_{species}_100Ksubset.vcf.gz'
  output:
    touch('vcf/stats/daphnia_{species}_100Ksubset.frq.done')
  log: 'log/stats/daphnia_{species}_100Ksubset.frq.log'
  threads: 12
  message: """--- Calculate allele frequency ---"""
  shell:
    """
    vcftools --gzvcf {input} --freq2 --out vcf/stats/daphnia_{wildcards.species}_100Ksubset --max-alleles 2 2> {log}
    """


rule mean_depth_individual:
  input:
    'vcf/subset/daphnia_{species}_100Ksubset.vcf.gz'
  output:
    touch('vcf/stats/daphnia_{species}_100Ksubset.idepth.done')
  log: 'log/stats/daphnia_{species}_100Ksubset.idepth.log'
  threads: 12
  message: """--- Calculate mean depth per individual ---"""
  shell:
    """
    vcftools --gzvcf {input} --depth --out vcf/stats/daphnia_{wildcards.species}_100Ksubset 2> {log}
    """


rule mean_depth_site:
  input:
    'vcf/subset/daphnia_{species}_100Ksubset.vcf.gz'
  output:
    touch('vcf/stats/daphnia_{species}_100Ksubset.ldepth.mean.done')
  log: 'log/stats/daphnia_{species}_100Ksubset.ldepth.mean.log'
  threads: 12
  message: """--- Calculate mean depth per site ---"""
  shell:
    """
    vcftools --gzvcf {input} --site-mean-depth --out vcf/stats/daphnia_{wildcards.species}_100Ksubset 2> {log}

    """

rule site_quality:
  input:
    'vcf/subset/daphnia_{species}_100Ksubset.vcf.gz'
  output:
    touch('vcf/stats/daphnia_{species}_100Ksubset.lqual.done')
  log: 'log/stats/daphnia_{species}_100Ksubset.lqual.log'
  threads: 12
  message: """--- Calculate site quality  ---"""
  shell:
    """
    vcftools --gzvcf {input} --site-quality --out vcf/stats/daphnia_{wildcards.species}_100Ksubset 2> {log}
    """


rule missingdata_individual:
  input:
    'vcf/subset/daphnia_{species}_100Ksubset.vcf.gz'
  output:
    touch('vcf/stats/daphnia_{species}_100Ksubset.imiss.done')
  log: 'log/stats/daphnia_{species}_100Ksubset.imiss.log'
  threads: 12
  message: """--- Calculate proportion of missing data per individual ---"""
  shell:
    """
    vcftools --gzvcf {input} --missing-indv --out vcf/stats/daphnia_{wildcards.species}_100Ksubset 2> {log}
    """



rule missingdata_site:
  input:
    'vcf/subset/daphnia_{species}_100Ksubset.vcf.gz'
  output:
    touch('vcf/stats/daphnia_{species}_100Ksubset.lmiss.done')
  log: 'log/stats/daphnia_{species}_100Ksubset.lmiss.log'
  threads: 12
  message: """--- Calculate proportion of missing data per site ---"""
  shell:
    """
    vcftools --gzvcf {input} --missing-site --out vcf/stats/daphnia_{wildcards.species}_100Ksubset 2> {log}
    """


rule print_MQ_DP:
  input:
    'vcf/subset/daphnia_{species}_100Ksubset.vcf.gz'
  output:
    'vcf/stats/daphnia_{species}_100Ksubset.MQ_DP.txt'
  log: 'log/stats/daphnia_{species}_100Ksubset.MQ_DP.log'
  threads: 12
  message: """--- For each site print MQ, DP values ---"""
  shell:
    """
    bcftools query {input} -f'%MQ\t%DP\n' > {output} 2> {log}
    """
