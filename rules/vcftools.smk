rule filter_sites:
  input: 
    vcf = 'bcftools/daphnia_init_{chrom}_{species}.vcf.gz'
  output:
    'bcftools/{species}/{species}_{chrom}_filtered.vcf.gz'
  log: 'log/{species}_{chrom}_filtered.log'
  threads: 24
  message: """--- Generate AllSites vcf using BCFtools (mpileup/call) ---"""
  shell:
    """
    vcftools --gzvcf my_vcf.vcf.gz --remove-indels --max-missing 0.8 --minQ 30 --min-meanDP 20 --max-meanDP 500 --recode --stdout | gzip -c > {output}
    --threads {threads} 2> {log}
    """
