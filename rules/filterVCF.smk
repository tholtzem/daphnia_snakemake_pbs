rule tabix:
  input:
    'vcf/concat/tags/daphnia_{species}_alltags_Srenamed.vcf.gz'
  output:
    'vcf/concat/tags/daphnia_{species}_alltags_Srenamed.vcf.gz.tbi'
  log:
    "log/tabix_daphnia_{species}_alltags_Srenamed.vcf.log"
  shell:
    """
    tabix -p vcf {input}
    """

rule HardFilterVCF_Allsites:
  input: 
    vcf = 'vcf/concat/tags/daphnia_{species}_alltags_Srenamed.vcf.gz',
    tabix = 'vcf/concat/tags/daphnia_{species}_alltags_Srenamed.vcf.gz.tbi'
  output:
    'vcf/filtered/daphnia_{species}_Allsites_hardfiltered_1.vcf.gz'
  log: 'log/hardfiltered_Allsites_{species}.log'
  threads: 2
  message: """--- Hard filter AllSites VCF (GATKs best practices) ---"""
  shell:
    """
    bcftools filter {input.vcf} -e 'FS>60.0 && MQ<40' --SnpGap 5 -Oz -o {output} 2> {log}
    """


rule further_filtering:
  input: 
    'vcf/filtered/daphnia_{species}_Allsites_hardfiltered_1.vcf.gz'
  output:
    'vcf/filtered/daphnia_{species}_Allsites_hardfiltered_2.vcf.gz'
  log: 'log/hardfiltered_Allsites_{species}.log'
  threads: 2
  message: """--- Further filtering of AllSites VCF ---"""
  shell:
    """
    vcftools --gzvcf {input} --remove-indels --max-missing 0.8 --min-meanDP 20 --max-meanDP 60 --recode --stdout | gzip -c > {output}
    """



rule filterVCF_Invariant:
  input: 
    'vcf/filtered/daphnia_{species}_Allsites_hardfiltered_2.vcf.gz'
  output:
    vcf = 'vcf/filtered/daphnia_{species}_invariant.vcf.gz',
    tabix = 'vcf/filtered/daphnia_{species}_invariant.vcf.gz.tbi'
  log: 'log/{species}_invariantVCF.log'
  threads: 2
  message: """--- Create a filtered VCF containing only invariant sites ---"""
  shell:
    """
    vcftools --gzvcf {input} --max-maf 0 --recode --stdout | gzip -c > {output.vcf} && tabix -p vcf {output.tabix} 2> {log}
    """


rule filterVCF_Variant:
  input: 
    'vcf/filtered/daphnia_{species}_Allsites_hardfiltered_2.vcf.gz'
  output:
    vcf = 'vcf/filtered/daphnia_{species}_variant.vcf.gz',
    tabix = 'vcf/filtered/daphnia_{species}_variant.vcf.gz.tbi'
  log: 'log/{species}_variantVCF.log'
  threads: 2
  message: """--- Create a filtered VCF containing only variant sites ---"""
  shell:
    """
    vcftools --gzvcf {input} --mac 1 --minQ 30 --minDP 20 --maxDP 60 --minGQ 30  --recode --stdout | gzip -c > {output.vcf} && tabix -p vcf {output.tabix} 2> {log}
    """


rule finalVCF:
  input:
    invariant = 'vcf/filtered/daphnia_{species}_invariant.vcf.gz',
    tb1 = 'vcf/filtered/daphnia_{species}_invariant.vcf.gz.tbi',
    variant = 'vcf/filtered/daphnia_{species}_variant.vcf.gz',
    tb2 = 'vcf/filtered/daphnia_{species}_variant.vcf.gz.tbi'
  output:
    'vcf/filtered/daphnia_{species}_final.vcf.gz'
  log: 'log/daphnia_{species}_final.vcf.log'
  threads: 12
  message: """ --- Combine variantVCF and invariantVCF --- """
  shell:
    """
    bcftools concat --allow-overlaps {input.invariant} {input.variant} -Oz -o {output} 2> {log}
    """




