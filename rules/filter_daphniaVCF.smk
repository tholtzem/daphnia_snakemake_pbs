rule bcftools_mpileup_call_all:
  input:
    ref = config['ref'],
    samples = 'list/daphnia.list'
  output:
    'vcf/all/daphnia_{chrom}.vcf.gz'
  log:
    'log/daphnia_{chrom}.vcf.log'
  threads: 36
  message: """--- Generate AllSites vcf using BCFtools (mpileup/call) ---"""
  shell:
    """
    bcftools mpileup -f {input.ref} -b {input.samples} -r {wildcards.chrom} -a AD,DP,SP | bcftools call -m -Oz -f GQ -o {output} --threads {threads} 2> {log}
    """

#rule bcftools_merge:
#  input:
#    vcf = 'vcf/concat/tags/daphnia_{species}_alltags.vcf.gz',
#    tbi = 'vcf/concat/tags/daphnia_{species}_alltags.vcf.gz.tbi'
#  output:
#    'vcf/allpops/daphnia.vcf.gz'
#  log:
#    "log/daphnia.vcf.log"
#  shell:
#    """
#    bcftools merge -p vcf {input}
#    """

#rule HardFilterVCF_Allsites:
#  input: 
#    vcf = 'vcf/concat/tags/daphnia_{species}_alltags.vcf.gz',
#    tabix = 'vcf/concat/tags/daphnia_{species}_alltags.vcf.gz.tbi'
#  output:
#   'vcf/filtered/daphnia_{species}_Allsites_hardfiltered.vcf.gz'
#  log: 'log/hardfiltered_Allsites_{species}.log'
#  threads: 2
#  message: """--- Hard filter AllSites VCF (GATKs best practices) ---"""
#  shell:
#    """
#    bcftools filter -e 'FS>60.0 || MQ<40' -Oz -o {output} {input} 2> {log}
#    """


#rule filterVCF_Invariant:
#  input: 
#    'vcf/filtered/daphnia_{species}_Allsites_hardfiltered.vcf.gz'
#  output:
#    vcf = 'vcf/filtered/daphnia_{species}_invariant.vcf.gz',
#    tabix = 'vcf/filtered/daphnia_{species}_invariant.vcf.gz.tbi'
#  log: 'log/{species}_invariantVCF.log'
#  threads: 2
#  message: """--- Create a filtered VCF containing only invariant sites ---"""
#  shell:
#    """
#    vcftools --gzvcf {input} --max-maf 0 --remove-indels --max-missing 0.8 --min-meanDP 20 --max-meanDP 60 --minDP 20 --maxDP 60  --recode --stdout | gzip -c > {output.vcf} && tabix -p vcf {output.tabix} 2> {log}
#    """


#rule filterVCF_Variant:
#  input: 
#    'vcf/filtered/daphnia_{species}_Allsites_hardfiltered.vcf.gz'
#  output:
#    vcf = 'vcf/filtered/daphnia_{species}_variant.vcf.gz',
#    tabix = 'vcf/filtered/daphnia_{species}_variant.vcf.gz.tbi'
#  log: 'log/{species}_variantVCF.log'
#  threads: 2
#  message: """--- Create a filtered VCF containing only variant sites ---"""
#  shell:
#    """
#    vcftools --gzvcf {input} --mac 1 --remove-indels --max-missing 0.8 --minQ 30 --min-meanDP 20 --max-meanDP 60 --minDP 20 --maxDP 60 --minGQ 30  --recode --stdout | gzip -c > {output.vcf} && tabix -p vcf {output.tabix} 2> {log}
#    """

#rule finalVCF:
#  input:
#    invariant = 'vcf/filtered/daphnia_{species}_invariant.vcf.gz',
#    tb1 = 'vcf/filtered/daphnia_{species}_invariant.vcf.gz.tbi',
#   variant = 'vcf/filtered/daphnia_{species}_variant.vcf.gz',
#    tb2 = 'vcf/filtered/daphnia_{species}_variant.vcf.gz.tbi'
#  output:
#    'vcf/filtered/daphnia_{species}_final.vcf.gz'
#  log: 'log/daphnia_{species}_final.vcf.log'
#  threads: 12
#  message: """ --- Combine variantVCF and invariantVCF --- """
#  shell:
#    """
#    bcftools concat --allow-overlaps {input.invariant} {input.variant} -Oz -o {output} 2> {log}
#    """
