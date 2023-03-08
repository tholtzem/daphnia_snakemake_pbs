rule reheader_VCF:
  input:
    vcf = 'vcf/daphnia_9pops.vcf.gz',
    samplelist = 'list/vcf_daphnia_9pops.list' 
  output:
    'vcf/daphnia_9pops.reheaded.vcf.gz'
  log: 'log/daphnia_9pops_reheaded.log'
  threads: 12
  message: """--- Change sample names with ambiguous names in header, change names to capital letters only --- """
  shell:
    """
    bcftools view {input.vcf} -Oz | bcftools reheader -s {input.samplelist} -o {output} 2> {log}
    """


rule sort_VCF:
  input:
    'vcf/daphnia_9pops.reheaded.vcf.gz'
  output:
    samples = 'list/vcf_daphnia_9pops_samples_sorted.list',
    vcf = 'vcf/daphnia_9pops.reheaded.sorted.vcf.gz'
  log: 'log/daphnia_9pops.reheaded.sorted.log'
  threads: 12
  message: """--- Sort the order of the samples in the VCF file alphabetically (handy for downstream-analysis) ---"""
  shell:
    """
    bcftools query -l {input} | sort > {output.samples} &&
    bcftools view -S {output.samples} {input} -Oz -o {output.vcf} 2> {log}
    """
