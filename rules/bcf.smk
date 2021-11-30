rule bcf_mpileup_call:
  input:
    ref = config['ref'],
    touched = 'list/makePOPlist.done'
  output:
    'bcf/daphnia_init_{chrom}_{species}.bcf'
  log:
    'log/daphnia_init_{chrom}_{species}_bcf.log'
  threads: 36
  message: """--- Generate AllSites vcf using BCFtools (mpileup/call) ---"""
  shell:
    """
    bcftools mpileup -f {input.ref} -b list/{wildcards.species}.list -r {wildcards.chrom} | bcftools call -m -Ou -f GQ -o {output} --threads {threads} 2> {log}
    """

rule index_bcf:
  input:
    'bcf/daphnia_init_{chrom}_{species}.bcf'
  output:
    'bcf/daphnia_init_{chrom}_{species}.bcf.csi'
  log:
    'log/daphnia_init_{chrom}_{species}_indexBCF.log'
  shell:
    """
    bcftools index {input} {output} 2> {log}
    """

rule bcf2vcf:
  input:
    bcf = 'bcf/daphnia_init_{chrom}_{species}.bcf',
    idx = 'bcf/daphnia_init_{chrom}_{species}.bcf.csi'
  output:
    'vcf/daphnia_init_{chrom}_{species}.vcf.gz'
  log:
    'log/daphnia_init_{chrom}_{species}_VCF.log'
  shell:
    """
    bcftools convert --threads {threads} -Oz -o {output} {input.bcf} 2> {log}
    """

