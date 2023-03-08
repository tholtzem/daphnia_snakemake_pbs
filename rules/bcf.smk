
#rule bcftools_prefilter:
#  input:
#    'vcf/concat/daphnia_{species}_alltags.vcf.gz'
#  output:
#    'vcf/concat/daphnia_{species}_prefiltered.vcf.gz'
#  log:
#    'log/daphnia_{species}_prefiltered.log'
#  message: """ --- Prefilter the allSites VCF using bcftools, set genotypes according to filtering criteria and tabix output vcf --- """
#  shell:
#    """
#    bcftools view --exclude-types indels --max-alleles 2 {input}  | bcftools +setGT - -- -n . -t q -e 'FORMAT/DP>=10&(GQ>=30|RGQ>=30)' | bgzip -c > {output} && tabix -p vcf {output}
#    """
