rule ls_dedupBAM:
  input:
    #'depth/{sample}.depth.gz'
  output:
    'list/dedupBAM.list'
  log: 'log/dedupBAM.list.log'
  message: """--- Creating a sample list of deduplicated bam files ---"""
  shell:
    """
    ls deDup/*dedup.bam > {output} 2> {log}
    """

rule makePop_lists:
  input:#
  output:
    touch('list/makePOPlist.done')
  log: 'log/makePOPlist.log'
  message: """--- Creating population lists of deduplicated bam files ---"""
  shell:
    """
    scripts/makePOP_list.sh
    """

rule index_bam:
  input:
    'deDup/{sample}.bb.RAPID.dedup.bam'
  output:
    'deDup/{sample}.bb.RAPID.dedup.bam.bai'
  log: 'log/{sample}.index_bam.log'
  message: """ --- Indexing bam files --- """
  shell:
    """
    samtools index -b {input} 2> {log}
    """


rule bcftools_mpileup_call:
  input:
    ref = config['ref'],
    touched = 'list/makePOPlist.done'#,
    #bai = 'deDup/{sample}.bb.RAPID.dedup.bam.bai'
  output:
    #'vcf/daphnia_init_{chrom}_{species}.vcf.gz'
    'bcf/daphnia_init_{chrom}_{species}.bcf'
  log:
    'log/daphnia_init_{chrom}_{species}_bcf.log'
  threads: 36
  message: """--- Generate AllSites vcf/bcf using BCFtools (mpileup/call) ---"""
  shell:
    """
    bcftools mpileup -Ou -f {input.ref} -b list/{wildcards.species}.list -r {wildcards.chrom} -a AD,DP --threads {threads} | bcftools call -m -Ob -f GQ -o {output} --threads {threads} 2> {log}
    """


rule bcftools_mpileup_all:
  input:
    ref = config['ref'],
    bam_list = 'list/daphnia118.list'
  output:
    mpileup = 'bcf/all/daphnia_init_{chrom}.mpileup.bcf'
  log:
    'log/daphnia_init_{chrom}_mpileup.log'
  threads: 24
  message: """--- Generate AllSites bcf containing GLs with BCFtools mpileup ---"""
  shell:
    """
    bcftools mpileup -Ou -o {output.mpileup} -f {input.ref} -b {input.bam_list} -r {wildcards.chrom} -a AD,DP --threads {threads} 2> {log}
    """


rule bcftools_call_all:
  input:
    mpileup = 'bcf/all/daphnia_init_{chrom}.mpileup.bcf'
  output:
    calls = 'bcf/all/daphnia_init_{chrom}.call.bcf'
  log:
    'log/daphnia_init_{chrom}_bcf.log'
  threads: 24
  message: """--- Generate AllSites bcf using BCFtools call ---"""
  shell:
    """
    bcftools call -m -Ob -f GQ -o {output.calls} {input.mpileup} --threads {threads} 2> {log}
    """


rule list_chromVCFs:
  input:
    #'vcf/daphnia_init_{chrom}_{species}.vcf.gz'
  output:
    'list/VCF_{species}.list'
  message: """ --- Create a list of chromosome VCFs for each species --- """
  shell:
    """
      ls -v vcf/daphnia_init_*_{wildcards.species}.vcf.gz > {output}
    """

rule bcftools_concat:
  input:
    'list/VCF_{species}.list'
    #'vcf/daphnia_init_{chrom}_{species}.vcf.gz'
  output:
    #touch('vcf/concat/daphnia_concat_{species}.done')
    'vcf/concat/daphnia_{species}.vcf.gz'
  log:
    'log/daphnia_ConcatVCF_{species}.log'
  threads: 12
  message: """--- Concatenate chromosome VCFs into one VCF ---"""
  shell:
    """
    bcftools concat -f {input} -Oz -o {output} --threads {threads} 2> {log}
    """

rule bcftools_filltags:
  input:
    'vcf/concat/daphnia_{species}.vcf.gz'
  output:
    'vcf/concat/tags/daphnia_{species}_alltags.vcf.gz'
  log:
    'log/daphnia_{species}_alltags_vcf.log'
  message: """--- Add all tags ---"""
  shell:
    """
    bcftools +fill-tags {input} -Oz -o {output} -- -t all 2> {log}
    """

#rule bcftools_merge:
#  input:
#    #'vcf/concat/tags/daphnia_{species}_alltags.vcf.gz'
#  output:
#    'vcf/concat/tags/daphnia.vcf.gz'
#  log:
#    'log/daphnia_{species}_alltags_vcf.log'
#  message: """--- Add all tags ---"""
#  shell:
#    """
#    bcftools merge daphnia_*_alltags.vcf.gz -Oz -o {output} 2> {log}
#    """


rule vcf_randomsample:
  input:
    'vcf/concat/tags/daphnia_{species}_alltags.vcf.gz'
  output:
    'vcf/concat/tags/daphnia_{species}_alltags_rndsubset.vcf.gz'
  log: 'log/daphnia_{species}_alltags_rndsubset.log'
  threads: 24
  message: """--- Randomly subsample vcf ---"""
  shell:
    """
    bcftools view {input} | vcfrandomsample -r 0.0007 | bgzip -c > {output} 2> {log}
    """


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


