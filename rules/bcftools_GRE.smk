rule bcftools_mpileup_call_GRE14:
  input:
    ref = config['ref'],
    touched = 'list/makePOPlist.done',
    #bam = 'deDup/GRE1408.bb.RAPID.dedup.bam',
    #bai = 'deDup/GRE1408.bb.RAPID.dedup.bam.bai'
  output:
    'bcf/GRE14/GRE1408_{chrom}.bcf'
  log:
    'log/GRE1408_{chrom}.log'
  threads: 24
  message: """--- Generate AllSites vcf/bcf using BCFtools (mpileup/call) ---"""
  shell:
    """
    bcftools mpileup -Ou -f {input.ref} deDup/GRE1408.bb.RAPID.dedup.bam -r {wildcards.chrom} -a AD,DP --threads {threads} | bcftools call -m -Ob -f GQ -o {output} --threads {threads} 2> {log}
    """


rule list_chromBCFs_GRE:
  input:
    #'bcf/GRE14/GRE1408_{chrom}.bcf'
  output:
    'list/BCF_GRE1408.list'
  log:
    'log/BCF_GRE1408.log'
  message: """ --- Create a list of chromosome VCFs for each species --- """
  shell:
    """
    ls -v bcf/GRE14/GRE1408_*.bcf > list/BCF_GRE1408.list
    """


rule bcftools_concat_GRE:
  input:
    'list/BCF_GRE1408.list'
  output:
    bcf = 'bcf/concat/GRE1408.bcf',
    csi = 'bcf/concat/GRE1408.bcf.csi'
  log:
    'log/GRE1408_ConcatBCF.log'
  threads: 12
  message: """--- Concatenate all chromosome-BCFs into one single BCF file and add all tags ---"""
  shell:
    """
    bcftools concat -f {input} | bcftools +fill-tags -Ob -o {output.bcf} --threads {threads} &&
    bcftools index {output.bcf} {output.csi} 2> {log}
    """

rule bcftools_merge_galeata:
  input:
    #'bcf/daphnia_{species}.vcf.gz',
  output:
    'bcf/concat/daphnia_galeata2.bcf'
  log:
    'log/daphnia_galeata2.log'
  threads: 12
  message: """--- Merge bcf galeata ---"""
  shell:
    """
    DIR=(bcf/concat)
    bcftools merge $DIR/daphnia_galeata.bcf $DIR/GRE1408.bcf | bcftools +fill-tags -Ob -o {output} -- -t all 2> {log}
    """
