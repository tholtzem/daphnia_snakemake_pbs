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
    'bcf/daphnia_init_{chrom}_{species}.bcf'
  log:
    'log/daphnia_init_{chrom}_{species}_bcf.log'
  threads: 24
  message: """--- Generate AllSites vcf/bcf using BCFtools (mpileup/call) ---"""
  shell:
    """
    bcftools mpileup -Ou -f {input.ref} -b list/{wildcards.species}.list -r {wildcards.chrom} -a AD,DP --threads {threads} | bcftools call -m -Ob -f GQ -o {output} --threads {threads} 2> {log}
    """



rule list_chromBCFs:
  input:
    #'bcf/daphnia_init_{chrom}_{species}.bcf'
  output:
    'list/BCF_{species}.list'
  log:
    'log/BCF_{species}.list.log'
  message: """ --- Create a list of chromosome VCFs for each species --- """
  shell:
    """
      ls -v bcf/daphnia_init_*_{wildcards.species}.bcf > {output}
    """


rule bcftools_concat:
  input:
    'list/BCF_{species}.list'
  output:
    bcf = 'bcf/concat/daphnia_{species}.bcf',
    csi = 'bcf/concat/daphnia_{species}.bcf.csi'
  log:
    'log/daphnia_ConcatBCF_{species}.log'
  threads: 12
  message: """--- Concatenate all chromosome-BCFs into one single BCF file and add all tags ---"""
  shell:
    """
    bcftools concat -f {input} | bcftools +fill-tags -Ob -o {output.bcf} --threads {threads} &&
    bcftools index {output.bcf} {output.csi} 2> {log}
    """


#rule bcftools_merge:
#  input:
#    #bcf = 'bcf/concat/daphnia_{species}.bcf',
#  output:
#    'bcf/merged/daphnia_9pops.bcf'
#  log:
#    'log/daphnia_9pops.log'
#  threads: 12
#  message: """--- Add all tags ---"""
#  shell:
#    """
#    DIR=(bcf/concat)
#    bcftools merge $DIR/daphnia_cucullata.bcf $DIR/daphnia_dentifera.bcf $DIR/daphnia_galeata.bcf $DIR/daphnia_lacustris.bcf $DIR/daphnia_longispina.bcf $DIR/daphnia_longispinaFIN.bcf $DIR/daphnia_mendotae.bcf $DIR/daphnia_umbra.bcf $DIR/daphnia_zschokkei.bcf -Ob -o {output} 2> {log}
#    """


rule bcftools_merge_10pops:
  input:
    #'bcf/daphnia_{species}.vcf.gz',
  output:
    'bcf/merged/daphnia_10pops.bcf'
  log:
    'log/daphnia_10pops.log'
  threads: 12
  message: """--- Merge ---"""
  shell:
    """
    DIR=(bcf/concat)
    bcftools merge $DIR/daphnia_cucullata.bcf $DIR/daphnia_curvirostris.bcf $DIR/daphnia_dentifera.bcf $DIR/daphnia_galeata.bcf $DIR/GRE1408.bcf $DIR/daphnia_lacustris.bcf $DIR/daphnia_longispina.bcf $DIR/daphnia_longispinaFIN.bcf $DIR/daphnia_mendotae.bcf $DIR/daphnia_umbra.bcf $DIR/daphnia_zschokkei.bcf | bcftools +fill-tags -Ob -o {output} -- -t all 2> {log}
    """


rule bcf2vcf_10pops:
  input:
    'bcf/merged/daphnia_10pops.bcf'
  output:
    'vcf/daphnia_10pops.vcf.gz'
  log:
    'log/daphnia_10pops_BCF2VCF.log'
  message: """ --- Convert bcf2vcf --- """
  threads: 12
  shell:
    """
    bcftools index -f --csi {input} --threads {threads} &&
    bcftools convert --threads {threads} -Oz -o {output} {input} 2> {log}
    """


rule bcf2vcf:
  input:
    'bcf/concat/daphnia_{species}.bcf'
  output:
    'vcf/daphnia_{species}.vcf.gz'
  log:
   'log/daphnia_{species}_BCF2VCF.log'
  message: """ --- Convert bcf2vcf --- """
  threads: 12
  shell:
    """
    bcftools index -f --csi {input} --threads {threads} &&
    bcftools convert --threads {threads} -Oz -o {output} {input} 2> {log}
    """
