localrules: list_chromBCFs

rule bamIndex:
  input:
    '/scratch/c7701178/mach2/DAPHNIA/daphnia_snakemake_pbs/realigned/HiC/minid95/{sample}.bb.HiC.realigned.bam'
  output:
    '/scratch/c7701178/mach2/DAPHNIA/daphnia_snakemake_pbs/realigned/HiC/minid95/{sample}.bb.HiC.realigned.bai'
  threads: 1
  log: "log/HiC/minid95/indexBam_{sample}.bb.HiC.realigned.log"
  message: """--- Indexing with samtools ---"""
  shell:
    """
    samtools index {input} {output}
    """


rule bcftools_mpileup_call:
  input:
    ref = config['ref_HiC'],
    bams = 'list/pops/realignedBAM_10pops.list'
    #touched = 'list/makePOPlist.done'#,
    #bai = 'deDup/{sample}.bb.RAPID.dedup.bam.bai'
  output:
    'bcf/daphnia_init_10pops_{chrom}_new.bcf'
  log:
    'log/daphnia_init_10pops10pops_{chrom}_new.bcf.log'
  threads: 1
  resources: mem_mb=800, walltime="72:00:00"
  message: """--- Generate AllSites vcf/bcf using BCFtools (mpileup/call) ---"""
  shell:
    """
    bcftools mpileup -Ou -f {input.ref} -b {input.bams} -r {wildcards.chrom} -a AD,DP,SP --threads {threads} | bcftools call -m -Ob -f GQ -o {output} --threads {threads} 2> {log}
    """


rule list_chromBCFs:
  input:
    #'bcf/daphnia_init_10pops_{chrom}.bcf'
  output:
    'list/daphnia_init_10pops_BCFchrom_new.list'
  log:
    'log/daphnia_init_10pops_BCFchrom_new.list.log'
  message: """ --- Create a list of chromosome VCFs for each species --- """
  shell:
    """
    ls -v bcf/daphnia_init_10pops_HiC_scaffold_*_new.bcf  > {output}
    """


rule bcftools_concat:
  input:
    'list/daphnia_init_10pops_BCFchrom_new.list'
  output:
    bcf = 'bcf/concat/daphnia_init_10pops.bcf',
    csi = 'bcf/concat/daphnia_init_10pops.bcf.csi'
  log:
    'log/daphnia_init_10pops_ConcatBCF.log'
  threads: 2
  resources: mem_mb=200, walltime="1:30:00"
  message: """--- Concatenate all chromosome-BCFs into one single BCF file and add all tags ---"""
  shell:
    """
    bcftools concat -f {input} | bcftools +fill-tags -Ob -o {output.bcf} --threads {threads} &&
    bcftools index {output.bcf} {output.csi} 2> {log}
    """


#rule bcftools_merge_10pops:
#  input:
#    #'bcf/daphnia_{species}.vcf.gz',
#  output:
#    'bcf/merged/daphnia_10pops.bcf'
#  log:
#    'log/daphnia_10pops.log'
#  threads: 12
#  message: """--- Merge ---"""
#  shell:
#    """
#    DIR=(bcf/concat)
#    bcftools merge $DIR/daphnia_cucullata.bcf $DIR/daphnia_curvirostris.bcf $DIR/daphnia_dentifera.bcf $DIR/daphnia_galeata.bcf $DIR/GRE1408.bcf $DIR/daphnia_lacustris.bcf $DIR/daphnia_longispina.bcf $DIR/daphnia_longispinaFIN.bcf $DIR/daphnia_mendotae.bcf $DIR/daphnia_umbra.bcf $DIR/daphnia_zschokkei.bcf | bcftools +fill-tags -Ob -o {output} -- -t all 2> {log}
#    """


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
