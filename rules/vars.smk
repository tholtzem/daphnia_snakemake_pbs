rule callVars1:
  input:
    #sam = lambda wildcards: getAllSams(wildcards.allTAs),
    ref = config['ref'],
    id_list = ""
  output:
    vcf = "vars/ants_initial.vcf"
  threads: 24
  message: """--- Calling variants (bbtools) for T. alpestre samples ---"""
  shell: 
    """
    callvariants.sh t={threads} list={input.id_list ref={input.ref} ploidy=2 multisample out={output.vcf}
    """


rule qualCalc:
  input:
    bam = "bbmap/{sample}_trm.bam",
    vcf = "vars/ants_initial.vcf"
    ref = config['ref']
  output:
    "placeholder_for_output"
  threads: 12
  message: """--- Calculating true quality with bbtools ---"""
  shell:
    """
    calctruequality.sh t={threads} vcf={input.vcf} ref={input.ref} id= [[a,b,c]] 
    """


rule bbdukQualRecal:
  input:
    bam = bbmap/{sample}.bam'

