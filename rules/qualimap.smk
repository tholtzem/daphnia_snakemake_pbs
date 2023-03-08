rule qualimap:
  input:
    deDup = "deDup/HiC/minid95/{sample}.bb.HiC.dedup.bam"
  output:
    directory('qmap/{sample}.bb.HiC.dedup')
  log: 'log/HiC/minid95/{sample}.qmap.log'
  message: """--- Qualimap bams ---"""
  threads: 12 
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 qualimap bamqc -bam {input.deDup} -outdir {output} --java-mem-size=100g 2> {log}
    """


rule qFilter_bams:
  input:
    bam = 'deDup/{sample}.bb.RAPID.dedup.bam'
  output:
    bam = 'qmap/{sample}.bb.RAPID.dedup.q{quality}.bam'
  log: 'log/{sample}.{quality}.bam.log'
  message: """--- Filter bams with quality treshold ---"""
  threads: 12 
  shell:
    """
    samtools view -F 4 -Shu -q {wildcards.quality} {input.bam} | samtools sort - -o {output.bam} 2> {log}
    """


rule qualimap_q_bams:
  input:
    bam = 'qmap/{sample}.bb.RAPID.dedup.q{quality}.bam'
  output:
    touch('qmap/{sample}.q{quality}qmap.done')
  log: 'log/{sample}.q{quality}.qmap.log'
  message: """--- Qualimap filtered bams ---"""
  threads: 12 
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 qualimap bamqc -bam {input.bam} -outdir qmap/{wildcards.sample}/q{wildcards.quality}/ --java-mem-size=100g 2> {log}
    """
