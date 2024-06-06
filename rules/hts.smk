ruleorder: bamIndex > bamIndex_mrgd
ruleorder: clip_overlap > ls_ClipBam > list_indels
localrules: ls_ClipBam, plot_summary

#rule fastqc_raw:
#  input:
#    r1 = lambda wildcards: getFqHome(wildcards.sample)[0],
#    r2 = lambda wildcards: getFqHome(wildcards.sample)[1]
#  output:
#    "raw/qc/fastqc/{sample}_R1_001_fastqc.html",
#    "raw/qc/fastqc/{sample}_R1_001_fastqc.zip",
#    "raw/qc/fastqc/{sample}_R2_001_fastqc.html",
#    "raw/qc/fastqc/{sample}_R2_001_fastqc.zip"
#  log: "log/{sample}.qc.log"
#  threads: 12
#  message: """ Quality check of raw data with FastQC before trimming. """
#  shell:
#    """
#    fastqc -o raw/qc/fastqc/ -f fastq {input.r1} & fastqc -o raw/qc/fastqc/ -f fastq {input.r2}
#    """
#
#
#rule multiqc_raw:
#  input:
#  output:
#    "raw/qc/multiqc/"
#  log: "log/raw_multiqc_report.log"
#  message: """ Multiqc of raw reads """
#  shell:
#    """
#    multiqc raw/qc/fastqc/ -o {output}
#    """
#
#
#rule bbduk_trm:
#  input:
#    r1 = lambda wildcards: getFqHome(wildcards.sample)[0],
#    r2 = lambda wildcards: getFqHome(wildcards.sample)[1],
#    adapters = config["adapters"]
#  output:
#    r1 = "trm/{sample}_R1.trmd.fq.gz",
#    r2 = "trm/{sample}_R2.trmd.fq.gz",
#  log:
#    'log/bbduk_trm_{sample}.log'
#  threads: 12
#  message: """ --- High sensitivity adapter trimming and polyG removal and minlength --- """
#  shell:
#    """
#    /apps/uibk/bin/sysconfcpus -n 12 bbduk.sh -Xms110g -Xmx110g in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} ref={input.adapters} trimpolygright=10 qtrim=rl trimq=10 ordered=t ktrim=r k=21 mink=10 hdist=2 hdist2=1 minlength=40 stats=trm/stats/{wildcards.sample}.stats1 refstats=trm/stats/{wildcards.sample}.refstats1 bhist=trm/hist/{wildcards.sample}.bhist1 qhist=trm/hist/{wildcards.sample}.qhist1 lhist=trm/hist/{wildcards.sample}.lhist1 tpe tbo 2> {log}
#    """
#
#
#rule bbduk_trmfilt:
#  input:
#    r1 = "trm/{sample}_R1.trmd.fq.gz",
#    r2 = "trm/{sample}_R2.trmd.fq.gz",
#    phiX = config["phiX"]
#  output:
#    r1 = "trm/{sample}_R1.trmdfilt.fq.gz",
#    r2 = "trm/{sample}_R2.trmdfilt.fq.gz"
#  log:
#    'log/bbduk_trmfilt_{sample}.log'
#  threads: 12
#  message: """ --- Additional quality and PhiX filtering and removal of reads with more than 1 N --- """
#  shell:
#    """
#    /apps/uibk/bin/sysconfcpus -n 12 bbduk.sh -Xms110g -Xmx110g in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} ref={input.phiX} maq=10 k=31 hdist=1 maxns=1 ordered=t stats=trm/stats/{wildcards.sample}.stats2 bhist=trm/hist/{wildcards.sample}.bhist2 qhist=trm/hist/{wildcards.sample}.qhist2 lhist=trm/hist/{wildcards.sample}.lhist2 tpe tbo 2> {log} """
#
#
#rule fastqc_trm:
#  input:
#    r1 = "trm/{sample}_R1.trmdfilt.fq.gz",
#    r2 = "trm/{sample}_R2.trmdfilt.fq.gz"
#  output:
#    "trm/qc/fastqc/{sample}_R1.trmdfilt_fastqc.html",
#    "trm/qc/fastqc/{sample}_R1.trmdfilt_fastqc.zip",
#    "trm/qc/fastqc/{sample}_R2.trmdfilt_fastqc.html",
#   "trm/qc/fastqc/{sample}_R2.trmdfilt_fastqc.zip"
# log: "log/{sample}_trmdfilt_qc.log"
#  threads: 1
#  message: """ Quality check with FastQC after trimming with bbduk """
#  shell:
#    """
#    fastqc -o trm/qc/fastqc/ -f fastq {input.r1} & fastqc -o trm/qc/fastqc/ -f fastq {input.r2} 2> {log}
#    """
#
#
#rule multiqc_trm:
#  output:
#    "trm/qc/multiqc/"
#  log: "log/trm_multiqc_report.log"
#  message: """ Multiqc of trimmed reads """
#  shell:
#    """
#    multiqc trm/qc/fastqc/ -o {output}
#    """
#
#
#rule kraken2:
#  input:
#    r1 = "trm/{sample}_R1.trmdfilt.fq.gz",
#    r2 = "trm/{sample}_R2.trmdfilt.fq.gz",
#    DB = config["KRAKEN2_DB"]
#  output:
#     touch("KRAKEN2_RESULTS/{sample}_kraken2_classified.done")
#  log: "log/{sample}_kraken2_classified.log"
#  shell:
#    """
#    module load ncbi-blast/2.7.1
#    /apps/uibk/bin/sysconfcpus -n 12 kraken2 --db {input.DB} --threads 12 --paired --confidence 0.9 --classified-out KRAKEN2_RESULTS/{wildcards.sample}_R#.trmdfilt.classified.fq --unclassified-out KRAKEN2_RESULTS/{wildcards.sample}_R#.trmdfilt.unclassified.fq --report KRAKEN2_RESULTS/{wildcards.sample}.kraken_conf_0.9.report --output KRAKEN2_RESULTS/{wildcards.sample}.kraken_conf_0.9.out {input.r1} {input.r2} 2> {log}
#    """
#
#
#rule process_kraken2:
#  input:
#    touched = "KRAKEN2_RESULTS/{sample}_kraken2_classified.done",
#    r1 = "KRAKEN2_RESULTS/{sample}_R_1.trmdfilt.classified.fq",
#    r2 = "KRAKEN2_RESULTS/{sample}_R_2.trmdfilt.classified.fq"
#  output:
#    r1 = "KRAKEN2_RESULTS/{sample}_R1.trmdfilt.classified_retain.fq",
#    r2 = "KRAKEN2_RESULTS/{sample}_R2.trmdfilt.classified_retain.fq"
#  log: "log/{sample}.classified_retained.log"
#  message:
#    """
#    Process kraken2 report: from the classified reads grep all reads that are classified as opisthokonta or higher tax ids as they could still contain Daphnia reads (the settings below are for the big database containing protozoa, plants, fungi and humans. Clades below opisthokonta do not need to be considered for this database (no reads until human stuff) but this needs to be adapted depending on database
#    """
#  shell:
#    """
#    grep -A3 --no-group-separator -e 'kraken:taxid|1$' -e 'kraken:taxid|131567$' -e 'kraken:taxid|2759$' -e 'kraken:taxid|33154$' {input.r1} > {output.r1} && grep -A3 --no-group-separator -e 'kraken:taxid|1$' -e 'kraken:taxid|131567$' -e 'kraken:taxid|2759$' -e 'kraken:taxid|33154$' {input.r2} > {output.r2} 2> {log}
#    """
#
#
#rule join_kraken_reads:
#  input:
#    touched = "KRAKEN2_RESULTS/{sample}_kraken2_classified.done",
#    unclassified_r1 = "KRAKEN2_RESULTS/{sample}_R_1.trmdfilt.unclassified.fq",
#    unclassified_r2 = "KRAKEN2_RESULTS/{sample}_R_2.trmdfilt.unclassified.fq",
#    retained_r1 = "KRAKEN2_RESULTS/{sample}_R1.trmdfilt.classified_retain.fq",
#    retained_r2 = "KRAKEN2_RESULTS/{sample}_R2.trmdfilt.classified_retain.fq"
#  output:
#    r1 = "KRAKEN2_RESULTS/{sample}_R1.trmdfilt.keep.fq.gz",
#    r2 = "KRAKEN2_RESULTS/{sample}_R2.trmdfilt.keep.fq.gz"
#  log: "log/{sample}.join_kraken_reads.log"
#  message: """ Join the unclassified reads and the retained classified reads, remove intermediate files if desired """
#  shell:
#    """
#    cat {input.unclassified_r1} {input.retained_r1} | gzip > {output.r1}  &&
#    cat {input.unclassified_r2} {input.retained_r2} | gzip > {output.r2} 2> {log} 
#    """
#
#
#rule fastqc_KRAKEN:
#  input:
#    r1 = "KRAKEN2_RESULTS/{sample}_R1.trmdfilt.keep.fq.gz",
#    r2 = "KRAKEN2_RESULTS/{sample}_R2.trmdfilt.keep.fq.gz"
#  output:
#    "KRAKEN2_RESULTS/qc/fastqc/{sample}_R1.trmdfilt.keep_fastqc.html",
#    "KRAKEN2_RESULTS/qc/fastqc/{sample}_R2.trmdfilt.keep_fastqc.html",
#    "KRAKEN2_RESULTS/qc/fastqc/{sample}_R1.trmdfilt.keep_fastqc.zip",
#    "KRAKEN2_RESULTS/qc/fastqc/{sample}_R2.trmdfilt.keep_fastqc.zip"
#  log: "log/{sample}.qcKRAKEN.log"
#  threads: 12
#  message: """ Fastqc quality check of kept reads after cleaning with KRAKEN """
#  shell:
#    """
#    fastqc -o KRAKEN2_RESULTS/qc/fastqc/ -f fastq {input.r1} & fastqc -o KRAKEN2_RESULTS/qc/fastqc/ -f fastq {input.r2} 2> {log}
#    """
#
#
#rule multiqc_KRAKEN:
#  input:
#    #"KRAKEN2_RESULTS/qc/fastqc/{sample}_R1_trmdfilt.keep_fastqc.html",
#    #"KRAKEN2_RESULTS/qc/fastqc/{sample}_R2_trmdfilt.keep_fastqc.html"
#  output:
#    "KRAKEN2_RESULTS/qc/multiqc/"
#  log: "log/KRAKEN_multiqc_report.log"
#  message: """ Multiqc of kept reads after cleaning with KRAKEN """
#  shell:
#    """
#    multiqc KRAKEN2_RESULTS/qc/fastqc/ -o {output}
#    """
#

rule bb_indexRef:
  input:
    ref = config['ref_HiC'],
  output:
    touch("ref/bb_indexRef.done")
  log: "log/bb_indexRef.log" 
  threads: 12
  message: """ --- Index reference genome for bbmap --- """
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 bbmap.sh -Xmx110g t={threads} ref={input.ref}
    """


rule bbmap:
  input:
    kraken_r1 = lambda wildcards: getKrakenHome(wildcards.sample)[0],
    kraken_r2 = lambda wildcards: getKrakenHome(wildcards.sample)[1],
    ref = config['ref_HiC'],
    idxref = "ref/bb_indexRef.done",
    SAMPLETABLE = "list/reference_clones_readgroups_sorted_uniq.tsv"
  output:
    #bam = 'bbmap/HiC/minid76/{sample}.bam',
    bam = "/home/uibk/c7701178/local/bbmap/HiC/minid76/{sample}.bam",
    #sam = 'bbmap/HiC/{sample}.sam',
    bhist = "bbmap/HiC/minid76/stats/bhist/{sample}.bhist.txt",
    qhist = "bbmap/HiC/minid76/qhist/{sample}.qhist.txt",
    lhist = "bbmap/HiC/minid76/lhist/{sample}.lhist.txt",
    covstats = "bbmap/HiC/minid76/cov/{sample}.covstats.txt",
    covhist = "bbmap/HiC/minid76/cov/{sample}.covhist.txt",
    basecov = "bbmap/HiC/minid76/cov/{sample}.basecov.txt",
    bincov = "bbmap/HiC/minid76/cov/{sample}.bincov.txt"
  log: "log/HiC/minid76/bbmap_HiC_{sample}.log"
  threads: 24
  message: """ --- Mapping reads to reference genome with minid 0.76 (default), convert 2 bam, exclude unmapped reads, only keep reads with minq => 20 --- """
  shell:
    """
    id=`echo {input.kraken_r1} | sed -e "s/_R_1.trmdfilt.keep.fastq.gz$//" | cut -f 9 -d '/'`
    echo $id
    
    RG_ID=`grep -P "${{id}}\t" {input.SAMPLETABLE} | cut -f 9`
    RG_LB=`grep -P "${{id}}\t" {input.SAMPLETABLE} | cut -f 13`
    RG_SM=`grep -P "${{id}}\t" {input.SAMPLETABLE} | cut -f 11`
    RG_PL=`grep -P "${{id}}\t" {input.SAMPLETABLE} | cut -f 12`
    RG_PU=`grep -P "${{id}}\t" {input.SAMPLETABLE} | cut -f 10`

    echo $RG_ID
    echo $RG_LB
    echo $RG_SM
    echo $RG_PL
    echo $RG_PU
    
    /apps/uibk/bin/sysconfcpus -n 24 bbmap.sh -Xmx200g t={threads} ref={input.ref} in1={input.kraken_r1} in2={input.kraken_r2} out=stdout.sam minid=0.76 k=13 bw=0 ordered=t rgid=$RG_ID rglb=$RG_LB rgsm=$RG_SM rgpl=$RG_PL rgpu=$RG_PU overwrite=f unpigz=t bhist={output.bhist} qhist={output.qhist} lhist={output.lhist} covstats={output.covstats} covhist={output.covhist} basecov={output.basecov} bincov={output.bincov} | samtools view -F 4 -Shu -q 20 | samtools sort - -o {output.bam} 2> {log}
    """


rule bbmap_minid95:
  input:
    kraken_r1 = lambda wildcards: getKrakenHome(wildcards.sample)[0],
    kraken_r2 = lambda wildcards: getKrakenHome(wildcards.sample)[1],
    ref = config['ref_HiC'],
    idxref = "ref/bb_indexRef.done",
    SAMPLETABLE = ancient('list/reference_clones_readgroups_sorted_uniq.tsv')
  output:
    bam = "bbmap/HiC/minid95/{sample}.bam",
    bhist = "bbmap/HiC/minid95/stats/bhist/{sample}.bhist.txt",
    qhist = "bbmap/HiC/minid95/stats/qhist/{sample}.qhist.txt",
    lhist = "bbmap/HiC/minid95/stats/lhist/{sample}.lhist.txt",
    covstats = "bbmap/HiC/minid95/stats/cov/{sample}.covstats.txt",
    covhist = "bbmap/HiC/minid95/stats/cov/{sample}.covhist.txt",
    basecov = "bbmap/HiC/minid95/stats/cov/{sample}.basecov.txt",
    bincov = "bbmap/HiC/minid95/stats/cov/{sample}.bincov.txt"
  log: "log/HiC/minid95/bbmap_HiC_{sample}.log"
  threads: 24
  message: """ --- Mapping reads to reference genome with minid 0.95, convert 2 bam, exclude unmapped reads, only keep reads with minq => 20 --- """
  shell:
    """
    id=`echo {input.kraken_r1} | sed -e "s/_R_1.trmdfilt.keep.fastq.gz$//" | cut -f 9 -d '/'`
    echo $id
    
    RG_ID=`grep -P "${{id}}\t" {input.SAMPLETABLE} | cut -f 9`
    RG_LB=`grep -P "${{id}}\t" {input.SAMPLETABLE} | cut -f 13`
    RG_SM=`grep -P "${{id}}\t" {input.SAMPLETABLE} | cut -f 11`
    RG_PL=`grep -P "${{id}}\t" {input.SAMPLETABLE} | cut -f 12`
    RG_PU=`grep -P "${{id}}\t" {input.SAMPLETABLE} | cut -f 10`

    echo $RG_ID
    echo $RG_LB
    echo $RG_SM
    echo $RG_PL
    echo $RG_PU
    
    /apps/uibk/bin/sysconfcpus -n 24 bbmap.sh -Xmx200g t={threads} ref={input.ref} in1={input.kraken_r1} in2={input.kraken_r2} out=stdout.sam minid=0.95 k=13 bw=0 ordered=t rgid=$RG_ID rglb=$RG_LB rgsm=$RG_SM rgpl=$RG_PL rgpu=$RG_PU overwrite=f unpigz=t bhist={output.bhist} qhist={output.qhist} lhist={output.lhist} covstats={output.covstats} covhist={output.covhist} basecov={output.basecov} bincov={output.bincov} | samtools view -F 4 -Shu -q 20 | samtools sort - -o {output.bam} 2> {log}
    """  


rule bamIndex:
  input:
    "bbmap/HiC/{sample}.bam"
  output:
    "bbmap/HiC/{sample}.bam.bai"
  threads: 2
  log: "log/HiC/minid95/indexBam_HiC_{sample}.log"
  message: """--- Indexing with samtools ---"""
  shell:
    """
    samtools index {input} {output}
    """


rule mergeBAM:
  output:
    mrgd = "bbmap/HiC/minid95/merged/{sample}.merged.bb.HiC.bam"
  log: "log/HiC/minid95/mergeBam_HiC_{sample}.log"
  threads: 12
  message: """--- Merging bam files with samtools merge ---"""
  shell:
    """
    samtools merge -f {output.mrgd} bbmap/HiC/minid95/{wildcards.sample}*.bam
    """


rule bamIndex_mrgd:
  input:
    mrg = "bbmap/HiC/minid95/merged/{sample}.merged.bb.HiC.bam"
  output:
    bai = "bbmap/HiC/minid95/merged/{sample}.merged.bb.HiC.bam.bai"
  log: "log/HiC/minid95/indexMrgd_HiC_{sample}.log"
  threads: 2
  message: """--- Indexing with samtools ---"""
  shell:
    """
    samtools index {input.mrgd} {output.bai}
    """


rule remove_duplicates:
  input:
    mrgd = "bbmap/HiC/minid95/merged/{sample}.merged.bb.HiC.bam",
    bai = "bbmap/HiC/minid95/merged/{sample}.merged.bb.HiC.bam.bai"
  output:
    deDup = "deDup/HiC/minid95/{sample}.bb.HiC.dedup.bam",
    metrics = "deDup/HiC/minid95/{sample}.bb.HiC.dedup.metrics.txt"
  log: "log/HiC/minid95/{sample}.dedup.bam.log"
  threads: 12
  message: """--- Removing duplicates of merged bam files with Picard ---"""
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar /home/uibk/c7701178/.conda/envs/da/share/picard-2.27.2-0/picard.jar MarkDuplicates TMP_DIR=./ MAX_RECORDS_IN_RAM=15000000 REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 CREATE_INDEX=true INPUT={input.mrgd} OUTPUT={output.deDup} METRICS_FILE={output.metrics} 2> {log}
    """


rule mark_duplicates:
  input:
    mrg = "bbmap/HiC/minid95/merged/{sample}.merged.bb.HiC.bam",
    bai = "bbmap/HiC/minid95/merged/{sample}.merged.bb.HiC.bam.bai"
  output:
    mrkDup = "deDup/HiC/minid95/{sample}.bb.HiC.mrkdup.bam",
    metrics = "deDup/HiC/minid95/{sample}.bb.HiC.mrkdup.metrics.txt"
  log: "log/HiC/minid95/{sample}.mrkdup.bam.log"
  threads: 12
  message: """--- Marking duplicates of merged bam files with Picard ---"""
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar /home/uibk/c7701178/.conda/envs/da/share/picard-2.27.2-0/picard.jar MarkDuplicates TMP_DIR=./ MAX_RECORDS_IN_RAM=15000000 REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 CREATE_INDEX=true INPUT={input.mrgd} OUTPUT={output.deDup} METRICS_FILE={output.metrics} 2> {log}
    """


rule clip_overlap:
  input:
    deDup = "deDup/HiC/minid95/{sample}.bb.HiC.dedup.bam",
  output:
    clip = "deDup/HiC/minid95/{sample}.bb.HiC.overlapclipped.bam" 
  log: "log/HiC/minid95/{sample}.bb.HiC.overlapclipped.log"
  threads: 12
  message:
    """ Clip overlapping paired end reads """
  shell:
    """
    bam clipOverlap --in {input.deDup} --out {output.clip} --stats 2> {log}
    """


rule Index_clippedBAM:
  input:
    clip = "deDup/HiC/minid95/{sample}.bb.HiC.overlapclipped.bam"
  output:
    idx = "deDup/HiC/minid95/{sample}.bb.HiC.overlapclipped.bam.bai"
  log: 'log/HiC/minid95/{sample}.bb.HiC.overlapclipped.bam.bai.log'
  threads: 2
  message: """--- Indexing clipped BAM files with samtools ---"""
  shell:
    """
    samtools index {input.clip} {output.idx} 2> {log}
    """


rule refIndex:
  input:
    ref = config['ref_HiC']
  output:
    "ref/Dgaleata_M5_PBasm.FINAL.fasta.fai"
  log: "log/Dgaleata_M5_PBasm.FINAL_refIndex.log"
  shell:
    """
    samtools faidx {input.ref} 2> {log}
    """


rule ref_Dict:
  input:
    ref = config['ref_HiC']
  output:
    "ref/Dgaleata_M5_PBasm.FINAL.dict"
  log: "log/Dgaleata_M5_PBasm.FINAL_refDict.log"
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar /home/uibk/c7701178/.conda/envs/da/share/picard-2.27.2-0/picard.jar CreateSequenceDictionary R={input.ref} O={output} 2> {log}
    """


rule ls_ClipBam:
  input:
    clip = "deDup/HiC/minid95/{sample}.bb.HiC.overlapclipped.bam"
  output:
    touch("deDup/HiC/minid95/{sample}.added2ClippedList.done")
  log: 'log/HiC/minid95/{sample}.added2ClippedList.log'
  message: """--- Creating a sample list of clipped bam files for indel-realignement ---"""
  shell:
    """
    ls {input.clip} >> list/HiC_overlapclippedBAM.list 2> {log}
    """


rule list_indels:
  input:
    ref = config['ref_HiC'],
    idx_ref = "ref/Dgaleata_M5_PBasm.FINAL.fasta.fai",
    dict_ref = "ref/Dgaleata_M5_PBasm.FINAL.dict",
    cliped_list = "list/HiC_overlapclippedBAM.list"
  output:
    indels = "list/HiC_indels.list"
  log: "log/HiC/minid95/HiC_listIndels.log"
  threads: 48
  message:
    """ Create list of potential indels """
  shell:
    """
    GATK=(/home/uibk/c7701178/.conda/envs/da/opt/gatk-3.8/GenomeAnalysisTK.jar)
    /apps/uibk/bin/sysconfcpus -n 48 java -Xmx440g -jar $GATK -T RealignerTargetCreator -R {input.ref} -I list/HiC_overlapclippedBAM.list -o {output.indels} -drf BadMate 2> {log}
    """


rule realign_indel:
  input:
    clip = "deDup/HiC/minid95/{sample}.bb.HiC.overlapclipped.bam",
    ref = config['ref_HiC'],
    indels = 'list/HiC_indels.list'
  output:
    realigned = 'realigned/HiC/minid95/{sample}.bb.HiC.realigned.bam'
  log: 'log/HiC/minid95/{sample}.bb.HiC.realigned.log'
  threads: 12
  message:
    """ Realign in-dels """
  shell:
    """
    GATK=(/home/uibk/c7701178/.conda/envs/da/opt/gatk-3.8/GenomeAnalysisTK.jar)
    /apps/uibk/bin/sysconfcpus -n 12 java -Xmx110g -jar $GATK -T IndelRealigner -R {input.ref} -I {input.clip} -targetIntervals {input.indels} -o {output.realigned} --consensusDeterminationModel USE_READS 2> {log}
    """


rule samtools_coverage:
  input:
    realigned = 'realigned/HiC/minid95/{sample}.bb.HiC.realigned.bam'
  output:
    coverage = 'depth/HiC/minid95/{sample}.bb.HiC.realigned.bam.coverage.hist'
  log: 'log/HiC/minid95/{sample}.bb.HiC.realigned.bam.coverage.hist.log'
  threads: 12
  message: """ Produce an ASCII-art histogram of coverage per chromosome   """
  shell:
    """
    samtools coverage -A -o {output.coverage} {input.realigned} 2> {log}
    """



rule samtools_depth:
  input:
    realigned = 'realigned/HiC/minid95/{sample}.bb.HiC.realigned.bam'
  output:
    depth = "depth/HiC/minid95/{sample}.bb.HiC.realigned.bam.depth.gz"
  log: "log/HiC/minid95/{sample}.bb.HiC.realigned.bam.depth.log"
  threads: 12
  message:
    """ Count per position depth per sample using samtools depth """
  shell:
    """
    samtools depth -aa {input.realigned} | cut -f3 | gzip > {output.depth} 2> {log}
    """


rule read_depth:
  input:
    #'depth/depth.list'
  output:
    args2 = "depth/HiC/minid95/stats/references_depth_statistics.txt",
  log: 'log/HiC/minid95/genome_stats.log'
  threads: 12
  message:
    """ --- Running Rscript to plot the genome-wide distribution of coverage --- """
  shell:
    """
    depthlist=(depth/HiC/minid95/depth.list)
    ls depth/HiC/minid95/*.bb.HiC.realigned.bam.depth.gz | cut -f4 -d '/' > $depthlist 
    Rscript scripts/read_depth.R $depthlist {output.args2} 2> {log} 
    """


rule plot_summary:
  input:
    args1 = "depth/HiC/minid95/stats/references_depth_statistics.txt",
  output:
    args2 = "depth/HiC/minid95/stats/references_depth_hist_dfAll.pdf",
    args3 = "depth/HiC/minid95/stats/references_depth_hist_df10.pdf",
    args4 = "depth/HiC/minid95/stats/references_depth_df1_boxplot.pdf",
    args5 = "depth/HiC/minid95/stats/references_depth_df10_boxplot.pdf",
    args6 = "depth/HiC/minid95/stats/references_depthFilter.list",
    args7 = "depth/HiC/minid95/stats/references_depth_NonZero_hist_perSample.pdf",
    args8 = "depth/HiC/minid95/stats/references_dedupBAM_depth1.list",
    args9 = "depth/HiC/minid95/stats/references_dedupBAM_depth10.list",
    args10 = "depth/HiC/minid95/stats/references_dedupBAM_depth_under10.list"
  log: 'log/HiC/minid95/plot_summary.log'
  threads: 12
  message:
    """ --- Running Rscript to plot depth summary per sample and output depth filters --- """
  shell:
    """
    Rscript scripts/plot_summary.R {input.args1} {output.args2} {output.args3} {output.args4} {output.args5} {output.args6} {output.args7} {output.args8} {output.args9} {output.args10} 2> {log}
    """


rule Qmap_deDup:
  input:
    deDup = "deDup/HiC/minid95/{sample}.bb.HiC.dedup.bam"
  output:
    directory("qmap/{sample}.bb.HiC.dedup")
  log: "log/HiC/minid95/{sample}.qmapdedup.log"
  message: """--- Qualimap deduplicated bams ---"""
  threads: 12 
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 qualimap bamqc -bam {input.deDup} -outdir {output} --java-mem-size=100g 2> {log}
    """


rule Qmap_real:
  input:
    realigned = "realigned/HiC/minid95/{sample}.bb.HiC.realigned.bam"
  output:
    directory("qmap/HiC/minid95/{sample}.bb.HiC.realigned")
  log: "log/HiC/minid95/{sample}.qmaprealigned.log"
  message: """--- Qualimap realigned bams ---"""
  threads: 12 
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 12 qualimap bamqc -bam {input.realigned} -outdir {output} --java-mem-size=100g 2> {log}
    """


