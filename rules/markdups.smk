
rule remove_duplicates:
  input:
    mrgd = 'bbmap/rapid/{ID}.merged.bb.RAPID.bam'
  output:
    deDup = 'deDup/{ID}.bb.RAPID.dedup.bam',
    metrics = 'deDup/{ID}.bb.RAPID.dedup.metrics.txt'
  threads: 12
  message: """--- Removing duplicates of merged bam files with Picard ---"""
  shell:
    """
    picard MarkDuplicates TMP_DIR=./ MAX_RECORDS_IN_RAM=15000000 REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 CREATE_INDEX=true INPUT={input.mrgd} OUTPUT={output.deDup} METRICS_FILE={output.metrics}
    """

rule mark_duplicates:
  input:
    mrgd =  'bbmap/rapid/{ID}.merged.bb.RAPID.bam'
  output:
    mrkDup = 'markDup/{ID}.bb.RAPID.mrkdup.bam',
    metrics = 'markDup/{ID}.bb.RAPID.mrkdup.metrics.txt'
  threads: 12
  message:
    """--- Marking duplicates of merged bam files with Picard ---"""
  shell:
    """
   picard MarkDuplicates TMP_DIR=./ MAX_RECORDS_IN_RAM=15000000 REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 CREATE_INDEX=true INPUT={input.mrgd} OUTPUT={output.mrkDup} METRICS_FILE={output.metrics}
   """

