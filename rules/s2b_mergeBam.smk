rule sam2bam:
  input:
    sam = '/home/uibk/c7701178/scratch/DAPHNIA/04_mapped/bbmap/sam/RAPID/{sample}.sam.gz'
  output:
    bam = 'bbmap/rapid/{sample}.bam'
  message: """--- Converting sam files to sorted bam ---"""
  threads: 12 
  shell:
    """
    samtools view -F 4 -Shu {input.sam} | samtools sort - -o {output.bam}
    """

rule bamIndex:
  input:
    aln = 'bbmap/rapid/{sample}'	
  output:
    idx = "bbmap/rapid/{sample}.bai"
  threads: 1
  message: """--- Indexing with samtools ---"""
  shell:
    """
    samtools index {input.aln} {output.idx}
    """

#rule replacegroups:
#  output:
#    'replaceReadGroup_umbra.done'
#  threads: 12
#  message: """ replace read groups umbra """
#  shell:
#    """
#    ./scripts/ReplaceReadGroup.sh 
#    """

rule mergeBAM:
  input:
    aln = 'new_d.csv'
  output:
    mrgd = touch('mytask.done')
  threads: 12
  message: """--- Merging bam files with samtools merge ---"""
  run:
    df = pd.read_csv(input.aln, index_col=0, keep_default_na=False)
    filecontents = df.to_dict("split")
    filecontents = dict(zip(filecontents["index"], filecontents["data"]))
    for i in filecontents.items():
      filter_object = filter(lambda x: x != "", i[1])
      without_empty_strings = list(filter_object)
      inputVar = " ".join(without_empty_strings)
      outputVar = i[0]
      final_string = str("samtools merge" + " " + "bbmap/rapid/" + outputVar + ".merged.bb.RAPID.bam" + " " + inpu)
      alt_string = str("ln -s" + " " + "/home/uibk/c7701178/scratch/DAPHNIA/daphnia_snakemake_pbs/" + inputVar + " " + "bbmap/rapid/" + outputVar + ".merged.bb.RAPID.bam")
      final_string = str("samtools merge -f" + " " + "bbmap/rapid/" + outputVar + ".merged.bb.RAPID.bam" + " " + inputVar)
      alt_string = str("ln -s" + " " + "/home/tania/github/daphnia_snakemake_pbs/" + inputVar + " " + "bbmap/rapid/" + outputVar + ".merged.bb.RAPID.bam")
      if len(without_empty_strings) > 1:
        os.system(final_string)
      else:
        os.system(alt_string)
