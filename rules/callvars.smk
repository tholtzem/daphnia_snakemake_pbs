rule callVars:
  input:
    ref = config['ref'],
    bamlist = 'list/daphnia118.list'
  output:
    'vars/daphnia_init_bb.vcf.gz'
  log:
    'log/daphnia_init_bb.log'
  threads: 36
  message: """--- Call variants of merged bam files using bbtools ---"""
  shell:
    """
    /apps/uibk/bin/sysconfcpus -n 36 callvariants2.sh -Xmx330g t={threads} list={input.bamlist} out={output} ref={input.ref} multisample=t extended=t ploidy=2 2> {log}
    """ 

