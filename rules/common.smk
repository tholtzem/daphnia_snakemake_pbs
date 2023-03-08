import os
import pandas as pd

###### Config file and sample sheets #####


configfile: 
  "config/config.yaml"

## load sample info ##
#samples_information = pd.read_csv("list/reference_clones_readgroups_sorted_uniq.tsv", sep='\t', index_col=False)

samples_information = pd.read_csv("list/reference_clones_readgroups_sorted_uniq.tsv", sep='\t').set_index("fastq_prefix", drop=False)
#samples_information = pd.read_csv("list/Da_ResteggsOnly_metadata_LC_LZ.csv", sep='\t').set_index("fastq_prefix", drop=False)
sample_dir = list(samples_information['path2fastq'])
fastq_prefix = list(samples_information['fastq_prefix'])
sample_ID = list(samples_information['sample_id'])


## more sample information
#metadata = pd.read_csv("list/sequenceIDs_update_20210701.csv", sep=',', index_col=False)
#sample = list(metadata['ID'])
#species = list(pd.unique(metadata['species']))
#region = list(metadata['region'])
#latitude = list(metadata['latitude'])
#longitude = list(metadata['longitude'])
#lake = list(metadata['lake'])


## load chromosom info
#chromosom_information = pd.read_csv("list/dgal_rapid_ChromInfo.csv", sep=',', index_col=False)
#chromosom_information = pd.read_csv("list/dgal_rapid_Chrom_dgal1_dgal10.csv", sep=',', index_col=False)
### chromosom_information
#chromosom_names = list(map(str, chromosom_information['chrom']))
#chromosom_length= list(chromosom_information['length'])
#chromosom_start = list(chromosom_information['start'])
#chromosom_end = list(chromosom_information['end'])

## load prefix list of VCF chromosoms
#chromosom_vcfPrefix = pd.read_csv("list/vcfChrom_prefix.list", sep=',', index_col=False)

###### helper functions ######

#def getFqHome(sample):
#  return(list(os.path.join(samples_dict[sample],"{0}_{1}_001.fastq.gz".format(sample,pair)) for pair in ['R1','R2']))

#def getTrmHome(sample):
#  return(list(os.path.join('trm', "{0}_{1}_trm.fq.gz".format(sample,pair)) for pair in ['R1','R2']))

#def getTrmFiltHome(sample):
#  return(list(os.path.join('trm', "{0}_{1}_trmfilt.fq.gz".format(sample,pair)) for pair in ['R1','R2']))

def get_fastq_prefix(sample):
  return samples_information.loc[index, "fastq_prefix"]

def getKrakenHome(sample):
    return(list(os.path.join('/home/uibk/c7701125/scratch/EXCHANGE/CLONE_DATA/KRAKEN2', "{0}_{1}.trmdfilt.keep.fastq.gz".format(sample,pair)) for pair in ['R_1','R_2']))

#new_d = {}
#for key, value in d.items():
#  if not d[key] in new_d:
#    new_d[d[key]] = [key]
#  else:
#    new_d[d[key]].append(key)
#
#df = pd.DataFrame.from_dict(new_d, orient="index")
#df.to_csv("list/new_d.csv")

#def get_input_for_list_ChromVCF(wildcards):
#  ck_output = checkpoints.bcftools_mpileup_call.get(**wildcards).output[0]
#  VCF, = glob_wildcards(os.path.join(ck_output, "{sample}.vcf.gz"))
#  return expand(os.path.join(ck_output, "{SAMPLE}.vcf.gz"), SAMPLE=VCF)
