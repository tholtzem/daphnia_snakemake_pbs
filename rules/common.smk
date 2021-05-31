import os
import pandas as pd

###### Config file and sample sheets #####


configfile: 
  "config/config.yaml"

# load sample info 
#samples_information = pd.read_csv("samples.txt", sep='\t', index_col=False)
samples_information = pd.read_csv("list/samRAPID.txt", sep='\t', index_col=False)
sample_names = list(samples_information['sample'])
sample_locations = list(samples_information['location'])
sample_ID = list(samples_information['id'])
sample_bam = list(samples_information['bams'])
samples_set = zip(sample_names, sample_locations)
samples_dict = dict(zip(sample_names, sample_locations))
d = dict(zip(sample_bam, sample_ID))

# load chromosom info
chromosom_information = pd.read_csv("list/dgal_rapid_ChromInfo.csv", sep=',', index_col=False)
#chromosom_information
chromosom_names = list(map(str, chromosom_information['chrom']))
chromosom_length= list(chromosom_information['length'])
chromosom_start = list(chromosom_information['start'])
chromosom_end = list(chromosom_information['end'])

###### helper functions ######

#def getFqHome(sample):
#  return(list(os.path.join(samples_dict[sample],"{0}_{1}_001.fastq.gz".format(sample,pair)) for pair in ['R1','R2']))

#def getTrmHome(sample):
#  return(list(os.path.join('trm', "{0}_{1}_trm.fq.gz".format(sample,pair)) for pair in ['R1','R2']))

#def getTrmFiltHome(sample):
#  return(list(os.path.join('trm', "{0}_{1}_trmfilt.fq.gz".format(sample,pair)) for pair in ['R1','R2']))

new_d = {}
for key, value in d.items():
  if not d[key] in new_d:
    new_d[d[key]] = [key]
  else:
    new_d[d[key]].append(key)

df = pd.DataFrame.from_dict(new_d, orient="index")
df.to_csv("list/new_d.csv")

