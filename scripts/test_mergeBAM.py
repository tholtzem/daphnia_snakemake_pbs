import os
import pandas as pd

## load sample info ##
#samples_information = pd.read_csv("samples.txt", sep='\t', index_col=False)
samples_information = pd.read_csv("samRAPID.txt", sep='\t', index_col=False)
sample_names = list(samples_information['sample'])
sample_locations = list(samples_information['location'])
sample_ID = list(samples_information['id'])
sample_unit = list(samples_information['units'])
samples_set = zip(sample_names, sample_locations)
samples_dict = dict(zip(sample_names, sample_locations))
d = dict(zip(sample_unit, sample_ID))

## Create new dictionary ##
new_d = {}
for key, value in d.items():
    if not d[key] in new_d:
        new_d[d[key]] = [key]
    else:
        new_d[d[key]].append(key)

for i in new_d.items():
    #print(i[0])
    #print(i[1])
    inputVar = " ".join(i[1])
    #print(inputVar)
    outputVar = i[0]
    #print(outputVar)
    ## define more strings ##
    final_string = str(" echo samtools merge" + " " + outputVar + ".merged.bb.RAPID.bam" + " " + inputVar)
    alt_string = str(" echo ln -s" + " " + inputVar + " " + outputVar + ".merged.bb.RAPIP.bam")
    #print(final_string)
    #print(alt_string)
    if len(i[1]) > 1:
        os.system(final_string)
    else:
        os.system(alt_string)


