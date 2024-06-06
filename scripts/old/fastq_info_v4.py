#!/home/uibk/c7701125/.conda/envs/daph/bin/python

## Reads all *fastq.gz files from directory or given as a list and generates a table with information for readgroups and builds readgroups. Can take a path to directory with fastq.gz files (-F) and, if only specific files should be used, a list of fastq.gz files (in the latter case use -I and -F). Additional sample information can be added (location, lat, lon, species....) with the -S option (Note that first column in this file must be named 'clone_ID' and match 'sample' names in the *fastq.gz file). This version of the script generates a pandas dataframe and outputs it also as csv file

## ToDO: for Pfäffikersee (and other samples with Umlaut) change names for readgroups!!!

import argparse, sys, os, glob, gzip, pandas as pd

parser = argparse.ArgumentParser(description="define input and output file")

parser.add_argument('-I', "--Input", type=argparse.FileType('r'), help = "Optional list of *fastq.gz files", required = False)
parser.add_argument('-O', "--Output", type=argparse.FileType('w'), help = "File information table", default=sys.stdout, required = True)
parser.add_argument('-F', "--Path_to_fastqc", type=str, help = "path to directory with  *fastq.gz files", required = False)
parser.add_argument('-S', "--Sample_Info", type=argparse.FileType('r'), help = "Additional sample information file. First column must be named 'clone_ID' and  match 'sample' name in fastq.gz file", required = False)

args = parser.parse_args()
OF = args.Output
PATH = args.Path_to_fastqc

if args.Input is not None:
	IF = args.Input
	fastqlist = IF.read().splitlines()
	#print(fastqlist)
else:
	#print("Read files from directory")
	#print(PATH + '*fastqc.gz')
	fastqlist = [os.path.basename(fastqfiles) for fastqfiles in glob.glob(PATH + '/*fastq.gz')]
	fastqlist.sort()
	#print(fastqlist)

if args.Sample_Info is not None:
	SI = args.Sample_Info
	info_df = pd.read_csv(SI) 
	#print(info_df)
else: 
	info_df=None

column_names=["file", "sample","library_info","instrument","run_id","flow_cell","lane","ID","PU","SM","PL","LB"]
#OF.write(','.join(column_names) + '\n')

linelist = []
for fastq in fastqlist:
	split_name = fastq.split("_")
	#print(fastq.split("_"))
	with gzip.open (PATH + '/' + fastq,'rt') as f:
		firstline = f.readline().strip().split(':')
	
	tags  = ['file', 'sample','library_info','instrument','run_id','flow_cell','lane']
	row_values = [fastq] + split_name[0:2] + firstline[0:4]
	d = dict(zip(tags,row_values))
	
	SM = d['sample']
	PL = 'illumina'
	if d['library_info'].startswith('2'):
		LB = 'LIB2of' + SM
	else:
		LB = 'LIB1of' + SM

	ID = '_'.join([d['flow_cell'], d['lane'], SM, LB]) 
	PU = '_'.join([d['flow_cell'], d['lane'], SM, LB])

	linelist.append([fastq] + split_name[0:2] + firstline[0:4] + [ID, PU, SM, PL, LB])
	fastq_df = pd.DataFrame(linelist, columns=column_names)

#print(fastq_df)
#print(info_df)

if info_df is not None: 
	merged_info = pd.merge(fastq_df, info_df, left_on='SM', right_on='clone_ID')
	#print(merged_info)
	merged_info.to_csv(OF, index=False)
else:
	fastq_df.to_csv(OF, index=False)	

OF.close()

if args.Input is not None:
	IF.close()
