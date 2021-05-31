# daphnia_snakemake_PBS

======================================================

Snakemake workflow for *Daphnia* DNA analysis on the mach2 HPC cluster with the PBS-torque batch job submission system. This repository was initially based on the [ta_dna_snakemake_pbs Tutorial](https://github.com/schimar/ta_dna_snakemake_pbs),the [Snakemake Cluster Tutorial](https://github.com/SchlossLab/snakemake_cluster_tutorial.git) and the [Software Carpentry lesson repository](https://hpc-carpentry.github.io/hpc-python/17-cluster/). For more information on
snakemake itself (https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html). 


======================================================

## conda and other [dependencies](https://github.com/schimar/ta_dna_snakemake_pbs/blob/main/envs/s21.yaml)   

create environment from yaml file (in envs/):
```
# run these two once, to create the environment:
conda init bash
conda env create -f envs/s21.yaml

# with this, you can activate the environment with all [dependencies](https://github.com/schimar/ta_dna_snakemake_pbs/blob/main/envs/s21.yaml):
conda activate ta

# (also, when ssh'ing onto mach2, you can activate the env and then do a dry-run of your workflow) 
## how to submit the main snakemake job:
qsub code/clusterSnakemake.pbs

# if you've added new software to install to the conda environment, then you can update:
conda env update --name ta --file envs/s21.yaml
```



