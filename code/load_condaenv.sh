!/bin/bash

# properties = {properties}

module load Anaconda3/2021.04/miniconda-base-2021.04
source $UIBK_CONDA_PROFILE
conda activate da

{exec_job}

