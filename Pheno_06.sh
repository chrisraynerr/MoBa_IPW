#!/bin/bash
#SBATCH --output=/tsd/p805/data/durable/projects/crayner/logs/%x.%j.out

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd /tsd/p805/data/durable/projects/crayner/GenPar/
module purge; module load R-bundle-CRAN/2023.12-foss-2023a

Rscript --vanilla Pheno_Step06_SsbRegistryTrioVars.R 

## *****************************************************************************
## USAGE:
# SCRIPT=Pheno_Step06.sh
# tr -d '\r' < ${SCRIPT} > ${SCRIPT}.tmp && mv ${SCRIPT}.tmp ${SCRIPT}; chmod +x ${SCRIPT}
# sbatch --account=p805 --time=6:00:00 --cpus-per-task=2 --mem-per-cpu=50G --job-name=Pheno_Step06 ${SCRIPT}
## *****************************************************************************
