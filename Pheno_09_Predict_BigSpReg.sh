#!/bin/bash
#SBATCH --output=/tsd/p805/data/durable/projects/crayner/logs/%x.%j.out

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
# export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd /tsd/p805/data/durable/projects/crayner/GenPar/
module purge; module load R-bundle-CRAN/2023.12-foss-2023a

Rscript --vanilla Pheno_Step09_Predict_BigSpReg.R --Y=${1} --S=${2} #--DfY=${3} --DfX=${4}  

## *****************************************************************************
## USAGE:
# cd /tsd/p805/data/durable/projects/crayner/GenPar/
# SCRIPT=Pheno_Step09_Predict_BigSpReg.sh
# tr -d '\r' < ${SCRIPT} > ${SCRIPT}.tmp && mv ${SCRIPT}.tmp ${SCRIPT}; chmod +x ${SCRIPT}
# sbatch --account=p805 --partition=bigmem --time=6:00:00 --cpus-per-task=4 --mem-per-cpu=50G --job-name=Predict_FaQF2Participation ${SCRIPT}  FaQF2Participation FaQ1Participation
# sbatch --account=p805 --partition=bigmem --time=6:00:00 --cpus-per-task=4 --mem-per-cpu=50G --job-name=Predict_FaBiobankParticipation ${SCRIPT}  FaBiobankParticipation FaQ1Participation
# sbatch --account=p805 --partition=bigmem --time=6:00:00 --cpus-per-task=4 --mem-per-cpu=50G --job-name=Predict_MoBiobankParticipation ${SCRIPT}  MoBiobankParticipation MoQ1Participation
# sbatch --account=p805 --partition=bigmem --time=6:00:00 --cpus-per-task=4 --mem-per-cpu=50G --job-name=Predict_ChBiobankParticipation ${SCRIPT}  ChBiobankParticipation MoQ1Participation


# sbatch --account=p805 --partition=bigmem --time=6:00:00 --cpus-per-task=4 --mem-per-cpu=50G --job-name=Predict_${OUTCOME} ${SCRIPT}  MoQ9Participation MoQ1Participation

# # MoQ9Participation MoQ1Participation #
# for OUTCOME in FaQ1Participation FaQF2Participation MoBiobankParticipation FaBiobankParticipation ChBiobankParticipation
# do
# sbatch --account=p805 --partition=bigmem --time=6:00:00 --cpus-per-task=4 --mem-per-cpu=50G --job-name=Predict_${OUTCOME} ${SCRIPT} ${OUTCOME}
# done
## *****************************************************************************
# SCRIPT=Pheno_Step09_Predict_BigSpReg.sh
# tr -d '\r' < ${SCRIPT} > ${SCRIPT}.tmp && mv ${SCRIPT}.tmp ${SCRIPT}; chmod +x ${SCRIPT}
# sbatch --account=p805 --partition=bigmem --time=6:00:00 --cpus-per-task=4 --mem-per-cpu=50G --job-name=Predict_${OUTCOME} \
# Pheno_Step09_Predict_BigSpReg.sh \
# "MoQ9Participation" "MoQ1Participation""ContinuedParticipation_Q1Vars_PreImputationRF_IDvars.rds" "ContinuedParticipation_Q1Vars_PreImputationRF.rds"
