#!/bin/bash

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# source ~/config.sh

cd /tsd/p805/data/durable/projects/crayner/GenPar

DATASET=InitialParticipation_NarrowEligibility_TrioSSBVars.rds
DATASET=ContinuedParticipation_Q1Vars.rds
SCRIPT1=Pheno_Step07_SingleImputation.R
SCRIPT2=Pheno_Step07_SingleImputation.sh

## TRANSFER DATASETS:
# sh DataTransferCluster.sh ${TSD}/GenPar/${DATASET} ${DATASET}
# sh DataTransferCluster.sh ${TSD}/GenPar/${PREDICTORS} ${PREDICTORS}

## TRANSFER SCRIPTS:
# sh DataTransferCluster.sh $TSD/GenPar/Pheno_Step01_SingleImputation.R Pheno_Step01_SingleImputation.R
# sh DataTransferCluster.sh ${TSD}/GenPar/Data/${DATASET} Data/${DATASET}

# tr -d '\r' < ${SCRIPT1} > ${SCRIPT1}.tmp && mv ${SCRIPT1}.tmp ${SCRIPT1}; chmod +x ${SCRIPT1}

# module purge; module load R/4.2.1-foss-2022a
# module purge; module load R/4.3.2-gfbf-2023a

cat > ${SCRIPT2} << EOT
#!/bin/bash
#SBATCH --output=/tsd/p805/data/durable/projects/crayner/logs/%x.%j.out

cd /tsd/p805/data/durable/projects/crayner/GenPar/
module purge; module load R-bundle-CRAN/2023.12-foss-2023a

Rscript --vanilla Pheno_Step07_SingleImputation.R --dataset ${DATASET}

EOT

SCRIPT2=Pheno_Step07_SingleImputation.sh
tr -d '\r' < ${SCRIPT2} > ${SCRIPT2}.tmp && mv ${SCRIPT2}.tmp ${SCRIPT2}; chmod +x ${SCRIPT2}

sbatch \
--account=p805 \
--partition=bigmem \
--time=5-00:00:00 \
--cpus-per-task=32 \
--mem-per-cpu=20G \
--job-name=MLIMimputation_${DATASET} \
${SCRIPT2} 

# --partition=bigmem \

# SCRIPT="Pheno_Step05_SingleImputation_sbatch.sh"
# tr -d '\r' < ${SCRIPT} > ${SCRIPT}.tmp && mv ${SCRIPT}.tmp ${SCRIPT}; chmod +x ${SCRIPT}
# sh ${SCRIPT}
