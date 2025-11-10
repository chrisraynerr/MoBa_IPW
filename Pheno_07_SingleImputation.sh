#!/bin/bash
#SBATCH --output=/tsd/p805/data/durable/projects/crayner/logs/%x.%j.out

cd /tsd/p805/data/durable/projects/crayner/GenPar/
module purge; module load R-bundle-CRAN/2023.12-foss-2023a

Rscript --vanilla Pheno_07_SingleImputation.R --dataset InitialParticipation_NarrowEligibility_TrioSSBVars.rds

