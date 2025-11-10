Scripts to replicate analyses described here: https://osf.io/preprints/osf/ymk37_v1 \n\n

Pheno_Functions: Custom R functions used in this project \n\n

Pheno_01:        Computes participation phenotypes in MoBa participants\n
Pheno_02-04:     Uses registry data to identify eligible non-participants from the population\n
Pheno_05-06:     Extracts registry variables to use as predictors of participation\n
Pheno_07:        Performs single imputation of missing data in the predictors of participation\n
Pheno_08:        Tidies up the dataset in preparation for the prediction model\n
Pheno_09:        Performs an Elastic net penalized regression model for the specificed participation outcome\n
Pheno_10:        Prepares a dataset to compare population and sample estimates\n
Pheno_11:        Tests for differences in sample characteristics between the target population and analytical samples (sample representativenesss)\n
Pheno_12:        Tests for differences in univariable associations between the target population and analytical samples (estimate representativenesss)\n
