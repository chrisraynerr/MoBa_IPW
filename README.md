Scripts to replicate analyses described here: https://osf.io/preprints/osf/ymk37_v1  
  
Pheno_Functions: Custom R functions used in this project  
  
Pheno_01:        Computes participation phenotypes in MoBa participants  
Pheno_02-04:     Uses registry data to identify eligible non-participants from the population  
Pheno_05-06:     Extracts registry variables to use as predictors of participation  
Pheno_07:        Performs single imputation of missing data in the predictors of participation  
Pheno_08:        Tidies up the dataset in preparation for the prediction model  
Pheno_09:        Performs an Elastic net penalized regression model for the specificed participation outcome  
Pheno_10:        Prepares a dataset to compare population and sample estimates  
Pheno_11:        Tests for differences in sample characteristics between the target population and analytical samples (sample representativenesss)  
Pheno_12:        Tests for differences in univariable associations between the target population and analytical samples (estimate representativenesss)  
