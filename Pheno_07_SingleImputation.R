#!/usr/bin/Rscript

cat("\n *** Script for:","\n *** \t - performing single imputation of missing predictors")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Preparing workspace \n")

remove(list = ls())
source("profile.R")

LIB<-"/tsd/p805/data/durable/projects/crayner/software/R/4.3.2"
.libPaths(c(LIB,.libPaths()))

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... loading packages \n")

packages <- c("optparse","data.table","dplyr","stringr","tidyr","rlang","magrittr","glmnet",
              "jomo","mice","missRanger","md.log","lme4","mitml","memuse","futile.logger")

if(!"RCurl" %in% rownames(installed.packages())){install.packages("RCurl",dependencies=T)}
if(!"h20" %in% rownames(installed.packages(lib.loc=LIB))){
  Sys.setenv(H2O_JAR_PATH="/tsd/p805/data/durable/projects/crayner/software/h2o-3.40.0.2/h2o.jar")
  install.packages("/tsd/p805/data/durable/projects/crayner/software/h2o-3.40.0.2/R/h2o_3.40.0.2.tar.gz",repos=NULL,type="source",dependencies=T)
  Sys.unsetenv("/tsd/p805/data/durable/projects/crayner/software/h2o-3.40.0.2/h2o.jar")
}
if(!"mlim" %in% rownames(installed.packages(lib.loc=LIB))){
  install.packages("/tsd/p805/data/durable/projects/crayner/software/mlim_0.3.0.tar.gz",repos=NULL,type="source",dependencies=T)
}
library(h2o);library(mlim)
lapply(packages, LoadPkg)

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... getting options \n")

opt_list=list(optparse::make_option("--dataset",action="store",default="InitialParticipation_NarrowEligibility_TrioSSBVars.rds",type="character",help="Scale label"))
opt=optparse::parse_args(optparse::OptionParser(option_list=opt_list))

DATASET  <- opt$dataset
OUTNAME1 <- stringr::str_replace(string=DATASET, pattern=".rds", replacement="_PreImputationRF")
OUTNAME2 <- stringr::str_replace(string=DATASET, pattern=".rds", replacement="_ImputedMLIM")

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... impute data \n")
set.seed(2024)

Df <- readRDS(paste0("Data/",DATASET))
IG <- names(Df %>% dplyr::select(matches("ID",ignore.case=F),matches("LNR",ignore.case=F),matches("Participation",ignore.case=F),matches("BaN",ignore.case=F)))

saveRDS(Df %>% dplyr::select(all_of(IG)),paste0("Data/",OUTNAME1,"_IDvars.rds")) 

if(!file.exists(paste0("Data/",OUTNAME1,".rds"))){
  cat("\n ... performing RF pre-imputation \n")
  preRF  <- mlim.preimpute(data=Df %>% select(!c(IG)),preimpute="RF",seed=2024)
  saveRDS(preRF,paste0("Data/",OUTNAME1,".rds"))
}

cat("\n ... performing primary Elastic-Net imputation procedure \n")

preRF <- readRDS(paste0("Data/",OUTNAME1,".rds"))
IG    <- c(IG, names(preRF%>%dplyr::select(where(~is.factor(.)&&n_distinct(.)>50),where(~is.character(.)&&n_distinct(.)>50))))

if(file.exists(paste0("Data/",OUTNAME2,".mlim"))){
  
  impEN <- mlim(
    load=paste0("Data/",OUTNAME2,".mlim"),
    data=Df,ignore=IG,preimputed.data=preRF,m=1,flush=T,
    algos="ELNET",seed=2024,tuning_time=36000, 
    save=paste0("Data/",OUTNAME2,".mlim"), 
    report=paste0("Data/",OUTNAME2,".md"), 
    verbosity="debug") 
  
} else {
  impEN <- mlim(
    data=Df,ignore=IG,preimputed.data=preRF,m=1,flush=T,
    algos="ELNET", seed=2024, tuning_time=36000, 
    save=paste0("Data/",OUTNAME2,".mlim"), 
    report=paste0("Data/",OUTNAME2,".md"), 
    verbosity="debug") 
}

saveRDS(impEN,paste0("Data/",OUTNAME2,".rds"))
sum<-mlim.summarize(impEN)
saveRDS(sum,paste0("Data/",OUTNAME2,"_summary.rds"))


rm(list=ls())
q(save="no")
