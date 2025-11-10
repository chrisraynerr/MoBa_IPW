#!/usr/bin/Rscript

cat("\n *** Script for:","\n *** \t - cleaning full dataset in preparation for participation prediction models")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Preparing workspace \n")

remove(list = ls())

source("profile.R")
pkgs <- c("data.table","dplyr","stringr","tidyr","bigstatsr","ggplot2")

for(p in pkgs){  LoadOrInstall(p) }

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
cat("\n ... getting options \n")

opt_list=list(
  optparse::make_option("--DfY",action="store",default="InitialParticipation_NarrowEligibility_TrioSSBVars.rds",type="character",help="Outcome dataset"),
  optparse::make_option("--DfX",action="store",default="InitialParticipation_NarrowEligibility_TrioSSBVars_PreImputationRF.rds",type="character",help="Outcome dataset"),
  optparse::make_option("--Y",action="store",default="MoQ1Participation",type="character",help="Outcome label"),
  optparse::make_option("--outdir",action="store",default=F,type="character",help="Output directory if previously computed FBMs")
)
opt=optparse::parse_args(optparse::OptionParser(option_list=opt_list))

DATA_Y   <- opt$DfY 
DATA_X   <- opt$DfX
OUTCOME  <- opt$Y
OUTDIR   <- opt$outdir
WHO      <- str_sub(OUTCOME,1,2)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
cat("\n ... loading data \n")

Df <-readRDS(paste0("Data/",DATA_Y))%>%select(LNR,ID=all_of(paste0(WHO,"LNR")),everything())
Y  <-Df%>%select(LNR,ID,matches("Participation",ignore.case=F))
DfY<-Df%>%select(LNR,ID,all_of(OUTCOME))
DfX<-readRDS(paste0("Data/",DATA_X))%>%select(-matches("Participation",ignore.case=F))
Df <-DfY%>%bind_cols(DfX)%>%na.omit()%>%arrange(ID,desc(!!sym(OUTCOME)),ChDOB)%>%distinct(ID,.keep_all=T)

Df <-Df%>%select(-ID,-matches("LNR",ignore.case=F))%>% 
  dplyr::mutate(across(where(~n_distinct(.)<10)&!matches("Particpation"),~as.character(.x)))%>% 
  dplyr::mutate(across(matches("Municipality"),~as.character(.x)))%>% 
  dplyr::mutate(across(where(is.character),function(x) factor(x,levels=names(sort(table(x),decreasing=T)))))

Mat <- Y %>% bind_cols(Df %>% select(!matches("Participation")))
saveRDS(Mat,"Data/NarrowEligibility_TrioSSBVars_AnalysisData.rds")


cat("\n ... Categorical variables ...")
  
DfF <- Df  %>% select(all_of(OUTCOME),where(is.factor)) %>% 
 mutate(across(where(is.factor),~as.character(.x))) %>% 
 mutate(across(where(is.factor),~str_to_title(.x,))) %>% 
 mutate(across(where(is.factor),~str_remove_all(.x," "))) %>% 
 mutate(across(where(is.factor),~str_replace_all(.x,"[^[:alnum:]]","_"))) %>% 
 mutate(across(where(is.factor),~str_sub(.x,1,20))) %>% 
 mutate(across(where(is.factor),~ifelse(is.na(.x),"Missing",.x))) %>% 
 select(where(~n_distinct(.)>1))

Fct <- paste0(names(DfF)[-1],collapse = " + ")
Fct <- model.matrix.lm(as.formula(paste0(OUTCOME,"~",Fct)),data=DfF,na.action = "na.pass")
Fct <- Fct[,-1]

# Identify binary variables
BIN <- colnames(Fct)[apply(Fct, 2, function(x) length(unique(x)) == 2)]
# Calculate minor group sizes and percentages
MIN <- sapply(BIN, function(var) min(table(Fct[, var])))
PER <- MIN/nrow(Fct)*100
# Filter binary variables based on minor group percentage threshold
DROP<- BIN[PER<2]
length(DROP)
# Drop the selected variables from the model matrix
Fct <- Fct[,!colnames(Fct) %in% DROP]

saveRDS(DROP,"AnalysisData_DroppedCategories.rds") 
tail(DROP)

MISSING <- apply(Fct,2, FUN = function(x){sum(is.na(x))})
max(MISSING)

Fct <- apply(Fct, 2, function(x) {  
 x[is.na(x)] <- mean(x, na.rm = TRUE)
 x
})

cat("\n ... Continuous variables ...")

winsorize <- function(x,lowerP=0.01,upperP=0.99) {
 lowerB<-quantile(x,lowerP,na.rm=T);upperB<-quantile(x,upperP,na.rm=T)
 x[x<lowerB]<-lowerB; x[x>upperB]<-upperB; return(x)
}

Num <- Df  %>% 
 select(-all_of(OUTCOME)) %>% select(!matches("Participation")) %>% 
 select(!where(is.factor)) %>% 
 mutate(across(everything(),as.numeric)) %>% 
 mutate(across(everything(), ~winsorize(.x))) %>% 
 # mutate(across(everything(), ~(.x-mean(.x,na.rm=T))/sd(.x,na.rm=T))) %>% 
 mutate(across(everything(), ~replace(., is.na(.), mean(., na.rm=T)))) %>% 
 dplyr::select(where(~mean(is.na(.))<1)) %>% 
 as.matrix()

MISSING <- apply(Num,2, FUN = function(x){sum(is.na(x))})
max(MISSING)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
cat("\n ... join predictor matrix \n")

Mat <- cbind.data.frame(Num,Fct)
Var <- names(as.data.frame(Mat))

cat("\n there are",nrow(Mat)," individuals\n")
cat("\n there are",length(Var)," predictors\n")
cat("\n predictors include: \n", Var)

Mat <- Y %>% bind_cols(Mat)

saveRDS(Mat,"NarrowEligibility_TrioSSBVars_AnalysisDataDummyCoded.rds")

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
remove(list = ls())
q(save="no")
exit()



