#!/usr/bin/Rscript
cat("\n *** Script for:","\n *** \t - preparing dataset for quantification and correction of bias")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Preparing workspace \n")

remove(list = ls())

source("profile.R")
source("Pheno_Functions.R")

pkgs <- c("data.table","dplyr","stringr","tidyr","ggplot2","survey", "flextable","officer")

for(p in pkgs){  LoadOrInstall(p) }


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... loading data \n")

DfY <- readRDS(paste0("Data/","InitialParticipation_NarrowEligibility_TrioSSBVars.rds")) %>% 
  select(LNR,matches("Participation")) %>% 
  select(LNR,matches("MoQ"))

DfX <- readRDS(paste0("Data/","InitialParticipation_NarrowEligibility_TrioSSBVars_PreImputationRF.rds")) %>% 
  select(-matches("Participation",ignore.case=F))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... merging and tidying data \n")

Df  <- DfY %>% bind_cols(DfX) %>% na.omit() %>% arrange(LNR,MoQ1Participation,ChDOB)%>%distinct(LNR,.keep_all=T)

Df  <- Df  %>% 
  dplyr::mutate(across(where(~n_distinct(.)<10)&!matches("Particpation"),~as.character(.x)))%>% 
  dplyr::mutate(across(matches("Municipality"),~as.character(.x)))%>% 
  dplyr::mutate(across(where(is.character),function(x) factor(x,levels=names(sort(table(x),decreasing=T)))))

Vars <- 
  c("LNR","ChDOB","FaAge","MoAge","ChSex",
    "MoHoursWorked","FaHoursWorked",
    "MoHouseholdIncome_IF_LOG10","FaHouseholdIncome_IF_LOG10",
    "MoEduYears","FaEduYears",
    "MoMoDead","MoFaDead","FaMoDead","FaFaDead",
    "ChNumberNorwegianParents","MoNumberNorwegianParents","FaNumberNorwegianParents",
    "MoEduCode","FaEduCode",
    "MoBirthPlace","FaBirthPlace",
    "MoCivilStatus_2","FaCivilStatus_2","MoCivilStatus","FaCivilStatus",
    "MoEmploymentStatus","FaEmploymentStatus",
    "MoResidenceSettlementUrbanicity","FaResidenceSettlementUrbanicity",
    "MoResidenceSettlementPopulationDensity","FaResidenceSettlementPopulationDensity",
    "ChEng05Zsc","ChEng08Zsc","ChEng09Zsc","ChMath05Zsc","ChMath08Zsc","ChMath09Zsc"
    )

X <- Df %>% 
  select(matches("MoQ1|MoQ9"),any_of(Vars)) %>% 
  mutate(MoUniversityDegree=ifelse(grepl(x=MoEduCode,"Level[6-8]"),1,0)) %>% 
  mutate(FaUniversityDegree=ifelse(grepl(x=FaEduCode,"Level[6-8]"),1,0)) %>% 
  mutate(MoBornInNorway=ifelse(grepl(x=MoBirthPlace,"Norway"),1,0)) %>% 
  mutate(FaBornInNorway=ifelse(grepl(x=FaBirthPlace,"Norway"),1,0)) %>% 
  mutate(across(matches("Dead|NumberNorw"),~as.integer(.x))) %>% 
  mutate(across(matches("Dead"),~(.x*-1)+2)) %>% 
  mutate(ChNumberOfLivingGrandparents=MoMoDead+FaMoDead+MoMoDead+FaMoDead) %>% 
  mutate(ChNumberNorwegianGrandparents=MoNumberNorwegianParents+FaNumberNorwegianParents)%>%  
  select(-matches("EduCode"),-matches("BirthPlace"),-matches("Dead"),-MoNumberNorwegianParents,-FaNumberNorwegianParents) %>% 
  mutate(ChYOB=as.integer(str_sub(ChDOB,1,4))) %>% 
  mutate(MoCivilStatus_Married=ifelse(grepl(x=MoCivilStatus_2,"1_Married_Cohabiting"),1,0)) %>% 
  mutate(FaCivilStatus_Married=ifelse(grepl(x=FaCivilStatus_2,"1_Married_Cohabiting"),1,0)) %>% 
  mutate(MoEmploymentStatus_WageEarner=ifelse(grepl(x=MoEmploymentStatus,"Wage"),1,0)) %>% 
  mutate(FaEmploymentStatus_WageEarner=ifelse(grepl(x=FaEmploymentStatus,"Wage"),1,0)) %>% 
  mutate(MoEmploymentStatus_Unemployed=ifelse(grepl(x=MoEmploymentStatus,"Outside"),1,0)) %>% 
  mutate(FaEmploymentStatus_Unemployed=ifelse(grepl(x=FaEmploymentStatus,"Outside"),1,0)) %>% 
  mutate(MoResidence_Urban=ifelse(grepl(x=MoResidenceSettlementUrbanicity,"T"),1,0)) %>% 
  mutate(FaResidence_Urban=ifelse(grepl(x=FaResidenceSettlementUrbanicity,"T"),1,0)) %>% 
  mutate(ChSex_Male=ifelse(ChSex=="1",1,0)) %>% 
  mutate(MoHouseholdIncome_IF=10^(MoHouseholdIncome_IF_LOG10)) %>% 
  mutate(FaHouseholdIncome_IF=10^(FaHouseholdIncome_IF_LOG10)) %>% 
  mutate(ChEngZsc=rowMeans(across(c(ChEng05Zsc,ChEng08Zsc,ChEng09Zsc)))) %>% 
  mutate(ChMathZsc=rowMeans(across(c(ChMath05Zsc,ChMath08Zsc,ChMath09Zsc)))) %>% 
  mutate(ChMeanNationalTestScore=rowMeans(across(c(ChMathZsc,ChEngZsc)))) %>% 
  select(
    LNR,
    MoQ1Participation,MoQ9Participation, 
    MoAge,
    MoBornInNorway,
    MoCivilStatus_Married,
    MoEmploymentStatus_WageEarner,
    MoEmploymentStatus_Unemployed,
    MoHoursWorked,
    MoHouseholdIncome_IF,
    MoUniversityDegree,
    MoEduYears,
    MoResidence_Urban,
    MoResidenceSettlementPopulationDensity,
    ChYOB,ChSex_Male,
    ChNumberNorwegianParents,
    ChMeanNationalTestScore
    )

X  <- X %>% select(!matches("^PP_"))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... saving data \n")

saveRDS(X,"Data/NarrowEligibility_BiasVariables.rds")

remove(list = ls())

