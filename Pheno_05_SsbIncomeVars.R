#!/usr/bin/Rscript

cat("\n *** Script for:", "\n *** \t - extracting income variables")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Preparing workspace \n")

remove(list = ls())

source("profile.R")
pkgs <- c("data.table","dplyr","stringr","tidyr", "foreign","haven","lubridate","purrr")
for(p in pkgs){  LoadOrInstall(p) }

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Load narrow eligibility criteria\n")

EliID  <- data.table::fread("Data/InitialParticipation_NarrowEligibility_CriteriaVars_UniqueMoLNR.csv")  %>% 
  mutate(across(everything(),~ifelse(.x==""|.x==" ",NA,.x)))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Identify income variables and correction for inflation\n")

IncNames    <- fread("DataDictionary/VarNames_Innekt.txt",header=F) %>% 
  mutate(V1=str_remove_all(str_to_title(str_replace_all(V1,"_"," "))," "))%>%distinct() %>%
  group_by(V1)%>%mutate(V1=ifelse(row_number()>1,paste0(V1,"_",row_number()),V1))%>%ungroup()

CPI    <- read.table("DataDictionary/ConsumerPriceIndex1931_2024.txt",h=T)
BASE   <- CPI[CPI$Year==min(CPI$Year),]$Average
CPI$IF <- CPI$Average/BASE
CPI    <- CPI[,c("Year","IF")] %>% na.omit()


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Load income data, filter to specified dates for each participant, and save required data\n")

PATHS  <- list.files(paste0(SSB_REGISTER,"SSB/01_data/data_v5.0/csv"),pattern="INNTEKT.csv",full.names=T)

freadFilter <- function(x){fread(x,header=T,data.table=F,na.strings=c(NA,"NA",""))%>%dplyr::select(LNR=w19_0634_lnr,Year=aargang,everything())%>% 
    filter(Year>1996&Year<2011)%>%filter(LNR%in%EliID$MoLNR|LNR%in%EliID$FaLNR) }

IncDf <- freadFilter(PATHS)

data.table::fwrite(IncDf,"Data/InitialParticipation_NarrowEligibility_IncomeVars.csv")
