#!/usr/bin/Rscript

cat("\n *** Script for:",
    "\n *** \t - identifying all mothers in the population who were BROADLY eligible to initially participate in MoBa",
    "\n *** \t - computing baseline participation phenotypes in MoBa",
    "\n *** \t - computing Continued participation phenotypes in MoBa"
)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ### Preparing workspace \n")

remove(list = ls())

source("profile.R")
pkgs <- c("data.table","dplyr","stringr","tidyr", "foreign","haven","lubridate")
for(p in pkgs){  LoadOrInstall(p) }

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ### Load linkage file and create baseline participation phenotype \n")

Link   <- data.table::fread("../MoBa_SSB_IDs_Trios20250415.csv", na.strings = c(NA,"NA",""))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ### Creating baseline participation phenotype \n")

Q1Part <- foreign::read.spss(paste0(MOBA_ORIGINAL,"sav/PDB2601_status_v12.sav"), to.data.frame=TRUE) %>%
  dplyr::mutate(across(everything(), ~ na_if(str_trim(.x), ""))) %>%
  dplyr::mutate(
    MoQ1Selection = case_when (status_skj1_v12 %in% c("x", "1") ~ 1,status_skj1_v12 %in% c("2", "3") ~ 0),
    MoQ1Invited = ifelse(!is.na(status_skj1_v12), 1, 0),
    FaQ1Selection = case_when (status_skjF_v12 %in% c("x", "1") ~ 1, status_skjF_v12 %in% c("2", "3") ~0),
    FaQ1Invited = ifelse(!is.na(status_skjF_v12), 1, 0)
  ) %>% 
  dplyr::mutate(
    MoQ1Participation = ifelse(MoQ1Invited==1,MoQ1Selection,NA),
    FaQ1Participation = ifelse(FaQ1Invited==1,FaQ1Selection,NA)
  ) %>% 
  dplyr::select(PrID=PREG_ID_2601, matches("Select|Invited|Partici")) %>% 
  dplyr::mutate(across(everything(), as.integer)) %>% 
  inner_join(Link,"PrID")%>% 
  dplyr::mutate(MoBaselineParticipation=1) %>%#ifelse((is.na(MoID)|MoID==""|MoID==" "),0,1)) %>% # this is actually ascertainment!!!
  dplyr::mutate(FaBaselineParticipation=1) %>%#ifelse((is.na(FaID)|FaID==""|FaID==" "),0,1)) %>%
  dplyr::mutate(MoBiobankParticipation=ifelse((is.na(MoIID)|MoIID==""|MoIID==" "),0,1)) %>%
  dplyr::mutate(FaBiobankParticipation=ifelse((is.na(FaIID)|FaIID==""|FaIID==" "),0,1)) %>% 
  dplyr::mutate(ChBiobankParticipation=ifelse((is.na(ChIID)|ChIID==""|ChIID==" "),0,1)) %>% 
  # dplyr::select(LNR,MoLNR,FaLNR,PrID,BaN,matches("Partici"))%>% 
  mutate(across(everything(),~ifelse(.x==""|.x==" ",NA,.x))) %>% 
  filter(!(is.na(LNR)&is.na(MoLNR)&is.na(FaLNR)))

table(Q1Part$MoQ1Invited,Q1Part$MoQ1Participation,useNA="a");
table(Q1Part$FaQ1Invited,Q1Part$FaQ1Participation,useNA="a")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ### Creating continued participation phenotype \n")

Cont <- read.spss(paste0(MOBA_ORIGINAL,"sav/PDB2601_status_v12.sav"), to.data.frame = TRUE) %>%
  mutate(across(everything(), ~ na_if(str_trim(.x), ""))) %>%
  ## Make variables for participation (1) vs no participation (0)
  mutate(
    MoQ1Part = factor(case_when (status_skj1_v12 %in% c("x", "1") ~ 1,status_skj1_v12 %in% c("2", "3") ~ 0)),
    MoQ2Part = factor(case_when (status_skj2_v12 %in% c("x", "1") ~ 1,status_skj2_v12 %in% c("2", "3") ~ 0)),
    MoQ3Part = factor(case_when (status_skj3_v12 %in% c("x", "1") ~ 1,status_skj3_v12 %in% c("2", "3") ~ 0)),
    MoQ4Part = factor(case_when (status_skj4_v12 %in% c("x", "1") ~ 1,status_skj4_v12 %in% c("2", "3") ~ 0)),
    MoQ5Part = factor(case_when (status_skj5_v12 %in% c("x", "1") ~ 1,status_skj5_v12 %in% c("2", "3") ~ 0)),
    MoQ6Part = factor(case_when (status_skj6_v12 %in% c("x", "1") ~ 1,status_skj6_v12 %in% c("2", "3") ~ 0)),
    MoQ7Part = factor(case_when (status_skj5y_v12 %in% c("x", "1") ~ 1,status_skj5y_v12 %in% c("2", "3") ~ 0)),
    MoQ8Part = factor(case_when (status_skj7y_v12 %in% c("x", "1") ~ 1,status_skj7y_v12 %in% c("2", "3") ~ 0)),
    MoQ9Part = factor(case_when (status_skj8y_v12 %in% c("x", "1") ~ 1,status_skj8y_v12 %in% c("2", "3") ~ 0)),
    ## Makes variable for invitation (1) vs not (0)
    MoQ1Inv = factor(ifelse(!is.na(status_skj1_v12), 1, 0)),
    MoQ2Inv = factor(ifelse(!is.na(status_skj2_v12), 1, 0)),
    MoQ3Inv = factor(ifelse(!is.na(status_skj3_v12), 1, 0)),
    MoQ4Inv = factor(ifelse(!is.na(status_skj4_v12), 1, 0)),
    MoQ5Inv = factor(ifelse(!is.na(status_skj5_v12), 1, 0)),
    MoQ6Inv = factor(ifelse(!is.na(status_skj6_v12), 1, 0)),
    MoQ7Inv = factor(ifelse(!is.na(status_skj5y_v12), 1, 0)),
    MoQ8Inv = factor(ifelse(!is.na(status_skj7y_v12), 1, 0)),
    MoQ9Inv = factor(ifelse(!is.na(status_skj8y_v12), 1, 0)),
    ## Make variables for fathers participation (1) vs no participation (0)
    FaQF1Part = factor(case_when (status_skjF_v12 %in% c("x", "1")  ~ 1, status_skjF_v12 %in% c("2", "3") ~0)),
    FaQF2Part = factor(case_when (status_skjF2_v12 %in% c("x", "1") ~ 1, status_skjF2_v12 %in% c("2", "3") ~0)),
    ## Make variables for fathers invitation (1) vs not (0)
    FaQF1Inv = factor(ifelse(!is.na(status_skjF_v12), 1, 0)),
    FaQF2Inv = factor(ifelse(!is.na(status_skjF2_v12), 1, 0))
  ) %>%
  ## Make variables for continued participation (1) vs no participation (0)
  mutate(
    MoQ2Participation  = case_when(MoQ2Inv  == 1 & MoQ2Part  == 1 ~ 1,MoQ2Inv == 1 & MoQ2Part == 0 ~ 0),
    MoQ3Participation  = case_when(MoQ3Inv  == 1 & MoQ3Part  == 1 ~ 1,MoQ3Inv == 1 & MoQ3Part == 0 ~ 0),
    MoQ4Participation  = case_when(MoQ4Inv  == 1 & MoQ4Part  == 1 ~ 1,MoQ4Inv == 1 & MoQ4Part == 0 ~ 0),
    MoQ5Participation  = case_when(MoQ5Inv  == 1 & MoQ5Part  == 1 ~ 1,MoQ5Inv == 1 & MoQ5Part == 0 ~ 0),
    MoQ6Participation  = case_when(MoQ6Inv  == 1 & MoQ6Part  == 1 ~ 1,MoQ6Inv == 1 & MoQ6Part == 0 ~ 0),
    MoQ7Participation  = case_when(MoQ7Inv  == 1 & MoQ7Part  == 1 ~ 1,MoQ7Inv == 1 & MoQ7Part == 0 ~ 0),
    MoQ8Participation  = case_when(MoQ8Inv  == 1 & MoQ8Part  == 1 ~ 1,MoQ8Inv == 1 & MoQ8Part == 0 ~ 0),
    MoQ9Participation  = case_when(MoQ9Inv  == 1 & MoQ9Part  == 1 ~ 1,MoQ9Inv == 1 & MoQ9Part == 0 ~ 0),
    FaQF2Participation = case_when(FaQF2Inv == 1 & FaQF2Part == 1 ~ 1,FaQF2Inv == 1 & FaQF2Part == 0 ~ 0)) %>%
  mutate(across(matches("MoQ\\dInv$|FaQF\\dInv$|MoQ\\dPart$|FaQF\\dPart$"), ~ as.numeric(.x) - 1)) %>%
  mutate(across(matches("MoQ\\dInv$|FaQF\\dInv$|MoQ\\dPart$|FaQF\\dPart$"), ~ ifelse(is.na(.x),0,.x))) %>%
  ## Make variables for proportion of possible particiFaion
  mutate(FaNumInvited = rowSums(across(matches("FaQF\\dInv$")))) %>%
  mutate(MoNumInvited = rowSums(across(matches("MoQ\\dInv$")))) %>%
  mutate(FaNumPart = rowSums(across(matches("FaQF\\dPart$")))) %>%
  mutate(MoNumPart = rowSums(across(matches("MoQ\\dPart$")))) %>%  
  mutate(FaContinuedParticipation = residuals(lm(FaNumPart~FaNumInvited))) %>%
  mutate(MoContinuedParticipation = residuals(lm(MoNumPart~MoNumInvited))) %>%    
  select(PrID=PREG_ID_2601, matches("Participation"), matches("Cont"), matches("Num")) %>% mutate(PrID=as.integer(PrID)) %>% 
  distinct(PrID,.keep_all = T) %>% 
  inner_join(Link,"PrID")
  
Cont<-Cont %>% select(LNR,MoLNR,FaLNR,PrID,matches("Participation")) #,FaQF2Participation,MoContinuedParticipation,FaContinuedParticipation)

MoQ10 <- read.spss(paste0(MOBA_ORIGINAL,"sav_0523/PDB2601_Q14yrs_v12.sav"), to.data.frame=T) %>% 
  dplyr::select(PrID=PREG_ID_2601,matches("UM"))%>%
  dplyr::distinct(PrID,.keep_all=T) %>% 
  dplyr::mutate(across(everything(), ~as.character(.x))) %>% 
  dplyr::mutate(across(everything(), ~ifelse(.x==""|.x==" ",NA,.x))) %>% as.tibble() %>%  
  dplyr::mutate(across(matches("UB"), ~ifelse(is.na(.x)|.x==0,1,0),.names="{.col}_NA")) %>% 
  dplyr::mutate(TotalItems = ncol(.)-2) %>% 
  rowwise()%>%dplyr::mutate(IncompleteParticipation=sum(c_across(matches("_NA")))) %>% ungroup() %>%
  dplyr::mutate(PropIncompleteItems =IncompleteParticipation/TotalItems) %>% 
  dplyr::filter(PropIncompleteItems<.9) %>% 
  dplyr::mutate(MoQ10Participation=1) %>% 
  dplyr::select(PrID,MoQ10Participation) %>% 
  mutate(PrID=as.numeric(PrID))

ChQ10 <- read.spss(paste0(MOBA_ORIGINAL,"sav_0523/PDB2601_Q14yrs_CHILD_v12.sav"), to.data.frame=T) %>% 
  dplyr::select(PrID=PREG_ID_2601,BaN=BARN_NR,matches("UB"))%>%
  dplyr::mutate(across(everything(), ~as.character(.x))) %>% 
  dplyr::mutate(across(everything(), ~ifelse(.x==""|.x==" ",NA,.x))) %>% as.tibble() %>%  
  dplyr::select(where(~mean(is.na(.))<.1)) %>% 
  dplyr::mutate(TotalItems = ncol(.)-2) %>% 
  dplyr::mutate(across(matches("UB"), ~ifelse(is.na(.x)|.x==0,1,0),.names="{.col}_NA")) %>% 
  rowwise()%>%dplyr::mutate(ChQ10IncompleteParticipation=sum(c_across(matches("_NA")))) %>% ungroup() %>% 
  dplyr::mutate(ChQ10Participation=1) %>% 
  dplyr::select(PrID,BaN,ChQ10Participation)

ChQ11 <- read.spss(paste0(MOBA_ORIGINAL,"sav_0523/PDB2601_Q16yrs_v12.sav"), to.data.frame=T) %>% 
  dplyr::select(PrID=PREG_ID_2601,BaN=BARN_NR,matches("YA"))%>%
  dplyr::mutate(across(everything(), ~as.character(.x))) %>% 
  dplyr::mutate(across(everything(), ~ifelse(.x==""|.x==" ",NA,.x))) %>% as.tibble() %>%  
  dplyr::select(where(~mean(is.na(.))<.1)) %>% 
  dplyr::mutate(TotalItems = ncol(.)-2) %>% 
  dplyr::mutate(across(matches("YA"), ~ifelse(is.na(.x)|.x==0,1,0),.names="{.col}_NA")) %>% 
  rowwise()%>%dplyr::mutate(ChQ11IncompleteParticipation=sum(c_across(matches("_NA")))) %>% ungroup() %>% 
  dplyr::mutate(ChQ11Participation=1) %>% 
  dplyr::select(PrID,BaN,ChQ11Participation)

ChCon <- ChQ10 %>% full_join(ChQ11,by=c("PrID","BaN")) %>% mutate(across(everything(),as.numeric))

Cont <- Cont %>% left_join(MoQ10) %>% left_join(ChCon)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ### Saving participation phenotype files \n")

Part <-
  Q1Part%>%select(PrID,matches("LNR"),matches("Participation")) %>% 
  full_join(Cont%>%select(LNR,BaN,matches("Participation")),"LNR") %>% distinct(LNR,.keep_all=T) %>% 
  select(LNR,MoLNR,FaLNR,PrID,BaN,matches("Participation"))
    
data.table::fwrite(Part,"Data/MoBa_ParticipationOutcomes.csv",na=NA)

remove(list = ls())



