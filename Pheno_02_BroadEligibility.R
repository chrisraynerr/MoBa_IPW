#!/usr/bin/Rscript

cat("\n *** Script for:",
    "\n *** \t - identifying all pregnancies in the population who were BROADLY eligible to participate in MoBa")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Preparing workspace \n")

remove(list = ls())

source("profile.R")
pkgs <- c("data.table","dplyr","stringr","tidyr", "foreign","haven","lubridate","countrycode")
for(p in pkgs){  LoadOrInstall(p) }

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### First, load this dataset for use with cleaning DEM variables and GEO variables\n")

Countries <- data.table::fread("DataDictionary/CountryCodes.txt",sep=";") %>% as.data.frame() %>% 
  filter(!(name1 %in% c("Bonaire", "British Indian Ocean Territory (the)", "Cocos (Keeling) Islands (the)", "Kosovo", "Stateless", "United States Minor Outlying Islands (the)", "Unspecified") ))

Countries$continent <- countrycode::countrycode(sourcevar=Countries[,"name1"],origin="country.name",destination="continent")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Load linkage file and participation phenotypes\n")

Part <- data.table::fread("Data/MoBa_ParticipationOutcomes.csv",na.strings=c(NA,"NA",""))
table(is.na(Part$LNR));table(is.na(Part$MoLNR));table(is.na(Part$FaLNR))
test <-Part %>% filter(MoInitialParticipation==1);length(unique(test$LNR));length(unique(test$MoLNR));length(unique(test$FaLNR));
test <-Part %>% filter(MoQ1Participation==1);length(unique(test$LNR));length(unique(test$MoLNR));length(unique(test$FaLNR));

# C=113792; M=94504; F=94719 Unique in initial participation
# C=103749; M=86586; F=86770 unique in Q1

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Load registry based demographic data and restrict to those borni in the recruitment window\n")

DEMS<-list.files(paste0(SSB_REGISTER,"SSB/01_data/data_v5.0/csv/"),pattern="POPULATION_FASTE_OPPLYSNINGER",full.names = T)

# First get min and max participants birthdays as broad criteria
DOB<-data.table::fread(DEMS,header=T,data.table=F,na.strings=c(NA,"NA","")) %>% 
  dplyr::select(LNR=w19_0634_lnr,DOB=foedsels_aar_mnd) %>% 
  dplyr::filter(LNR%in%Part$LNR) %>% inner_join(Part)

MinDOB <- min(DOB$DOB,na.rm=T) # 199909
MaxDOB <- max(DOB$DOB,na.rm=T) # 200908

rm(DOB)

DemDf  <- data.table::fread(DEMS,header=T,data.table=F,na.strings=c(NA,"NA","")) %>% 
  dplyr::select(
    LNR=w19_0634_lnr,MoLNR=lopenr_mor,FaLNR=lopenr_far,
    ChDOB=foedsels_aar_mnd,ChSex=kjoenn,ChDOD=doeds_aar_mnd,
    ChBirthCountry=fodeland,ChNationalHeritage=landbak3gen,ChImmigrationStatus=invkat) %>% 
  dplyr::filter(ChDOB>=MinDOB&ChDOB<=MaxDOB) # Filter data for only children born between 1999 - 2009 # from july 1999

MoDemDf  <- data.table::fread(DEMS,header=T,data.table=F,na.strings=c(NA,"NA","")) %>% 
  dplyr::select(
    MoLNR=w19_0634_lnr,MoMoLNR=lopenr_mor,MoFaLNR=lopenr_far,
    MoDOB=foedsels_aar_mnd,MoDOD=doeds_aar_mnd,
    MoBirthCountry=fodeland,MoNationalHeritage=landbak3gen,MoImmigrationStatus=invkat) %>% 
  dplyr::filter(MoLNR%in%na.omit(DemDf$MoLNR))

FaDemDf  <- data.table::fread(DEMS,header=T,data.table=F,na.strings=c(NA,"NA","")) %>% 
  dplyr::select(
    FaLNR=w19_0634_lnr,FaMoLNR=lopenr_mor,FaFaLNR=lopenr_far,
    FaDOB=foedsels_aar_mnd,FaDOD=doeds_aar_mnd,
    FaBirthCountry=fodeland,FaNationalHeritage=landbak3gen,FaImmigrationStatus=invkat) %>% 
  dplyr::filter(FaLNR%in%na.omit(DemDf$FaLNR)) 

# First get childs data (by filtering on birthday)
# Then join parental data

TrDemDf <- DemDf %>% 
  dplyr::left_join(MoDemDf %>% filter(!is.na(MoLNR)),"MoLNR") %>% 
  dplyr::left_join(FaDemDf %>% filter(!is.na(FaLNR)),"FaLNR") %>% 
  mutate(across(matches("Heritage"), function(.x) factor(.x, levels=c(Countries$code),labels=c(Countries$name2)))) %>% 
  mutate(across(matches("Country"), function(.x) factor(.x, levels=c(Countries$code),labels=c(Countries$name2)))) %>% 
  mutate(across(matches("Country"), function(.x) factor(.x, levels=c(Countries$name2),labels=c(Countries$continent)),.names="{stringr::str_replace(.col, 'Country', 'Continent')}")) %>%
  mutate(across(matches("ImmigrationStatus"), function(.x) factor(.x, levels=c("A","F","C","B","E","G"), 
         labels=c( "0_Norwegian_NorwegianParents", "Norwegian_OneImmigrantParent",  "Norwegian_BothImmigrantParents", "Immigrant",  "Immigrant_OneNorwegianParent", "Immigrant_BothNorwegianParent")))) %>% 
  mutate(across(matches("ImmigrationStatus"), function(.x) factor(.x, levels=c( "0_Norwegian_NorwegianParents", "Norwegian_OneImmigrantParent",  "Norwegian_BothImmigrantParents", "Immigrant",  "Immigrant_OneNorwegianParent", "Immigrant_BothNorwegianParent"),
         labels=c(0,0,0,1,1,1)),.names="{.col}01")) %>% 
  mutate(across(ends_with("ImmigrationStatus"), function(.x) factor(.x, levels=c( "0_Norwegian_NorwegianParents", "Norwegian_OneImmigrantParent",  "Norwegian_BothImmigrantParents", "Immigrant",  "Immigrant_OneNorwegianParent", "Immigrant_BothNorwegianParent"),
         labels=c(2,1,0,0,1,2)),.names="{stringr::str_replace(.col, 'ImmigrationStatus', 'NumberNorwegianParents')}")) %>%
  mutate(across(where(is.factor),~as.character(.x))) %>% 
  mutate(across(matches("BirthCountry$"), function(.x) str_replace_all(.x,"Norway","0_Norway" ))) %>% 
  mutate(ChBirthPlace=ifelse(ChBirthCountry=="0_Norway","0_Norway",ifelse(ChBirthCountry=="Sweden"|ChBirthCountry=="Finland"|ChBirthCountry=="Denmark","Scandinavia",ChBirthContinent))) %>% 
  mutate(MoBirthPlace=ifelse(MoBirthCountry=="0_Norway","0_Norway",ifelse(MoBirthCountry=="Sweden"|MoBirthCountry=="Finland"|MoBirthCountry=="Denmark","Scandinavia",MoBirthContinent))) %>% 
  mutate(FaBirthPlace=ifelse(FaBirthCountry=="0_Norway","0_Norway",ifelse(FaBirthCountry=="Sweden"|FaBirthCountry=="Finland"|FaBirthCountry=="Denmark","Scandinavia",FaBirthContinent))) %>% 
  mutate(across(matches("LNR"),~na_if(.,"")))%>% 
  distinct(LNR,MoLNR,FaLNR,.keep_all=T)
  # matches("DOB"),matches("DOD"),matches("BirthPlace"),matches("ImmigrationStatus"))

rm(DemDf,MoDemDf,FaDemDf)

TrDemDf <-TrDemDf %>% full_join(Part %>% select(LNR,matches("Part")),by=c("LNR")) %>% 
  mutate(across(matches("Part"),~ifelse(is.na(.x),0,.x)))

table(TrDemDf$ChBirthPlace,TrDemDf$ChImmigrationStatus)

test <- TrDemDf %>% filter(MoQ1Participation==1)
table(test$ChBirthPlace,test$ChImmigrationStatus)

prop.table(table(TrDemDf$MoQ1Participation))
prop.table(table(TrDemDf$MoInitialParticipation))

test <- TrDemDf %>% filter(ChBirthPlace=="0_Norway"|MoInitialParticipation==1)

prop.table(table(test$MoQ1Participation))
prop.table(table(test$MoInitialParticipation))

test2 <- TrDemDf %>% filter(ChBirthPlace!="0_Norway"&MoInitialParticipation==1)
# 142 children not born in Norway who participate
table(test$ChBirthPlace)
rm(test);rm(test2)

TrDemDf <- TrDemDf %>% filter(ChBirthPlace=="0_Norway"|MoInitialParticipation==1)

AllBirthsLNR <- TrDemDf %>% select(LNR,MoLNR,FaLNR,ChDOB)
data.table::fwrite(AllBirthsLNR,"Data/InitialParticipation_BroadEligibility_TrioLNR.csv")
sum(is.na(AllBirthsLNR$LNR))
sum(is.na(AllBirthsLNR$MoLNR))
sum(is.na(AllBirthsLNR$FaLNR))


test1 <- data.table::fread("Data/InitialParticipation_BroadEligibility_TrioLNR.csv")
test2 <- data.table::fread("Data/MoBa_ParticipationOutcomes.csv")
test3 <- test1 %>% filter(LNR %in% test2$LNR) %>% na.omit()


## THERE WERE LOTS OF MISSING LNR -- TRY TO RECOVER SOME MISSING IDS USING SLEKT

slekt <- fread(paste0(SSB_REGISTER,"SSB/01_data/data_v5.0/csv/POPULATION_SLEKT.csv"),na.strings=c(NA,"NA","")) %>% 
  select(LNR=w19_0634_lnr,MoLNR=lopenr_mor,FaLNR=lopenr_far) %>% 
  dplyr::filter(LNR%in%AllBirthsLNR$LNR) %>% mutate(across(matches("LNR"), ~na_if(.,"")))

Df0    <- AllBirthsLNR %>% filter(!(is.na(LNR)|is.na(MoLNR)|is.na(FaLNR)))
DfChNA <- AllBirthsLNR %>% filter(is.na(LNR))
DfMoNA <- AllBirthsLNR %>% filter(is.na(MoLNR)) %>% left_join(slekt,"LNR") %>% mutate(MoLNR=ifelse(!is.na(MoLNR.x),MoLNR.x,MoLNR.y))%>% mutate(FaLNR=ifelse(!is.na(FaLNR.x),FaLNR.x,FaLNR.y)) %>% select(LNR,MoLNR,FaLNR)
DfFaNA <- AllBirthsLNR %>% filter(is.na(FaLNR)) %>% left_join(slekt,"LNR") %>% mutate(MoLNR=ifelse(!is.na(MoLNR.x),MoLNR.x,MoLNR.y))%>% mutate(FaLNR=ifelse(!is.na(FaLNR.x),FaLNR.x,FaLNR.y)) %>% select(LNR,MoLNR,FaLNR)

AllBirthsLNR2 <- Df0 %>% bind_rows(DfMoNA) %>% bind_rows(DfFaNA) %>% distinct()
sum(is.na(AllBirthsLNR2$LNR))
sum(is.na(AllBirthsLNR2$MoLNR))
sum(is.na(AllBirthsLNR2$FaLNR))
sum(is.na(AllBirthsLNR2$MoLNR)&is.na(AllBirthsLNR2$FaLNR))

AllBirthsLNR2 <- AllBirthsLNR2 %>% 
  mutate(LNR=ifelse(is.na(LNR),paste0("MissChLNR_MoLNR",MoLNR),LNR)) %>% 
  mutate(FaLNR=ifelse(is.na(FaLNR),paste0("MissFaLNR_ChLNR",LNR),FaLNR))%>% 
  mutate(MoLNR=ifelse(is.na(MoLNR),paste0("MissMoLNR_ChLNR",LNR),MoLNR)) %>% 
  filter(!grepl("_NA",FaLNR))

sum(is.na(AllBirthsLNR2$LNR))
sum(is.na(AllBirthsLNR2$MoLNR))
sum(is.na(AllBirthsLNR2$FaLNR))
sum(is.na(AllBirthsLNR2$MoLNR)&is.na(AllBirthsLNR2$FaLNR))


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ### Saving linkage file to all eligible pregnancies in population \n")

data.table::fwrite(AllBirthsLNR2,"Data/InitialParticipation_BroadEligibility_TrioLNR.csv")

remove(list = ls())