#!/usr/bin/Rscript

cat("\n *** Script for:",
    "\n *** \t - creating a file with social, demographic and geographical predictor variables at the time of participation")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Preparing workspace \n")

remove(list = ls())

source("profile.R")
pkgs <- c("data.table","dplyr","stringr","tidyr", "foreign","haven","lubridate","purrr")
for(p in pkgs){  LoadOrInstall(p) }

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### First, load this dataset for use with cleaning DEM variables and GEO variables\n")

Countries <- data.table::fread("DataDictionary/CountryCodes.txt",sep=";") %>% as.data.frame() %>% 
  filter(!(name1 %in% c("Bonaire", "British Indian Ocean Territory (the)", "Cocos (Keeling) Islands (the)", "Kosovo", "Stateless", "United States Minor Outlying Islands (the)", "Unspecified") ))
Countries$continent <- countrycode::countrycode(sourcevar=Countries[,"name1"],origin="country.name",destination="continent")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Load narrow eligibility criteria\n")

EliID  <- data.table::fread("Data/InitialParticipation_NarrowEligibility_CriteriaVars.csv")  %>% 
  mutate(across(everything(),~ifelse(.x==""|.x==" ",NA,.x))) %>% mutate(ChYOB=as.integer(str_sub(ChDOB,1,4)))

table(is.na(EliID$ChDOB))  

Window<-EliID%>%select(LNR,MoLNR,FaLNR,ChDOB)%>%mutate(Year_1=as.integer(str_sub(ChDOB,1,4)),Year_0=Year_1-1,Year_2=Year_1+1)%>%pivot_longer(cols=matches("Year"),values_to="Year",names_to="Occasion")%>%mutate(Occasion=str_remove_all(Occasion,"Year_"))
WindowLong<-Window %>% select(LNR=MoLNR,Year) %>% bind_rows(Window %>% select(LNR=FaLNR,Year)) %>% na.omit()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Load lifetime/time-invariant demographic variables\n")

# One row per individual - in the eligible subset
DEMS   <- list.files(paste0(PATH,"registers/SSB/01_data/data_v5.0/csv/"),pattern="POPULATION_FASTE_OPPLYSNINGER",full.names = T)

DemDf  <- data.table::fread(DEMS,header=T,data.table=F,na.strings=c(NA,"NA","")) %>% 
  dplyr::select(
    LNR=w19_0634_lnr,MoLNR=lopenr_mor,FaLNR=lopenr_far,
    ChDOB=foedsels_aar_mnd,ChSex=kjoenn,ChDOD=doeds_aar_mnd,
    ChBirthCountry=fodeland,ChNationalHeritage=landbak3gen,ChImmigrationStatus=invkat) %>% 
  mutate(across(everything(),~ifelse(.x==""|.x==" ",NA,.x)))%>% 
  dplyr::filter(na.omit(LNR%in%EliID$LNR)) # Filter data for only children born between 1999 - 2009 # from july 1999

MoDemDf  <- data.table::fread(DEMS,header=T,data.table=F,na.strings=c(NA,"NA","")) %>% 
  dplyr::select(
    MoLNR=w19_0634_lnr,MoMoLNR=lopenr_mor,MoFaLNR=lopenr_far,
    MoDOB=foedsels_aar_mnd,MoDOD=doeds_aar_mnd,
    MoBirthCountry=fodeland,MoNationalHeritage=landbak3gen,MoImmigrationStatus=invkat) %>% 
  mutate(across(everything(),~ifelse(.x==""|.x==" ",NA,.x)))%>% 
  dplyr::filter(MoLNR%in%na.omit(EliID$MoLNR))

FaDemDf  <- data.table::fread(DEMS,header=T,data.table=F,na.strings=c(NA,"NA","")) %>% 
  dplyr::select(
    FaLNR=w19_0634_lnr,FaMoLNR=lopenr_mor,FaFaLNR=lopenr_far,
    FaDOB=foedsels_aar_mnd,FaDOD=doeds_aar_mnd,
    FaBirthCountry=fodeland,FaNationalHeritage=landbak3gen,FaImmigrationStatus=invkat) %>% 
  mutate(across(everything(),~ifelse(.x==""|.x==" ",NA,.x)))%>% 
  dplyr::filter(FaLNR%in%na.omit(EliID$FaLNR)) 

GrID <- na.omit(unique(c(MoDemDf$MoMoLNR,MoDemDf$MoFaLNR,FaDemDf$FaMoLNR,FaDemDf$FaFaLNR)))

GrDemDf  <- data.table::fread(DEMS,header=T,data.table=F,na.strings=c(NA,"NA","")) %>% 
  dplyr::select(LNR=w19_0634_lnr,DOD=doeds_aar_mnd,BirthCountry=fodeland,) %>% na.omit() %>% 
  dplyr::filter(LNR%in%GrID) %>% 
  mutate(Norwegian=ifelse(BirthCountry==0,1,0))

table(is.na(GrDemDf$DOD))

MoDemDf  <- MoDemDf %>% 
  left_join(GrDemDf %>% select(MoMoLNR=LNR,MoMoDOD=DOD,MoMoNorwegian=Norwegian))%>% 
  left_join(GrDemDf %>% select(MoFaLNR=LNR,MoFaDOD=DOD,MoFaNorwegian=Norwegian))

FaDemDf  <- FaDemDf %>% 
  left_join(GrDemDf %>% select(FaMoLNR=LNR,FaMoDOD=DOD,FaMoNorwegian=Norwegian))%>% 
  left_join(GrDemDf %>% select(FaFaLNR=LNR,FaFaDOD=DOD,FaFaNorwegian=Norwegian))

# First get childs data (by filtering on birthday)
# Then join parental data

TrDemDf <- DemDf %>% 
  dplyr::left_join(MoDemDf %>% filter(!is.na(MoLNR)),"MoLNR") %>% 
  dplyr::left_join(FaDemDf %>% filter(!is.na(FaLNR)),"FaLNR") %>% 
  dplyr::mutate(across(matches("DOB"), function(.x) lubridate::ymd(paste0(substr(.x, 1, 4), "-", substr(.x, 5, 6), "-01")), .names = "{.col}_2")) %>%
  dplyr::mutate(FaAge = lubridate::interval(FaDOB_2,ChDOB_2)%/%months(1)/12) %>% dplyr::mutate(MoAge = lubridate::interval(MoDOB_2,ChDOB_2)%/%months(1)/12) %>% 
  dplyr::mutate(ParentalAgeGap = FaAge - MoAge) %>% dplyr::select(-matches("DOB_2")) %>% 
  mutate(ChMOB=as.integer(str_sub(ChDOB,5,6)))%>% mutate(ChYOB=as.integer(str_sub(ChDOB,1,4)))%>%mutate(FaYOB=as.integer(str_sub(FaDOB,1,4)))%>%mutate(MoYOB=as.integer(str_sub(MoDOB,1,4)))%>% 
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
  distinct(LNR,MoLNR,FaLNR,.keep_all=T) %>% 
  mutate(across(matches("DOD",ignore.case=F),~ifelse(!is.na(.x),1,0),.names="{stringr::str_replace(.col, 'DOD', 'Dead')}"))
# matches("DOB"),matches("DOD"),matches("BirthPlace"),matches("ImmigrationStatus"))

table(TrDemDf$MoFaDead)
table(TrDemDf$MoMoDead)
table(TrDemDf$FaFaDead)
table(TrDemDf$FaMoDead)

rm(DemDf,MoDemDf,FaDemDf)

table(is.na(TrDemDf$LNR));table(is.na(TrDemDf$MoLNR));table(is.na(TrDemDf$FaLNR))

data.table::fwrite(TrDemDf %>% select(-matches("Partici")),"Data/InitialParticipation_NarrowEligibility_TrioDemVars.csv")

rm(DemDf);rm(Part);rm(TrDemDf);rm(Countries)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Load longitudinal/time-varying residential area characteristics variables\n")

Tspath <- list.files(paste0(PATH,"registers/SSB/01_data/data_v5.0/csv/"),pattern="TS_AREAL_1997|TS_AREAL_1998|TS_AREAL_1999|TS_AREAL_200|TS_AREAL_2011",full.names=T)
Names  <- str_remove_all(basename(Tspath),"POPULATION_TS_AREAL_|.csv")

freadFilter <- function(x, i){
  Df  <- fread(x,header=T,data.table=F,na.strings=c(NA,"NA",""))%>% 
    mutate(across(everything(),~ifelse(.x=="",NA,.x)))%>% 
    mutate(across(everything(),~ifelse(.x==" ",NA,.x)))%>%
    select(LNR=all_of(1),everything()) %>% 
    rename_with(~str_replace_all(str_to_title(.),"_",""),.cols=2:ncol(.)) %>% 
    filter(LNR %in% EliID$MoLNR |LNR %in% EliID$FaLNR) %>%
    purrr::set_names("LNR",paste0(names(.)[-1],"_",i))%>%
    mutate(Year=as.integer(paste0(i))) %>% 
    inner_join(WindowLong) %>% 
    rename_with(~str_to_title(.),.cols=2:ncol(.)) %>% 
    mutate(across(matches("Tskode"), function(x) str_to_upper(x))) %>% 
    mutate(across(matches("Areal"), function(x) str_replace_all(x,",","\\."))) %>% 
    mutate(across(everything(),as.character)) %>% select(-Year)
}
getmode <- function(v) {uniqv<-unique(v);uniqv[which.max(tabulate(match(v,uniqv)))]}

TsDf<-
  lapply(seq_along(Tspath), function(z) freadFilter(Tspath[z],Names[z]))
TsDf<-
  plyr::join_all(TsDf,"LNR","full")
TsDf<-
  TsDf%>%tidyr::pivot_longer(cols=!c(LNR),names_to=c("Var","Year"),names_sep="_",values_to="Value")%>%na.omit()
Tskode<-
  TsDf%>%filter(Var=="Tskode")%>%select(-Year,-Var)%>%group_by(LNR)%>%summarize(ResidenceSettlementUrbanicity=getmode(Value))
Tsstor<-
  TsDf%>%filter(Var=="Tsstor")%>%select(-Year,-Var)%>%
  mutate(Value=as.integer(str_pad(Value,width=2,side="left",pad="1")))%>%
  mutate(Value=ifelse(Value==19|Value==99,NA,Value))%>%na.omit()%>%
  group_by(LNR)%>%summarize(ResidenceSettlementPopulationDensity=median(Value,na.rm=T))
Areal<-
  TsDf%>%filter(Var=="Areal")%>%select(-Year,-Var)%>%
  mutate(Value=as.numeric(Value))%>%na.omit()%>%
  group_by(LNR)%>%summarize(ResidenceSettlementSizeKm2=mean(Value,na.rm=T))
TsDf <-
  Tskode%>%full_join(Tsstor)%>%full_join(Areal)
MoTsDf<-
  TsDf %>% filter(LNR%in%EliID$MoLNR) %>% rename_with(~paste0("Mo",.))
FaTsDf<-
  TsDf %>% filter(LNR%in%EliID$FaLNR) %>% rename_with(~paste0("Fa",.))
TsDf<-
  EliID%>%select(LNR,MoLNR,FaLNR,MoResidenceDistrict_Rec,FaResidenceDistrict_Rec)%>%left_join(MoTsDf)%>%left_join(FaTsDf)%>%
  group_by(MoResidenceDistrict_Rec) %>% 
  tidyr::fill(all_of(c("MoResidenceSettlementUrbanicity","MoResidenceSettlementPopulationDensity","MoResidenceSettlementSizeKm2")),.direction="downup") %>% ungroup() %>% 
  group_by(FaResidenceDistrict_Rec) %>% 
  tidyr::fill(all_of(c("FaResidenceSettlementUrbanicity","FaResidenceSettlementPopulationDensity","FaResidenceSettlementSizeKm2")),.direction="downup") %>% ungroup()

data.table::fwrite(TsDf %>% select(!matches("Rec")),"Data/InitialParticipation_NarrowEligibility_TrioTsVars.csv")
rm(TsDf);rm(FaTsDf);rm(MoTsDf);rm(Tskode);rm(Tsstor);rm(Areal)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Load longitudinal/time-varying civil status variables\n")

MsPath <- list.files(paste0(PATH,"registers/SSB/01_data/data_v5.0/csv/"),pattern="SIVILSTAND_",full.names = T)

MsDf   <- 
  fread(MsPath,header=T,data.table=F,na.strings=c(NA,"NA","")) %>% select(LNR=all_of(1),everything()) %>% 
  filter(LNR%in%EliID$MoLNR|LNR %in% EliID$FaLNR)%>% 
  pivot_longer(-LNR,values_to="CivilStatus",names_to="Year") %>% 
  mutate(Year=as.integer(str_remove_all(Year,"sivilstand_"))) %>% 
  group_by(LNR) %>% arrange(Year) %>% fill(CivilStatus, .direction="downup") %>% ungroup() %>% 
  filter(Year>1996&Year<2011)

MoMs <- MsDf %>% mutate(MoLNR=LNR)%>% inner_join(EliID%>% select(MoLNR,ChYOB),relationship="many-to-many")%>% 
  filter(Year==ChYOB)%>% select(MoLNR,ChYOB,MoCivilStatus=CivilStatus) %>% distinct()
FaMs <- MsDf %>% mutate(FaLNR=LNR)%>% inner_join(EliID%>% select(FaLNR,ChYOB),relationship="many-to-many")%>% 
  filter(Year==ChYOB)%>% select(FaLNR,ChYOB,FaCivilStatus=CivilStatus)  %>% distinct()

MsDf <- EliID %>%select(LNR,MoLNR,FaLNR,ChYOB) %>% 
  left_join(MoMs,c("MoLNR","ChYOB")) %>% 
  left_join(FaMs,c("FaLNR","ChYOB")) %>% 
  mutate(across(matches("CivilStatus"), function(.x) factor(.x,levels=c(1,2,3,4,5,6,7,8,9),
                                                            labels=c("2_Unmarried","1_Married","3_Widowed","4_Divorced","5_Separated","6_Registered_partner","7_Separated_partner","8_Divorced_partner","9_Surviving_partner")),
                .names="{.col}_1"))%>% 
  mutate(across(matches("CivilStatus$"), function(.x) factor(.x,levels=c(1,2,3,4,5,6,7,8,9),
                                                             labels=c("Unmarried","1_Married_Cohabiting","Widowed_Other","Separated_Divorced","Separated_Divorced","1_Married_Cohabiting","Separated_Divorced","Separated_Divorced","Widowed_Other")),
                .names="{.col}_2"))

data.table::fwrite(MsDf,"Data/InitialParticipation_NarrowEligibility_TrioMsVars.csv") 

table(MsDf$MoCivilStatus_1)
table(MsDf$MoCivilStatus_2)

rm(MoMs);rm(FaMs);rm(MsDf)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Load longitudinal/time-varying income variables\n")

IncNames<-
  fread("DataDictionary/VarNames_Innekt.txt",header=F) %>% 
  mutate(V1=str_remove_all(str_to_title(str_replace_all(V1,"_"," "))," ")) %>% distinct() %>%
  group_by(V1) %>%
  mutate(V1=ifelse(row_number() > 1, paste0(V1, "_", row_number()), V1)) %>%
  ungroup()

CPI    <- read.table("DataDictionary/ConsumerPriceIndex1931_2024.txt",h=T)
BASE   <- CPI[CPI$Year==min(CPI$Year),]$Average
CPI$IF <- CPI$Average/BASE
CPI    <- CPI[,c("Year","IF")] %>% na.omit()
rint <- function(x){n<-length(x);ranks<-rank(x);scores<-qnorm((ranks-0.5)/n);return(scores)}

IncDf<- 
  data.table::fread("Data/InitialParticipation_NarrowEligibility_IncomeVars.csv")  %>%
  rename(any_of(setNames(IncNames$V2, IncNames$V1))) %>% 
  mutate(across(everything(),~ifelse(.x==""|.x==" ",NA,.x)))%>% 
  distinct() %>% inner_join(CPI,"Year") %>% 
  mutate(across(!c(matches("LNR"),matches("Year"),matches("lopenr")), function(.x) round(.x/IF,0),.names="{.col}_IF")) %>% select(-IF) %>% 
  mutate(across(matches("_IF"), function(.x) log(.x+abs(min(.x,na.rm=T))+1,base=10),.names="{.col}_LOG10")) %>% 
  select(LNR,Year,matches("LOG10")) %>% select(-matches("IF_IF"))

MoInc<-
  IncDf%>%select(LNR,Year,matches("LOG10"))%>%filter(LNR%in%Window$MoLNR)%>%
  inner_join(Window%>%select(LNR=MoLNR,Year))%>%distinct()%>%select(-Year)%>%
  group_by(LNR)%>%summarize(across(everything(),mean,na.rm=T))%>%
  rename_with(~paste0("Mo",.))

FaInc<-
  IncDf%>%select(LNR,Year,matches("LOG10"))%>%filter(LNR%in%Window$FaLNR)%>%
  inner_join(Window%>%select(LNR=FaLNR,Year))%>%distinct()%>%select(-Year)%>%
  group_by(LNR)%>%summarize(across(everything(),mean,na.rm=T))%>%
  rename_with(~paste0("Fa",.))

IncDf <- EliID %>% select(LNR,MoLNR,FaLNR) %>% left_join(MoInc,"MoLNR")%>% left_join(FaInc,"FaLNR")

data.table::fwrite(IncDf,"Data/InitialParticipation_NarrowEligibility_TrioIncVars.csv") 

rm(MoInc);rm(FaInc);rm(IncDf)
# IncDf<-data.table::fread("Data/InitialParticipation_NarrowEligibility_TrioIncVars.csv") 


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Load longitudinal/time-varying occupational status variables\n")

OccCodes <- fread("DataDictionary/OccupationalCodes.txt") %>% 
  mutate(code2=ifelse(level==2, str_pad(code,width=2,side="left",pad="0"),
                      ifelse(level==3, str_pad(code,width=3,side="left",pad="0"),
                             ifelse(level==4, str_pad(code,width=4,side="left",pad="0"),code))))

SecCodes <- fread("DataDictionary/SectorCodes.txt") %>% filter(coding==1987&level==3) 

GeoCodes <- 
  readxl::read_excel("DataDictionary/eventualchanges_norgeo.xlsx") %>% 
  drop_na() %>% filter(oldName != "Nothing from API") %>% 
  transform(oldMunicipality=substr(oldCode,1,4)) %>% 
  transform(newMunicipality=substr(currentCode,1,4)) %>%  
  select(oldMunicipality,newMunicipality,changeOccurred) |> 
  unique() |> filter(oldMunicipality!=newMunicipality) %>% 
  arrange(oldMunicipality,desc(changeOccurred)) %>% 
  distinct(oldMunicipality,.keep_all = T)

EmpPath <- list.files(paste0(PATH,"registers/SSB/01_data/data_v5.0/csv/"),pattern="REGSYS_200|REGSYS_2010|REGSYS_2011",full.names = T)
Names <- str_remove_all(basename(EmpPath),"EMPLOYMENT_REGSYS_|.csv")

freadFilter <- function(x, i){
  Df  <- fread(x,data.table=F,colClasses="character",na.strings=c("NA","")) %>%
    select(LNR=all_of(1),CompanyID=all_of(2),EnterpriseID=all_of(3),everything()) %>% 
    filter(LNR%in%EliID$MoLNR|LNR%in%EliID$FaLNR) %>%
    mutate(across(everything(),~ifelse(.x=="",NA,.x)))%>% 
    mutate(across(everything(),~ifelse(.x==" ",NA,.x)))%>%
    rename_with(~str_replace_all(str_to_title(.),"_",""),.cols=4:ncol(.)) %>% 
    rename_with(~str_replace_all(.,"Arbarbtim","HoursWorked")) %>% 
    rename_with(~str_replace_all(.,"Arbkommnr","WorkLocationMunicipality")) %>% 
    rename_with(~str_replace_all(.,"Frtksektor","IndustrySector")) %>% 
    rename_with(~str_replace_all(.,"Arbyrke","OccupationalCode")) %>% 
    rename_with(~str_replace_all(.,'Arbhovedarbeid', 'PrimaryOrSecondaryJob')) %>%   
    rename_with(~str_replace_all(.,'Arbarbmarkstatus', 'EmploymentStatus')) %>%   
    rename_with(~str_replace_all(.,'Arbarbtid', 'HoursWorkedCategory')) %>% 
    rename_with(~str_replace_all(.,'Virknace1', 'CompanyMainIndustry'))%>% 
    select(-matches("Virkkild")) %>% #-matches("CompanyMainIndustry"),
    select(-matches("persontype"),-matches("2014"),-matches("imp"),-matches("publ"),-matches("ID",ignore.case=F))%>% 
    filter(PrimaryOrSecondaryJob!=2)%>%
    # filter(PrimaryOrSecondaryJob!=0)%>%
    arrange(LNR,desc(HoursWorked)) %>% 
    distinct(LNR,.keep_all=T) %>% 
    purrr::set_names("LNR",paste0(names(.)[-1],"_",i))%>% 
    tidyr::pivot_longer(cols=!c(LNR),names_to=c("Var","Year"),names_sep="_",values_to="Value")%>%na.omit()%>% 
    mutate(Year=as.integer(Year)) %>% inner_join(WindowLong,by=c("LNR","Year"))
}

EmpDf <- lapply(seq_along(EmpPath), function(z) freadFilter(EmpPath[z],Names[z]))
EmpDf <- bind_rows(EmpDf)
table(EmpDf$Var)

getmode <- function(v) {uniqv<-unique(v);uniqv[which.max(tabulate(match(v,uniqv)))]}

EmpDf<-EmpDf%>%distinct(LNR,Var,Year,.keep_all=T)%>%select(-Year)%>%group_by(LNR,Var)%>%summarize(Value=getmode(Value))
EmpDf<-EmpDf%>%tidyr::pivot_wider(id_cols=c(LNR),names_from=c(Var),values_from=Value,values_fn=list(Value=getmode))%>%ungroup()

AllCodes     <- EmpDf %>% select(matches("Municipality")) %>% apply(FUN=unique,2,na.rm=T) %>% unlist() %>% as.data.frame()
NoChange     <- unique(AllCodes$. [!(AllCodes$. %in% GeoCodes$oldMunicipality)])
same_District<-data.frame(oldMunicipality=c(NoChange),newDistrict=c(NoChange))

EmpDf<-EmpDf%>%
  mutate(PrimaryOrSecondaryJob=factor(PrimaryOrSecondaryJob,levels=c(0,1,2),labels=c("NoEmployment","PrimaryEmployer","SecondaryEmployer")))%>% 
  mutate(EmploymentStatus=factor(EmploymentStatus,levels=c(0,1,2,3,4),labels=c("OutsideWorkforce","WageEarner","Independent","Free","Qualification")))%>% 
  mutate(HoursWorkedCategory=factor(HoursWorkedCategory,levels=c("1","2","3"),labels=c("0_19hrs","20_29hrs","30plus"))) %>% 
  mutate(HoursWorked=as.numeric(HoursWorked))%>% 
  mutate(across(matches("OccupationalCode"),function(x) factor(
    str_sub(x,end=1), levels=c(OccCodes[which(OccCodes$level==1),]$code),labels=c(OccCodes[which(OccCodes$level==1),]$Name)),
    .names="{stringr::str_replace(.col, 'OccupationalCode', 'OccupationalCodeLevel1')}"))%>% 
  mutate(across(matches("OccupationalCode$"),function(x) factor(
    str_sub(x,end=2), levels=c(OccCodes[which(OccCodes$level==2),]$code),labels=c(OccCodes[which(OccCodes$level==2),]$Name)),
    .names="{stringr::str_replace(.col, 'OccupationalCode', 'OccupationalCodeLevel2')}"))%>% 
  mutate(across(matches("OccupationalCode$"),function(x) factor(
    str_sub(x,end=3), levels=c(OccCodes[which(OccCodes$level==3),]$code),labels=c(OccCodes[which(OccCodes$level==3),]$Name)),
    .names="{stringr::str_replace(.col, 'OccupationalCode', 'OccupationalCodeLevel3')}"))%>% 
  mutate(across(matches("OccupationalCode$"),function(x) factor(
    str_sub(x,end=4), levels=c(OccCodes[which(OccCodes$level==4),]$code),labels=c(OccCodes[which(OccCodes$level==4),]$Name)),
    .names="{stringr::str_replace(.col, 'OccupationalCode', 'OccupationalCodeLevel4')}"))%>% 
  mutate(across(matches("IndustrySector"),function(x) factor(x,levels=c(SecCodes$code),labels=c(SecCodes$name))))%>%
  mutate(across(matches("Municipality"),function(.x) as.character(factor(.x,levels=c(GeoCodes$oldMunicipality),labels=c(GeoCodes$newMunicipality))))) %>%  
  select(-matches("HoursWorked_"),-matches("OccupationalCode_"),-matches("OccupationalCodeLevel4_"),-matches("OccupationalCodeLevel3_"),-matches("Occasion"),-matches("publ"),-matches("Primary"))

MoEmpDf<-EmpDf%>%filter(LNR%in%EliID$MoLNR)%>%rename_with(~paste0("Mo",.))%>%distinct()
FaEmpDf<-EmpDf%>%filter(LNR%in%EliID$FaLNR)%>%rename_with(~paste0("Fa",.))%>%distinct()

EmpDf<-EliID %>% select(LNR,MoLNR,FaLNR)%>%left_join(MoEmpDf)%>%left_join(FaEmpDf)%>%distinct(LNR,MoLNR,FaLNR,.keep_all=T)

EmpDf <-  EmpDf %>%
  mutate(across(where(is.factor),as.character)) %>% 
  mutate(across(c(where(is.character)&!matches("LNR")),function(x)ifelse(is.na(x)&MoEmploymentStatus=="OutsideWorkforce","OutsideWorkforce",x)))

data.table::fwrite(EmpDf,"Data/InitialParticipation_NarrowEligibility_TrioEmpVars.csv") 
# EmpDf <- data.table::fread("Data/InitialParticipation_NarrowEligibility_TrioEmpVars.csv",data.table=F,colClasses="character",na.strings=c("NA","")) 

rm(GeoCodes);rm(OccCodes);rm(SecCodes);rm(same_District);rm(AllCodes)
rm(MoEmpDf);rm(FaEmpDf);rm(EmpDf) 


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Load longitudinal/time-varying educational attainment variables\n")

InPath <- list.files(paste0(PATH,"registers/SSB/01_data/data_v5.0/csv/"),pattern="EDUCATION_BU_UTD",full.names=T)

EduDf  <- data.table::fread(InPath,header=T,data.table=F,na.strings=c(NA,"NA","")) %>% 
  select(LNR=all_of(1),matches("nivaa")) %>% select(LNR,matches(paste0(1996:2011))) %>% 
  filter(LNR %in% EliID$MoLNR |LNR %in% EliID$FaLNR)%>% 
  mutate(across(everything(),~ifelse(.x=="",NA,.x)))%>% 
  mutate(across(everything(),~ifelse(.x==" ",NA,.x)))%>% 
  rename_with(~str_replace_all(.,"bu_nivaa","EduCode"),.cols=2:ncol(.)) %>% 
  tidyr::pivot_longer(cols=!c(LNR),names_to=c("Var","Year"),names_sep="_",values_to="EduCode")%>%na.omit()%>%select(-Var)

MoEduDf<-EduDf%>%inner_join(Window%>%select(LNR=MoLNR,Year)%>%mutate(Year=as.character(Year)))%>%
  rename_with(~paste0("Mo",.))%>%select(-MoYear)%>%
  group_by(MoLNR)%>%summarize(MoEduCode=max(MoEduCode))%>%
  mutate(MoEduYears=as.numeric(as.character(factor(MoEduCode,levels=c(0,1,2,3,4,5,6,7,8),labels=c(0,7,10,13,14,15,17,19,22)))))%>%
  distinct()

FaEduDf<-EduDf%>%inner_join(Window%>%select(LNR=FaLNR,Year)%>%mutate(Year=as.character(Year)))%>%
  rename_with(~paste0("Fa",.))%>%select(-FaYear)%>%
  group_by(FaLNR)%>%summarize(FaEduCode=max(FaEduCode))%>%
  mutate(FaEduYears=as.numeric(as.character(factor(FaEduCode,levels=c(0,1,2,3,4,5,6,7,8),labels=c(0,7,10,13,14,15,17,19,22)))))%>%
  distinct()

EduDf <- EliID %>%select(LNR,MoLNR,FaLNR) %>% 
  left_join(MoEduDf,c("MoLNR"))%>% 
  left_join(FaEduDf,c("FaLNR"))

data.table::fwrite(EduDf,"Data/InitialParticipation_NarrowEligibility_TrioEduVars.csv") 

rm(MoEduDf);rm(FaEduDf);rm(EduDf)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Load educational test score variables\n")

NTsc <- 
  fread(paste0(PATH,"registers/SSB/01_data/data_v5.0/csv/EDUCATION_NASJONALE_PROVER.csv"),header=T,data.table=F,na.strings=c(NA,"NA","")) %>% 
  select(LNR=all_of(1),everything()) %>% 
  filter(LNR %in% EliID$MoLNR |LNR %in% EliID$FaLNR|LNR %in% EliID$LNR)

## Z-SCORE READING (NPLES050809)
ReSc <- NTsc %>% 
  filter(grepl("NPLES", x=PROVE)) %>% 
  mutate(POENG=ifelse(is.na(POENG),0,POENG)) %>% 
  filter(DELTATTSTATUS=="D"&POENG>0) %>%
  distinct(LNR,PROVE,POENG,.keep_all=TRUE) %>%
  select(LNR,PROVE,POENG,AARGANG) %>%
  group_by(LNR,PROVE) %>% 
  pivot_wider(id_cols=LNR,names_from=PROVE,values_from=POENG,values_fn=mean,unused_fn=max)  %>%
  ungroup() %>% group_by(AARGANG) %>% 
  mutate(across(matches("NPLES"), function(x) ((x-mean(x,na.rm=T))/sd(x,na.rm=T)),.names="{.col}Zsc")) %>% 
  ungroup() %>% select(LNR,matches("Zsc"))

## Z-SCORE MATH (NPREG050809)
MaSc <- NTsc %>% 
  filter(grepl("NPREG", x=PROVE)) %>% 
  mutate(POENG=ifelse(is.na(POENG),0,POENG)) %>% 
  filter(DELTATTSTATUS=="D"&POENG>0) %>% 
  distinct(LNR,PROVE,POENG,.keep_all=TRUE) %>%  
  select(LNR,PROVE,POENG,AARGANG) %>%
  group_by(LNR,PROVE) %>%
  pivot_wider(id_cols=LNR,names_from=PROVE,values_from=POENG,values_fn=mean,unused_fn=max)  %>%
  ungroup() %>% group_by(AARGANG) %>% 
  mutate(across(matches("NPREG"), function(x) ((x-mean(x,na.rm=T))/sd(x,na.rm=T)),.names="{.col}Zsc")) %>% 
  ungroup() %>% select(LNR,matches("Zsc"))

NTsc <- ReSc %>% full_join(MaSc,"LNR")
rm(MaSc);rm(ReSc)
names(NTsc) <- str_replace_all(str_replace_all(names(NTsc),"NPREG","ChMath"),"NPLES","ChEng")
NTsc <- EliID %>%select(LNR,MoLNR,FaLNR) %>% left_join(NTsc,"LNR")
data.table::fwrite(NTsc%>%select(LNR,MoLNR,FaLNR,matches("Zsc")),"Data/InitialParticipation_NarrowEligibility_TrioNTScVars.csv") 
rm(NTsc)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Merge all datasets \n")

remove(list = ls())

Df <- 
  data.table::fread("Data/InitialParticipation_NarrowEligibility_CriteriaVars_UniqueMoLNR.csv",data.table=F,na.strings=c(NA,"")) %>% 
  dplyr::select(LNR,MoLNR,FaLNR,matches("Partici"),matches("District_Rec")) %>% 
    dplyr::left_join(
      data.table::fread("Data/InitialParticipation_NarrowEligibility_TrioDemVars.csv",data.table=F,na.strings=c(NA,""))  %>% 
        dplyr::select(-matches("Partici")) %>% 
        dplyr::distinct(LNR,MoLNR,FaLNR,.keep_all=T),c("LNR","MoLNR","FaLNR")) %>% 
  dplyr::left_join(
      data.table::fread("Data/InitialParticipation_NarrowEligibility_TrioMsVars.csv",data.table=F,na.strings=c(NA,""))   %>% 
        dplyr::select(-matches("ChYOB"),-matches("Partici")) %>% 
        dplyr::distinct(LNR,MoLNR,FaLNR,.keep_all=T),c("LNR","MoLNR","FaLNR")) %>%
  dplyr::left_join(
      data.table::fread("Data/InitialParticipation_NarrowEligibility_TrioTsVars.csv",data.table=F,na.strings=c(NA,""))   %>% 
        dplyr::select(-matches("ChYOB"),-matches("Partici")) %>% 
        dplyr::distinct(LNR,MoLNR,FaLNR,.keep_all=T),c("LNR","MoLNR","FaLNR")) %>%
  dplyr::left_join(
      data.table::fread("Data/InitialParticipation_NarrowEligibility_TrioEmpVars.csv",data.table=F,na.strings=c(NA,""))  %>% 
        dplyr::select(-matches("ChYOB"),-matches("Partici")) %>% 
        dplyr::distinct(LNR,MoLNR,FaLNR,.keep_all=T),c("LNR","MoLNR","FaLNR")) %>%
  dplyr::left_join(
      data.table::fread("Data/InitialParticipation_NarrowEligibility_TrioIncVars.csv",data.table=F,na.strings=c(NA,""))  %>% 
        dplyr::select(-matches("ChYOB"),-matches("Partici")) %>% 
        dplyr::distinct(LNR,MoLNR,FaLNR,.keep_all=T),c("LNR","MoLNR","FaLNR")) %>%
  dplyr::left_join(
      data.table::fread("Data/InitialParticipation_NarrowEligibility_TrioNTScVars.csv",data.table=F,na.strings=c(NA,"")) %>% 
        dplyr::select(-matches("ChYOB"),-matches("Partici")) %>% 
        dplyr::distinct(LNR,MoLNR,FaLNR,.keep_all=T),c("LNR","MoLNR","FaLNR")) %>%
  dplyr::left_join(
      data.table::fread("Data/InitialParticipation_NarrowEligibility_TrioEduVars.csv",data.table=F,na.strings=c(NA,""))  %>% 
        dplyr::select(-matches("ChYOB"),-matches("Partici")) %>% 
        dplyr::distinct(LNR,MoLNR,FaLNR,.keep_all=T),c("LNR","MoLNR","FaLNR"))%>% 
  dplyr::mutate(across(everything(),~ifelse(.x==""|.x==" ",NA,.x)))

data.table::fwrite(Df,"Data/InitialParticipation_NarrowEligibility_TrioSSBVars.csv") 

remove(list = ls())



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Tidying main dataset \n")

Df<-data.table::fread("Data/InitialParticipation_NarrowEligibility_TrioSSBVars.csv") 
library(dplyr)
  
NAMES <- c(
  "LNR","MoLNR","FaLNR",
  #"MoMoLNR","MoFaLNR","FaMoLNR","FaFaLNR",
  "MoQ1Participation","FaQ1Participation","MoInitialParticipation","FaInitialParticipation",
  "MoBiobankParticipation","FaBiobankParticipation","ChBiobankParticipation","MoQ9Participation","FaQF2Participation",
  "MoContinuedParticipation","FaContinuedParticipation",
  "MoResidenceDistrict_Rec","FaResidenceDistrict_Rec",
  "ChSex","ChDOB","FaAge","MoAge","ParentalAgeGap",
  # "ChMOB","ChYOB","ChDOD","FaYOB","MoYOB","MoDOB","MoDOD","FaDOB","FaDOD",
  "ChDead","MoDead","MoMoDead","MoFaDead","FaDead","FaMoDead","FaFaDead",
  # "ChBirthCountry","ChNationalHeritage","ChImmigrationStatus",
  # "MoBirthCountry","MoNationalHeritage","MoImmigrationStatus",
  # "FaBirthCountry","FaNationalHeritage","FaImmigrationStatus",
  # "MoMoDOD","MoFaDOD","FaMoDOD","FaFaDOD",
  # "ChBirthContinent","MoBirthContinent","FaBirthContinent",
  # "ChImmigrationStatus01","MoImmigrationStatus01","FaImmigrationStatus01",
  "ChNumberNorwegianParents","MoNumberNorwegianParents","FaNumberNorwegianParents",
  "ChBirthPlace","MoBirthPlace","FaBirthPlace",
  "MoCivilStatus_2","FaCivilStatus_2",
  # "MoCivilStatus","FaCivilStatus","MoCivilStatus_1","FaCivilStatus_1",
  "MoEmploymentStatus","MoHoursWorked","MoHoursWorkedCategory",
  "MoOccupationalCode","MoWorkLocationMunicipality","MoIndustrySector1987",
  "MoOccupationalCodeLevel1","MoOccupationalCodeLevel2","MoOccupationalCodeLevel3",# "MoOccupationalCodeLevel4",
  #"MoCompanyMainIndustrysn94","MoCompanyMainIndustrysn07",#"MoCompanyMainIndustrysn02",
  "FaEmploymentStatus","FaHoursWorked","FaHoursWorkedCategory",
  "FaOccupationalCode","FaWorkLocationMunicipality","FaIndustrySector1987",
  "FaOccupationalCodeLevel1","FaOccupationalCodeLevel2","FaOccupationalCodeLevel3",# "FaOccupationalCodeLevel4",
  #"FaCompanyMainIndustrysn94","FaCompanyMainIndustrysn07",#"FaCompanyMainIndustrysn02",
  
  "MoProfessionalIncome_IF_LOG10","MoSalaryIncome_IF_LOG10","MoNetBusinessIncome_IF_LOG10",
  "MoCapitalIncome_IF_LOG10","MoInterestIncome_IF_LOG10","MoStockDividends_IF_LOG10","MoRealizationGains_IF_LOG10",
  "MoRealizationLoss_IF_LOG10","MoOtherCapitalIncome_IF_LOG10","MoTransfers_IF_LOG10","MoTaxableTransfers_IF_LOG10",
  "MoPensionsFromTheNationalInsurance_IF_LOG10","MoRetirementPensions_IF_LOG10","MoDisabilityPension_IF_LOG10",
  "MoWorkSettlementAllowance_IF_LOG10","MoServicePensions_IF_LOG10","MoUnemploymentBenefit_IF_LOG10","MoSickPay_IF_LOG10",
  "MoParentalAllowance_IF_LOG10","MoTaxFreeTransfers_IF_LOG10","MoChildBenefit_IF_LOG10","MoHousingBenefit_2_IF_LOG10",
  "MoStudiestipend_IF_LOG10","MoDependentDeduction_IF_LOG10","MoSocialAssistance_IF_LOG10","MoBasicAndAuxiliaryAllowance_IF_LOG10",
  "MoCashSupport_IF_LOG10","MoTotalIncome_IF_LOG10","MoEqualizedTaxAndNegativeTransfers_IF_LOG10","MoEqualizedTax_IF_LOG10",
  "MoIncomeAfterTax_IF_LOG10","MoInterestExpenses_IF_LOG10","MoProfessionalIncome_2_IF_LOG10","MoSalaryIncome_2_IF_LOG10",
  "MoNetBusinessIncome_2_IF_LOG10","MoTransfers_2_IF_LOG10","Moskpl_overf_IF_LOG10","MoServicePensions_2_IF_LOG10",
  
  "FaProfessionalIncome_IF_LOG10","FaSalaryIncome_IF_LOG10","FaNetBusinessIncome_IF_LOG10","FaCapitalIncome_IF_LOG10",
  "FaInterestIncome_IF_LOG10","FaStockDividends_IF_LOG10","FaRealizationGains_IF_LOG10","FaRealizationLoss_IF_LOG10",
  "FaOtherCapitalIncome_IF_LOG10","FaTransfers_IF_LOG10","FaTaxableTransfers_IF_LOG10",
  "FaPensionsFromTheNationalInsurance_IF_LOG10","FaRetirementPensions_IF_LOG10","FaDisabilityPension_IF_LOG10",
  "FaWorkSettlementAllowance_IF_LOG10","FaServicePensions_IF_LOG10","FaUnemploymentBenefit_IF_LOG10","FaSickPay_IF_LOG10",
  "FaParentalAllowance_IF_LOG10","FaTaxFreeTransfers_IF_LOG10","FaChildBenefit_IF_LOG10","FaHousingBenefit_2_IF_LOG10",
  "FaStudiestipend_IF_LOG10","FaDependentDeduction_IF_LOG10","FaSocialAssistance_IF_LOG10","FaBasicAndAuxiliaryAllowance_IF_LOG10",
  "FaCashSupport_IF_LOG10","FaTotalIncome_IF_LOG10","FaEqualizedTaxAndNegativeTransfers_IF_LOG10","FaEqualizedTax_IF_LOG10",
  "FaIncomeAfterTax_IF_LOG10","FaInterestExpenses_IF_LOG10","FaProfessionalIncome_2_IF_LOG10","FaSalaryIncome_2_IF_LOG10",
  "FaNetBusinessIncome_2_IF_LOG10","FaTransfers_2_IF_LOG10","Faskpl_overf_IF_LOG10","FaServicePensions_2_IF_LOG10",
  
  "MoHouseholdIncome_IF_LOG10","FaHouseholdIncome_IF_LOG10",
  "ChEng09Zsc","ChEng05Zsc","ChEng08Zsc","ChMath05Zsc","ChMath09Zsc","ChMath08Zsc",
  "MoResidenceSettlementUrbanicity","MoResidenceSettlementPopulationDensity","MoResidenceSettlementSizeKm2",
  "FaResidenceSettlementUrbanicity","FaResidenceSettlementPopulationDensity","FaResidenceSettlementSizeKm2",
  "MoEduCode","MoEduYears","FaEduCode","FaEduYears"
  )

Df1<-Df%>%
  dplyr::select(NAMES)%>%
  dplyr::mutate(across(matches("Municipality"),function(x) stringr::str_pad(x,width=4,side="left",pad="0"))) %>% 
  dplyr::mutate(across(matches("District"),function(x) stringr::str_pad(x,width=6,side="left",pad="0"))) %>% 
  dplyr::mutate(across(matches("District"),~stringr::str_sub(.x,1,4),.names="{stringr::str_replace(.col, 'District', 'Municipality')}"))%>% 
  dplyr::select(matches("LNR"),matches("Municipality"),where(~is.integer(.)||is.numeric(.)),where(~is.character(.)&&n_distinct(.)<100))%>% 
  dplyr::mutate(across(everything(),~ifelse(.x==""|.x==" ",NA,.x)))%>%
  dplyr::mutate(across(everything(),~ifelse(.x=="0000"|.x=="000000",NA,.x)))%>%
  dplyr::mutate(across(matches("EduCode"),~ifelse(!is.na(.x),paste0("Level",.x),.x)))%>% 
  dplyr::mutate(across(where(~n_distinct(.)<10),~as.character(.x)))%>% 
  dplyr::mutate(across(matches("Municipality"),~as.character(.x)))%>% 
  dplyr::mutate(across(where(is.character),function(x) factor(x,levels=names(sort(table(x),decreasing=T)))))

table(Df1$FaResidenceSettlementUrbanicity)
table(Df1$FaWorkLocationMunicipality)
hist(Df1$FaParentalAllowance_IF_LOG10)
table(is.na(Df1$FaParentalAllowance_IF_LOG10))
table(is.na(Df1$MoQ9Participation))
table((Df1$FaEduCode))
str((Df1$MoMoDead))

saveRDS(Df1,"Data/InitialParticipation_NarrowEligibility_TrioSSBVars.rds")





# What is the highest level of education that you have achieved to date? 
# 1="8th grade or less"
# 2="Some high school"
# 3="High school graduate"
# 4="Some vocational/technical training"
# 5="Completed vocational/technical training (after high school)"
# 6="Some college"
# 7="Completed college (bachelor's degree)"
# 8="Some graduate school"
# 9="Completed a master's degree"
# 10="Some graduate training beyond a master's degree"
# 11="Completed a doctoral degree"
# 12="Some post baccalaureate professional education (e.g., law school, med school, nurse)"
# 13="Completed post baccalaureate professional education"
                                                                         
# 1 = EduYears = 7
# 2 = EduYears = 10
# 3 = EduYears = 13
# 4 = EduYears = 13
# 5 = EduYears = 15
# 6 = EduYears = 13 
# 7 = EduYears = 19
# 8 = EduYears = 19
# 9 = EduYears = 19
# 10 = EduYears = 19
# 12 = EduYears = 19
# 11 = EduYears = 22
# 13 = EduYears = 22

# NUS CODES 
# 0 - No education and pre-school education
# 1 - Primary education
# 2 - Lower secondary education
# 3 - Upper secondary education, basic education
# 4 - Upper secondary, final year
# 5 - Post-secondary non-tertiary education
# 6 - First stage of tertiary education, undergraduate level
# 7 - First stage of tertiary education, graduate level
# 8 - Second stage of tertiary education (postgraduate education)
# 9 - Unspecified

# 0=0,1=7,2=10,3=13,4=14,5=15,6=17,7=19,8=22
