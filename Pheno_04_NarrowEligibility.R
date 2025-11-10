#!/usr/bin/Rscript


cat("\n *** Script for:",
    "\n *** \t - creating a NARROW subset of the BROADLY eligible participants in MoBa (to get closer to the 40% reported")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ### Preparing workspace \n")
remove(list = ls())

source("profile.R")
pkgs <- c("data.table","dplyr","stringr","tidyr", "readr","haven","psych","purrr")
for(p in pkgs){  LoadOrInstall(p) }

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... Loading and merging datasets \n")

Df<-data.table::fread("Data/InitialParticipation_BroadEligibility_TrioLNR.csv",na.strings=c(NA,""))%>% 
  mutate(across(everything(),~ifelse(.x==""|.x==" ",NA,.x)))

Part<-data.table::fread("Data/MoBa_ParticipationOutcomes.csv",na.strings=c(NA,"NA","")) %>% select(LNR,matches("Part"))%>% 
  mutate(across(everything(),~ifelse(.x==""|.x==" ",NA,.x)))

Df<-Df%>%left_join(Part,"LNR")%>%mutate(across(matches("Participation"),~ifelse(is.na(.x),0,.x)))

Geo<-data.table::fread("Data/InitialParticipation_BroadEligibility_TrioGeoVars.csv",na.strings=c(NA,""))%>% 
  mutate(across(everything(),~ifelse(.x==""|.x==" ",NA,.x)))

Df<-Df%>%full_join(Geo %>% select(-MoLNR,-FaLNR),by="LNR") %>%
  mutate(MoResidenceDistrict_Rec=ifelse(is.na(MoResidenceDistrict_Rec),MoResidenceDistrict_Mode,MoResidenceDistrict_Rec))%>%
  mutate(FaResidenceDistrict_Rec=ifelse(is.na(FaResidenceDistrict_Rec),FaResidenceDistrict_Mode,FaResidenceDistrict_Rec))

rm(Geo,Part)

length(unique(Df$MoResidenceDistrict_Rec))
table(is.na(Df$MoResidenceDistrict_Rec),is.na(Df$FaResidenceDistrict_Rec))
table(is.na(Df$MoResidenceDistrict_Mode),is.na(Df$MoResidenceDistrict_Mode))

cat(paste0(
  "# N=",nrow(Df),
  "; 0=",table(Df$MoQ1Participation,useNA="a")[1]," (",round(prop.table(table(Df$MoQ1Participation,useNA="a")),2)[1],
  "); 1=",table(Df$MoQ1Participation,useNA="a")[2]," (",round(prop.table(table(Df$MoQ1Participation,useNA="a")),2)[2],")"
   ))
## N=581526; 0=477778 (0.82); 1=103748 (0.18)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... Filtering dataset \n")

## FILTER OUT INDIVIUALS WITHOUT KONWN RESIDENCE IN 3 YEARS AROUND  BIRTH 

Df0  <- Df %>% dplyr::filter(is.na(MoResidenceDistrict_Rec) & is.na(FaResidenceDistrict_Rec) & MoInitialParticipation==0) 
Df1  <- Df %>% dplyr::filter(is.na(MoResidenceDistrict_Rec) & is.na(FaResidenceDistrict_Rec) & MoInitialParticipation==1) 
Df   <- Df %>% dplyr::filter(!LNR %in% Df0$LNR) %>% dplyr::distinct(LNR,.keep_all=T)

cat(paste0(
  "# N=",nrow(Df),
  "; 0=",table(Df$MoQ1Participation,useNA="a")[1]," (",round(prop.table(table(Df$MoQ1Participation,useNA="a")),2)[1],
  "); 1=",table(Df$MoQ1Participation,useNA="a")[2]," (",round(prop.table(table(Df$MoQ1Participation,useNA="a")),2)[2],")"
))
## N=575151; 0=471403 (0.82); 1=103748 (0.18)

## FILTER OUT INDIVIUALS WHO WERE BORN IN DISTRICTS WHERE NO MOBA PARTICIPANTS WERE BORN

Df<-Df %>% mutate(ChDOB=as.integer(ChDOB))

## table(Df$ChMOB,Df$MoInitialParticipation,useNA="a");table(Df$ChYOB,Df$MoInitialParticipation,useNA="a")

sum(is.na(Df$ChDOB))
min(Df[which(Df$MoInitialParticipation==1),]$ChDOB,na.rm=T)
max(Df[which(Df$MoInitialParticipation==1),]$ChDOB,na.rm=T)

RecrDates_ResidenceDistrict <- 
  Df %>% 
  dplyr::filter(MoInitialParticipation==1) %>% 
  dplyr::group_by(MoResidenceDistrict_Rec ) %>% #MoDelomrade_T0) %>%
  dplyr::mutate(MinDob = min(ChDOB)) %>% # Create variable ("min_dob_mun") identifying first MoBa birth in each municipality to indicate start date of recruitment for each municipality
  dplyr::mutate(MaxDob = max(ChDOB)) %>% # Create variable ("min_dob_mun") identifying first MoBa birth in each municipality to indicate start date of recruitment for each municipality
  tidyr::fill(MinDob,.direction = "downup") %>% # Only MoBa children have a municipality recruitment start date value - impute for rest of children in population 
  tidyr::fill(MaxDob,.direction = "downup") %>% # Only MoBa children have a municipality recruitment start date value - impute for rest of children in population 
  dplyr::ungroup()%>%
  dplyr::select(MoResidenceDistrict_Rec,MinDob,MaxDob) %>% 
  dplyr::distinct()

RecrDates0  <- 
  RecrDates_ResidenceDistrict %>% 
  dplyr::select(MoResidenceDistrict_Rec,MinDob,MaxDob) %>% dplyr::distinct() %>% 
  dplyr::mutate(MinDob_date = as.character(MinDob)) %>% 
  dplyr::mutate(MinDob_date = as.Date(paste0(substr(MinDob_date,1,4),"-",substr(MinDob_date,5,6),"-01"))) %>% 
  dplyr::mutate(MaxDob_date = as.character(MaxDob)) %>% 
  dplyr::mutate(MaxDob_date = as.Date(paste0(substr(MaxDob_date,1,4),"-",substr(MaxDob_date,5,6),"-01"))) %>% 
  dplyr::mutate(Window = as.integer(str_remove_all(difftime(MaxDob_date,MinDob_date,units="weeks")," weeks"))/4.345) %>% 
  dplyr::mutate(Window = ifelse(Window==0,1,round(Window,1)))

data.table::fwrite(RecrDates0 %>% select(MoResidenceDistrict_Rec,MinDob,MaxDob,Window) %>% distinct(), 
                   "Data/InitialParticipation_NarrowEligibility_RecruitmentWindow.csv")

Df1 <- Df %>% dplyr::filter(MoInitialParticipation==1) 

Df0 <- Df %>% dplyr::filter(MoInitialParticipation==0) %>% 
  dplyr::full_join(RecrDates0)  

Df0 <- Df0 %>% 
  filter(!is.na(Window)) %>% 
  dplyr::group_by(MoResidenceDistrict_Rec)  %>% 
  dplyr::filter(ChDOB >= MinDob) %>% 
  dplyr::filter(ChDOB <= MaxDob) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(names(Df1))

Df   <- Df1 %>% bind_rows(Df0) 
rm(Df0);rm(Df1)

cat(paste0(
  "# N=",nrow(Df),
  "; 0=",table(Df$MoQ1Participation,useNA="a")[1]," (",round(prop.table(table(Df$MoQ1Participation,useNA="a")),2)[1],
  "); 1=",table(Df$MoQ1Participation,useNA="a")[2]," (",round(prop.table(table(Df$MoQ1Participation,useNA="a")),2)[2],")"
))
## N=417175; 0=313436 (0.75); 1=103739 (0.25)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... Saving dataset \n")

data.table::fwrite(Df, "Data/InitialParticipation_NarrowEligibility_CriteriaVars_MultiPregnancy.csv") 


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... Applying additional filtering \n")

## FILTER INDIVIUALS WITH MULTIPLE PREGNANCIES PRIORITISING PARTICIPATION / CHILDS AGE

Df<-data.table::fread("Data/InitialParticipation_NarrowEligibility_CriteriaVars_1Pregnancy.csv") 

Na <- rowMeans(is.na(Df)); Df$PropMissing <- Na; hist(Na)
Df <- Df %>% dplyr::arrange(MoLNR,-MoInitialParticipation,-MoContinuedParticipation,ChDOB) %>% dplyr::distinct(MoLNR,.keep_all=T)
hist(Df$PropMissing)

cat(paste0(
  "# N=",nrow(Df),
  "; 0=",table(Df$MoQ1Participation,useNA="a")[1]," (",round(prop.table(table(Df$MoQ1Participation,useNA="a")),2)[1],
  "); 1=",table(Df$MoQ1Participation,useNA="a")[2]," (",round(prop.table(table(Df$MoQ1Participation,useNA="a")),2)[2],")"
))

# N=296987; 0=210523 (0.71); 1=86464 (0.29)

data.table::fwrite(Df %>% dplyr::select(-PropMissing), "Data/InitialParticipation_NarrowEligibility_CriteriaVars_UniqueMoLNR.csv") 




## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... Plotting nation-wide patterns of participation \n")

# BY DISTRICT

Window <- fread("Data/InitialParticipation_NarrowEligibility_RecruitmentWindow.csv")  %>% 
  dplyr::mutate(across(matches("District"),function(x) stringr::str_pad(x,width=6,side="left",pad="0"))) %>% 
  dplyr::mutate(across(matches("District"),function(x) ifelse(x==000000,NA,x)))

test   <- fread("Data/InitialParticipation_NarrowEligibility_CriteriaVars_UniqueMoLNR.csv")  %>% 
  dplyr::mutate(across(matches("District"),function(x) stringr::str_pad(x,width=6,side="left",pad="0"))) %>% 
  dplyr::mutate(across(matches("District"),function(x) ifelse(x==000000,NA,x))) %>% 
  dplyr::select(LNR,MoQ1Participation,MoResidenceDistrict_Rec) %>% inner_join(Window) %>% 
  group_by(MoResidenceDistrict_Rec) %>% mutate(ProportionParticipated = mean(MoQ1Participation)) %>% 
  ungroup() %>% select(MoResidenceDistrict_Rec,MinDob,MaxDob,Window,ProportionParticipated) %>% distinct()

data.table::fwrite(test, "Data/InitialParticipation_ProportionParticipatedByDistrict.csv") 

# BY MUNICIPALITY

Window <- fread("Data/InitialParticipation_NarrowEligibility_RecruitmentWindow.csv")  %>% 
  dplyr::mutate(across(matches("District"),function(x) stringr::str_pad(x,width=6,side="left",pad="0"))) %>% 
  dplyr::mutate(across(matches("District"),function(x) ifelse(x==000000,NA,x))) %>% 
  dplyr::mutate(across(matches("District"),function(x) stringr::str_sub(x,1,4),
                       .names="{stringr::str_replace(.col,'District','Municipality')}"))

test   <- fread("Data/InitialParticipation_NarrowEligibility_CriteriaVars_UniqueMoLNR.csv")  %>% 
  dplyr::mutate(across(matches("District"),function(x) stringr::str_pad(x,width=6,side="left",pad="0"))) %>% 
  dplyr::mutate(across(matches("District"),function(x) ifelse(x==000000,NA,x))) %>% 
  dplyr::mutate(across(matches("District"),function(x) stringr::str_sub(x,1,4),
                       .names="{stringr::str_replace(.col,'District','Municipality')}")) %>% 
  dplyr::select(LNR,MoQ1Participation,MoResidenceMunicipality_Rec) %>% 
  group_by(MoResidenceMunicipality_Rec) %>% 
  mutate(ProportionParticipated = mean(MoQ1Participation,na.rm=T)) %>% 
  mutate(NumberParticipated = sum(MoQ1Participation,na.rm=T)) %>% 
  ungroup() %>%
  inner_join(
    Window %>% select(-MoResidenceDistrict_Rec)%>%group_by(MoResidenceMunicipality_Rec) %>% 
      mutate(MinDob = min(MinDob,na.rm=T)) %>% mutate(MaxDob = max(MaxDob,na.rm=T)) %>% 
      mutate(Window = MaxDob-MinDob) %>% ungroup() %>% distinct()
    ) %>% 
  select(MoResidenceMunicipality_Rec,MinDob,MaxDob,Window,ProportionParticipated,NumberParticipated) %>% 
  distinct()

N <- sum(test$NumberParticipated)
test$SampleProportion <- test$NumberParticipated/N
data.table::fwrite(test, "Data/InitialParticipation_ByMunicipality.csv") 


