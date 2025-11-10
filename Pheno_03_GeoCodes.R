#!/usr/bin/Rscript

cat("\n *** Script for:","\n *** \t - re-coding geographical variables prior to GEO-code linked recruitment windows")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Preparing workspace \n")

remove(list = ls())

source("profile.R")
pkgs <- c("data.table","dplyr","stringr","tidyr", "foreign","haven","lubridate","purrr")
for(p in pkgs){  LoadOrInstall(p) }

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Load linkage file and participation phenotypes\n")

EliID <- data.table::fread("Data/InitialParticipation_BroadEligibility_TrioLNR.csv")

Window<-EliID%>%select(LNR,MoLNR,FaLNR,ChDOB)%>%mutate(Year_1=as.integer(str_sub(ChDOB,1,4)),Year_0=Year_1-1,Year_2=Year_1+1)%>%pivot_longer(cols=matches("Year"),values_to="Year",names_to="Occasion")%>%mutate(Occasion=str_remove_all(Occasion,"Year_"))
WindowLong<-Window %>% select(ChLNR=LNR,LNR=MoLNR,Year) %>% bind_rows(Window %>% select(ChLNR=LNR,LNR=FaLNR,Year))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Obtain longitudinal geographical variables\n")

GeoPath <- list.files(paste0(SSB_REGISTER,"SSB/01_data/data_v5.0/csv/"),pattern="POPULATION_GRUNNKRETS_",full.names = T)
GeoPath <- GeoPath[which(grepl(x=GeoPath,"1998|1999|200|2010|2011"))]
Names   <- str_remove_all(basename(GeoPath),"POPULATION_GRUNNKRETS_|.csv")

freadFilter <- function(x, i){
  Df  <- fread(x,header=T,data.table=F,colClasses="character",na.strings=c(NA,"NA","")) %>% na.omit() %>%
    purrr::set_names("LNR","Kommune","Grunnkretser","Delomraade") %>%
    mutate(Year=as.integer(paste0(i))) %>% 
    inner_join(WindowLong) %>% 
    mutate(Municipality=str_pad(Kommune,width=4,side="left",pad="0")) %>%
    mutate(Delomraade=str_pad(Delomraade,width=2,side="left",pad="0")) %>%
    mutate(District = paste0(Municipality,Delomraade)) %>%
    select(ChLNR,LNR,District) %>% 
    mutate(across(everything(),~ifelse(.x=="000000",NA,.x))) %>% 
    filter(!is.na(District))
}

GeoLS <- lapply(seq_along(GeoPath), function(z) freadFilter(GeoPath[z],Names[z]))
GeoLS <- GeoLS[sapply(GeoLS, nrow) > 0]
GeoDf<-purrr::reduce(GeoLS,bind_rows)%>%filter(!is.na(District))%>%group_by(ChLNR)%>%mutate(N=row_number())%>%ungroup()%>%pivot_wider(id_cols=c(ChLNR,LNR),names_from=N,names_prefix="District",values_from=District)
get3<-function(x){x[!is.na(x)][1:3]}
GeoDf2<-GeoDf%>%select(matches("LNR"))%>%bind_cols(GeoDf%>%select(-matches("LNR"))%>%apply(1,get3)%>%t()%>%as.data.frame()%>%setNames(paste0("District",1:3)))
GeoDf<-GeoDf2
rm(GeoLS,GeoDf2)
GeoDf<-GeoDf%>%mutate(across(everything(),~ifelse(.x=="000000",NA,.x)))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n *** Importing data on how Geo codes have changed over time")

GeoCodes<-readxl::read_excel("DataDictionary/MunicipalityCodes_Changes.xlsx") %>% 
  tidyr::drop_na()%>%dplyr::filter(oldName!="Nothing from API")%>% 
  dplyr::mutate(oldMunicipality=substr(oldCode,1,4))%>%
  dplyr::mutate(newMunicipality=substr(currentCode,1,4))%>% 
  dplyr::mutate(oldDistrict=substr(oldCode,1,6))%>%
  dplyr::mutate(newDistrict=substr(currentCode,1,6))%>% 
  dplyr::select(oldMunicipality,newMunicipality,oldDistrict,newDistrict,changeOccurred) |> 
  dplyr::distinct()|>dplyr::filter(oldDistrict!=newDistrict) %>% 
  dplyr::arrange(oldDistrict,desc(changeOccurred)) %>% 
  dplyr::distinct(oldDistrict,.keep_all=T)

DisCodes<-GeoCodes%>%group_by(newDistrict)%>%summarize(oldDistrictCodes=list(oldDistrict))
MunCodes<-GeoCodes%>%group_by(newMunicipality)%>%summarize(oldMunicipalityCodes=list(oldMunicipality))

UpdDistrict<-function(x){for(i in 1:nrow(GeoCodes)){
  old<-GeoCodes$oldDistrictCodes[[i]];new<-GeoCodes$newDistrict[i];x[x%in%old]<-new};return(x)
}
UpdMunicipality<-function(x){for(i in 1:nrow(GeoCodes)){
  old<-GeoCodes$oldMunicipalityCodes[[i]];new<-GeoCodes$newMunicipality[i];x[x%in%old]<-new};return(x)
}

DisDf<-GeoDf%>%mutate(across(3:ncol(.),UpdDistrict))
MunDf<-GeoDf%>%mutate(across(3:ncol(.),~substr(.x,1,4)))%>%mutate(across(3:ncol(.),UpdMunicipality))

length(unique(unlist(DisDf%>%select(matches("District")))))
length(unique(unlist(DisDf%>%select(matches("District"))%>%mutate(across(everything(),~substr(.x,1,4))))))

length(unique(unlist(MunDf%>%select(matches("District")))))
length(unique(unlist(MunDf%>%select(matches("District"))%>%mutate(across(everything(),~substr(.x,1,4))))))

MunDf <- NULL

getmode <- function(v) {uniqv<-unique(v);uniqv[which.max(tabulate(match(v,uniqv)))]}

MoGeo<-
  EliID %>% select(LNR,MoLNR,ChDOB)%>%inner_join(DisDf%>%select(MoLNR=LNR,LNR=ChLNR,matches("District")))%>%
  rename_with(~str_replace_all(.,"District","MoResidenceDistrict_"))%>%
  mutate(MoResidenceDistrict_Mode=apply(.[,3:5],1,getmode))

FaGeo<-
  EliID %>% select(LNR,FaLNR)%>%inner_join(DisDf%>%select(FaLNR=LNR,LNR=ChLNR,matches("District")))%>%
  rename_with(~str_replace_all(.,"District","FaResidenceDistrict_"))%>%
  mutate(FaResidenceDistrict_Mode=apply(.[,3:5],1,getmode))

GeoDf<-MoGeo%>%full_join(FaGeo,c("LNR"))%>%distinct(LNR,.keep_all=T)

GeoDf<-GeoDf%>%
  mutate(ChMOB=as.integer(str_sub(ChDOB,5,6)))%>%
  mutate(MoResidenceDistrict_Rec=ifelse(ChMOB>5,MoResidenceDistrict_2,MoResidenceDistrict_1))%>% 
  mutate(FaResidenceDistrict_Rec=ifelse(ChMOB>5,FaResidenceDistrict_2,FaResidenceDistrict_1))%>%
  mutate(MoResidenceDistrict_Rec=ifelse(is.na(MoResidenceDistrict_Rec),MoResidenceDistrict_Mode,MoResidenceDistrict_Rec))%>%
  mutate(FaResidenceDistrict_Rec=ifelse(is.na(FaResidenceDistrict_Rec),FaResidenceDistrict_Mode,FaResidenceDistrict_Rec))

GeoDf<-GeoDf%>%select(LNR,MoLNR,FaLNR,MoResidenceDistrict_Rec,MoResidenceDistrict_Mode,FaResidenceDistrict_Rec,FaResidenceDistrict_Mode)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ### Saving geographical codes for all eligible pregnancies in population \n")

data.table::fwrite(GeoDf,"Data/InitialParticipation_BroadEligibility_TrioGeoVars.csv") 

remove(list = ls())


