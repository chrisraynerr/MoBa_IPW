#!/usr/bin/Rscript

cat("\n *** Script for:", "\n *** \t - extracting kuhr variables")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Preparing workspace \n")

remove(list = ls())

source("profile.R")
pkgs <- c("data.table","dplyr","stringr","tidyr", "foreign","haven","lubridate","purrr")
for(p in pkgs){  LoadOrInstall(p) }

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Identifying directories\n")

inDir   <- paste0(KUHR_REGISTER,"data")
outDir1 <- paste0(KUHR_REGISTER,"prepared/crayner")
outDir0 <- "Data/"
dir.create(outDir0) 
dir.create(outDir1)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Creating code lists and defintions of required outcomes\n")

HypADHD  <-c("P81","F90.0","F90.1","F90.8","F90.9") 
Anxiety   <-c("P74","F41.0","F41.1","F41.3","F41.8","F41.9"); 
PhobiaOCD<-c("P79","F40.0","F40.1","F40.2","F40.8","F40.9","F42.0","F42.1","F42.2","F42.8","F42.9")
Depression<-c("P76","F32.0","F32.1","F32.2","F32.3","F32.8","F32.9","F33.0","F33.1","F33.2","F33.3",
              "F33.4","F33.8","F33.9","F34.1","F34.8","F34.9","F38.0","F38.1","F38.8","F39","F41.2","F53.0")

CODES  <-c(HypADHD,Anxiety,PhobiaOCD,Depression)
CODES2 <-paste0(",",CODES,",")
LABELS <-c("Hyperkinetic disorder", "Disturbance of activity and attention", "Hyperkinetic conduct disorder", "Other hyperkinetic disorders", "Hyperkinetic disorder, unspecified", "Anxiety disorder", "Panic disorder", "Generalized anxiety disorder", "Other mixed anxiety disorders", "Other specified anxiety disorders", "Anxiety disorder, unspecified", "Other neurotic disorders", "Agoraphobia", "Social phobias", "Specific (isolated) phobias", "Other phobic anxiety disorders", "Phobic anxiety disorder, unspecified", "Obsessive-compulsive disorder", "Predominantly obsessional thoughts or ruminations", "Predominantly compulsive acts", "Other obsessive-compulsive disorders", "Obsessive-compulsive disorder, unspecified", "Depressive disorder", "Mild depressive episode", "Moderate depressive episode", "Severe depressive episode without psychotic symptoms", "Severe depressive episode with psychotic symptoms", "Other depressive episodes", "Depressive episode, unspecified", "Recurrent depressive disorder, current episode mild", "Recurrent depressive disorder, current episode moderate", "Recurrent depressive disorder, current episode severe without psychotic symptoms", "Recurrent depressive disorder, current episode severe with psychotic symptoms", "Recurrent depressive disorder, current episode unspecified", "Other recurrent depressive disorders", "Recurrent depressive disorder, unspecified", "Dysthymia", "Other persistent mood disorders", "Persistent mood disorder, unspecified", "Other single mood disorders", "Other recurrent mood disorders", "Other specified mood disorders", "Unspecified mood disorder", "Mixed anxiety and depressive disorder", "Mental and behavioural disorders associated with the puerperium, not elsewhere classified")
LABELS <-stringr::str_to_title(LABELS)
LABELS <-stringr::str_remove_all(LABELS,"[^[:alnum:]]")
LABELS <-data.frame(Diagnosis=CODES,DiagnosisLabels=LABELS) %>% 
  mutate(DiagnosisCategory=
      ifelse(Diagnosis%in%HypADHD,"HypADHD",
        ifelse(Diagnosis%in%Anxiety,"Anxiety",
          ifelse(Diagnosis%in%Depression,"Depression",
            ifelse(Diagnosis%in%PhobiaOCD,"PhobiaOCD",NA)))))


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Loading and creating ID lists for participants\n")

idList1 <- data.table::fread(paste0(outDir0,"/AllEligiblePregnancies.csv"))
idList2 <- unique(c(idList1$LNR,idList1$MoLNR,idList1$FaLNR))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Retrieving primary diagnoses per year\n")

log   <- file(paste0(outDir0,"/Records.txt"),open="w")
PATHS <- list.files(inDir,pattern="Data fra KUHR 23-6795",full.names=T)

for (p in PATHS) {
  Y  <- stringr::str_sub(basename(p),1,4) # get the year
  
  df <- data.table::fread(p) %>% 
    dplyr::mutate(across(everything(),~ifelse(.x=="",NA,.x)))%>% 
    dplyr::mutate(across(everything(),~ifelse(.x==" ",NA,.x)))%>%
    dplyr::select(LNR=PASIENTLOPENUMMER,Diagnosis=DIAGNOSER,HealthProf=PRAKSISTYPE)%>%
    dplyr::filter(LNR%in%idList2)%>%
    dplyr::mutate(Diagnosis=stringr::str_replace_all(Diagnosis,",.*$",""))%>%
    dplyr::filter(Diagnosis%in%CODES)%>%
    dplyr::left_join(LABELS,"Diagnosis")%>% 
    dplyr::filter(!is.na(Diagnosis)) %>% 
    dplyr::select(LNR,Diagnosis=DiagnosisCategory) %>% na.omit()
  
  df  <- df %>% dplyr::filter(LNR%in%idList2)%>%dplyr::group_by(LNR,Diagnosis)%>%
    dplyr::summarise(Count=n())%>%dplyr::mutate(Year=paste0(Y))
  
  Ch <- df %>% dplyr::filter(LNR%in%idList1$LNR)
  Mo <- df %>% dplyr::filter(LNR%in%idList1$MoLNR)
  Fa <- df %>% dplyr::filter(LNR%in%idList1$FaLNR)
  
  cat(paste0("\nIn ",Y,":"),file=log,append=T)
  cat(paste0("\nThere are ", nrow(Ch)," recorded diagnoses across ",nrow(Ch%>%distinct(LNR))," children\n"),file=log,append=T);
  cat(paste0("\nThere are ", nrow(Mo)," recorded diagnoses across ",nrow(Mo%>%distinct(LNR))," mothers\n"),file=log,append=T);
  cat(paste0("\nThere are ", nrow(Fa)," recorded diagnoses across ",nrow(Fa%>%distinct(LNR))," fathers\n"),file=log,append=T);
  
  data.table::fwrite(df,paste0(outDir1,"/KuhrDiagnoses_",Y,".csv"),col.names=T,na=NA)
}
close(log)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Creating trio files per diagnosis\n")

log   <- file(paste0(outDir0,"/Prevalences.txt"),open="w")

PATHS <- list.files(outDir0,pattern="KuhrDiagnoses_[1-2]",full.names=T)

for( d in c("Depression","Anxiety","PhobiaOCD","HypADHD")){
  df0 <- 
    lapply(PATHS,fread) %>% do.call("bind_rows",.)%>%
    dplyr::group_by(LNR,Diagnosis)%>% dplyr::mutate(EarliestRecord=min(Year),TotalCount=sum(Count)) %>% 
    dplyr::ungroup() %>% dplyr::select(LNR,Diagnosis,EarliestRecord,TotalCount)%>% dplyr::distinct()%>% 
    filter(Diagnosis==paste0(d))%>%dplyr::mutate(Diagnosis=ifelse(!is.na(Diagnosis),1,NA))
  
  df1 <- idList1%>%
    dplyr::left_join(df0 %>% dplyr::select(LNR,ChDiagnosis=Diagnosis,ChEarliestRecord=EarliestRecord,ChTotalCount=TotalCount))%>%
    dplyr::left_join(df0 %>% dplyr::select(MoLNR=LNR,MoDiagnosis=Diagnosis,MoEarliestRecord=EarliestRecord,MoTotalCount=TotalCount))%>%
    dplyr::left_join(df0 %>% dplyr::select(FaLNR=LNR,FaDiagnosis=Diagnosis,FaEarliestRecord=EarliestRecord,FaTotalCount=TotalCount)) %>% 
    dplyr::mutate(across(everything(),~ifelse(is.na(.x),0,.x)))
  
  data.table::fwrite(df1, paste0(outDir1, "/KuhrDiagnoses_", d, ".csv"),col.names=T,na=NA)
  
  cat(paste0("\nPrevalences for ",d,":"),file=log,append=T)
  cat(paste0("\nChildren ", mean(df1$ChDiagnosis,na.rm=T)," \n"),file=log,append=T)
  cat(paste0("\nMothers ",  mean(df1$MoDiagnosis,na.rm=T)," \n"),file=log,append=T)
  cat(paste0("\nFathers ",  mean(df1$FaDiagnosis,na.rm=T)," \n"),file=log,append=T)
  
  cat(paste0("\nPeak incidences for ",d,":"),file=log,append=T)
  cat(paste0("\nChildren ", median(df1%>%filter(ChEarliestRecord!=0)%>%.$ChEarliestRecord,na.rm=T)," \n"),file=log,append=T)
  cat(paste0("\nMothers ",  median(df1%>%filter(MoEarliestRecord!=0)%>%.$MoEarliestRecord,na.rm=T)," \n"),file=log,append=T)
  cat(paste0("\nFathers ",  median(df1%>%filter(FaEarliestRecord!=0)%>%.$FaEarliestRecord,na.rm=T)," \n"),file=log,append=T)
}
close(log)

rm(df0);rm(df1)



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Creating 1 trio file containing all diagnoses\n")

PATHS     <- list.files(outDir1,pattern="KuhrDiagnoses_[A-Z]",full.names=T);PATHS
NAMES     <- str_remove_all(str_remove_all(str_remove_all(basename(PATHS),"KuhrDiagnoses_"),".csv"),"_");NAMES
RenameCol <- function(df,name){colnames(df)[4:12]<-paste0(colnames(df)[4:12],"_",name);return(df)}

dfL       <- lapply(PATHS,fread)
dfL       <- mapply(FUN=RenameCol,dfL,NAMES,SIMPLIFY=F)
df0       <- dfL%>%purrr::reduce(.,inner_join,by=c("LNR","MoLNR","FaLNR"))

data.table::fwrite(df0,paste0(outDir1, "/KuhrDiagnoses_Analysis.csv"),col.names=T,na=NA)



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Merging with main data \n")

X  <-  readRDS(paste0(outDir0,"/NarrowEligibility_BiasVariables.rds")) %>% 

D  <-  data.table::fread(paste0(outDir1,"/KuhrDiagnoses_Analysis.csv")) %>%
  mutate(ChDiagnosis_Internalising=pmax(ChDiagnosis_Depression,ChDiagnosis_Anxiety,ChDiagnosis_PhobiaOCD)) %>%
  mutate(MoDiagnosis_Internalising=pmax(MoDiagnosis_Depression,MoDiagnosis_Anxiety,MoDiagnosis_PhobiaOCD)) %>%
  as.data.frame() #%>%
  # select(LNR,matches("HypADHD"),matches("Internalising"))

X <- X %>% inner_join(D,by=c("LNR"))
Y <- X %>% filter(MoQ1Participation==1)
Z <- X %>% filter(MoQ9Participation==1)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Reporting prevalences \n")

log <- file(paste0(outDir0,"/DiagnosisTables.txt"),open="w")

cat(paste0("\nChildren \n"),file=log,append=T)
cat(paste0("\nFull Population \n"),file=log,append=T)
cat(paste0("\nDepression:\n"),file=log,append=T)
cat(table(X$ChDiagnosis_Depression),file=log,append=T)
#cat(table(X$MoQ1Participation),file=log,append=T)
cat(paste0("\nAnxiety:\n"),file=log,append=T)
cat(table(X$ChDiagnosis_Anxiety),file=log,append=T)
cat(paste0("\nPhobia:\n"),file=log,append=T)
cat(table(X$ChDiagnosis_PhobiaOCD),file=log,append=T)
cat(paste0("\nInternalising:\n"),file=log,append=T)
cat(table(X$ChDiagnosis_Internalising),file=log,append=T)
cat(paste0("\nADHD:\n"),file=log,append=T)
cat(table(X$ChDiagnosis_HypADHD),file=log,append=T)
cat(paste0("\nQ1 sample \n"),file=log,append=T)
cat(paste0("\nDepression:\n"),file=log,append=T)
cat(table(Y$ChDiagnosis_Depression),file=log,append=T)
cat(paste0("\nAnxiety:\n"),file=log,append=T)
cat(table(Y$ChDiagnosis_Anxiety),file=log,append=T)
cat(paste0("\nPhobia:\n"),file=log,append=T)
cat(table(Y$ChDiagnosis_PhobiaOCD),file=log,append=T)
cat(paste0("\nInternalising:\n"),file=log,append=T)
cat(table(Y$ChDiagnosis_Internalising),file=log,append=T)
cat(paste0("\nADHD:\n"),file=log,append=T)
cat(table(Y$ChDiagnosis_HypADHD),file=log,append=T)
cat(paste0("\nQ9 sample \n"),file=log,append=T)
cat(paste0("\nDepression:\n"),file=log,append=T)
cat(table(Z$ChDiagnosis_Depression),file=log,append=T)
cat(paste0("\nAnxiety:\n"),file=log,append=T)
cat(table(Z$ChDiagnosis_Anxiety),file=log,append=T)
cat(paste0("\nPhobia:\n"),file=log,append=T)
cat(table(Z$ChDiagnosis_PhobiaOCD),file=log,append=T)
cat(paste0("\nInternalising:\n"),file=log,append=T)
cat(table(Z$ChDiagnosis_Internalising),file=log,append=T)
cat(paste0("\nADHD:\n"),file=log,append=T)
cat(table(Z$ChDiagnosis_HypADHD),file=log,append=T)

cat(paste0("\nMothers\n"),file=log,append=T)
cat(paste0("\nFull Population \n"),file=log,append=T)
cat(paste0("\nDepression:\n"),file=log,append=T)
cat(table(X$MoDiagnosis_Depression),file=log,append=T)
cat(paste0("\nAnxiety:\n"),file=log,append=T)
cat(table(X$MoDiagnosis_Anxiety),file=log,append=T)
cat(paste0("\nPhobia:\n"),file=log,append=T)
cat(table(X$MoDiagnosis_PhobiaOCD),file=log,append=T)
cat(paste0("\nInternalising:\n"),file=log,append=T)
cat(table(X$MoDiagnosis_Internalising),file=log,append=T)
cat(paste0("\nADHD:\n"),file=log,append=T)
cat(table(X$MoDiagnosis_HypADHD),file=log,append=T)
cat(paste0("\nQ1 sample \n"),file=log,append=T)
cat(paste0("\nDepression:\n"),file=log,append=T)
cat(table(Y$MoDiagnosis_Depression),file=log,append=T)
cat(paste0("\nAnxiety:\n"),file=log,append=T)
cat(table(Y$MoDiagnosis_Anxiety),file=log,append=T)
cat(paste0("\nPhobia:\n"),file=log,append=T)
cat(table(Y$MoDiagnosis_PhobiaOCD),file=log,append=T)
cat(paste0("\nInternalising:\n"),file=log,append=T)
cat(table(Y$MoDiagnosis_Internalising),file=log,append=T)
cat(paste0("\nADHD:\n"),file=log,append=T)
cat(table(Y$MoDiagnosis_HypADHD),file=log,append=T)
cat(paste0("\nQ9 sample \n"),file=log,append=T)
cat(paste0("\nDepression:\n"),file=log,append=T)
cat(table(Z$MoDiagnosis_Depression),file=log,append=T)
cat(paste0("\nAnxiety:\n"),file=log,append=T)
cat(table(Z$MoDiagnosis_Anxiety),file=log,append=T)
cat(paste0("\nPhobia:\n"),file=log,append=T)
cat(table(Z$MoDiagnosis_PhobiaOCD),file=log,append=T)
cat(paste0("\nInternalising:\n"),file=log,append=T)
cat(table(Z$MoDiagnosis_Internalising),file=log,append=T)
cat(paste0("\nADHD:\n"),file=log,append=T)
cat(table(Z$MoDiagnosis_HypADHD),file=log,append=T)

close(log)

remove(list = ls())
