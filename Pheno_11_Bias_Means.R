#!/usr/bin/Rscript
cat("\n *** Script for:","\n *** \t - quantifying bias in means and prevalences")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Preparing workspace \n")

remove(list = ls())

source("profile.R")
source("Pheno_Functions.R")

pkgs <- c("data.table","dplyr","stringr","tidyr","ggplot2","survey", "flextable","officer")

for(p in pkgs){  LoadOrInstall(p) }

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... loading data \n")

X<-readRDS("./Data/NarrowEligibility_BiasVariables.rds")

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... Calculating Bonferroni \n")

pca<-prcomp(X%>%select(-LNR,-matches("Participation"),-matches("IPW")),scale=TRUE)
var<-cumsum(pca$sdev^2/sum(pca$sdev^2))
nPC<-min(which(var>=0.995))
BONF<-0.05/nPC
BONF<-0.05/(nPC+2)
BONF<-as.numeric(formatC(BONF,format="e",digits=2))

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... loading weights \n")

PP_PopQ1<-fread("Predict/MoQ1Participation_202412/MoQ1Participation_FullPred.csv")%>%
  select(LNR,PP=pred)%>%filter(LNR%in%X$LNR)%>%
  mutate(IPW=1/PP)

PP_PopQ9<-fread("Predict/MoQ9Participation_202412/MoQ9Participation_FullPred.csv")%>%
  select(LNR,PP=pred)%>%filter(LNR%in%X$LNR)%>% 
  mutate(IPW=1/PP)

PP_Q1Q9<-fread("Predict/MoQ9ParticipationInMoQ1Participation_202412/MoQ9Participation_FullPred.csv")%>%
  select(LNR,PP=pred)%>%filter(LNR%in%X$LNR)%>% 
  mutate(IPW=1/PP)

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... merging data & weights \n")

df_Pop<- X %>% select(-matches("Participation"))%>%as.data.frame()
df_Q1 <- X %>% filter(MoQ1Participation==1)%>%select(-matches("Participation"))%>%as.data.frame()
df_Q9 <- X %>% filter(MoQ9Participation==1)%>%select(-matches("Participation"))%>%as.data.frame()

df_PopQ1<-df_Q1%>%inner_join(PP_PopQ1,"LNR")
df_PopQ9<-df_Q9%>%inner_join(PP_PopQ9,"LNR")
df_Q1Q9<-df_Q9%>%inner_join(PP_Q1Q9,"LNR")


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... estimating unweighted means \n")

VARS<-names(df_Pop %>% select(-LNR))

Pop<-bind_rows(lapply(VARS,function(i){X<-df_Pop[,i];GetUnwMean(X=X,I=paste0(i))}))
BL<-bind_rows(lapply(VARS,function(i){X<-df_Q1[,i];GetUnwMean(X=X,I=paste0(i))}))
FU<-bind_rows(lapply(VARS,function(i){X<-df_Q9[,i];GetUnwMean(X=X,I=paste0(i))}))

Df1 <- 
  Pop %>%
    inner_join(BL,"Variable")%>%inner_join(FU,"Variable") %>%
    select(Variable,
           Pop_Mean=all_of(2),Pop_SD=all_of(3),Pop_N=all_of(4),
           Q1_Mean=all_of(6),Q1_SD=all_of(7),Q1_N=all_of(8),
           Q9_Mean=all_of(10),Q9_SD=all_of(11),Q9_N=all_of(12)
    )

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... comparing unweighted means \n")

Pop.BL<-PairwiseComparisons(Pop,BL) %>% mutate(Sample="BL",Ref="Pop",Weight="Unweighted")
Pop.FU<-PairwiseComparisons(Pop,FU)%>% mutate(Sample="FU",Ref="Pop",Weight="Unweighted")
BL.FU<-PairwiseComparisons(BL,FU)%>% mutate(Sample="FU",Ref="BL",Weight="Unweighted")

Df2 <- 
  Pop.BL%>%select(Variable,PopQ1_Zp=Z_plog10,PopQ1_Fp=F_plog10) %>% 
  inner_join(
    Pop.FU%>%select(Variable,PopQ9_Zp=Z_plog10,PopQ9_Fp=F_plog10)) %>% 
  inner_join(
    BL.FU%>%select(Variable,Q1Q9_Zp=Z_plog10,Q1Q9_Fp=F_plog10)
)

ssb <- Df1 %>% inner_join(Df2,"Variable")


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Because KUHR is restricted, this was done by the data amanger which is why there is a separate dataset
##  First create a dataframe in the same format as ssb and then join

kuhr<- fread("N:/durable/projects/crayner/GenPar/Bias/KuhrVars_UnweightedMeans.csv")
Pop_kuhr<- fread("N:/durable/projects/crayner/GenPar/Bias/KuhrVarsPop_UnweightedMeans.csv")
BL_kuhr<- fread("N:/durable/projects/crayner/GenPar/Bias/KuhrVarsQ1_UnweightedMeans.csv")
FU_kuhr<- fread("N:/durable/projects/crayner/GenPar/Bias/KuhrVarsQ9_UnweightedMeans.csv")

Df1 <- Pop_kuhr%>%
  inner_join(BL_kuhr,"Variable")%>%inner_join(FU_kuhr,"Variable") %>%
  select(Variable,
         Pop_Mean=all_of(2),Pop_SD=all_of(3),Pop_N=all_of(4),
         Q1_Mean=all_of(6),Q1_SD=all_of(7),Q1_N=all_of(8),
         Q9_Mean=all_of(10),Q9_SD=all_of(11),Q9_N=all_of(12)
  )

Pop.BL_kuhr<-PairwiseComparisons(Pop_kuhr,BL_kuhr) %>% mutate(Sample="BL",Ref="Pop",Weight="Unweighted")
Pop.FU_kuhr<-PairwiseComparisons(Pop_kuhr,FU_kuhr)%>% mutate(Sample="FU",Ref="Pop",Weight="Unweighted")
BL.FU_kuhr<-PairwiseComparisons(BL_kuhr,FU_kuhr)%>% mutate(Sample="FU",Ref="BL",Weight="Unweighted")

Df2 <- 
  Pop.BL_kuhr%>%select(Variable,PopQ1_Zp=Z_plog10,PopQ1_Fp=F_plog10) %>% 
  inner_join(
    Pop.FU_kuhr%>%select(Variable,PopQ9_Zp=Z_plog10,PopQ9_Fp=F_plog10)) %>% 
  inner_join(
    BL.FU_kuhr%>%select(Variable,Q1Q9_Zp=Z_plog10,Q1Q9_Fp=F_plog10)
  )

kuhr <- Df1 %>% inner_join(Df2,"Variable")



### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... tabulating unweighted means: T \n")

UnweightedMeans_T <-
  ssb %>% 
  bind_rows(kuhr) %>%
  mutate(Who = ifelse(grepl(x=Variable,"^Mo"),"Mothers",ifelse(grepl(x=Variable,"^Fa"),"Fathers","Children"))) %>%
  mutate(Who = factor(Who, levels=c("Mothers","Fathers","Children"))) %>% 
  filter(Who!="Fathers") %>% 
  mutate(Variable = str_remove_all(Variable,"Mo|Fa|Ch")) %>%
  mutate(Variable = str_replace_all(Variable,"ResidenceSettlement","Hometown ")) %>%
  mutate(Variable = str_replace_all(Variable,"Residence","Hometown ")) %>%
  mutate(Variable = str_replace_all(Variable,"HouseholdIncome_IF","Household Income")) %>%
  mutate(Variable = str_replace_all(Variable,"Edu","Education")) %>%
  mutate(Variable = str_replace_all(Variable,"Math","Math test scores")) %>%
  mutate(Variable = str_replace_all(Variable,"Eng","English test scores")) %>%
  mutate(Variable = str_replace_all(Variable,"Zsc","")) %>%
  mutate(Variable = str_replace_all(Variable,"Status","Status (")) %>%
  mutate(Variable = str_replace_all(Variable,"Male","(Male")) %>%
  mutate(Variable = str_replace_all(Variable,"([a-z])([A-Z])","\\1 \\2")) %>%
  mutate(Variable = str_replace_all(Variable," IF",""))%>%
  mutate(Variable = str_replace_all(Variable,"Hyp",""))%>%
  mutate(Variable = str_replace_all(Variable,"Diagnosis",""))%>%
  mutate(Variable = str_replace_all(Variable,"_"," "))%>%
  mutate(Variable = str_replace_all(Variable,"^ ",""))%>%
  mutate(Variable = str_replace_all(Variable,"YOB","Year of birth"))%>%
  mutate(Variable = str_replace_all(Variable,"(\\()(.*)([a-z,2]$)","\\1\\2\\3\\)")) %>%
  mutate(Variable = str_replace_all(Variable,"\\( ","\\(")) %>%
  filter(!grepl(x=Variable,"Mental|Urban|Earner")) %>%
  mutate(Order =
    ifelse(grepl(x=Variable,"Age|Sex|birth"),"A", ifelse(grepl(x=Variable,"Born"),"B",
    ifelse(grepl(x=Variable,"Parents",fixed=T),"C",ifelse(grepl(x=Variable,"Grandparents"),"D",
    ifelse(grepl(x=Variable,"Civil"),"E", ifelse(grepl(x=Variable,"Education|University|test"),"F",
    ifelse(grepl(x=Variable,"Income"),"G", ifelse(grepl(x=Variable,"Employment"),"H",
    ifelse(grepl(x=Variable,"Hours"),"I", ifelse(grepl(x=Variable,"Worked"),"J",
    ifelse(grepl(x=Variable,"Hometown"),"K","L"
    )))))))))))) %>% 
  arrange(Who,Order) %>% 
  select(-Order) %>% 
  select(Who,everything())

names(UnweightedMeans_T)<-str_replace_all(names(UnweightedMeans_T),"Q1","BL")
names(UnweightedMeans_T)<-str_replace_all(names(UnweightedMeans_T),"Q9","FU")

UnweightedMeans_T <- UnweightedMeans_T %>% mutate(Variable=factor(Variable,levels=unique(UnweightedMeans_T$Variable)))
 
names(UnweightedMeans_T)<-str_replace_all(names(UnweightedMeans_T),"Q1","BL")
names(UnweightedMeans_T)<-str_replace_all(names(UnweightedMeans_T),"Q9","FU")

saveRDS(UnweightedMeans_T,"Bias/SsbVars_UnweightedMeans.rds")

library(dplyr)

UnweightedMeans_T2 <- 
  UnweightedMeans_T %>% 
  mutate(across(matches("Mean|SD|Zp|Fp",~as.numeric))) %>% 
  mutate(across(matches("Mean|SD",~round(.x,digits=3)))) %>% 
  mutate(across(matches("Zp|Fp",~round(.x,digits=2))))


WordTable(
  UnweightedMeans_T %>% 
    mutate(across(matches("Mean|SD",function(x) round(x,3)))) %>% 
    mutate(across(matches("Zp|Fp",function(x) format(x, digits = 2, scientific = TRUE))))
    ,
  "Bias/SsbVars_UnweightedMeans.docx"
  )


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... tabulating unweighted means: FOR PLOTS \n")

UnweightedMeans_P <- 
  bind_rows(list(
    Pop.BL,Pop.FU,BL.FU,
    Pop.BL_kuhr,Pop.FU_kuhr,BL.FU_kuhr)) %>% 
  select(Variable,Cohens,Cohens_2.5,Cohens_97.5,Z_stat,P=Z_p,Sample,Ref)  %>% 
  mutate(Who = ifelse(grepl(x=Variable,"^Mo"),"Mothers",ifelse(grepl(x=Variable,"^Fa"),"Fathers","Children"))) %>%
  mutate(Who = factor(Who, levels=c("Mothers","Fathers","Children"))) %>% 
  filter(Who!="Fathers") %>% 
  mutate(Variable = str_remove_all(Variable,"Mo|Fa|Ch")) %>%
  mutate(Variable = str_replace_all(Variable,"ResidenceSettlement","Hometown ")) %>%
  mutate(Variable = str_replace_all(Variable,"Residence","Hometown ")) %>%
  mutate(Variable = str_replace_all(Variable,"HouseholdIncome_IF","Household Income")) %>%
  mutate(Variable = str_replace_all(Variable,"Edu","Education")) %>%
  mutate(Variable = str_replace_all(Variable,"Math","Math test scores")) %>%
  mutate(Variable = str_replace_all(Variable,"Eng","English test scores")) %>%
  mutate(Variable = str_replace_all(Variable,"Zsc","")) %>%
  mutate(Variable = str_replace_all(Variable,"Status","Status (")) %>%
  mutate(Variable = str_replace_all(Variable,"Male","(Male")) %>%
  mutate(Variable = str_replace_all(Variable,"([a-z])([A-Z])","\\1 \\2")) %>%
  mutate(Variable = str_replace_all(Variable," IF",""))%>%
  mutate(Variable = str_replace_all(Variable,"Hyp",""))%>%
  mutate(Variable = str_replace_all(Variable,"Diagnosis",""))%>%
  mutate(Variable = str_replace_all(Variable,"_"," "))%>%
  mutate(Variable = str_replace_all(Variable,"YOB","Year of birth"))%>%
  mutate(Variable = str_replace_all(Variable,"(\\()(.*)([a-z,2]$)","\\1\\2\\3\\)")) %>%
  mutate(Variable = str_replace_all(Variable,"\\( ","\\("))%>%
  mutate(Variable = str_replace_all(Variable,"_"," "))%>%
  filter(!grepl(x=Variable,"Mental|Urban|Earner")) %>%
  mutate(Order =
    ifelse(grepl(x=Variable,"Age|Sex|birth"),"A", ifelse(grepl(x=Variable,"Born"),"B",
    ifelse(grepl(x=Variable,"Parents",fixed=T),"C",ifelse(grepl(x=Variable,"Grandparents"),"D",
    ifelse(grepl(x=Variable,"Civil"),"E", ifelse(grepl(x=Variable,"Education|University|test"),"F",
    ifelse(grepl(x=Variable,"Income"),"G", ifelse(grepl(x=Variable,"Employment"),"H",
    ifelse(grepl(x=Variable,"Hours"),"I", ifelse(grepl(x=Variable,"Worked"),"J",
    ifelse(grepl(x=Variable,"Hometown"),"K","L"
    )))))))))))) %>% arrange(Order) %>% select(-Order) %>% select(Who,everything())%>% #,Q1Q9))
  # mutate(`P-value`=ifelse(P<as.numeric(BONF),paste0("<",BONF),paste0(">",BONF))) %>% 
  # mutate(`P-value`=factor(`P-value`,levels=c(paste0("<",BONF),paste0(">",BONF))))%>% #,Q1Q9))
    mutate(`P-value`=ifelse(P<as.numeric(BONF),paste0("P<",BONF),paste0("P>",BONF))) %>% 
  mutate(`P-value`=factor(`P-value`,levels=c(paste0("P<",BONF),paste0("P>",BONF))))

UnweightedMeans_P <- UnweightedMeans_P %>% mutate(Variable=factor(Variable,levels=unique(UnweightedMeans_P$Variable)))

saveRDS(UnweightedMeans_P,"Bias/SsbVars_UnweightedMeansForPlots.rds")






### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... estimating weighted means \n")

VARS<-names(df_PopQ1%>%select(-IPW,-LNR))

w<-df_PopQ1$IPW
Q1w<-bind_rows(lapply(VARS,function(i){X<-df_PopQ1[,i];GetWeiMean(X=X,I=paste0(i),W=w)}))

w<-df_PopQ9$IPW
Q9w1<-bind_rows(lapply(VARS,function(i){X<-df_PopQ9[,i];GetWeiMean(X=X,I=paste0(i),W=w)}))

w<-df_Q1Q9$IPW
Q9w2<-bind_rows(lapply(VARS,function(i){X<-df_Q1Q9[,i];GetWeiMean(X=X,I=paste0(i),W=w)}))

Df1 <- 
  Pop%>%
  inner_join(Q1w,"Variable")%>%
  inner_join(Q9w1,"Variable")%>%
  inner_join(Q9w2,"Variable")%>%
  select(Variable,
         Pop_Mean=all_of(2),Pop_SD=all_of(3),Pop_N=all_of(4),
         Q1w_Mean=all_of(6),Q1w_SD=all_of(7),Q1w_N=all_of(8),
         Q9w1_Mean=all_of(10),Q9w1_SD=all_of(11),Q9w1_N=all_of(12),
         Q9w2_Mean=all_of(14),Q9w2_SD=all_of(15),Q9w2_N=all_of(16) 
  ) 

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... comparing weighted means \n")

Pop.BLw  <- PairwiseComparisons(Pop,Q1w) %>% mutate(Sample="BL",Ref="Pop",Weight="BL|Pop")
Pop.FUw1 <- PairwiseComparisons(Pop,Q9w1)%>% mutate(Sample="FU",Ref="Pop",Weight="FU|Pop")
Pop.FUw2 <- PairwiseComparisons(Pop,Q9w2)%>% mutate(Sample="FU",Ref="Pop",Weight="FU|BL")
BL.FUw2  <- PairwiseComparisons(BL,Q9w2)%>% mutate(Sample="FU",Ref="BL",Weight="FU|BL")

Df2 <-
  Pop.BLw%>%select(Variable,PopQ1w_Zp=Z_plog10,PopQ1w_Fp=F_plog10) %>%
  inner_join(
    Pop.FUw1%>%select(Variable,PopQ9w1_Zp=Z_plog10,PopQ9w1_Fp=F_plog10)) %>%
  inner_join(
    Pop.FUw2%>%select(Variable,PopQ9w2_Zp=Z_plog10,PopQ9w2_Fp=F_plog10)) %>%
  inner_join(
    BL.FUw2%>%select(Variable,Q1Q9w2_Zp=Z_plog10,Q1Q9w2_Fp=F_plog10)
  )


ssb <- Df1 %>% inner_join(Df2,"Variable")


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... KUHR variables \n")
kuhr     <- fread("N:/durable/projects/crayner/GenPar/Bias/KuhrVars_WeightedMeans.csv")
Pop_kuhr<- fread("N:/durable/projects/crayner/GenPar/Bias/KuhrVarsPop_UnweightedMeans.csv")
BLw_kuhr <- fread("N:/durable/projects/crayner/GenPar/Bias/KuhrVarsQ1_WeightedMeans.csv")
FUw1_kuhr<- fread("N:/durable/projects/crayner/GenPar/Bias/KuhrVarsQ9_WeightedMeans.csv")
FUw2_kuhr<- fread("N:/durable/projects/crayner/GenPar/Bias/KuhrVarsQ9m_WeightedMeans.csv")
 
Df1 <- 
  Pop_kuhr%>%
  inner_join(BLw_kuhr,"Variable")%>%
  inner_join(FUw1_kuhr,"Variable")%>%
  inner_join(FUw2_kuhr,"Variable")%>%
  select(Variable,
         Pop_Mean=all_of(2),Pop_SD=all_of(3),Pop_N=all_of(4),
         Q1w_Mean=all_of(6),Q1w_SD=all_of(7),Q1w_N=all_of(8),
         Q9w1_Mean=all_of(10),Q9w1_SD=all_of(11),Q9w1_N=all_of(12),
         Q9w2_Mean=all_of(14),Q9w2_SD=all_of(15),Q9w2_N=all_of(16) 
  ) 

Pop.BLw_kuhr <-PairwiseComparisons(Pop_kuhr,BLw_kuhr) %>% mutate(Sample="BL",Ref="Pop",Weight="BL|Pop")
Pop.FUw1_kuhr<-PairwiseComparisons(Pop_kuhr,FUw1_kuhr)%>% mutate(Sample="FU",Ref="Pop",Weight="FU|Pop")
Pop.FUw2_kuhr<-PairwiseComparisons(Pop_kuhr,FUw2_kuhr)%>% mutate(Sample="FU",Ref="Pop",Weight="FU|BL")
BL.FUw2_kuhr <-PairwiseComparisons(BL_kuhr,FUw2_kuhr)%>% mutate(Sample="FU",Ref="BL",Weight="FU|BL")

Df2 <-
  Pop.BLw_kuhr%>%select(Variable,PopQ1w_Zp=Z_plog10,PopQ1w_Fp=F_plog10) %>%
  inner_join(
    Pop.FUw1_kuhr%>%select(Variable,PopQ9w1_Zp=Z_plog10,PopQ9w1_Fp=F_plog10)) %>%
  inner_join(
    Pop.FUw2_kuhr%>%select(Variable,PopQ9w2_Zp=Z_plog10,PopQ9w2_Fp=F_plog10)) %>%
  inner_join(
    BL.FUw2_kuhr%>%select(Variable,Q1Q9w2_Zp=Z_plog10,Q1Q9w2_Fp=F_plog10)
  )

kuhr <- Df1 %>% inner_join(Df2,"Variable")

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... tabulating weighted means: T \n")

WeightedMeans_T1 <-
  ssb %>% 
  bind_rows(kuhr) %>%
  mutate(Who = ifelse(grepl(x=Variable,"^Mo"),"Mothers",ifelse(grepl(x=Variable,"^Fa"),"Fathers","Children"))) %>%
  mutate(Who = factor(Who, levels=c("Mothers","Fathers","Children"))) %>% 
  filter(Who!="Fathers") %>% 
  mutate(Variable = str_remove_all(Variable,"Mo|Fa|Ch")) %>%
  mutate(Variable = str_replace_all(Variable,"ResidenceSettlement","Hometown ")) %>%
  mutate(Variable = str_replace_all(Variable,"Residence","Hometown ")) %>%
  mutate(Variable = str_replace_all(Variable,"HouseholdIncome_IF","Household Income")) %>%
  mutate(Variable = str_replace_all(Variable,"Edu","Education")) %>%
  mutate(Variable = str_replace_all(Variable,"Math","Math test scores")) %>%
  mutate(Variable = str_replace_all(Variable,"Eng","English test scores")) %>%
  mutate(Variable = str_replace_all(Variable,"Zsc","")) %>%
  mutate(Variable = str_replace_all(Variable,"Status","Status (")) %>%
  mutate(Variable = str_replace_all(Variable,"Male","(Male")) %>%
  mutate(Variable = str_replace_all(Variable,"([a-z])([A-Z])","\\1 \\2")) %>%
  mutate(Variable = str_replace_all(Variable," IF",""))%>%
  mutate(Variable = str_replace_all(Variable,"Hyp",""))%>%
  mutate(Variable = str_replace_all(Variable,"Diagnosis",""))%>%
  mutate(Variable = str_replace_all(Variable,"_"," "))%>%
  mutate(Variable = str_replace_all(Variable,"YOB","Year of birth"))%>%
  mutate(Variable = str_replace_all(Variable,"(\\()(.*)([a-z,2]$)","\\1\\2\\3\\)")) %>%
  mutate(Variable = str_replace_all(Variable,"\\( ","\\(")) %>%
  filter(!grepl(x=Variable,"Mental|Urban|Earner")) %>%
  mutate(Order =
    ifelse(grepl(x=Variable,"Age|Sex|birth"),"A", ifelse(grepl(x=Variable,"Born"),"B",
    ifelse(grepl(x=Variable,"Parents",fixed=T),"C",ifelse(grepl(x=Variable,"Grandparents"),"D",
    ifelse(grepl(x=Variable,"Civil"),"E", ifelse(grepl(x=Variable,"Education|University|test"),"F",
    ifelse(grepl(x=Variable,"Income"),"G", ifelse(grepl(x=Variable,"Employment"),"H",
    ifelse(grepl(x=Variable,"Hours"),"I", ifelse(grepl(x=Variable,"Worked"),"J",
    ifelse(grepl(x=Variable,"Hometown"),"K","L"
    )))))))))))) %>% arrange(Who,Order) %>% select(-Order) %>% select(Who,everything())

WeightedMeans_T1 <- WeightedMeans_T1 %>% mutate(Variable=factor(Variable,levels=unique(WeightedMeans_T1$Variable)))

names(WeightedMeans_T1)<-str_replace_all(names(WeightedMeans_T1),"Q1","BL")
names(WeightedMeans_T1)<-str_replace_all(names(WeightedMeans_T1),"Q9","FU")

saveRDS(WeightedMeans_T1,"Bias/SsbVars_WeightedMeans_T1.rds")

WordTable(
  WeightedMeans_T1 %>% 
    mutate(across(matches("Mean|SD",function(x) round(x,2)))) %>% 
    mutate(across(matches("p$",function(x) format(x, digits = 2, scientific = TRUE))))
    ,
  "Bias/SsbVars_WeightedMeans_T1.docx")

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... tabulating weighted means: T2 \n")

WeightedMeans_T2 <-
  BL%>%inner_join(Q9w2,"Variable") %>% 
  select(Variable,
         Q1_Mean=all_of(2),Q1_SD=all_of(3),Q1_N=all_of(4), #Q1_Type=all_of(5),
         Q9w_Mean=all_of(6),Q9w_SD=all_of(7),Q9w_N=all_of(8) #,Q1w_Type=all_of(9)
         ) %>% 
  inner_join(
    BL.FUw2%>%select(Variable,Q1Q9w_Zp=Z_p,Q1Q9w_Fp=F_p)
  ) %>% 
  mutate(Who = ifelse(grepl(x=Variable,"^Mo"),"Mothers",ifelse(grepl(x=Variable,"^Fa"),"Fathers","Children"))) %>%
  mutate(Who = factor(Who, levels=c("Mothers","Fathers","Children"))) %>% 
  filter(Who!="Fathers") %>% 
  mutate(Variable = str_remove_all(Variable,"Mo|Fa|Ch")) %>%
  mutate(Variable = str_replace_all(Variable,"ResidenceSettlement","Hometown ")) %>%
  mutate(Variable = str_replace_all(Variable,"Residence","Hometown ")) %>%
  mutate(Variable = str_replace_all(Variable,"HouseholdIncome_IF","Household Income")) %>%
  mutate(Variable = str_replace_all(Variable,"Edu","Education")) %>%
  mutate(Variable = str_replace_all(Variable,"Math","Math test scores")) %>%
  mutate(Variable = str_replace_all(Variable,"Eng","English test scores")) %>%
  mutate(Variable = str_replace_all(Variable,"Zsc","")) %>%
  mutate(Variable = str_replace_all(Variable,"Status","Status (")) %>%
  mutate(Variable = str_replace_all(Variable,"Male","(Male")) %>%
  mutate(Variable = str_replace_all(Variable,"([a-z])([A-Z])","\\1 \\2")) %>%
  mutate(Variable = str_replace_all(Variable," IF",""))%>%
  mutate(Variable = str_replace_all(Variable,"Hyp",""))%>%
  mutate(Variable = str_replace_all(Variable,"Diagnosis",""))%>%
  mutate(Variable = str_replace_all(Variable,"_"," "))%>%
  mutate(Variable = str_replace_all(Variable,"YOB","Year of birth"))%>%
  mutate(Variable = str_replace_all(Variable,"(\\()(.*)([a-z,2]$)","\\1\\2\\3\\)")) %>%
  mutate(Variable = str_replace_all(Variable,"\\( ","\\(")) %>%
  filter(!grepl(x=Variable,"Mental|Urban|Earner")) %>%
  mutate(Order =
    ifelse(grepl(x=Variable,"Age|Sex|birth"),"A", ifelse(grepl(x=Variable,"Born"),"B",
    ifelse(grepl(x=Variable,"Parents",fixed=T),"C",ifelse(grepl(x=Variable,"Grandparents"),"D",
    ifelse(grepl(x=Variable,"Civil"),"E", ifelse(grepl(x=Variable,"Education|University|test"),"F",
    ifelse(grepl(x=Variable,"Income"),"G", ifelse(grepl(x=Variable,"Employment"),"H",
    ifelse(grepl(x=Variable,"Hours"),"I", ifelse(grepl(x=Variable,"Worked"),"J",
    ifelse(grepl(x=Variable,"Hometown"),"K","L"
    )))))))))))) %>% arrange(Who,Order) %>% select(-Order) %>% select(Who,everything())

WeightedMeans_T2 <- WeightedMeans_T2 %>% mutate(Variable=factor(Variable,levels=unique(WeightedMeans_T2$Variable)))

names(WeightedMeans_T2)<-str_replace_all(names(WeightedMeans_T2),"Q1","BL")
names(WeightedMeans_T2)<-str_replace_all(names(WeightedMeans_T2),"Q9","FU")

saveRDS(WeightedMeans_T2,"Bias/SsbVars_WeightedMeans_T2.rds")


WordTable(
  WeightedMeans_T2 %>% 
    mutate(across(matches("Mean|SD",function(x) round(x,2)))) %>% 
    mutate(across(matches("p$",function(x) format(x, digits = 2, scientific = TRUE))))
    ,
  "Bias/SsbVars_WeightedMeans_T2.docx")



### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... tabulating weighted means: FOR PLOTS \n")

WeightedMeans_P <- 
  bind_rows(list(
    Pop.BLw,Pop.FUw1,Pop.FUw2,
    Pop.BLw_kuhr,Pop.FUw1_kuhr,Pop.FUw2_kuhr,
    BL.FUw2 ,BL.FUw2_kuhr
    )) %>% 
  select(Variable,Cohens,Cohens_2.5,Cohens_97.5,Z_stat,P=Z_p,Sample,Ref,Weight)  %>% 
  mutate(Who = ifelse(grepl(x=Variable,"^Mo"),"Mothers",ifelse(grepl(x=Variable,"^Fa"),"Fathers","Children"))) %>%
  mutate(Who = factor(Who, levels=c("Mothers","Fathers","Children"))) %>% 
  filter(Who!="Fathers") %>% 
  mutate(Variable = str_remove_all(Variable,"Mo|Fa|Ch")) %>%
  mutate(Variable = str_replace_all(Variable,"ResidenceSettlement","Hometown ")) %>%
  mutate(Variable = str_replace_all(Variable,"Residence","Hometown ")) %>%
  mutate(Variable = str_replace_all(Variable,"HouseholdIncome_IF","Household Income")) %>%
  mutate(Variable = str_replace_all(Variable,"Edu","Education")) %>%
  mutate(Variable = str_replace_all(Variable,"Math","Math test scores")) %>%
  mutate(Variable = str_replace_all(Variable,"Eng","English test scores")) %>%
  mutate(Variable = str_replace_all(Variable,"Zsc","")) %>%
  mutate(Variable = str_replace_all(Variable,"Status","Status (")) %>%
  mutate(Variable = str_replace_all(Variable,"Male","(Male")) %>%
  mutate(Variable = str_replace_all(Variable,"([a-z])([A-Z])","\\1 \\2")) %>%
  mutate(Variable = str_replace_all(Variable," IF",""))%>%
  mutate(Variable = str_replace_all(Variable,"Hyp",""))%>%
  mutate(Variable = str_replace_all(Variable,"Diagnosis",""))%>%
  mutate(Variable = str_replace_all(Variable,"_"," "))%>%
  mutate(Variable = str_replace_all(Variable,"YOB","Year of birth"))%>%
  mutate(Variable = str_replace_all(Variable,"(\\()(.*)([a-z,2]$)","\\1\\2\\3\\)")) %>%
  mutate(Variable = str_replace_all(Variable,"\\( ","\\(")) %>%
  filter(!grepl(x=Variable,"Mental|Urban|Earner")) %>%
  mutate(Order =
    ifelse(grepl(x=Variable,"Age|Sex|birth"),"A", ifelse(grepl(x=Variable,"Born"),"B",
    ifelse(grepl(x=Variable,"Parents",fixed=T),"C",ifelse(grepl(x=Variable,"Grandparents"),"D",
    ifelse(grepl(x=Variable,"Civil"),"E", ifelse(grepl(x=Variable,"Education|University|test"),"F",
    ifelse(grepl(x=Variable,"Income"),"G", ifelse(grepl(x=Variable,"Employment"),"H",
    ifelse(grepl(x=Variable,"Hours"),"I", ifelse(grepl(x=Variable,"Worked"),"J",
    ifelse(grepl(x=Variable,"Hometown"),"K","L"
    )))))))))))) %>% arrange(Order) %>% select(-Order) %>% select(Who,everything()) %>%
    mutate(`P-value`=ifelse(P<as.numeric(BONF),paste0("P<",BONF),paste0("P>",BONF))) %>% 
  mutate(`P-value`=factor(`P-value`,levels=c(paste0("P<",BONF),paste0("P>",BONF))))

WeightedMeans_P <- WeightedMeans_P %>% mutate(Variable=factor(Variable,levels=unique(WeightedMeans_P$Variable)))

names(WeightedMeans_P)<-str_replace_all(names(WeightedMeans_P),"Q1","BL")
names(WeightedMeans_P)<-str_replace_all(names(WeightedMeans_P),"Q9","FU")

saveRDS(WeightedMeans_P,"Bias/SsbVars_WeightedMeansForPlot1.rds")

remove(list = ls())

