#!/usr/bin/Rscript
cat("\n *** Script for:","\n *** \t - quantifying bias in univariable regression analysis")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Preparing workspace \n")

remove(list = ls())

source("profile.R")
source("Pheno_Functions.R")

pkgs <- c("data.table","dplyr","stringr","tidyr","ggplot2","survey", "flextable","officer")

for(p in pkgs){  LoadOrInstall(p) }

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Load data

X<-readRDS("./Data/NarrowEligibility_BiasVariables.rds")

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Calculate Bonferroni

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
## outcome exposure variable pairs

YVARS<-c("MoHouseholdIncome_IF","MoHouseholdIncome_IF","MoHouseholdIncome_IF","MoHouseholdIncome_IF","MoHouseholdIncome_IF","MoHouseholdIncome_IF",
         "ChMeanNationalTestScore","ChMeanNationalTestScore","ChMeanNationalTestScore","ChMeanNationalTestScore","ChMeanNationalTestScore"
         )

XVARS<-c("MoAge","MoEduYears","MoHoursWorked","MoBornInNorway","MoCivilStatus_Married","MoResidence_Urban",
         "MoHouseholdIncome_IF","MoEduYears","MoCivilStatus_Married","ChSex_Male","ChNumberNorwegianParents"
)



df_Pop<- X %>% select(LNR,any_of(unique(YVARS)),any_of(unique(XVARS)))%>%as.data.frame()%>%
  mutate(across(!LNR&where(~n_distinct(.)>2),~winsorize(.x)))%>%
  mutate(across(!LNR&where(~n_distinct(.)>2),~(.x-mean(.x,na.rm=T))/sd(.x,na.rm=T)))

df_Q1 <- X %>% filter(MoQ1Participation==1)%>% select(LNR,any_of(unique(YVARS)),any_of(unique(XVARS)))%>%as.data.frame()%>%
  mutate(across(!LNR&where(~n_distinct(.)>2),~winsorize(.x)))%>%
  mutate(across(!LNR&where(~n_distinct(.)>2),~(.x-mean(.x,na.rm=T))/sd(.x,na.rm=T)))

df_Q9 <- X %>% filter(MoQ9Participation==1)%>% select(LNR,any_of(unique(YVARS)),any_of(unique(XVARS)))%>%as.data.frame()%>%
  mutate(across(!LNR&where(~n_distinct(.)>2),~winsorize(.x)))%>%
  mutate(across(!LNR&where(~n_distinct(.)>2),~(.x-mean(.x,na.rm=T))/sd(.x,na.rm=T)))

bv_Pop<-lapply(seq_along(YVARS),function(i){BvUnweightedModel(dat=df_Pop,x=XVARS[i],y=YVARS[i])})
bv_Q1<-lapply(seq_along(YVARS),function(i){BvUnweightedModel(dat=df_Q1,x=XVARS[i],y=YVARS[i])})
bv_Q9<-lapply(seq_along(YVARS),function(i){BvUnweightedModel(dat=df_Q9,x=XVARS[i],y=YVARS[i])})

tab_Pop.Q1<-bind_rows(lapply(seq_along(YVARS),function(i){GetResultsTab(X=XVARS[i],Y=YVARS[i],Model1=bv_Pop[[i]],Model2=bv_Q1[[i]])}))
tab_Pop.Q1$Ref<-"Pop";tab_Pop.Q1$Sample<-"Q1";tab_Pop.Q1$Weight<-"Unweighted"

tab_Pop.Q9<-bind_rows(lapply(seq_along(YVARS),function(i){GetResultsTab(X=XVARS[i],Y=YVARS[i],Model1=bv_Pop[[i]],Model2=bv_Q9[[i]])}))
tab_Pop.Q9$Ref<-"Pop";tab_Pop.Q9$Sample<-"Q9";tab_Pop.Q9$Weight<-"Unweighted"

tab_Q1.Q9<-bind_rows(lapply(seq_along(YVARS),function(i){GetResultsTab(X=XVARS[i],Y=YVARS[i],Model1=bv_Q1[[i]],Model2=bv_Q9[[i]])}))
tab_Q1.Q9$Ref<-"Q1";tab_Q1.Q9$Sample<-"Q9";tab_Q1.Q9$Weight<-"Unweighted"


## Weighted assoc

df_PopQ1<-df_Q1%>%inner_join(PP_PopQ1,"LNR")%>%select(-LNR,-PP)
df_PopQ9<-df_Q9%>%inner_join(PP_PopQ9,"LNR")%>%select(-LNR,-PP)
df_Q1Q9<-df_Q9%>%inner_join(PP_Q1Q9,"LNR")%>%select(-LNR,-PP)

library(survey)
bv_Q1w<-lapply(seq_along(YVARS),function(i){BvWeightedModel(dat=df_PopQ1,x=XVARS[i],y=YVARS[i],w="IPW")})
bv_Q9w1<-lapply(seq_along(YVARS),function(i){ BvWeightedModel(dat=df_PopQ9,x=XVARS[i],y=YVARS[i],w="IPW")})
bv_Q9w2<-lapply(seq_along(YVARS),function(i){ BvWeightedModel(dat=df_Q1Q9,x=XVARS[i],y=YVARS[i],w="IPW")})

tab_Pop.Q1w<-bind_rows(lapply(seq_along(YVARS),function(i){GetResultsTab(X=XVARS[i],Y=YVARS[i],Model1=bv_Pop[[i]],Model2=bv_Q1w[[i]])}))
tab_Pop.Q1w$Ref<-"Pop";tab_Pop.Q1w$Sample<-"Q1";tab_Pop.Q1w$Weight<-"Q1"

tab_Pop.Q9w1<-bind_rows(lapply(seq_along(YVARS),function(i){GetResultsTab(X=XVARS[i],Y=YVARS[i],Model1=bv_Pop[[i]],Model2=bv_Q9w1[[i]])}))
tab_Pop.Q9w1$Ref<-"Pop";tab_Pop.Q9w1$Sample<-"Q9";tab_Pop.Q9w1$Weight<-"Q9"

tab_Pop.Q9w2<-bind_rows(lapply(seq_along(YVARS),function(i){GetResultsTab(X=XVARS[i],Y=YVARS[i],Model1=bv_Pop[[i]],Model2=bv_Q9w2[[i]])}))
tab_Pop.Q9w2$Ref<-"Pop";tab_Pop.Q9w2$Sample<-"Q9";tab_Pop.Q9w2$Weight<-"Q9|Q1"

tab_Q1.Q9w2<-bind_rows(lapply(seq_along(YVARS),function(i){GetResultsTab(X=XVARS[i],Y=YVARS[i],Model1=bv_Q1[[i]],Model2=bv_Q9w2[[i]])}))
tab_Q1.Q9w2$Ref<-"Q1";tab_Q1.Q9w2$Sample<-"Q9";tab_Q1.Q9w2$Weight<-"Q9|Q1"

Coef_Dem <- 
  bind_rows(
    tab_Pop.Q1,tab_Pop.Q9,tab_Q1.Q9,
    tab_Pop.Q1w,tab_Pop.Q9w1,tab_Pop.Q9w2,tab_Q1.Q9w2
    )
  
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... LOAD KUHR ANALYSIS RESULTS \n")

Coef_kuhr <- 
  fread("N:/durable/projects/crayner/KuhrDiagnoses/Bias/KuhrVars_AllCoefficients.csv") %>% 
  filter(grepl(x=Y,"MoHousehold")|grepl(x=Y,"ChMath")) %>% 
  filter(!grepl(x=X,"Burden")) %>% 
  filter(!grepl(x=X,"^Fa")) %>% 
  mutate(Weight=str_remove_all(Weight,"IPP ")) %>% filter(Weight!="Q9x") %>% 
  mutate(Sample=ifelse(Ref=="Q1","Q9",Sample)) %>% 
  mutate(Weight=ifelse(Ref=="Q1"&Sample=="Q9"&Weight!="Unweighted","Q9|Q1",Weight))%>% 
  distinct(X,Y,Ref,Sample,Weight,.keep_all=T)

Coef_kuhr <-
  Coef_kuhr %>% 
  bind_rows(
    Coef_kuhr %>% filter(Ref=="Pop") %>% select(X,Y,Coeff1,SE1,Ref) %>% distinct() %>% 
      inner_join(
        Coef_kuhr %>% filter(Sample=="Q9"&Weight=="Unweighted") %>% select(X,Y,Coeff2,SE2,Sample,Weight) %>% distinct(),
        by=c("X","Y")) %>% 
      mutate(Zdiff=(Coeff1-Coeff2)/sqrt(SE1^2+SE2^2)) %>% 
      mutate(Zdiff_P=2*pnorm(abs(Zdiff),lower.tail=F))
  ) %>% 
  bind_rows(
    Coef_kuhr %>% filter(Ref=="Pop") %>% select(X,Y,Coeff1,SE1,Ref) %>% distinct() %>% 
      inner_join(
        Coef_kuhr %>% filter(Sample=="Q9"&Weight=="Q9|Q1") %>% select(X,Y,Coeff2,SE2,Sample,Weight) %>% distinct(),
        by=c("X","Y"))%>%
      mutate(Zdiff=(Coeff1-Coeff2)/sqrt(SE1^2+SE2^2)) %>% 
      mutate(Zdiff_P=2*pnorm(abs(Zdiff),lower.tail=F))
  )

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... BIND ROWS AND WRITE TO FILE \n")

COEF_TAB <-
  bind_rows(Coef_Dem,Coef_kuhr)%>% 
  mutate(Who.Y=ifelse(grepl(x=Y,"^Mo"),"Mothers",ifelse(grepl(x=Y,"^Fa"),"Fathers","Children"))) %>%
  mutate(Who.Y=factor(Who.Y,levels=c("Mothers","Fathers","Children"))) %>% 
  mutate(Who.X=ifelse(grepl(x=X,"^Mo"),"Mothers",ifelse(grepl(x=X,"^Fa"),"Fathers","Children"))) %>%
  mutate(Who.X=factor(Who.X,levels=c("Mothers","Fathers","Children"))) %>% 
  filter(Who.X!="Fathers") %>% filter(Who.Y!="Fathers") %>% 
  mutate(across(c("X","Y"),~str_remove_all(.x,"Mo|Fa|Ch"))) %>%
  mutate(across(c("X","Y"),~str_replace_all(.x,"ResidenceSettlement","Hometown "))) %>%
  mutate(across(c("X","Y"),~str_replace_all(.x,"Residence","Hometown "))) %>%
  mutate(across(c("X","Y"),~str_replace_all(.x,"HouseholdIncome_IF","Household Income")) )%>%
  mutate(across(c("X","Y"),~str_replace_all(.x,"Edu","Education"))) %>%
  mutate(across(c("X","Y"),~str_replace_all(.x,"Math","National Test Score"))) %>%
  mutate(across(c("X","Y"),~str_replace_all(.x,"Mean",""))) %>%
  mutate(across(c("X","Y"),~str_replace_all(.x,"Eng","English test scores"))) %>%
  mutate(across(c("X","Y"),~str_replace_all(.x,"Zsc",""))) %>%
  mutate(across(c("X","Y"),~str_replace_all(.x,"Hyp","")))%>%
  mutate(across(c("X","Y"),~str_replace_all(.x,"Status","Status ("))) %>%
  mutate(across(c("X","Y"),~str_replace_all(.x,"Male","(Male"))) %>%
  mutate(across(c("X","Y"),~str_replace_all(.x,"([a-z])([A-Z])","\\1 \\2"))) %>%
  mutate(across(c("X","Y"),~str_replace_all(.x," IF","")))%>%
  mutate(across(c("X","Y"),~str_replace_all(.x,"_"," ")))%>%
  mutate(across(c("X","Y"),~str_replace_all(.x,"YOB","Year of birth")))%>%
  mutate(across(c("X","Y"),~str_replace_all(.x,"(\\()(.*)([a-z,2]$)","\\1\\2\\3\\)")))%>%
  mutate(across(c("X","Y"),~str_replace_all(.x,"\\( ","\\("))) %>% 
  mutate(across(c("X","Y"),~str_replace_all(.x,"Hyp","")))%>%
  mutate(across(c("X","Y"),~str_replace_all(.x,"Diagnosis","")))%>%
  mutate(across(c("X","Y"),~str_replace_all(.x,"^ ","")))%>%
  filter(!grepl(x=X,"Mental|Earner")) %>%
  filter(!grepl(x=Y,"Mental|Earner")) %>%
  filter(Y!="English test scores") %>% 
  mutate(X=ifelse(Who.Y=="Children" & Who.X=="Mothers", paste0(Who.X," ",X),
    ifelse(Who.Y=="Mothers" & Who.X=="Children", paste0(Who.X,"s ",X), paste0(X)))) %>% 
  mutate(X=factor(X,
    levels=c(
      "Sex (Male)","Age","Born In Norway","Number Norwegian Parents","Number Norwegian Grandparents",
      "Civil Status (Married)","Mothers Civil Status (Married)",
      "Education Years","Mothers Education Years",
      "Household Income","Mothers Household Income",
      "Hours Worked",
      "Hometown  Urban","Mothers Hometown  Urban",
      "Mothers ADHD",
      "Mothers Internalising",
      "ADHD",
      "Internalising",
      "Childrens ADHD",
      "Childrens Internalising")
    )) %>%
  distinct() %>% 
  mutate(`P-value`=ifelse(Zdiff_P<as.numeric(BONF),paste0("<",BONF),paste0(">",BONF))) %>% 
  mutate(`P-value`=factor(`P-value`,levels=c(paste0("<",BONF),paste0(">",BONF))))%>%
  mutate(Magnitude = ifelse(
    abs(round(Coeff1,3))>abs(round(Coeff2,3)), "Estimate is smaller in magnitude",
        ifelse(abs(round(Coeff1,3))<abs(round(Coeff2,3)), "Estimate is larger in magnitude","Same at 3 d.p"))) %>% 
  mutate(Direction = ifelse( 
             (round(Coeff1,3)<0&round(Coeff2,3)>0)|(round(Coeff1,3)>0&round(Coeff2,3)<0), 
         "Estimates have different direction of effect","Estimates are in same direction")) %>% 
  filter((Y=="Household Income"|Y=="National Test Score")) %>% 
  mutate(Analysis=ifelse(Weight=="Unweighted","Unweighted","Weighted")) %>% 
  mutate(Who.Y=paste0(Who.Y,"\n",Y)) %>% 
  mutate(Who.Y=factor(Who.Y,levels=sort(unique(Who.Y),decreasing=T))) 

fwrite(COEF_TAB,"Bias/FullWeightedandUnweightedCoeffTab.csv")











### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... KUHR diagnosese \n")


YVARS<-c("MoHouseholdIncome_IF","MoHouseholdIncome_IF","MoHouseholdIncome_IF",
         "MoHouseholdIncome_IF","MoHouseholdIncome_IF","MoHouseholdIncome_IF",
         "FaHouseholdIncome_IF","FaHouseholdIncome_IF","FaHouseholdIncome_IF",
         
         "ChMeanNationalTestScore","ChMeanNationalTestScore","ChMeanNationalTestScore",
         "ChMeanNationalTestScore","ChMeanNationalTestScore","ChMeanNationalTestScore",
         "ChMeanNationalTestScore","ChMeanNationalTestScore","ChMeanNationalTestScore",
         
        "ChDiagnosis_HypADHD","ChMentalHealthBurden","ChDiagnosis_Internalising",
        "MoDiagnosis_HypADHD","MoDiagnosis_Internalising","MoMentalHealthBurden",
        "FaDiagnosis_HypADHD","FaDiagnosis_Internalising","FaMentalHealthBurden",

        "ChDiagnosis_HypADHD","ChMentalHealthBurden","ChDiagnosis_Internalising",
        "MoDiagnosis_HypADHD","MoDiagnosis_Internalising","MoMentalHealthBurden",
        "FaDiagnosis_HypADHD","FaDiagnosis_Internalising","FaMentalHealthBurden",
        
        "ChDiagnosis_HypADHD","ChMentalHealthBurden","ChDiagnosis_Internalising",
        "MoDiagnosis_HypADHD","MoDiagnosis_Internalising","MoMentalHealthBurden",
        "FaDiagnosis_HypADHD","FaDiagnosis_Internalising","FaMentalHealthBurden"    
         )

XVARS<-c(
        "ChDiagnosis_HypADHD","ChMentalHealthBurden","ChDiagnosis_Internalising",
        "MoDiagnosis_HypADHD","MoDiagnosis_Internalising","MoMentalHealthBurden",
        "FaDiagnosis_HypADHD","FaDiagnosis_Internalising","FaMentalHealthBurden",
  
        "ChDiagnosis_HypADHD","ChMentalHealthBurden","ChDiagnosis_Internalising",
        "MoDiagnosis_HypADHD","MoDiagnosis_Internalising","MoMentalHealthBurden",
        "FaDiagnosis_HypADHD","FaDiagnosis_Internalising","FaMentalHealthBurden",
        
        "MoHouseholdIncome_IF","MoHouseholdIncome_IF","MoHouseholdIncome_IF",
         "MoHouseholdIncome_IF","MoHouseholdIncome_IF","MoHouseholdIncome_IF",
         "FaHouseholdIncome_IF","FaHouseholdIncome_IF","FaHouseholdIncome_IF",
        
        "ChNumberNorwegianGrandparents","ChNumberNorwegianGrandparents","ChNumberNorwegianGrandparents",
        "MoBornInNorway","MoBornInNorway","MoBornInNorway",
        "FaBornInNorway","FaBornInNorway","FaBornInNorway",
        
        "MoResidence_Urban","MoResidence_Urban","MoResidence_Urban",
        "MoResidence_Urban","MoResidence_Urban","MoResidence_Urban",
        "FaResidence_Urban","FaResidence_Urban","FaResidence_Urban"
)


df_Pop<- X %>% select(LNR,any_of(unique(YVARS)),any_of(unique(XVARS)))%>%as.data.frame()%>%
  mutate(across(!LNR&where(~n_distinct(.)>2),~winsorize(.x)))%>%
  mutate(across(!LNR&where(~n_distinct(.)>2),~(.x-mean(.x,na.rm=T))/sd(.x,na.rm=T)))

df_Q1 <- X %>% filter(MoQ1Participation==1)%>% select(LNR,any_of(unique(YVARS)),any_of(unique(XVARS)))%>%as.data.frame()%>%
  mutate(across(!LNR&where(~n_distinct(.)>2),~winsorize(.x)))%>%
  mutate(across(!LNR&where(~n_distinct(.)>2),~(.x-mean(.x,na.rm=T))/sd(.x,na.rm=T)))

df_Q9 <- X %>% filter(MoQ9Participation==1)%>% select(LNR,any_of(unique(YVARS)),any_of(unique(XVARS)))%>%as.data.frame()%>%
  mutate(across(!LNR&where(~n_distinct(.)>2),~winsorize(.x)))%>%
  mutate(across(!LNR&where(~n_distinct(.)>2),~(.x-mean(.x,na.rm=T))/sd(.x,na.rm=T)))

bv_Pop<-lapply(seq_along(YVARS),function(i){BvUnweightedModel(dat=df_Pop,x=XVARS[i],y=YVARS[i])})
bv_Q1<-lapply(seq_along(YVARS),function(i){BvUnweightedModel(dat=df_Q1,x=XVARS[i],y=YVARS[i])})
bv_Q9<-lapply(seq_along(YVARS),function(i){BvUnweightedModel(dat=df_Q9,x=XVARS[i],y=YVARS[i])})

tab_Pop.Q1<-bind_rows(lapply(seq_along(YVARS),function(i){GetResultsTab(X=XVARS[i],Y=YVARS[i],Model1=bv_Pop[[i]],Model2=bv_Q1[[i]])}))
tab_Pop.Q1$Ref<-"Pop";tab_Pop.Q1$Sample<-"Q1";tab_Pop.Q1$Weight<-"Unweighted"

tab_Pop.Q9<-bind_rows(lapply(seq_along(YVARS),function(i){GetResultsTab(X=XVARS[i],Y=YVARS[i],Model1=bv_Pop[[i]],Model2=bv_Q9[[i]])}))
tab_Pop.Q9$Ref<-"Pop";tab_Pop.Q9$Sample<-"Q9";tab_Pop.Q9$Weight<-"Unweighted"

tab_Q1.Q9<-bind_rows(lapply(seq_along(YVARS),function(i){GetResultsTab(X=XVARS[i],Y=YVARS[i],Model1=bv_Q1[[i]],Model2=bv_Q9[[i]])}))
tab_Q1.Q9$Ref<-"Q1";tab_Q1.Q9$Sample<-"Q9";tab_Q1.Q9$Weight<-"Unweighted"

## Weighted kuhr assocs

df_Pop.Q1<-df_Q1%>%inner_join(PP_PopQ1,"LNR")%>%select(-LNR,-PP)
df_Pop.Q9<-df_Q9%>%inner_join(PP_PopQ9,"LNR")%>%select(-LNR,-PP)
df_Pop.Q9x<-df_Q9%>%inner_join(PP_PopQ9x,"LNR")%>%select(-LNR,-PP)
df_Q1.Q9<-df_Q9%>%inner_join(PP_Q1Q9,"LNR")%>%select(-LNR,-PP)

library(survey)
bv_Pop.Q1<-lapply(seq_along(YVARS),function(i){BvWeightedModel(dat=df_Pop.Q1,x=XVARS[i],y=YVARS[i],w="IPW")})
bv_Pop.Q9<-lapply(seq_along(YVARS),function(i){ BvWeightedModel(dat=df_Pop.Q9,x=XVARS[i],y=YVARS[i],w="IPW")})
bv_Pop.Q9x<-lapply(seq_along(YVARS),function(i){ BvWeightedModel(dat=df_Pop.Q9x,x=XVARS[i],y=YVARS[i],w="IPW")})
bv_Q1.Q9<-lapply(seq_along(YVARS),function(i){ BvWeightedModel(dat=df_Q1.Q9,x=XVARS[i],y=YVARS[i],w="IPW")})

tab_Pop.Q1w<-bind_rows(lapply(seq_along(YVARS),function(i){GetResultsTab(X=XVARS[i],Y=YVARS[i],Model1=bv_Pop[[i]],Model2=bv_Pop.Q1[[i]])}))
tab_Pop.Q1w$Ref<-"Pop";tab_Pop.Q1w$Sample<-"Q1";tab_Pop.Q1w$Weight<-"IPP Q1"

tab_Pop.Q9w<-bind_rows(lapply(seq_along(YVARS),function(i){GetResultsTab(X=XVARS[i],Y=YVARS[i],Model1=bv_Pop[[i]],Model2=bv_Pop.Q9[[i]])}))
tab_Pop.Q9w$Ref<-"Pop";tab_Pop.Q9w$Sample<-"Q9";tab_Pop.Q9w$Weight<-"IPP Q9"

tab_Pop.Q9wx<-bind_rows(lapply(seq_along(YVARS),function(i){GetResultsTab(X=XVARS[i],Y=YVARS[i],Model1=bv_Pop[[i]],Model2=bv_Pop.Q9x[[i]])}))
tab_Pop.Q9wx$Ref<-"Pop";tab_Pop.Q9wx$Sample<-"Q9";tab_Pop.Q9wx$Weight<-"IPP Q9x"

tab_Q1.Q9w<-bind_rows(lapply(seq_along(YVARS),function(i){GetResultsTab(X=XVARS[i],Y=YVARS[i],Model1=bv_Q1[[i]],Model2=bv_Q1.Q9[[i]])}))
tab_Q1.Q9w$Ref<-"Q1";tab_Q1.Q9w$Sample<-"Q9";tab_Q1.Q9w$Weight<-"IPP Q9"


tab <- 
  bind_rows(
    tab_Pop.Q1,tab_Pop.Q9,tab_Q1.Q9,
    tab_Pop.Q1w,tab_Pop.Q9w,tab_Pop.Q9wx,tab_Q1.Q9w
    ) 

fwrite(tab,"N:/durable/projects/crayner/KuhrDiagnoses/Bias/KuhrVars_AllCoefficients.csv")
