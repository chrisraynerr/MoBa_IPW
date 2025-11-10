#!/usr/bin/Rscript

cat("\n *** Script for:","\n *** \t - predicting participation in MoBa via Elastic-net regression")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n### Preparing workspace \n")

remove(list = ls())

source("profile.R")
pkgs <- c("data.table","dplyr","stringr","tidyr","bigstatsr","ggplot2")

for(p in pkgs){  LoadOrInstall(p) }

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
cat("\n ... getting options \n")

opt_list = list(
  optparse::make_option("--DfY",action="store",default="InitialParticipation_NarrowEligibility_TrioSSBVars.rds",type="character",help="Outcome dataset"),
  optparse::make_option("--DfX",action="store",default="InitialParticipation_NarrowEligibility_TrioSSBVars_PreImputationRF.rds",type="character",help="Outcome dataset"),
  optparse::make_option("--Y",action="store",default=F,type="character",help="Outcome label"),
  optparse::make_option("--S",action="store",default=F,type="character",help="If filtering the sample"),
  optparse::make_option("--outdir",action="store",default=F,type="character",help="Output directory if previously computed FBMs")
)

opt = optparse::parse_args(optparse::OptionParser(option_list=opt_list))

DATA_Y   <- opt$DfY 
DATA_X   <- opt$DfX
OUTCOME  <- opt$Y
FILTER   <- opt$S
OUTDIR   <- opt$outdir
WHO      <- str_sub(OUTCOME,1,2)
WHO      <- ifelse(WHO=="Ch","",WHO)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
cat("\n ... create results directory \n")

if(OUTDIR==F){
  if(FILTER==F){
  tmp <- paste0("Predict/",OUTCOME,"_",format(Sys.time(), '%Y%m'),'/')
  dir.create(tmp,recursive=T)
} else {
  tmp <- paste0("Predict/",OUTCOME,"In",FILTER,"_",format(Sys.time(), '%Y%m'),'/')
  dir.create(tmp,recursive=T)
}} else {
  tmp <- paste0("Predict/",OUTDIR)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
cat("\n ... load data \n")

if(FILTER==F){
  DfY<-readRDS(paste0("Data/",DATA_Y))%>%select(LNR,ID=all_of(paste0(WHO,"LNR")),all_of(OUTCOME)) %>% 
    mutate(!!OUTCOME:=as.integer(!!sym(OUTCOME))) #%>% 
    # mutate(!!OUTCOME:=ifelse(max(!!sym(OUTCOME))==2,!!sym(OUTCOME)-1,!!sym(OUTCOME)))
  if(max(DfY[[paste0(OUTCOME)]],na.rm=T)==2){    DfY<-DfY %>% mutate(!!OUTCOME:=!!sym(OUTCOME)-1) }

  } else{
  
  DfY<-readRDS(paste0("Data/",DATA_Y))%>%select(LNR,ID=all_of(paste0(WHO,"LNR")),all_of(OUTCOME),all_of(FILTER))%>% 
    mutate(!!OUTCOME:=as.integer(!!sym(OUTCOME)))%>% 
    mutate(!!FILTER:=as.integer(!!sym(FILTER)))
  
  if(max(DfY[[paste0(OUTCOME)]],na.rm=T)==2){    DfY<-DfY %>% mutate(!!OUTCOME:=!!sym(OUTCOME)-1) }
  if(max(DfY[[paste0(FILTER)]],na.rm=T)==2){    DfY<-DfY %>% mutate(!!FILTER:=!!sym(FILTER)-1) }
}

table(DfY[[paste0(OUTCOME)]])

DfX  <- readRDS(paste0("Data/",DATA_X))%>%select(-matches("Participation",ignore.case=F))

if("ChDOB" %in% names(DfX)){
  DfX$YEAR<-DfX$ChDOB
  }else if("Year_GiveTheDateYouFilledInTheQuestionnaire" %in% names(DfX)){
    DfX$YEAR<-DfX$Year_GiveTheDateYouFilledInTheQuestionnaire
    }

Df<-
  DfY%>%bind_cols(DfX)%>%na.omit()%>%arrange(ID,desc(!!sym(OUTCOME)),YEAR)%>%distinct(ID,.keep_all=T) %>% select(-YEAR)

if(FILTER!=F){Df<-Df%>%filter(get(FILTER)==1)}

ID   <- Df%>%select(any_of(c("ID","LNR")))

Df   <- Df%>%select(-ID,-matches("LNR",ignore.case=F))%>% 
  dplyr::mutate(across(where(~n_distinct(.)<10)&!matches("Particpation"),~as.character(.x)))%>% 
  # dplyr::mutate(across(matches("Municipality"),~as.character(.x)))%>% 
  dplyr::mutate(across(where(is.character),function(x) factor(x,levels=names(sort(table(x),decreasing=T)))))

table(Df[[paste0(OUTCOME)]])

if(length(unique(Df[[paste0(OUTCOME)]]))>2){
  REGRESSION<-"LINEAR"
} else {
  REGRESSION<-"LOGISTIC"
}

table(Df[[paste0(OUTCOME)]])

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
cat("\n ... separate numerical data from factors (to dummy code) \n")

if(!file.exists(paste0(tmp,OUTCOME,"_PredMat.rds"))){
  
  cat("\n ... Categorical variables ...")
  
  DfF <- Df  %>% select(all_of(OUTCOME),where(is.factor)) %>% 
    mutate(across(where(is.factor),~as.character(.x))) %>% 
    mutate(across(where(is.factor),~str_to_title(.x,))) %>% 
    mutate(across(where(is.factor),~str_remove_all(.x," "))) %>% 
    mutate(across(where(is.factor),~str_replace_all(.x,"[^[:alnum:]]", "_"))) %>% 
    mutate(across(where(is.factor),~str_sub(.x,1,20))) %>% 
    mutate(across(where(is.factor),~ifelse(is.na(.x),"Missing",.x))) %>% 
    select(where(~n_distinct(.)>1))
  
  Fct <- paste0(names(DfF)[-1],collapse = " + ")
  Fct <- model.matrix.lm(as.formula(paste0(OUTCOME,"~",Fct)),data=DfF,na.action="na.pass")
  Fct <- Fct[,-1]
  
  # Identify binary variables
  BIN <- colnames(Fct)[apply(Fct, 2, function(x) length(unique(x)) == 2)]
  # Calculate minor group sizes and percentages
  MIN <- sapply(BIN, function(var) min(table(Fct[, var])))
  PER <- MIN/nrow(Fct)*100
  # Filter binary variables based on minor group percentage threshold
  DROP<- BIN[PER<2]
  # Drop the selected variables from the model matrix
  Fct <- Fct[,!colnames(Fct) %in% DROP]
  
  MISSING <- apply(Fct,2, FUN = function(x){sum(is.na(x))})
  max(MISSING)
  
  Fct <- apply(Fct, 2, function(x) {  
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x
  })
  
  cat("\n ... Continuous variables ...")
  
  winsorize <- function(x,lowerP=0.01,upperP=0.99) {
    lowerB<-quantile(x,lowerP,na.rm=T);upperB<-quantile(x,upperP,na.rm=T)
    x[x<lowerB]<-lowerB; x[x>upperB]<-upperB; return(x)
  }
  
  Num <- Df  %>% 
    select(-all_of(OUTCOME)) %>% select(!matches("Participation")) %>% 
    select(!where(is.factor)) %>% 
    mutate(across(everything(),as.numeric)) %>% 
    mutate(across(everything(), ~winsorize(.x))) %>% 
    mutate(across(everything(), ~(.x-mean(.x,na.rm=T))/sd(.x,na.rm=T))) %>% 
    mutate(across(everything(), ~replace(., is.na(.), mean(., na.rm=T)))) %>% 
    dplyr::select(where(~mean(is.na(.))<1)) %>% 
    as.matrix()
  
  MISSING <- apply(Num,2, FUN = function(x){sum(is.na(x))})
  max(MISSING)
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n ... join predictor matrix \n")

  Mat <- cbind(Num,Fct)
  Var <- names(as.data.frame(Mat))
  
  cat("\n there are",nrow(Mat)," individuals\n")
  cat("\n there are",length(Var)," predictors\n")
  cat("\n predictors include: \n", Var)
  
  Var <- t(as.matrix(Var))
  fwrite(Var,paste0(tmp,OUTCOME,"_PredVarNames.csv"),row.names=F,col.names=F) 
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n ... creating FBM \n")

  if(file.exists(paste0(tmp,OUTCOME,"_PredMat.rds"))){
    file.remove(paste0(tmp,OUTCOME,"_PredMat.rds"));file.remove(paste0(tmp,OUTCOME,"_PredMat.bk"))
  }
  PredMat <- as_FBM(Mat,type="float",is_read_only=F,backingfile=paste0(tmp,OUTCOME,"_PredMat"))
  PredMat$save()

}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
cat("\n ... runing the model in the full dataset \n")

if(!file.exists(paste0(tmp,OUTCOME,"_FullModel.rds"))){
  
  NCORES1  <- as.numeric(Sys.getenv("OMP_NUM_THREADS"))
  NCORES2  <- nb_cores()
  NCORES   <- max(NCORES1,NCORES2,na.rm=T)
  
  cat("\n*** USING", NCORES, "CORES ***\n")
  
  X   <- big_attach(paste0(tmp,OUTCOME,"_PredMat.rds"))
  N1  <- nrow(X)
  # N2  <- nrow(DfF)
  # N1==N2
  
  if(REGRESSION=="LINEAR"){
    Y<-Df%>%select(Y=all_of(OUTCOME))%>%mutate(Y=as.numeric(Y))%>%as.matrix() # THIS CANNOT BE OF CLASS DATA.FRAME # Error in xtfrm.data.frame(x) : cannot xtfrm data frames
    table(Y);str(Y)
    M<-big_spLinReg(X,Y,alphas=c(10^(-(0:4))),K=10,warn=F,ncores=NCORES)
  } else {
    Y<-Df%>%select(Y=all_of(OUTCOME))%>%mutate(Y=as.integer(Y)-1)%>%as.matrix() # THIS CANNOT BE OF CLASS DATA.FRAME # Error in xtfrm.data.frame(x) : cannot xtfrm data frames
    table(Y);str(Y)
    M<-big_spLogReg(X,Y,alphas=c(10^(-(0:4))),K=10,warn=F,ncores=NCORES)
  }
  
  saveRDS(M, paste0(tmp,OUTCOME,"_FullModel.rds"))

}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
cat("\n ... plotting the model performance \n")

M <- readRDS(paste0(tmp,OUTCOME,"_FullModel.rds"))
X <- big_attach(paste0(tmp,OUTCOME,"_PredMat.rds"))
if(REGRESSION=="LINEAR"){
  Y<-Df%>%select(Y=all_of(OUTCOME))%>%mutate(Y=as.numeric(Y))%>%as.matrix() # THIS CANNOT BE OF CLASS DATA.FRAME # Error in xtfrm.data.frame(x) : cannot xtfrm data frames
} else {
  Y<-Df%>%select(Y=all_of(OUTCOME))%>%mutate(Y=as.integer(Y)-1)%>%as.matrix() # THIS CANNOT BE OF CLASS DATA.FRAME # Error in xtfrm.data.frame(x) : cannot xtfrm data frames
}
table(Y)

ggsave(
  plot(M)+theme_classic()+scale_fill_viridis_d(option="A")+scale_color_viridis_d(option="A"),
  file=paste0(tmp,OUTCOME,"_FullModelCrossValidation.png"),
  device="png",dpi=320,width=3600,height=1200,units="px"
)

pred <- predict(M, X)
auc  <- AUC(pred, Y)
pred <- cbind(ID,pred) 
pred <- pred%>%as.data.frame()%>%mutate(Y=Y)%>%mutate(P=factor(Y,levels=c(0,1),labels=c("Non-participant","Participant")))
fwrite(pred, paste0(tmp,OUTCOME,"_FullPred.csv"))
cat("\n*** AUC  =", auc, "***\n")
ggsave(
  ggplot(pred,aes(x=pred,fill=P))+geom_density(alpha=0.4)+
  labs(x="Predicted probability of partcipation",subtitle=paste0("AUC = ",round(auc,2)),fill="") +
  theme_classic()+scale_fill_viridis_d(option="B")+theme(legend.position=c(0.9,0.2))+
  scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0)),
  file=paste0(tmp,OUTCOME,"_FullModelClassification.png"),
  device="png",dpi=320,width=1800,height=1200,units="px"
)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
cat("\n ... plotting the coefficients \n")

M       <- readRDS(paste0(tmp,OUTCOME,"_FullModel.rds"))
pred    <- fread(paste0(tmp,OUTCOME,"_FullPred.csv"))
pVars   <- fread(paste0(tmp,OUTCOME,"_PredVarNames.csv"),header=T)
pVars   <- colnames(pVars)
pVars1  <- str_sub(pVars,1,15)
pVars2  <- str_sub(pVars,-15,-1)
pVars   <- ifelse(nchar(pVars)>30, paste0(pVars1,"..",pVars2),pVars)
pVars   <- make.unique(pVars)

BetaLs  <- sapply(M[[1]], function(x) x[2])
BetaDf  <- do.call("rbind",BetaLs) 
avBeta  <- colMeans(BetaDf)

PredImp <- data.frame(pVars,avBeta)%>%mutate(abBeta=abs(avBeta))%>%arrange(-abBeta) %>% 
  mutate(Direction=ifelse(avBeta>0,"Increased probability","Decreased probability"))

fwrite(PredImp, paste0(tmp,OUTCOME,"_FullPred_Coefficients.csv"))

NPRED <- nrow(PredImp)

ggsave(
  ggplot(PredImp,aes(x=reorder(pVars,abBeta),y=abBeta,fill=Direction)) +
    geom_bar(stat="identity")+
    labs(x="Predictor",y="Absolute Beta",title=paste0("All ",NPRED," Predictors (ranked by importance)"),fill="Direction of effect") +
    coord_flip()+theme_classic()+scale_fill_viridis_d(option="C")+theme(legend.position=c(0.8,0.2)),
  file=paste0(tmp,OUTCOME,"_FullPred_AllPredictors.png"),device="png",dpi=320,width=3600,height=5600,units="px"
)

ggsave(
  ggplot(PredImp[1:20,],aes(x=reorder(pVars,abBeta),y=abBeta,fill=Direction)) +
    geom_bar(stat="identity")+
    labs(x="Predictor",y="Absolute Beta",title="Top 20 Predictors",fill="Direction of effect") +
    coord_flip()+theme_classic()+scale_fill_viridis_d(option="C")+theme(legend.position=c(0.8,0.2)),
  file=paste0(tmp,OUTCOME,"_FullPred_Top20Predictors.png"),device="png",dpi=320,width=3600,height=1200,units="px"
)

EliID  <- data.table::fread("Data/InitialParticipation_NarrowEligibility_CriteriaVars_UniqueMoLNR.csv")%>% 
  select(all_of(paste0(WHO,"LNR")),LNR)

Link   <- data.table::fread("/tsd/p805/data/durable/projects/crayner/MoBa_SSB_IDs_Trios20250415.csv") %>% 
  select(all_of(paste0(WHO,"LNR")),LNR,IID=all_of(paste0(WHO,"IID")))%>%filter(LNR%in%EliID$LNR) 

pred   <- pred %>% purrr::set_names(paste0(WHO,"LNR"),"LNR",paste0(OUTCOME,"PP"),paste0(OUTCOME),"P") %>% inner_join(Link,by=c("LNR",paste0(WHO,"LNR")))

fwrite(pred, paste0(tmp,OUTCOME,"_FullPredLinked.csv"))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
cat("\n ... FINISHED \n")

remove(list = ls())
q(save="no")
exit()



