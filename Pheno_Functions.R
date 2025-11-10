### WRITING FUNCTIONS

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
winsorize <- function(x,lowerP=0.01,upperP=0.99) {
  lowerB<-quantile(x,lowerP,na.rm=T);upperB<-quantile(x,upperP,na.rm=T)
  x[x<lowerB]<-lowerB; x[x>upperB]<-upperB; return(x)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GetUnwMean <- function(X,I){
  
  VARIABLE<-paste0(I);N<-sum(!is.na(X));L<-length(unique(na.omit(X)))
  
  if(L==2){
    EST <- mean(x=X,na.rm=T)
    VAR   <- sqrt((EST*(1-EST)/N))
    TAB  <- data.frame("Variable"=VARIABLE,"Estimate"=EST,"Variability"=VAR,"N"=N,"Type"="ProportionSE")   
    # } else if(is.integer(X)){
    #   EST <- median(x=X,na.rm=T)
    #   VAR   <- IQR(x=X,na.rm=T)
    #   TAB  <- data.frame("Variable"=VARIABLE,"Estimate"=EST,"Variability"=VAR,"N"=N,"Type"="MedianIQR")   
    } else {
      EST <- mean(X,na.rm=T)
      VAR   <- sd(X,na.rm=T)
      TAB  <- data.frame("Variable"=VARIABLE,"Estimate"=EST,"Variability"=VAR,"N"=N,"Type"="MeanSD")   
    }
  return(TAB)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wtd.med <- function(x,w) {
  w<-w[order(x,na.last=NA)];x<-x[order(x,na.last=NA)]
  prob<-cumsum(w)/sum(w)
  ps<-which(abs(prob-.5)==min(abs(prob-.5)))
  return(x[ps])
}
wtd.iqr <- function(x,w=F,na.rm=T) {
  Hmisc::wtd.quantile(x,probs=0.75,weight=w,na.rm=na.rm)-Hmisc::wtd.quantile(x,probs=0.25,weight=w,na.rm=na.rm)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GetWeiMean <- function(X,W,I){

  VARIABLE<-paste0(I);N<-sum(!is.na(X)&!is.na(W));L<-length(unique(na.omit(X)))
  
  if(L==2){
    # EST <- Hmisc::wtd.mean(x=X,w=W,na.rm=T)
    # VAR   <- sqrt((EST*(1-EST)/N))
    # VAR   <- sqrt(Hmisc::wtd.var(x=X,w=W,na.rm=T))
    dat  <- data.frame(var=X,w=W)
    des  <- svydesign(ids=~1,weights=~w,data=dat)
    est  <- svymean(~var,des,na.rm=TRUE)
    EST  <- as.numeric(est)
    VAR  <- as.numeric(SE(est))
    TAB  <- data.frame("Variable"=VARIABLE,"Estimate"=EST,"Variability"=VAR,"N"=N,"Type"="ProportionSE")   
    
  # } else if(is.integer(X)){
  #   EST <- wtd.med(x=X,w=W)
  #   VAR   <- wtd.iqr(x=X,w=W,na.rm=T)
  #   TAB  <- data.frame("Variable"=VARIABLE,"Estimate"=EST,"Variability"=VAR,"N"=N,"Type"="MedianIQR")   
    
  } else {
    EST <- Hmisc::wtd.mean(x=X,w=W,na.rm=T)
    VAR   <- sqrt(Hmisc::wtd.var(x=X,w=W,na.rm=T))
    TAB  <- data.frame("Variable"=VARIABLE,"Estimate"=EST,"Variability"=VAR,"N"=N,"Type"="MeanSD")   
    
  }
  return(TAB)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PairwiseComparisons <- function(A,B){
  DF <- A %>% full_join(B,"Variable") 
  
  DF_m <- DF %>% filter(Type.x=="MeanSD") %>% 
    rowwise() %>% 
    mutate(Z_stat  = (Estimate.x-Estimate.y)/sqrt(Variability.x^2/N.x+ Variability.y^2/N.y) ) %>% 
    mutate(Z_df    = (((Variability.x^2/N.x)+(Variability.y^2/N.y))^2/(( (Variability.x^2/N.x)^2/(N.x-1))+(((Variability.y^2/N.y)^2)/(N.y-1))))) %>% 
    mutate(Z_p     = 2*pt(-abs(Z_stat),df=Z_df)) %>% # WELCHES (MORE ROBUST TO NON-NORMALITY)
    mutate(Z_plog10= -(log(2)+pt(-abs(Z_stat),df=Z_df,log.p=TRUE))/log(10)) %>% 
    # mutate(Z_2.5   = uniroot(function(ncp) {pt(Z_stat,Z_df,ncp) - 0.975} , c(-1000, 1000))$root) %>% 
    # mutate(Z_97.5  = uniroot(function(ncp) {pt(Z_stat,Z_df,ncp) - 0.025} , c(-1000, 1000))$root) %>% 
    
    mutate(Cohens  = abs((Estimate.x-Estimate.y)/sqrt(((N.x-1)*Variability.x^2+(N.y-1)*Variability.y^2)/(N.x+N.y-2))))%>% # COHENS D
    mutate(df_p    = N.x + N.y - 2)%>%
    mutate(t_p     = Cohens / sqrt(1/N.x + 1/N.y))%>%
    mutate(ncp_lower = uniroot(function(ncp) { pt(t_p, df_p, ncp) - 0.975 }, c(-1000, 1000))$root)%>%
    mutate(ncp_upper = uniroot(function(ncp) { pt(t_p, df_p, ncp) - 0.025 }, c(-1000, 1000))$root)%>%
    mutate(Cohens_2.5  = ncp_lower * sqrt(1/N.x + 1/N.y))%>%
    mutate(Cohens_97.5 = ncp_upper * sqrt(1/N.x + 1/N.y))%>%
    # mutate(Cohens_2.5   = Z_2.5  * (sqrt(1/N.x + 1/N.y))) %>%
    # mutate(Cohens_97.5  = Z_97.5 * (sqrt(1/N.x + 1/N.y))) %>%
    
    mutate(F_stat=Variability.x^2/Variability.y^2) %>% # F-TEST (NOT AS ROBUST AS LEVENES BUT CAN USE SUMMARY DATA)
    mutate(F_p=pf(q=F_stat,df1=N.x-1,df2=N.y-1))%>%
    mutate(F_p=2*(min(F_p,(1-F_p))))%>%
    mutate(F_plog_L=pf(q=F_stat,df1=N.x-1,df2=N.y-1,log.p=T,lower.tail=T),
           F_plog_R=pf(q=F_stat,df1=N.x-1,df2=N.y-1,log.p=T,lower.tail=F)) %>%
    mutate(F_plog10=-(log(2)+pmin(F_plog_L,F_plog_R))/log(10)) %>% 
    ungroup()
  
  DF_p <- DF %>% filter(Type.x=="ProportionSE") %>% 
    rowwise() %>% 
    mutate(Z_stat=(Estimate.x-Estimate.y)^2/sqrt(Variability.x^2+Variability.y^2)^2) %>% # WALD statistic 
    mutate(Z_p=2*pchisq(Z_stat,df=1,lower.tail=F)) %>% 
    mutate(Z_plog10=-pchisq(Z_stat,df=1,lower.tail=FALSE,log.p=TRUE)/log(10)) %>%
    
    mutate(Cohens      = abs(2*asin(sqrt(Estimate.x))-2*asin(sqrt(Estimate.y)))) %>% # COHENS H
    mutate(Cohens_2.5  = Cohens - 1.96 * sqrt(1/N.x + 1/N.y)) %>% 
    mutate(Cohens_97.5 = Cohens + 1.96 * sqrt(1/N.x + 1/N.y)) %>% 
    ungroup()
  
  DF <- DF_m %>% bind_rows(DF_p)
  
  return(DF)
}


round2<-function(x){if(abs(x)>0){return(round(x,digits=2))}else {return(signif(x,digits=2))}}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BvUnweightedModel  <- function(dat,y,x){
  cat("\n*****************************************\n")
  cat("Performing unweighted regression analyses\n\n")
  cat("Outcome:\t\t",y,"\n")
  cat("Exposure:\t\t",x,"\n")

  Formula <- as.formula(paste0(y,"~",x))
  
  cat("\nFormula:\n")
  print(Formula)
  
  Nval <- length(unique(na.omit(dat[[paste0(y)]])))
  cat("\nOutcome has",Nval,"unique values\n")
  
  if(Nval>2){    
    cat("\nRunning linear regression\n")
    Model <- lm(formula=Formula,data=dat)
  }else{
    cat("\nRunning logistic regression\n")
    Model <- glm(formula=Formula,data=dat,family="binomial")
  }
  cat("\nModel completed\n")
  cat("\n*****************************************\n")
  return(Model)
}


UnweightedModel  <- function(dat,y,x,z){
  cat("\n*****************************************\n")
  cat("Performing unweighted regression analyses\n\n")
  cat("Outcome:\t\t",y,"\n")
  cat("Exposure:\t\t",x,"\n")
  cat(length(z),"Covariates:\t\t",z,"\n")
  
  if(length(z)!=0){
    Formula <- as.formula(paste0(y,"~",x,"+",paste0(z,collapse="+")))
  }else{
    Formula <- as.formula(paste0(y,"~",x))
  }
  
  cat("\nFormula:\n")
  print(Formula)
  
  Nval <- length(unique(na.omit(dat[[paste0(y)]])))
  cat("\nOutcome has",Nval,"unique values\n")
  
  if(Nval>2){    
    cat("\nRunning linear regression\n")
    Model <- lm(formula=Formula,data=dat)
  }else{
    cat("\nRunning logistic regression\n")
    Model <- glm(formula=Formula,data=dat,family="binomial")
  }
  cat("\nModel completed\n")
  cat("\n*****************************************\n")
  return(Model)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BvWeightedModel  <- function(dat,y,x,w){
  cat("\n*****************************************\n")
  cat("Performing *weighted* regression analyses\n\n")
  cat("Outcome:\t\t",y,"\n")
  cat("Exposure:\t\t",x,"\n")
  Formula<-
    as.formula(paste0(y,"~",x))
  dat<-
    dat%>% select(all_of(x),all_of(y),all_of(w)) %>% na.omit()
  cat("Weights:\t\t",w,"\n")
  dat$w <- dat[[paste0(w)]]
  
  cat("\nFormula:\n")
  print(Formula)
  
  design <- svydesign(ids=~1,data=dat,weights=~w)
  
  Nval <- length(unique(na.omit(dat[[paste0(y)]])))
  cat("\nOutcome has",Nval,"unique values\n")
  
  if(Nval>2){    
    cat("\nRunning linear regression\n")
    Model  <- svyglm(formula=Formula,design=design,family=gaussian,deff=T)
  }else{
    cat("\nRunning logistic regression\n")
    Model  <- svyglm(formula=Formula,design=design,family=binomial,deff=T)
  }
  cat("\nModel completed\n")
  cat("\n*****************************************\n")
  return(Model)
} 


WeightedModel  <- function(dat,y,x,w,z){
  cat("\n*****************************************\n")
  cat("Performing *weighted* regression analyses\n\n")
  cat("Outcome:\t\t",y,"\n")
  cat("Exposure:\t\t",x,"\n")
  cat(length(z),"Covariates:\t\t",z,"\n")
  
  if(length(z)!=0){
    Formula <- as.formula(paste0(y,"~",x,"+",paste0(z,collapse="+")))
  }else{
    Formula <- as.formula(paste0(y,"~",x))
  }
  
  dat <- dat %>% select(all_of(x),all_of(y),all_of(z),all_of(w)) %>% na.omit()
  
  cat("Weights:\t\t",w,"\n")
  
  dat$w <- dat[[paste0(w)]]
  
  cat("\nFormula:\n")
  print(Formula)
  
  design <- svydesign(ids=~1,data=dat,weights=~w)
  
  Nval <- length(unique(na.omit(dat[[paste0(y)]])))
  cat("\nOutcome has",Nval,"unique values\n")
  
  if(Nval>2){    
    cat("\nRunning linear regression\n")
    Model  <- svyglm(formula=Formula,design=design,family=gaussian,deff=T)
  }else{
    cat("\nRunning logistic regression\n")
    Model  <- svyglm(formula=Formula,design=design,family=binomial,deff=T)
  }
  cat("\nModel completed\n")
  cat("\n*****************************************\n")
  return(Model)
} 


GetResultsTab <- function(X,Y,Model1,Model2){
  X<-X; Y<-Y
  Coef1<-coef(summary(Model1))[2,1]
  SE1<-coef(summary(Model1))[2,2]
  Coef2<-coef(summary(Model2))[2,1]
  SE2<-coef(summary(Model2))[2,2]
  Z_SE<-sqrt(SE1^2+SE2^2)
  Z<-(Coef1-Coef2)/Z_SE
  P<-2*pnorm(-abs(Z))
  TAB<-data.frame("X"=X,"Y"=Y,"Coeff1"=Coef1,"SE1"=SE1,"Coeff2"=Coef2,"SE2"=SE2,"Zdiff"=Z,"Zdiff_P"=P)
  return(TAB)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CompareEffects <- function(Pmod,Umod,Wmod){
  Pcoef <- summary(Pmod)$coefficients %>% as.data.frame() %>% select(Pop_Estimate=Estimate)
  Ucoef <- summary(Umod)$coefficients %>% as.data.frame() %>% select(U_Estimate=Estimate)
  Wcoef <- summary(Wmod)$coefficients %>% as.data.frame() %>% select(W_Estimate=Estimate)
  COEF  <- bind_cols(Pcoef,Ucoef,Wcoef) %>% 
    mutate(UnweightedBias = Pop_Estimate - U_Estimate,
           WeightedBias   = Pop_Estimate - W_Estimate
    ) %>% 
    mutate(Improvement = ifelse(abs(WeightedBias)<abs(UnweightedBias), "Y","N"))
  return(COEF)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RunAllModels <- function(DATA,S1,S2=NULL,Y,X,Z,W1,W2=NULL,W3=NULL){
  
  # Population model
  PopModel    <- UnweightedModel(dat=DATA,y=Y,x=X,z=Z)
  PopCoef     <- coef(summary(PopModel))[2,1]
  PopSE       <- coef(summary(PopModel))[2,2]
  Pop         <- data.frame("Coeff"=PopCoef,"SE"=PopSE)
  
  # Initial participation unweighted model
  uModelS1    <- UnweightedModel(dat=DATA%>%dplyr::filter(.data[[S1]]==1),y=Y,x=X,z=Z)
  S1Coef      <- coef(summary(uModelS1))[2,1]
  S1SE        <- coef(summary(uModelS1))[2,2]
  S1zd        <- (PopCoef - S1Coef) / sqrt(PopSE^2 + S1SE^2)
  S1zdp       <- 2*pnorm(abs(S1zd), lower.tail = FALSE)
  S1u         <- data.frame("Coeff"=S1Coef,"SE"=S1SE,"Zdiff"=S1zd,"Zdiff_P"=S1zdp)
  TABLE       <- Pop %>% bind_rows(S1u)
  
  # Initial participation weighted model (Q1w)
  wModelS1W1  <- WeightedModel(dat=DATA%>%dplyr::filter(.data[[S1]]==1),y=Y,x=X,z=Z,w=W1)
  S1W1Coef    <- coef(summary(wModelS1W1))[2,1]
  S1W1SE      <- coef(summary(wModelS1W1))[2,2]
  S1W1zd      <- (PopCoef - S1W1Coef) / sqrt(PopSE^2 + S1W1SE^2)
  S1W1zdp     <- 2*pnorm(abs(S1W1zd), lower.tail = FALSE)
  S1W1        <- data.frame("Coeff"=S1W1Coef,"SE"=S1W1SE,"Zdiff"=S1W1zd,"Zdiff_P"=S1W1zdp)
  TABLE       <- TABLE %>% bind_rows(S1W1)
  
  if(!is.null(S2)){
    # Continued participation unweighted model
    uModelS2    <- UnweightedModel(dat=DATA%>%dplyr::filter(.data[[S2]]==1),y=Y,x=X,z=Z)
    S2Coef      <- coef(summary(uModelS2))[2,1]
    S2SE        <- coef(summary(uModelS2))[2,2]
    S2zd        <- (PopCoef - S2Coef) / sqrt(PopSE^2 + S2SE^2)
    S2zdp       <- 2*pnorm(abs(S2zd), lower.tail = FALSE)
    S2u         <- data.frame("Coeff"=S2Coef,"SE"=S2SE,"Zdiff"=S2zd,"Zdiff_P"=S2zdp)
    TABLE       <- TABLE %>% bind_rows(S2u)
    
    # Continued participation weighted model (Q1w)
    wModelS2W1  <- WeightedModel(dat=DATA%>%dplyr::filter(.data[[S2]]==1),y=Y,x=X,z=Z,w=W1)
    S2W1Coef    <- coef(summary(wModelS2W1))[2,1]
    S2W1SE      <- coef(summary(wModelS2W1))[2,2]
    S2W1zd      <- (PopCoef - S2W1Coef) / sqrt(PopSE^2 + S2W1SE^2)
    S2W1zdp     <- 2*pnorm(abs(S2W1zd), lower.tail = FALSE)
    S2W1        <- data.frame("Coeff"=S2W1Coef,"SE"=S2W1SE,"Zdiff"=S2W1zd,"Zdiff_P"=S2W1zdp)
    TABLE       <- TABLE %>% bind_rows(S2W1)
    
    if(!is.null(W2)){
      # Continued participation weighted model (Q9w)
      wModelS2W2  <- WeightedModel(dat=DATA%>%dplyr::filter(.data[[S2]]==1),y=Y,x=X,z=Z,w=W2)
      S2W2Coef    <- coef(summary(wModelS2W2))[2,1]
      S2W2SE      <- coef(summary(wModelS2W2))[2,2]
      S2W2zd      <- (PopCoef - S2W2Coef) / sqrt(PopSE^2 + S2W2SE^2)
      S2W2zdp     <- 2*pnorm(abs(S2W2zd), lower.tail = FALSE)
      S2W2        <- data.frame("Coeff"=S2W2Coef,"SE"=S2W2SE,"Zdiff"=S2W2zd,"Zdiff_P"=S2W2zdp)
      TABLE       <- TABLE %>% bind_rows(S2W2)
    }
    if(!is.null(W3)){
      # Continued participation unweighted model (Q9w)
      wModelS2W3  <- WeightedModel(dat=DATA%>%dplyr::filter(.data[[S2]]==1),y=Y,x=X,z=Z,w=W3)
      S2W3Coef    <- coef(summary(wModelS2W3))[2,1]
      S2W3SE      <- coef(summary(wModelS2W3))[2,2]
      S2W3zd      <- (PopCoef - S2W3Coef) / sqrt(PopSE^2 + S2W3SE^2)
      S2W3zdp     <- 2*pnorm(abs(S2W3zd), lower.tail = FALSE)
      S2W3        <- data.frame("Coeff"=S2W3Coef,"SE"=S2W3SE,"Zdiff"=S2W3zd,"Zdiff_P"=S2W3zdp)
      TABLE       <- TABLE %>% bind_rows(S2W3)
    }
  }
  
  cat(paste0("the population effect of ",X," on ",Y," is ",round(PopCoef,3)," (s.e=",round(PopSE,3),";N=",nobs(PopModel),")\n"))
  
  cat(paste0("the sample effect (S=1) of ",X," on ",Y," is ",round(S1Coef,3),"(s.e=",round(S1SE,3),";N=",nobs(uModelS1),")\n"))
  if(S1zdp <.05){
    cat(paste0("this effect is significantly different from the population effect (z=",round(S1zd,3),"; p=",round(S1zdp,3),")\n"))
  } else {cat(paste0("this effect does not differ from the population effect (z=",round(S1zd,3),"; p=",round(S1zdp,3),")\n"))}
  
  cat(paste0("the weighted sample effect (S=1; W=1) of",X," on ",Y," is ",round(S1W1Coef,3),"(s.e=",round(S1W1SE,3),")\n"))
  if(S1W1zdp <.05){
    cat(paste0("this effect is significantly different from the population effect (z=",round(S1W1zd,3),"; p=",round(S1W1zdp,3),")\n"))
  } else {cat(paste0("this effect does not differ from the population effect (z=",round(S1W1zd,3),"; p=",round(S1W1zdp,3),")\n"))}
  
  if(!is.null(S2)){
    cat(paste0("the sample effect (S=2) of ",X," on ",Y," is ",round(S2Coef,3),"(s.e=",round(S2SE,3),";N=",nobs(uModelS2),")\n"))
    if(S2zdp <.05){
      cat(paste0("this effect is significantly different from the population effect (z=",round(S2zd,3),"; p=",round(S2zdp,3),")\n"))
    } else {cat(paste0("this effect does not differ from the population effect (z=",round(S2zd,3),"; p=",round(S2zdp,3),")\n"))}
    
    cat(paste0("the weighted sample effect (S=1; W=1) of ",X," on ",Y," is ",round(S2W1Coef,3),"(s.e=",round(S2W1SE,3),")\n"))
    if(S2W1zdp <.05){
      cat(paste0("this effect is significantly different from the population effect (z=",round(S2W1zd,3),"; p=",round(S2W1zdp,3),")\n"))
    } else {cat(paste0("this effect does not differ from the population effect (z=",round(S2W1zd,3),"; p=",round(S2W1zdp,3),")\n"))}
  
    if(!is.null(W2)){
      cat(paste0("the weighted sample effect (S=2; W=2) of ",X," on ",Y," is ",round(S2W2Coef,3),"(s.e=",round(S2W2SE,3),")\n"))
      if(S2W2zdp <.05){
        cat(paste0("this effect is significantly different from the population effect (z=",round(S2W2zd,3),"; p=",round(S2W2zdp,3),")\n"))
      } else {cat(paste0("this effect does not differ from the population effect (z=",round(S2W2zd,3),"; p=",round(S2W2zdp,3),")\n"))}
    }
    
    if(!is.null(W3)){
      cat(paste0("the weighted sample effect (S=2; W=3) of ",X," on ",Y," is ",round(S2W3Coef,3),"(s.e=",round(S2W3SE,3),")\n"))
      if(S2W3zdp <.05){
        cat(paste0("this effect is significantly different from the population effect (z=",round(S2W3zd,3),"; p=",round(S2W3zdp,3),")\n"))
      } else {cat(paste0("this effect does not differ from the population effect (z=",round(S2W3zd,3),"; p=",round(S2W3zdp,3),")\n"))}
    }
  }
  return(TABLE)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
WordTable <- function(dt,fPath){
  # Create a flextable object
  ft <- flextable::flextable(dt)
  # Basic APA-like styling
  ft <- flextable::bold(ft,part="header") # Bold header
  ft <- flextable::font(ft,font="Times New Roman",part="all") # Set font to Times New Roman
  ft <- flextable::fontsize(ft,size=9,part="all") # Set font to Times New Roman
  ft <- flextable::autofit(ft) # Autofit columns
  ft <- flextable::align(ft,align="left",part="all") # Autofit columns
  
  # Create a Word document
  doc <- officer::read_docx()
  # Add the flextable to the document
  doc <- doc %>% officer::body_end_section_landscape()
  doc <- doc %>% flextable::body_add_flextable(value=ft,pos="on")
  # Save the Word document
  ls  <- officer::prop_section(page_size=page_size(orient="landscape"))
  flextable::save_as_docx(ft,path=fPath,pr_section=ls)
  # print(doc,target=paste0(fPath))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MeanFig_fn<-function(df){
  ggplot(
    df,
    aes(x=abs(Cohens),y=Variable,shape=`P-value`,colour=Analysis))+
    geom_vline(xintercept=0,linetype="dashed",colour="red")+
    geom_linerange(aes(xmin=Cohens_2.5,xmax=Cohens_97.5,colour=Analysis),linewidth=4,alpha = 0.3) +
    geom_point(aes(colour=Analysis,shape=`P-value`),size=2,alpha=.8)+
    # geom_pointrange(aes(xmin=abs(Cohens_2.5),xmax=abs(Cohens_97.5)),alpha=.8)+    
    # geom_errorbar(aes(xmin=abs(Cohens_2.5),xmax=abs(Cohens_97.5)))+    
    # geom_crossbar(aes(xmin=abs(Cohens_2.5),xmax=abs(Cohens_97.5)))+    
    facet_grid(cols=vars(Facet),rows=vars(Who),scales="free_y",space="free_y") + 
    scale_y_discrete(limits=rev)+
    scale_shape_manual(values=c(16,1))+
    scale_colour_manual(values=c("Unweighted"="#404688FF","Weighted"="#8FD744FF"))+
    labs(
      x="|Standardised mean differences|  (95% confidence intervals)",y="Variable",shape="P-value")+
    theme_classic()+
    # theme_tufte()+
    # theme_bw()+
    theme(
      legend.position="bottom",legend.justification="left",legend.title.position="top",
      legend.justification.bottom=c(0,1),
      legend.margin=margin(0,0,0,0),legend.key.height=unit(.2,"cm"),legend.key.width=unit(.2,"cm"),
      legend.key.size=unit(0.1,"lines"),
      legend.text=element_text(size=8,margin=margin(r=10,unit="pt")),
      legend.title=element_text(size=8,face="bold"),
      axis.title.x=element_text(hjust=0,size=10),axis.title.y=element_blank(),
      strip.background=element_rect(fill="#FFFFC5",colour="#FFFFC5")
      )+
    guides(
      color=guide_legend(order=1),
      shape=guide_legend(
        title = "Bonferoni corrected p-value")
      )+
    theme(strip.text.x=element_markdown(lineheight=1.2,hjust=0))
}



CoefFig_fn<-function(df){
  ggplot(
    df,
    aes(x=abs(Zdiff),y=X,shape=`P-value`,colour=Analysis))+
    geom_vline(xintercept=0,linetype="dashed",colour="red")+
    geom_linerange(aes(xmin=abs(Z_2.5),xmax=abs(Z_97.5),colour=Analysis),linewidth=4,alpha = 0.3)+
    geom_point(aes(colour=Analysis,shape=`P-value`),size=2,alpha=.8)+
    # geom_point(aes(colour=Analysis),size=3,alpha=.8)+
    facet_grid(cols=vars(Facet),rows=vars(Who.Y),scales="free_y",space="free_y") + 
    scale_y_discrete(limits=rev)+
    # scale_shape_manual(values=c("Yes (direction)"=8, "Yes (magnitude)"=16,"No"=1))+
    scale_shape_manual(values=c(1,16,8))+
    scale_colour_manual(values=c("Unweighted"="#404688FF","Weighted"="#8FD744FF"))+
    labs(x="|z-statistics|",y="Variable",shape="P-value")+
    theme_classic()+
    theme(
      legend.position="bottom",legend.justification="left",legend.title.position="top",legend.justification.bottom=c(0,1),
      legend.margin=margin(0,0,0,0),legend.key.height=unit(.2,"cm"),legend.key.width=unit(.2,"cm"),legend.key.size=unit(0.1,"lines"),
      legend.text=element_text(size=8,margin=margin(r=10,unit="pt")),legend.title=element_text(size=8,face="bold"),
      axis.title.x=element_text(hjust=0,size=10),axis.title.y=element_blank(),
      strip.background=element_rect(fill="#FFFFC5",colour="#FFFFC5")
      )+
    guides(
      color=guide_legend(order=1),
      shape=guide_legend(
        title = "Bonferoni corrected p-value")
    )+
    theme(strip.text.x=element_markdown(lineheight=1.2,hjust=0))
}


CoefFig_fn<-function(df){
  # MIN <- floor(min(df[["Coeffdiff_2.5"]]))
  # MAX <- floor(min(df[["Coeffdiff_97.5"]]))
  ggplot(
    df,
    aes(x=abs(Coeffdiff),y=X,shape=`P-value`,colour=Analysis))+
    geom_vline(xintercept=0,linetype="dashed",colour="red")+
    geom_linerange(aes(xmin=Coeffdiff_2.5,xmax=Coeffdiff_97.5,colour=Analysis),linewidth=4,alpha = 0.3)+
    geom_point(aes(colour=Analysis,shape=`P-value`),size=2,alpha=.8)+
    # geom_point(aes(colour=Analysis),size=3,alpha=.8)+
    facet_grid(cols=vars(Facet),rows=vars(Who.Y),scales="free_y",space="free_y") + 
    scale_y_discrete(limits=rev)+
    # scale_shape_manual(values=c("Yes (direction)"=8, "Yes (magnitude)"=16,"No"=1))+
    scale_shape_manual(values=c(1,16,8))+
    scale_colour_manual(values=c("Unweighted"="#404688FF","Weighted"="#8FD744FF"))+
    labs(
      x = expression(paste("|" * beta["diff"] * "| (95% confidence intervals)")),
      y="Variable",shape="P-value")+
    theme_classic()+
    theme(
      legend.position="bottom",legend.justification="left",legend.title.position="top",legend.justification.bottom=c(0,1),
      legend.margin=margin(0,0,0,0),legend.key.height=unit(.2,"cm"),legend.key.width=unit(.2,"cm"),legend.key.size=unit(0.1,"lines"),
      legend.text=element_text(size=8,margin=margin(r=10,unit="pt")),legend.title=element_text(size=8,face="bold"),
      axis.title.x=element_text(hjust=0,size=10),axis.title.y=element_blank(),
      strip.background=element_rect(fill="#FFFFC5",colour="#FFFFC5")
    )+
    guides(
      color=guide_legend(order=1),
      shape=guide_legend(
        title = "Bonferoni corrected p-value (type of effect)")
    )+
    theme(strip.text.x=element_markdown(lineheight=1.2,hjust=0))
}

