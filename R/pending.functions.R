

# ============= VARIANCE HETEROGENEITY FUNCTIONS ====================
# Function takes 2 methods. "Classic-Levene" and "Brown-Forsythe". The latter is more
# robust to skewed distributions. The default is "classic.Levene". Note: Selecting
# "Brown-Forsythe" gives the exact same result as levene.test function in lawstat package




# ============= META-ANALYSIS FUNCTIONS ====================
original.Stouffer.test <- function(p,Weights,log_p.value=FALSE) {
  #p is a vector of p-values
  #Weights is typically sample size
  #p=datatable[,RIF.MetaTau_p.value]; Weights=datatable[,N];log_p.value=FALSE
  #Weights=NA
  if(missing(Weights)) {
    w <- rep(1, length(p))/length(p)
  } else if(length(Weights)==length(p)){
    w <- sqrt(Weights)/sum(sqrt(Weights))
  } else {
    stop("Length of p and w must equal!")
  }
  Zi <- qnorm(p/2,lower=FALSE) # note /2 is modification for two-sided pvalues
  Z  <- sum(w*Zi)/sqrt(sum(w^2))
  if(log_p.value){
    p.val <- -1*pnorm(Z, lower.tail=FALSE, log=TRUE)
  } else{
    p.val <- pnorm(Z, lower.tail=FALSE)
  }
  return(c(Z,p.val))
  #return(as.vector(p.val))
}
original.Fischer.test <- function(p,log_p.value=FALSE) {
  #p is a vector of p-values
  #Weights is typically sample size
  #p=datatable[,RIF.MetaTau_p.value]
  T <- sum(-2*log(p))
  if(log_p.value){
    p.val <- pchisq(T, df=length(p)*2, lower.tail=FALSE, log=TRUE)
  } else{
    p.val <- pchisq(T, df=length(p)*2, lower.tail=FALSE)
  }
  return(c(T,p.val))
  #return(as.vector(p.val))
}
Stouffer.test <- function(x,Statistic="z",Weights,log_p.value=FALSE) {
  #z is a vector of z-values
  #Weights is typically sample size
  #z=datatable[,RIF.MetaTau_tval]; Weights=datatable[,N];log_p.value=FALSE
  #Weights=NA
  if(missing(Weights)) {
    w <- rep(1, length(x))/length(x)
  } else if(length(Weights)==length(x)){
    w <- sqrt(Weights)/sum(sqrt(Weights))
  } else {
    stop("Length of x and w must equal!")
  }
  if(Statistic=="z"){
    Zi <- abs(x)
    Z  <- sum(w*Zi)/sqrt(sum(w^2))
  } else if(Statistic=="p"){
    Zi <- qnorm(x/2,lower=FALSE) # note /2 is modification for two-sided pvalues
    Z  <- sum(w*Zi)/sqrt(sum(w^2))
  }
  if(log_p.value){
    p.val <- -1*pnorm(Z, lower.tail=FALSE, log=TRUE)
  } else{
    p.val <- pnorm(Z, lower.tail=FALSE)
  }
  return(c(Z,p.val))
  #return(as.vector(p.val))
}
Fischer.test <- function(x,Statistic="z",log_p.value=FALSE) {
  #z is a vector of z-values
  #Weights is typically sample size
  #z=datatable[,RIF.MetaTau_tval]
  if(Statistic=="z"){
    T <- sum(-2*pnorm(x,lower.tail=TRUE,log.p=TRUE))
  } else if(Statistic=="p"){
    T <- sum(-2*log(x))
  }
  if(log_p.value){
    p.val <- pchisq(T, df=length(x)*2, lower.tail=FALSE, log=TRUE)
  } else{
    p.val <- pchisq(T, df=length(x)*2, lower.tail=FALSE)
  }
  return(c(T,p.val))
  #return(as.vector(p.val))
}

MCQR.metafor.function <- function(datatable,fun1=Stouffer.test,fun2=Fischer.test,
                                  log_p.value=FALSE){
  #datatable=i.Results[SNP==i.SNPs[l]];fun1=Stouffer.test;fun2=Fischer.test;log_p.value=FALSE
  #colnames(datatable)
  # SUMMARIZE SNP INFO
  Results <- c(SNP=unique(datatable[,SNP]),
               Chr=unique(datatable[complete.cases(Chr)][,Chr]),
               N=sum(datatable[,N],na.rm=TRUE),
               avr.EAF=mean(datatable[,EAF],na.rm=TRUE),
               min.EAF=min(datatable[,EAF],na.rm=TRUE),
               max.EAF=max(datatable[,EAF],na.rm=TRUE))
  Note <- NULL
  # META-ANALYSIS OF MCQR.Median ESTIMATES
  meta <- tryCatch({rma(yi=MCQR.Median_Beta, sei=MCQR.Median_SE,
                        data=datatable[complete.cases(MCQR.Median_Beta,MCQR.Median_SE)],
                        method="REML",control=list(maxiter=1000, stepadj=.5))},
                   error=function(e) e)
  if(!class(meta)[1]=="rma.uni"){
    Results <- c(Results,MCQR.Median_Beta=NA,MCQR.Median_SE=NA,MCQR.Median_lower=NA,
                 MCQR.Median_upper=NA,MCQR.Median_tval=NA,
                 MCQR.Median_p.value=NA,MCQR.Median_I2=NA,MCQR.Median_QEp=NA)
    Note <- c(Note,"MCQR.Median Meta-Analysis Failed")
  } else {
    Results <- c(Results,MCQR.Median_Beta=meta$b,MCQR.Median_SE=meta$se,
                 MCQR.Median_lower=meta$ci.lb,MCQR.Median_upper=meta$ci.ub,
                 MCQR.Median_tval=meta$zval,MCQR.Median_p.value=meta$pval,
                 MCQR.Median_I2=meta$I2,MCQR.Median_QEp=meta$QEp)
  }
  # META-ANALYSIS OF abs_MCQR.Median ESTIMATES
  meta <- tryCatch({rma(yi=abs_MCQR.Median_Beta, sei=MCQR.Median_SE,
                        data=datatable[complete.cases(abs_MCQR.Median_Beta,MCQR.Median_SE)],
                        method="REML",control=list(maxiter=1000, stepadj=.5))},
                   error=function(e) e)
  if(!class(meta)[1]=="rma.uni"){
    Results <- c(Results,abs_MCQR.Median_Beta=NA,abs_MCQR.Median_SE=NA,abs_MCQR.Median_lower=NA,
                 abs_MCQR.Median_upper=NA,abs_MCQR.Median_tval=NA,
                 abs_MCQR.Median_p.value=NA,abs_MCQR.Median_I2=NA,abs_MCQR.Median_QEp=NA)
    Note <- c(Note,"abs_MCQR.Median Meta-Analysis Failed")
  } else {
    Results <- c(Results,abs_MCQR.Median_Beta=meta$b,abs_MCQR.Median_SE=meta$se,
                 abs_MCQR.Median_lower=meta$ci.lb,abs_MCQR.Median_upper=meta$ci.ub,
                 abs_MCQR.Median_tval=meta$zval,abs_MCQR.Median_p.value=meta$pval,
                 abs_MCQR.Median_I2=meta$I2,abs_MCQR.Median_QEp=meta$QEp)
  }

  # META-ANALYSIS OF MCQR.MetaTau ESTIMATES
  meta <- tryCatch({rma(yi=MCQR.MetaTau_Beta, sei=MCQR.MetaTau_SE,
                        data=datatable[complete.cases(MCQR.MetaTau_Beta,MCQR.MetaTau_SE)],
                        method="REML",control=list(maxiter=1000, stepadj=.5))},
                   error=function(e) e)
  if(!class(meta)[1]=="rma.uni"){
    Results <- c(Results,MCQR.MetaTau_Beta=NA,MCQR.MetaTau_SE=NA,MCQR.MetaTau_lower=NA,
                 MCQR.MetaTau_upper=NA,MCQR.MetaTau_tval=NA,
                 MCQR.MetaTau_p.value=NA,MCQR.MetaTau_I2=NA,MCQR.MetaTau_QEp=NA)
    Note <- c(Note,"MCQR.MetaTau Meta-Analysis Failed")
  } else {
    Results <- c(Results,MCQR.MetaTau_Beta=meta$b,MCQR.MetaTau_SE=meta$se,
                 MCQR.MetaTau_lower=meta$ci.lb,MCQR.MetaTau_upper=meta$ci.ub,
                 MCQR.MetaTau_tval=meta$zval,MCQR.MetaTau_p.value=meta$pval,
                 MCQR.MetaTau_I2=meta$I2,MCQR.MetaTau_QEp=meta$QEp)
  }

  # META-ANALYSIS OF abs_MCQR.MetaTau ESTIMATES
  meta <- tryCatch({rma(yi=abs_MCQR.MetaTau_Beta, sei=MCQR.MetaTau_SE,
                        data=datatable[complete.cases(abs_MCQR.MetaTau_Beta,MCQR.MetaTau_SE)],
                        method="REML",control=list(maxiter=1000, stepadj=.5))},
                   error=function(e) e)
  if(!class(meta)[1]=="rma.uni"){
    Results <- c(Results,abs_MCQR.MetaTau_Beta=NA,abs_MCQR.MetaTau_SE=NA,abs_MCQR.MetaTau_lower=NA,
                 abs_MCQR.MetaTau_upper=NA,abs_MCQR.MetaTau_tval=NA,
                 abs_MCQR.MetaTau_p.value=NA,abs_MCQR.MetaTau_I2=NA,abs_MCQR.MetaTau_QEp=NA)
    Note <- c(Note,"abs_MCQR.MetaTau Meta-Analysis Failed")
  } else {
    Results <- c(Results,abs_MCQR.MetaTau_Beta=meta$b,abs_MCQR.MetaTau_SE=meta$se,
                 abs_MCQR.MetaTau_lower=meta$ci.lb,abs_MCQR.MetaTau_upper=meta$ci.ub,
                 abs_MCQR.MetaTau_tval=meta$zval,abs_MCQR.MetaTau_p.value=meta$pval,
                 abs_MCQR.MetaTau_I2=meta$I2,abs_MCQR.MetaTau_QEp=meta$QEp)
  }

  if(is.null(Note)){
    Results <- c(Results,Notes=NA)
  } else {
    Note <- paste(Note, sep="", collapse=", ")
    Results <- c(Results,Notes=Note)
  }
  return(Results)
}
MUQR.metafor.function <- function(datatable,fun1=Stouffer.test,fun2=Fischer.test,
                                  log_p.value=FALSE){
  # datatable=i.Results[SNP==i.SNPs[l]];fun1=Stouffer.test;fun2=Fischer.test;log_p.value=FALSE
  # colnames(datatable)

  # SUMMARIZE SNP INFO
  Results <- c(SNP=unique(datatable[,SNP]),
               Chr=unique(datatable[complete.cases(Chr)][,Chr]),
               N=sum(datatable[,N],na.rm=TRUE),
               avr.EAF=mean(datatable[,EAF],na.rm=TRUE),
               min.EAF=min(datatable[,EAF],na.rm=TRUE),
               max.EAF=max(datatable[,EAF],na.rm=TRUE))
  Note <- NULL
  # META-ANALYSIS OF LN ESTIMATES
  meta <- tryCatch({rma(yi=LN_Beta, sei=LN_SE,
                        data=datatable[complete.cases(LN_Beta,LN_SE)],method="REML",
                        control=list(maxiter=1000, stepadj=.5))},
                   error=function(e) e)
  if(!class(meta)[1]=="rma.uni"){
    Results <- c(Results,LN_Beta=NA,LN_SE=NA,LN_lower=NA,LN_upper=NA,LN_tval=NA,
                 LN_p.value=NA,LN_I2=NA,LN_QEp=NA)
    Note <- c(Note,"LN Meta-Analysis Failed")
  } else {
    Results <- c(Results,LN_Beta=meta$b,LN_SE=meta$se,LN_lower=meta$ci.lb,
                 LN_upper=meta$ci.ub,LN_tval=meta$zval,LN_p.value=meta$pval,
                 LN_I2=meta$I2,LN_QEp=meta$QEp)
  }

  # META-ANALYSIS OF rt_LN ESTIMATES
  meta <- tryCatch({rma(yi=rt_LN_Beta, sei=rt_LN_SE,
                        data=datatable[complete.cases(rt_LN_Beta,rt_LN_SE)],method="REML",
                        control=list(maxiter=1000, stepadj=.5))},
                   error=function(e) e)
  if(!class(meta)[1]=="rma.uni"){
    Results <- c(Results,rt_LN_Beta=NA,rt_LN_SE=NA,rt_LN_lower=NA,rt_LN_upper=NA,rt_LN_tval=NA,
                 rt_LN_p.value=NA,rt_LN_I2=NA,rt_LN_QEp=NA)
    Note <- c(Note,"rt_LN Meta-Analysis Failed")
  } else {
    Results <- c(Results,rt_LN_Beta=meta$b,rt_LN_SE=meta$se,rt_LN_lower=meta$ci.lb,
                 rt_LN_upper=meta$ci.ub,rt_LN_tval=meta$zval,rt_LN_p.value=meta$pval,
                 rt_LN_I2=meta$I2,rt_LN_QEp=meta$QEp)
  }

  # META-ANALYSIS OF sc_LN ESTIMATES
  meta <- tryCatch({rma(yi=sc_LN_Beta, sei=sc_LN_SE,
                        data=datatable[complete.cases(sc_LN_Beta,sc_LN_SE)],
                        method="REML",control=list(maxiter=1000, stepadj=.5))},
                   error=function(e) e)
  if(!class(meta)[1]=="rma.uni"){
    Results <- c(Results,sc_LN_Beta=NA,sc_LN_SE=NA,
                 sc_LN_lower=NA,sc_LN_upper=NA,
                 sc_LN_tval=NA,sc_LN_p.value=NA,
                 sc_LN_I2=NA,sc_LN_QEp=NA)
    Note <- c(Note,"sc_LN Meta-Analysis Failed")
  } else {
    Results <- c(Results,sc_LN_Beta=meta$b,sc_LN_SE=meta$se,
                 sc_LN_lower=meta$ci.lb,sc_LN_upper=meta$ci.ub,
                 sc_LN_tval=meta$zval,sc_LN_p.value=meta$pval,
                 sc_LN_I2=meta$I2,sc_LN_QEp=meta$QEp)
  }

  # META-ANALYSIS OF sc_LN ESTIMATES
  meta <- tryCatch({rma(yi=log_LN_Beta, sei=log_LN_SE,
                        data=datatable[complete.cases(log_LN_Beta,log_LN_SE)],
                        method="REML",control=list(maxiter=1000, stepadj=.5))},
                   error=function(e) e)
  if(!class(meta)[1]=="rma.uni"){
    Results <- c(Results,log_LN_Beta=NA,log_LN_SE=NA,
                 log_LN_lower=NA,log_LN_upper=NA,
                 log_LN_tval=NA,log_LN_p.value=NA,
                 log_LN_I2=NA,log_LN_QEp=NA)
    Note <- c(Note,"log_LN Meta-Analysis Failed")
  } else {
    Results <- c(Results,log_LN_Beta=meta$b,log_LN_SE=meta$se,
                 log_LN_lower=meta$ci.lb,log_LN_upper=meta$ci.ub,
                 log_LN_tval=meta$zval,log_LN_p.value=meta$pval,
                 log_LN_I2=meta$I2,log_LN_QEp=meta$QEp)
  }

  # META-ANALYSIS OF MUQR.Median ESTIMATES
  meta <- tryCatch({rma(yi=MUQR.Median_Beta, sei=MUQR.Median_SE,
                        data=datatable[complete.cases(MUQR.Median_Beta,MUQR.Median_SE)],
                        method="REML",control=list(maxiter=1000, stepadj=.5))},
                   error=function(e) e)
  if(!class(meta)[1]=="rma.uni"){
    Results <- c(Results,MUQR.Median_Beta=NA,MUQR.Median_SE=NA,MUQR.Median_lower=NA,
                 MUQR.Median_upper=NA,MUQR.Median_tval=NA,
                 MUQR.Median_p.value=NA,MUQR.Median_I2=NA,MUQR.Median_QEp=NA)
    Note <- c(Note,"MUQR.Median Meta-Analysis Failed")
  } else {
    Results <- c(Results,MUQR.Median_Beta=meta$b,MUQR.Median_SE=meta$se,
                 MUQR.Median_lower=meta$ci.lb,MUQR.Median_upper=meta$ci.ub,
                 MUQR.Median_tval=meta$zval,MUQR.Median_p.value=meta$pval,
                 MUQR.Median_I2=meta$I2,MUQR.Median_QEp=meta$QEp)
  }

  # META-ANALYSIS OF rt_MUQR.Median ESTIMATES
  meta <- tryCatch({rma(yi=rt_MUQR.Median_Beta, sei=rt_MUQR.Median_SE,
                        data=datatable[complete.cases(rt_MUQR.Median_Beta,rt_MUQR.Median_SE)],
                        method="REML",control=list(maxiter=1000, stepadj=.5))},
                   error=function(e) e)
  if(!class(meta)[1]=="rma.uni"){
    Results <- c(Results,rt_MUQR.Median_Beta=NA,rt_MUQR.Median_SE=NA,rt_MUQR.Median_lower=NA,
                 rt_MUQR.Median_upper=NA,rt_MUQR.Median_tval=NA,
                 rt_MUQR.Median_p.value=NA,rt_MUQR.Median_I2=NA,rt_MUQR.Median_QEp=NA)
    Note <- c(Note,"rt_MUQR.Median Meta-Analysis Failed")
  } else {
    Results <- c(Results,rt_MUQR.Median_Beta=meta$b,rt_MUQR.Median_SE=meta$se,
                 rt_MUQR.Median_lower=meta$ci.lb,rt_MUQR.Median_upper=meta$ci.ub,
                 rt_MUQR.Median_tval=meta$zval,rt_MUQR.Median_p.value=meta$pval,
                 rt_MUQR.Median_I2=meta$I2,rt_MUQR.Median_QEp=meta$QEp)
  }

  # META-ANALYSIS OF sc_MUQR.Median ESTIMATES
  meta <- tryCatch({rma(yi=sc_MUQR.Median_Beta, sei=sc_MUQR.Median_SE,
                        data=datatable[complete.cases(sc_MUQR.Median_Beta,sc_MUQR.Median_SE)],
                        method="REML",control=list(maxiter=1000, stepadj=.5))},
                   error=function(e) e)
  if(!class(meta)[1]=="rma.uni"){
    Results <- c(Results,sc_MUQR.Median_Beta=NA,sc_MUQR.Median_SE=NA,sc_MUQR.Median_lower=NA,
                 sc_MUQR.Median_upper=NA,sc_MUQR.Median_tval=NA,
                 sc_MUQR.Median_p.value=NA,sc_MUQR.Median_I2=NA,sc_MUQR.Median_QEp=NA)
    Note <- c(Note,"sc_MUQR.Median Meta-Analysis Failed")
  } else {
    Results <- c(Results,sc_MUQR.Median_Beta=meta$b,sc_MUQR.Median_SE=meta$se,
                 sc_MUQR.Median_lower=meta$ci.lb,sc_MUQR.Median_upper=meta$ci.ub,
                 sc_MUQR.Median_tval=meta$zval,sc_MUQR.Median_p.value=meta$pval,
                 sc_MUQR.Median_I2=meta$I2,sc_MUQR.Median_QEp=meta$QEp)
  }
  # META-ANALYSIS OF log_MUQR.Median ESTIMATES
  meta <- tryCatch({rma(yi=log_MUQR.Median_Beta, sei=log_MUQR.Median_SE,
                        data=datatable[complete.cases(log_MUQR.Median_Beta,log_MUQR.Median_SE)],
                        method="REML",control=list(maxiter=1000, stepadj=.5))},
                   error=function(e) e)
  if(!class(meta)[1]=="rma.uni"){
    Results <- c(Results,log_MUQR.Median_Beta=NA,log_MUQR.Median_SE=NA,log_MUQR.Median_lower=NA,
                 log_MUQR.Median_upper=NA,log_MUQR.Median_tval=NA,
                 log_MUQR.Median_p.value=NA,log_MUQR.Median_I2=NA,log_MUQR.Median_QEp=NA)
    Note <- c(Note,"log_MUQR.Median Meta-Analysis Failed")
  } else {
    Results <- c(Results,log_MUQR.Median_Beta=meta$b,log_MUQR.Median_SE=meta$se,
                 log_MUQR.Median_lower=meta$ci.lb,log_MUQR.Median_upper=meta$ci.ub,
                 log_MUQR.Median_tval=meta$zval,log_MUQR.Median_p.value=meta$pval,
                 log_MUQR.Median_I2=meta$I2,log_MUQR.Median_QEp=meta$QEp)
  }

  # META-ANALYSIS OF abs_MUQR.Median ESTIMATES
  meta <- tryCatch({rma(yi=abs_MUQR.Median_Beta, sei=MUQR.Median_SE,
                        data=datatable[complete.cases(abs_MUQR.Median_Beta,MUQR.Median_SE)],
                        method="REML",control=list(maxiter=1000, stepadj=.5))},
                   error=function(e) e)
  if(!class(meta)[1]=="rma.uni"){
    Results <- c(Results,abs_MUQR.Median_Beta=NA,abs_MUQR.Median_SE=NA,abs_MUQR.Median_lower=NA,
                 abs_MUQR.Median_upper=NA,abs_MUQR.Median_tval=NA,
                 abs_MUQR.Median_p.value=NA,abs_MUQR.Median_I2=NA,abs_MUQR.Median_QEp=NA)
    Note <- c(Note,"abs_MUQR.Median Meta-Analysis Failed")
  } else {
    Results <- c(Results,abs_MUQR.Median_Beta=meta$b,abs_MUQR.Median_SE=meta$se,
                 abs_MUQR.Median_lower=meta$ci.lb,abs_MUQR.Median_upper=meta$ci.ub,
                 abs_MUQR.Median_tval=meta$zval,abs_MUQR.Median_p.value=meta$pval,
                 abs_MUQR.Median_I2=meta$I2,abs_MUQR.Median_QEp=meta$QEp)
  }

  # META-ANALYSIS OF MUQR.MetaTau ESTIMATES
  meta <- tryCatch({rma(yi=MUQR.MetaTau_Beta, sei=MUQR.MetaTau_SE,
                        data=datatable[complete.cases(MUQR.MetaTau_Beta,MUQR.MetaTau_SE)],
                        method="REML",control=list(maxiter=1000, stepadj=.5))},
                   error=function(e) e)
  if(!class(meta)[1]=="rma.uni"){
    Results <- c(Results,MUQR.MetaTau_Beta=NA,MUQR.MetaTau_SE=NA,MUQR.MetaTau_lower=NA,
                 MUQR.MetaTau_upper=NA,MUQR.MetaTau_tval=NA,
                 MUQR.MetaTau_p.value=NA,MUQR.MetaTau_I2=NA,MUQR.MetaTau_QEp=NA)
    Note <- c(Note,"MUQR.MetaTau Meta-Analysis Failed")
  } else {
    Results <- c(Results,MUQR.MetaTau_Beta=meta$b,MUQR.MetaTau_SE=meta$se,
                 MUQR.MetaTau_lower=meta$ci.lb,MUQR.MetaTau_upper=meta$ci.ub,
                 MUQR.MetaTau_tval=meta$zval,MUQR.MetaTau_p.value=meta$pval,
                 MUQR.MetaTau_I2=meta$I2,MUQR.MetaTau_QEp=meta$QEp)
  }

  # META-ANALYSIS OF rt_MUQR.MetaTau ESTIMATES
  meta <- tryCatch({rma(yi=rt_MUQR.MetaTau_Beta, sei=rt_MUQR.MetaTau_SE,
                        data=datatable[complete.cases(rt_MUQR.MetaTau_Beta,rt_MUQR.MetaTau_SE)],
                        method="REML",control=list(maxiter=1000, stepadj=.5))},
                   error=function(e) e)
  if(!class(meta)[1]=="rma.uni"){
    Results <- c(Results,rt_MUQR.MetaTau_Beta=NA,rt_MUQR.MetaTau_SE=NA,rt_MUQR.MetaTau_lower=NA,
                 rt_MUQR.MetaTau_upper=NA,rt_MUQR.MetaTau_tval=NA,
                 rt_MUQR.MetaTau_p.value=NA,rt_MUQR.MetaTau_I2=NA,rt_MUQR.MetaTau_QEp=NA)
    Note <- c(Note,"rt_MUQR.MetaTau Meta-Analysis Failed")
  } else {
    Results <- c(Results,rt_MUQR.MetaTau_Beta=meta$b,rt_MUQR.MetaTau_SE=meta$se,
                 rt_MUQR.MetaTau_lower=meta$ci.lb,rt_MUQR.MetaTau_upper=meta$ci.ub,
                 rt_MUQR.MetaTau_tval=meta$zval,rt_MUQR.MetaTau_p.value=meta$pval,
                 rt_MUQR.MetaTau_I2=meta$I2,rt_MUQR.MetaTau_QEp=meta$QEp)
  }

  # META-ANALYSIS OF sc_MUQR.MetaTau ESTIMATES
  meta <- tryCatch({rma(yi=sc_MUQR.MetaTau_Beta, sei=sc_MUQR.MetaTau_SE,
                        data=datatable[complete.cases(sc_MUQR.MetaTau_Beta,sc_MUQR.MetaTau_SE)],
                        method="REML",control=list(maxiter=1000, stepadj=.5))},
                   error=function(e) e)
  if(!class(meta)[1]=="rma.uni"){
    Results <- c(Results,sc_MUQR.MetaTau_Beta=NA,sc_MUQR.MetaTau_SE=NA,sc_MUQR.MetaTau_lower=NA,
                 sc_MUQR.MetaTau_upper=NA,sc_MUQR.MetaTau_tval=NA,
                 sc_MUQR.MetaTau_p.value=NA,sc_MUQR.MetaTau_I2=NA,sc_MUQR.MetaTau_QEp=NA)
    Note <- c(Note,"sc_MUQR.MetaTau Meta-Analysis Failed")
  } else {
    Results <- c(Results,sc_MUQR.MetaTau_Beta=meta$b,sc_MUQR.MetaTau_SE=meta$se,
                 sc_MUQR.MetaTau_lower=meta$ci.lb,sc_MUQR.MetaTau_upper=meta$ci.ub,
                 sc_MUQR.MetaTau_tval=meta$zval,sc_MUQR.MetaTau_p.value=meta$pval,
                 sc_MUQR.MetaTau_I2=meta$I2,sc_MUQR.MetaTau_QEp=meta$QEp)
  }

  # META-ANALYSIS OF log_MUQR.MetaTau ESTIMATES
  meta <- tryCatch({rma(yi=log_MUQR.MetaTau_Beta, sei=log_MUQR.MetaTau_SE,
                        data=datatable[complete.cases(log_MUQR.MetaTau_Beta,log_MUQR.MetaTau_SE)],
                        method="REML",control=list(maxiter=1000, stepadj=.5))},
                   error=function(e) e)
  if(!class(meta)[1]=="rma.uni"){
    Results <- c(Results,log_MUQR.MetaTau_Beta=NA,log_MUQR.MetaTau_SE=NA,log_MUQR.MetaTau_lower=NA,
                 log_MUQR.MetaTau_upper=NA,log_MUQR.MetaTau_tval=NA,
                 log_MUQR.MetaTau_p.value=NA,log_MUQR.MetaTau_I2=NA,log_MUQR.MetaTau_QEp=NA)
    Note <- c(Note,"log_MUQR.MetaTau Meta-Analysis Failed")
  } else {
    Results <- c(Results,log_MUQR.MetaTau_Beta=meta$b,log_MUQR.MetaTau_SE=meta$se,
                 log_MUQR.MetaTau_lower=meta$ci.lb,log_MUQR.MetaTau_upper=meta$ci.ub,
                 log_MUQR.MetaTau_tval=meta$zval,log_MUQR.MetaTau_p.value=meta$pval,
                 log_MUQR.MetaTau_I2=meta$I2,log_MUQR.MetaTau_QEp=meta$QEp)
  }

  # META-ANALYSIS OF abs_MUQR.MetaTau ESTIMATES
  meta <- tryCatch({rma(yi=abs_MUQR.MetaTau_Beta, sei=MUQR.MetaTau_SE,
                        data=datatable[complete.cases(abs_MUQR.MetaTau_Beta,MUQR.MetaTau_SE)],
                        method="REML",control=list(maxiter=1000, stepadj=.5))},
                   error=function(e) e)
  if(!class(meta)[1]=="rma.uni"){
    Results <- c(Results,abs_MUQR.MetaTau_Beta=NA,abs_MUQR.MetaTau_SE=NA,abs_MUQR.MetaTau_lower=NA,
                 abs_MUQR.MetaTau_upper=NA,abs_MUQR.MetaTau_tval=NA,
                 abs_MUQR.MetaTau_p.value=NA,abs_MUQR.MetaTau_I2=NA,abs_MUQR.MetaTau_QEp=NA)
    Note <- c(Note,"abs_MUQR.MetaTau Meta-Analysis Failed")
  } else {
    Results <- c(Results,abs_MUQR.MetaTau_Beta=meta$b,abs_MUQR.MetaTau_SE=meta$se,
                 abs_MUQR.MetaTau_lower=meta$ci.lb,abs_MUQR.MetaTau_upper=meta$ci.ub,
                 abs_MUQR.MetaTau_tval=meta$zval,abs_MUQR.MetaTau_p.value=meta$pval,
                 abs_MUQR.MetaTau_I2=meta$I2,abs_MUQR.MetaTau_QEp=meta$QEp)
  }

  # META-ANALYSIS OF MUQR.MetaTau p.values
  # meta <- fun1(x=datatable[complete.cases(MUQR.MetaTau_tval)][,MUQR.MetaTau_tval],
  #              Statistic="z",Weights=datatable[complete.cases(MUQR.MetaTau_tval)][,N],
  #              log_p.value=log_p.value)
  # Results <- c(Results,Stouf_MUQR.MetaTau_zval=meta[1],
  #              Stouf_MUQR.MetaTau_p.value=meta[2])
  # meta <- fun2(x=datatable[complete.cases(MUQR.MetaTau_tval)][,MUQR.MetaTau_tval],
  #              Statistic="z",log_p.value=log_p.value)
  # Results <- c(Results,Fisch_MUQR.MetaTau_Tval=meta[1],
  #              Fisch_MUQR.MetaTau_p.value=meta[2])

  # META-ANALYSIS OF LEVENE p.values
  meta <- fun1(x=datatable[complete.cases(Levene_p.value)][,Levene_p.value],
               Statistic="p",Weights=datatable[complete.cases(Levene_p.value)][,N],
               log_p.value=log_p.value)
  Results <- c(Results,Stouf_Levene_zval=meta[1],Stouf_Levene_p.value=meta[2])
  meta <- fun2(x=datatable[complete.cases(Levene_p.value)][,Levene_p.value],
               Statistic="p",log_p.value=log_p.value)
  Results <- c(Results,Fisch_Levene_Tval=meta[1],Fisch_Levene_p.value=meta[2])

  # META-ANALYSIS OF z2 ESTIMATES
  meta <- tryCatch({rma(yi=z2_Beta, sei=z2_SE,
                        data=datatable[complete.cases(z2_Beta,z2_SE)],method="REML",
                        control=list(maxiter=1000, stepadj=.5))},
                   error=function(e) e)
  if(!class(meta)[1]=="rma.uni"){
    Results <- c(Results,z2_Beta=NA,z2_SE=NA,z2_lower=NA,z2_upper=NA,z2_tval=NA,
                 z2_p.value=NA,z2_I2=NA,z2_QEp=NA)
    Note <- c(Note,"z2 Meta-Analysis Failed")
  } else {
    Results <- c(Results,z2_Beta=meta$b,z2_SE=meta$se,z2_lower=meta$ci.lb,
                 z2_upper=meta$ci.ub,z2_tval=meta$zval,z2_p.value=meta$pval,
                 z2_I2=meta$I2,z2_QEp=meta$QEp)
  }
  if(is.null(Note)){
    Results <- c(Results,Notes=NA)
  } else {
    Note <- paste(Note, sep="", collapse=", ")
    Results <- c(Results,Notes=Note)
  }
  return(Results)
}

# NEED A GENERIC FUNCTION



