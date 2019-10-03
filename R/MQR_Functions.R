
# ============= RECQUIRED PACKAGES ====================
require(data.table); require(quantreg); require(lubridate); require(foreach);
require(doParallel); require(MASS); library(metafor)

# ============= MQR FUNCTIONS ====================
MetaReg.MCQR.function <- function(Beta,boot.Bs,test="z"){
  to.rm <- which(colSums(is.na(boot.Bs))==nrow(boot.Bs)) # check for taus where bootstrapping failed.
  if(length(to.rm)>0){
    boot.Bs <- boot.Bs[,-to.rm] # remove boot.Bs columns where bootstrapping failed.
    Beta <- Beta[-to.rm]  # remove Betas where bootstrapping failed.
  }
  COV <- cov(boot.Bs)
  taus <- as.numeric(colnames(boot.Bs))
  m <- length(taus)
  A <- cbind(1, taus)
  k <- ncol(A)
  SigmaInv <- solve(COV)
  Beta.hat <- solve(t(A)%*%SigmaInv%*%A)%*%t(A)%*%SigmaInv%*%as.matrix(Beta)
  Beta.hat.var <- solve(t(A)%*%solve(COV)%*%A)
  Beta.hat.SE <- sqrt(diag(Beta.hat.var))
  tval <- Beta.hat/Beta.hat.SE
  #Joint_Chi <- pchisq(t(Beta.hat)%*%solve(Beta.hat.var)%*%Beta.hat, df=nrow(Beta.hat),
  #                    lower=FALSE)
  if(test=="z"){
    pval <- sapply(tval, function(x) 2*pnorm(abs(x), lower=FALSE))
  } else if(test=="t"){
    pval <- sapply(tval, function(x) 2*pt(abs(x), df=m-k, lower=FALSE))
  }
  #Result <- cbind(Beta.hat, Beta.hat.SE, tval, pval,rep(Joint_Chi, nrow(Beta.hat)))
  #dimnames(Result) <- list(c("(Intercept)","taus"),
  #                         c("MetaReg_Beta","MetaReg_SE","MetaReg_tval","MetaReg_p.value",
  #                           "MetaReg_Joint"))
  Result <- cbind(Beta.hat, Beta.hat.SE, tval, pval)
  dimnames(Result) <- list(c("(Intercept)","taus"),
                           c("MetaReg_Beta","MetaReg_SE","MetaReg_tval",
                             "MetaReg_p.value"))
  #rm(m,A,SigmaInv,Beta.hat.var,Beta.hat.SE,tval,pval)
  return(Result)
}

MCQR.function <- function(datatable, tau=seq(0.05, 0.95, by=0.05),y,g,covariates=NULL,
                          boot.m="mcmb",boot.R=200,seed=31371,enable.dither=TRUE,
                          cfun1=MetaReg.MCQR.function){
  # datatable=Batch.Data;tau=Taus;y=outcome;g=z.SNPs[i];covariates=covariates;seed=i.s
  # datatable=z.Data;tau=Taus;y=outcome;g="MCQR_geno";
  # covariates=NULL;seed=Seed
  # boot.m="mcmb";boot.R=200
  # cfun1=MetaReg.MCQR.function
  # datatable=DT;tau=Taus;y="pheno"; g="geno"
  g.ptm <- proc.time()
  if(is.null(covariates)){
    FML  <- as.formula(paste(y,g,sep=" ~ "))
  } else {
    covs <- paste(covariates,sep="",collapse=" + ")
    FML  <- as.formula(paste(paste(y,g,sep=" ~ "),covs,sep=" + "))
  }
  DT <- data.table(model.frame(FML, datatable))

  ifelse(nrow(DT)<2000,{rq_meth <- "br"},
         ifelse(nrow(DT)<100000,{rq_meth <- "fn"},
                {rq_meth <- "pfn"}))

  # CALCULATE QUANTILE ESTIMATES
  boot.Bs <- matrix(NA,nrow=boot.R,ncol=length(tau),dimnames=list(1:boot.R,tau))
  Beta <- Notes <- rep(NA,length(tau))

  # BOOTSTRAP ERROR ESTIMATES UISNG boot.m
  for(i in 1:length(tau)){
    Q.mod <- tryCatch({rq(FML, tau[i], DT, method=rq_meth)}, error=function(err){
      # print(paste0("rq model fitting at tau=",tau[i]," failed because : ",err))
      return(c("ERROR",paste0("rq(tau=",tau[i],") failed because : ",err[[1]])))})
    if(Q.mod[[1]]=="ERROR"){
      Notes[i] <- Q.mod[[2]]
      next
    }
    Beta[i] <- coef(Q.mod)[g]
    i.Bs  <- tryCatch({
      bs.M <- boot.m
      set.seed(seed)
      summary.rq(Q.mod, se="boot", covariance=TRUE, R=boot.R,bsmethod=bs.M)$B},
      error=function(err){
        return(c("ERROR",paste0("summary.rq(",tau[i],") ",boot.m, " boot failed because : ",
                                err[[1]])))})
    if(i.Bs[[1]]=="ERROR"){
      Notes[i] <- gsub("\n","",i.Bs[[2]],fixed=TRUE)
    } else {
      colnames(i.Bs) <- names(coef(Q.mod))
      boot.Bs[,i] <- i.Bs[,g]
    }
  }

  fail.taus <- which(colSums(is.na(boot.Bs))==nrow(boot.Bs))
  if(length(fail.taus)>0){
    successful.taus <- colnames(boot.Bs[,-fail.taus])
  } else {
    successful.taus <- colnames(boot.Bs)
  }
  # META-REGRESSION ANALYSIS
  Results <- tryCatch({cfun1(Beta,boot.Bs)}, error=function(err){
    # print(paste0("MR-CQR failed because : ",err))
    return(c("ERROR",paste0("MR-CQR failed because : ",err[[1]])))})
  # META-REGRESSION ANALYSIS CHECK : For some studies MR has singular matrix design, this is
  # fixed by using dither(), which is 'quantreg' version of jitter(), to introduce right-sided
  # noise to the response.
  if(Results[[1]]=="ERROR"){
    Notes <- c(Notes, Results[[2]])
    if(enable.dither){
      set.seed(seed)
      DT[,Dithered_Outcome:=dither(DT[[y]],type="right")]
      FML <- as.formula(paste("Dithered_Outcome",FML[[1]],FML[[3]]))

      # CALCULATE QUANTILE ESTIMATES
      boot.Bs <- matrix(NA,nrow=boot.R,ncol=length(tau),dimnames=list(1:boot.R,tau))
      Beta <- Notes.2 <- rep(NA,length(tau))

      # BOOTSTRAP ERROR ESTIMATES UISNG boot.m
      for(i in 1:length(tau)){
        Q.mod <- tryCatch({rq(FML, tau[i], DT, method=rq_meth)}, error=function(err){
          # print(paste0("rq model fitting at tau=",tau[i]," failed because : ",err))
          return(c("ERROR",paste0("rq(tau=",tau[i],") dithered failed because : ",err[[1]])))})
        if(Q.mod[[1]]=="ERROR"){
          Notes.2[i] <- Q.mod[[2]]
          next
        }
        Beta[i] <- coef(Q.mod)[g]
        i.Bs  <- tryCatch({
          bs.M <- boot.m
          set.seed(seed)
          summary.rq(Q.mod, se="boot", covariance=TRUE, R=boot.R,bsmethod=bs.M)$B},
          error=function(err){
            return(c("ERROR",paste0("summary.rq(",tau[i],") ",boot.m,
                                    " dithered boot failed because : ",err[[1]])))})
        if(i.Bs[[1]]=="ERROR"){
          Notes.2[i] <- i.Bs[[2]]
        } else {
          colnames(i.Bs) <- names(coef(Q.mod))
          boot.Bs[,i] <- i.Bs[,g]
        }
      }

      fail.taus <- which(colSums(is.na(boot.Bs))==nrow(boot.Bs))
      if(length(fail.taus)>0){
        successful.taus <- colnames(boot.Bs[,-fail.taus])
      } else {
        successful.taus <- colnames(boot.Bs)
      }
      # META-REGRESSION ANALYSIS
      Results <- tryCatch({cfun1(Beta,boot.Bs)}, error=function(err){
        # print(paste0("MR-CQR failed because : ",err))
        return(c("ERROR",paste0("MR-CQR dithered failed because : ",err[[1]])))})
      if(Results[[1]]=="ERROR"){
        Results <- c(MCQR.MetaTau_Beta=NA,MCQR.MetaTau_SE=NA,MCQR.MetaTau_tval=NA,
                     MCQR.MetaTau_p.value=NA)
        Notes <- c(Notes,Notes.2,Results[[2]],"MCQR failed")
        FML <- as.formula(paste(y,FML[[1]],FML[[3]]))
      } else {
        Notes <- c(Notes, Notes.2, "MCQR required dither to work")
        Results <- Results["taus",]
        names(Results) <- c("MCQR.MetaTau_Beta","MCQR.MetaTau_SE","MCQR.MetaTau_tval",
                            "MCQR.MetaTau_p.value")
      }
    } else {
      Notes <- c(Notes, "MCQR failed (dither not enabled)")
      Results <- c(MCQR.MetaTau_Beta=NA,MCQR.MetaTau_SE=NA,MCQR.MetaTau_tval=NA,
                   MCQR.MetaTau_p.value=NA)
    }
  } else {
    Results <- Results["taus",]
    names(Results) <- c("MCQR.MetaTau_Beta","MCQR.MetaTau_SE","MCQR.MetaTau_tval",
                        "MCQR.MetaTau_p.value")
  }
  g.etm <- proc.time()-g.ptm

  # Add MEAN Regression Estimates on Raw Scale.
  LN.sum <- coef(summary(lm(FML,DT)))
  LN.Results <- c(LN_Beta=LN.sum[g,1],LN_SE=LN.sum[g,2],
                  LN_tval=LN.sum[g,3],LN_p.value=LN.sum[g,4])

  # Add MEDIAN Regression Estimates
  Q.mod <- tryCatch({rq(FML, tau=0.5, DT,method=rq_meth)},error=function(err){
    return(c("ERROR",paste0("rq(Median) using ",rq_meth, " , failed because : ",err[[1]])))})
  if(Q.mod[[1]]=="ERROR"){
    Notes <- c(Notes,Q.mod[[2]])
    t.meth <- "br"
    Q.mod <- tryCatch({rq(FML, tau=0.5, DT,method=t.meth)},error=function(err){
      return(c("ERROR",paste0("rq(Median) using ",t.meth, " , failed because : ",err[[1]])))})
  }
  if(Q.mod[[1]]=="ERROR"){
    Notes <- c(Notes,Q.mod[[2]],"MCQR.Median fit failed")
    MEDIAN.Results <- c(MCQR.Median_Beta=NA,MCQR.Median_SE=NA,MCQR.Median_tval=NA,
                        MCQR.Median_p.value=NA)
  } else {
    Q.sum  <- tryCatch({
      bs.M <- boot.m
      set.seed(seed)
      summary.rq(Q.mod, se="boot", covariance=TRUE, R=boot.R, bsmethod=bs.M)},
      error=function(err){
        return(c("ERROR",paste0("summary.rq(Median) ",boot.m, " boot failed because : ",err[[1]])))})
    if(Q.sum[[1]]=="ERROR"){
      Notes <- c(Notes,Q.sum[[1]],"MCQR.Median boot failed")
      MEDIAN.Results <- c(MCQR.Median_Beta=NA,MCQR.Median_SE=NA,MCQR.Median_tval=NA,
                          MCQR.Median_p.value=NA)
    } else {
      boot.Bs <- Q.sum$B
      Med.Coef <- coefficients(Q.mod)
      Med.VNames <- names(Med.Coef)
      Med.P <- length(Med.Coef)
      Med.N <- length(Q.mod$resid)
      Med.Rdf <- Med.N - Med.P
      COV <- cov(boot.Bs)
      Med.Serr <- sqrt(diag(COV))
      Med.Coef <- array(Med.Coef, c(Med.P, 4))
      dimnames(Med.Coef) <- list(Med.VNames, c("Value", "Std. Error", "t value", "Pr(>|t|)"))
      Med.Coef[,2] <- Med.Serr
      Med.Coef[,3] <- Med.Coef[, 1]/Med.Coef[,2]
      Med.Coef[,4] <- if (Med.Rdf > 0)
        2 * (pt(abs(Med.Coef[, 3]), Med.Rdf,lower.tail=FALSE))
      MEDIAN.Results <- c(MCQR.Median_Beta=Med.Coef[g,1],MCQR.Median_SE=Med.Coef[g,2],
                          MCQR.Median_tval=Med.Coef[g,3],
                          MCQR.Median_p.value=Med.Coef[g,4])
    }
  }
  Notes <- Notes[!is.na(Notes)]
  if(length(Notes)>0){
    Notes <- paste(Notes,sep="",collapse=". ")
  } else {
    Notes <- NA
  }

  Results <- c(N=nrow(DT),EAF=sum(DT[[g]])/(2*nrow(DT)),LN.Results,MEDIAN.Results,
               Results,Successful.Taus=paste(successful.taus,sep="",collapse=","),
               No.Successful.Taus=length(successful.taus),rq.method=rq_meth,boot.method=bs.M,
               MCQR.TtC=g.etm[[3]],Notes=Notes)
  return(Results)
}

# Function to extract covariance matrix from multivariate OLS
cov.RIF <- function(object, pred){
  # object is class of multivariate lm
  # pred is a character for name of predictor
  #object=Models; pred="g"
  x <- object$qr$qr
  p <- object$rank
  m <- ncol(object$coefficients)

  xxinv <- diag(p)
  xxinv <- backsolve(qr(x)$qr[1:p, 1:p, drop = FALSE], xxinv)
  D.inv <- xxinv %*% t(xxinv)

  Sigma_hat <- var(object$residuals)
  COV <- kronecker(Sigma_hat,D.inv)

  vnames <- dimnames(x)[[2]]
  t <- which(vnames %in% pred)
  take <- seq(t,m*p,by=p)
  #rm(x, p, m, xxinv, D.inv, Sigma_hat,vnames, t)
  return(COV[take,take])
}

# First order RIF Transformation
RIF.Transformation <- function(y,taus){
  n <- length(y)
  qhat <- quantile(y,taus)
  bw0 <- bw.nrd0(y)
  Pf1 <- sapply(qhat,function(x) mean(dnorm((y-x),0,bw0)))
  RIF <- matrix(NA, ncol=length(taus),nrow=length(y))
  for(i in seq_along(taus)){
    dqdt <- sapply(y,function(x,tau_i,qhat_i,Pf1_i)
      (tau_i-!x>qhat_i)/Pf1_i,tau_i=taus[i], qhat_i=qhat[i], Pf1_i=Pf1[i])
    RIF[,i] <- qhat[i]+dqdt
  }
  return(RIF)
}

# Function to fit meta-regression models
MetaReg.MUQR.function <- function(Beta, COV, taus,test="z"){
  #taus=tau; test="z"
  m <- length(taus)
  A <- cbind(1, taus)
  k <- ncol(A)
  SigmaInv <- solve(COV)
  Beta.hat <- solve(t(A)%*%SigmaInv%*%A)%*%t(A)%*%SigmaInv%*%as.matrix(Beta)
  Beta.hat.var <- solve(t(A)%*%solve(COV)%*%A)
  Beta.hat.SE <- sqrt(diag(Beta.hat.var))
  tval <- Beta.hat/Beta.hat.SE
  #Joint_Chi <- pchisq(t(Beta.hat)%*%solve(Beta.hat.var)%*%Beta.hat, df=nrow(Beta.hat),
  #                    lower=FALSE)
  if(test=="z"){
    pval <- sapply(tval, function(x) 2*pnorm(abs(x), lower=FALSE))
  } else if(test=="t"){
    pval <- sapply(tval, function(x) 2*pt(abs(x), df=m-k, lower=FALSE))
  }
  #Result <- cbind(Beta.hat, Beta.hat.SE, tval, pval,rep(Joint_Chi, nrow(Beta.hat)))
  #dimnames(Result) <- list(c("(Intercept)","taus"),
  #                         c("MetaReg_Beta","MetaReg_SE","MetaReg_tval","MetaReg_p.value",
  #                           "MetaReg_Joint"))
  Result <- cbind(Beta.hat, Beta.hat.SE, tval, pval)
  dimnames(Result) <- list(c("(Intercept)","taus"),
                           c("MetaReg_Beta","MetaReg_SE","MetaReg_tval",
                             "MetaReg_p.value"))
  #rm(m,A,SigmaInv,Beta.hat.var,Beta.hat.SE,tval,pval)
  return(Result)
}

MUQR.function <- function(datatable, tau=seq(0.05, 0.95, by=0.05), y, g,covariates=NULL,
                          Univariable=TRUE,seed=31371,enable.dither=TRUE,
                          ufun1=RIF.Transformation,ufun2=cov.RIF,
                          ufun3=MetaReg.MUQR.function){
  #d atatable=Batch.Data; tau=Taus; y=outcome; g=l.SNPs[j];
  # datatable=Batch.Data;tau=Taus;y=outcome;g=z.SNPs[i];covariates=covariates;
  # Univariable=TRUE;
  # covariates=c(covariates,"STUDY");
  # ufun1=RIF.Transformation; ufun2=cov.RIF;ufun3=MetaReg.MUQR.function

  g.ptm <- proc.time()
  if(is.null(covariates)){
    FML  <- as.formula(paste(y,g,sep=" ~ "))
  } else {
    covs <- paste(covariates,sep="",collapse=" + ")
    FML  <- as.formula(paste(paste(y,g,sep=" ~ "),covs,sep=" + "))
  }
  DT <- data.table(model.frame(FML, datatable))
  DT[,y:=DT[[y]]]

  # Adjust for Median effects of g (i.e. residual scale)
  RIF <- ufun1(y=DT[,y], taus=0.5)
  Beta_0.5 <- coef(lm(as.formula(paste("RIF",g,sep=" ~ ")),data=DT))[g]
  DT[,y:= y - Beta_0.5*DT[[g]]]

  RIF <- ufun1(y=DT[,y], taus=tau) # compute RIF at taus of interest

  if(Univariable){
    FML <- as.formula(paste("RIF",g,sep=" ~ "))
  } else {
    FML <- as.formula(paste("RIF",as.character(FML)[1],as.character(FML)[3],
                            sep="",collapse=""))
  }

  Models <- lm(FML, DT)
  Beta <- Models$coefficients[2,]
  COV <- ufun2(Models, pred=g)

  # re-centred tau
  taus <- tau-0.5
  Results <- tryCatch({ufun3(Beta=Beta, COV=COV, taus=taus)}, error=function(err){
    # print(paste0("MR-UQR failed because : ",err))
    return(c("ERROR",paste0("MR-UQR failed because : ",err[[1]])))
  })
  # META-REGRESSION ANALYSIS CHECK : For some studies MR has singular matrix design, this is
  # fixed by using dither(), which is 'quantreg' version of jitter(), to introduce right-sided
  # noise to the response.
  if(Results[[1]]=="ERROR"){
    Notes <- Results[[2]]
    if(enable.dither){
      set.seed(seed)
      DT[,Dithered_Response:=dither(y,type="right")]
      # FML <- as.formula(paste("Dithered_Response",FML[[1]],FML[[3]]))
      # Adjust for Median effects of g (i.e. residual scale)
      RIF <- ufun1(y=DT[,Dithered_Response], taus=0.5)
      Beta_0.5 <- coef(lm(as.formula(paste("RIF",g,sep=" ~ ")),data=DT))[g]
      DT[,Residual_Dithered_Response:=Dithered_Response - Beta_0.5*DT[[g]]]

      RIF <- ufun1(y=DT[,Residual_Dithered_Response], taus=tau) # compute RIF at taus of interest
      Models <- lm(FML, DT)
      Beta <- Models$coefficients[2,]
      COV <- ufun2(Models, pred=g)
      Results <- tryCatch({ufun3(Beta=Beta, COV=COV, taus=taus)}, error=function(err){
        # print(paste0("MR-UQR failed because : ",err))
        return(c("ERROR",paste0("MR-UQR failed because : ",err[[1]])))})
      if(Results[[1]]=="ERROR"){
        Notes <- c(Notes,Results[[2]],"MUQR failed")
        Results <- c(MUQR.MetaTau_Beta=NA,MUQR.MetaTau_SE=NA,MUQR.MetaTau_tval=NA,
                     MUQR.MetaTau_p.value=NA)
      } else {
        Notes <- c(Notes, "MUQR required dither to work")
        Results <- Results["taus",]
        names(Results) <- c("MUQR.MetaTau_Beta","MUQR.MetaTau_SE","MUQR.MetaTau_tval",
                            "MUQR.MetaTau_p.value")
        # Add MEDIAN Regression Estimates
        RIF <- ufun1(y=DT[,Dithered_Response], taus=0.5)
        RIF.sum <- coef(summary(lm(FML,DT)))
        MEDIAN.Result <- c(MUQR.Median_Beta=RIF.sum[g,1],MUQR.Median_SE=RIF.sum[g,2],
                           MUQR.Median_tval=RIF.sum[g,3],MUQR.Median_p.value=RIF.sum[g,4])

        # Add MEAN Regression Estimates.
        FML <- as.formula(paste("Dithered_Response",FML[[1]],FML[[3]]))
        LN.sum <- coef(summary(lm(FML,DT)))
        LN.Results <- c(LN_Beta=LN.sum[g,1],LN_SE=LN.sum[g,2],
                        LN_tval=LN.sum[g,3],LN_p.value=LN.sum[g,4])
      }
    } else {
      Results <- c(MUQR.MetaTau_Beta=NA,MUQR.MetaTau_SE=NA,MUQR.MetaTau_tval=NA,
                   MUQR.MetaTau_p.value=NA)
      Notes <- c(Notes, "MUQR failed (dither not enabled)")
    }
  } else {
    Notes <- NA
    Results <- Results["taus",]
    names(Results) <- c("MUQR.MetaTau_Beta","MUQR.MetaTau_SE","MUQR.MetaTau_tval",
                        "MUQR.MetaTau_p.value")
    # Add MEDIAN Regression Estimates
    RIF <- ufun1(y=DT[[y]], taus=0.5)
    RIF.sum <- coef(summary(lm(FML,DT)))
    MEDIAN.Result <- c(MUQR.Median_Beta=RIF.sum[g,1],MUQR.Median_SE=RIF.sum[g,2],
                       MUQR.Median_tval=RIF.sum[g,3],MUQR.Median_p.value=RIF.sum[g,4])

    # Add MEAN Regression Estimates.
    FML <- as.formula(paste(y,FML[[1]],FML[[3]]))
    LN.sum <- coef(summary(lm(FML,DT)))
    LN.Results <- c(LN_Beta=LN.sum[g,1],LN_SE=LN.sum[g,2],
                    LN_tval=LN.sum[g,3],LN_p.value=LN.sum[g,4])
    # Result <- c(LN.Results,Result,MUQR.TtC=g.etm[[3]])
  }
  g.etm <- proc.time()-g.ptm
  if(length(Notes)>1){
    Notes <- paste(Notes,sep="",collapse=". ")
  }
  Results <- c(N=nrow(DT),EAF=sum(DT[[g]])/(2*nrow(DT)),LN.Results,MEDIAN.Result,Results,
               MUQR.TtC=g.etm[[3]],Notes=Notes)
  return(Results)
}

Scaling.Adjustor <- function(FML,data,weights,model,response,predictor,
                             method="Iteratively_ReWeighted"){
  # data=Results
  # browser()
  DT <- copy(data)
  if(missing(model)){
    DT[,w:=1/(DT[[weights]]^2)]
    if(method=="Iteratively_ReWeighted"){
      Model <- rlm(FML,data=DT,weights=w,method="M",wt.method="inv.var")
    } else if(method=="Median_rq"){
      Model <- rq(FML,tau=0.5,data=DT,weights=w,method="fn")
    } else if(method=="ols"){
      Model <- lm(FML,data=DT)
    }
    y.name <- colnames(Model$model)[1]
    y.se <- gsub("Beta","SE",y.name)
    x.name <- colnames(Model$model)[2]
    x.se <- gsub("Beta","SE",x.name)

    Pred <- Model$model[,x.name]*Model$coefficients[[2]] + Model$coefficients[[1]]
    Pred.se <- DT[[x.se]]*Model$coefficients[[2]] + Model$coefficients[[1]]
    Beta <- Model$residuals
    se <- sqrt(DT[[y.se]]^2 + Pred.se^2)
    tval <- Beta/se
    p.value <- 2*pnorm(abs(tval),lower=FALSE)
    Results <- cbind(Pred,Pred.se,Beta,se,tval,p.value)
    r.name <- gsub("Beta","",y.name)
    colnames(Results) <- c(paste0("Pred.",r.name,"Beta"),paste0("Pred.",r.name,"SE"),
                           paste0("Adj.",r.name,"Beta"),paste0("Adj.",r.name,"SE"),
                           paste0("Adj.",r.name,"tval"),paste0("Adj.",r.name,"p.value"))
    Results <- data.table(Results)
    # Model$model <- Model$residuals <- Model$weights <- NULL
    return(list(Model,Results))
  } else {
    Model <- copy(model)
    x.se <- gsub("Beta","SE",predictor)
    y.se <- gsub("Beta","SE",response)
    Pred <- DT[[predictor]]*Model$coefficients[[2]] + Model$coefficients[[1]]
    Pred.se <- DT[[x.se]]*Model$coefficients[[2]] + Model$coefficients[[1]]
    Beta <- DT[[response]]-Pred
    se <- sqrt(DT[[y.se]]^2 + Pred.se^2)
    tval <- Beta/se
    p.value <- 2*pnorm(abs(tval),lower=FALSE)
    Results <- cbind(Pred,Pred.se,Beta,se,tval,p.value)
    r.name <- gsub("Beta","",response)
    colnames(Results) <- c(paste0("Pred.",r.name,"Beta"),paste0("Pred.",r.name,"SE"),
                           paste0("Adj.",r.name,"Beta"),paste0("Adj.",r.name,"SE"),
                           paste0("Adj.",r.name,"tval"),paste0("Adj.",r.name,"p.value"))
    Results <- data.table(Results)
    return(Results)
  }
}

# # ============= VARIANCE HETEROGENEITY FUNCTIONS ====================
# # Function takes 2 methods. "Classic-Levene" and "Brown-Forsythe". The latter is more
# # robust to skewed distributions. The default is "classic.Levene". Note: Selecting
# # "Brown-Forsythe" gives the exact same result as levene.test function in lawstat package
# Levene.function <- function(datatable, y, g, covariates=NULL, method="Brown-Forsythe",
#                             log.p=FALSE){
#   #datatable=Data; method="Brown-Forsythe"
#   if(is.null(covariates)){
#     FML  <- as.formula(paste(y,g,sep=" ~ "))
#   } else {
#     covs <- paste(covariates,sep="",collapse=" + ")
#     FML  <- as.formula(paste(paste(y,g,sep=" ~ "),covs,sep=" + "))
#   }
#
#   DT <- data.table(model.frame(FML, datatable))
#   DT[,g:=ifelse(DT[[g]]<=0.5,0,
#                 ifelse(DT[[g]]<=1.5,1,
#                        ifelse(DT[[g]]>1.5,2, NA)))]
#   if(is.null(covariates)){
#     DT[,y:=DT[[y]]]
#   } else {
#     FML.r <- paste(y, covs, sep=" ~ ")
#     LMOD <- lm(FML,DT)
#     DT[,y:=resid(LMOD)]
#   }
#
#   DT <- DT[,.(y,g)]
#
#   # first calculate Zij (which is key) for each subject and add it to
#   # data.table.
#   # Use Median for "Brown–Forsythe" and Mean for "Classic-Levene". Otherwise tests
#   # are the same.
#   if(method=="Brown-Forsythe"){
#     DT[,Zij:=abs(y-median(y)),by=g]
#   } else if(method=="Classic-Levene"){
#     DT[,Zij:=abs(y-mean(y)),by=g]
#   }
#
#   # calculate Z..
#   Z.. <- mean(DT[,Zij])
#   # find the group number
#   k <- length(unique(DT[,g]))
#   # find the sample size
#   N <- nrow(DT)
#   # calculate Zij_Zi, which is the difference between Zij (for each
#   # subject) and the mean of Zij for each genotype, squared
#   DT[,Zij_Zi:=(Zij-mean(Zij))^2,by=g]
#   W <- (N-k)*sum(DT[,length(y)*(mean(Zij)-Z..)^2,by=g][,V1])/
#     ((k-1)*sum(DT[,Zij_Zi]))
#   if(log.p){
#     p <- pf(W,df1=k-1,df2=N-k, lower.tail=FALSE, log.p=TRUE)
#   } else {
#     p <- pf(W,df1=k-1,df2=N-k, lower.tail=FALSE)
#   }
#   Result <- c(Levene_Fval=W,Levene_df1=k-1,Levene_df2=N-k,Levene_p.value=p)
#   return(Result)
# }
#
# z2.function <- function(datatable, y, g,covariates=NULL){
#   #datatable=Data; tau=Taus; y="BMI"; g=SNPs[j]; Adjusted.for=NA
#   #covariates=c(covariates,"STUDY");
#   #fun1=RIF.Transformation; fun2=cov.RIF;fun3=MetaReg.function
#
#   if(is.null(covariates)){
#     FML  <- as.formula(paste(y,g,sep=" ~ "))
#   } else {
#     covs <- paste(covariates,sep="",collapse=" + ")
#     FML  <- as.formula(paste(paste(y,g,sep=" ~ "),covs,sep=" + "))
#   }
#   DT <- data.table(model.frame(FML, datatable))
#
#   DT[,z2:=(scale(DT[[y]]))^2]
#   FML <- as.formula(paste("z2",as.character(FML)[1],as.character(FML)[3],
#                           sep="",collapse=""))
#
#   # Add z2 Regression Estimates
#   LN.sum <- coef(summary(lm(FML,DT)))
#   Result <- c(z2_Beta=LN.sum[g,1],z2_SE=LN.sum[g,2],
#               z2_tval=LN.sum[g,3],z2_p.value=LN.sum[g,4])
#   return(Result)
# }
#
#
# # ============= META-ANALYSIS FUNCTIONS ====================
# original.Stouffer.test <- function(p,Weights,log_p.value=FALSE) {
#   #p is a vector of p-values
#   #Weights is typically sample size
#   #p=datatable[,RIF.MetaTau_p.value]; Weights=datatable[,N];log_p.value=FALSE
#   #Weights=NA
#   if(missing(Weights)) {
#     w <- rep(1, length(p))/length(p)
#   } else if(length(Weights)==length(p)){
#     w <- sqrt(Weights)/sum(sqrt(Weights))
#   } else {
#     stop("Length of p and w must equal!")
#   }
#   Zi <- qnorm(p/2,lower=FALSE) # note /2 is modification for two-sided pvalues
#   Z  <- sum(w*Zi)/sqrt(sum(w^2))
#   if(log_p.value){
#     p.val <- -1*pnorm(Z, lower.tail=FALSE, log=TRUE)
#   } else{
#     p.val <- pnorm(Z, lower.tail=FALSE)
#   }
#   return(c(Z,p.val))
#   #return(as.vector(p.val))
# }
# original.Fischer.test <- function(p,log_p.value=FALSE) {
#   #p is a vector of p-values
#   #Weights is typically sample size
#   #p=datatable[,RIF.MetaTau_p.value]
#   T <- sum(-2*log(p))
#   if(log_p.value){
#     p.val <- pchisq(T, df=length(p)*2, lower.tail=FALSE, log=TRUE)
#   } else{
#     p.val <- pchisq(T, df=length(p)*2, lower.tail=FALSE)
#   }
#   return(c(T,p.val))
#   #return(as.vector(p.val))
# }
# Stouffer.test <- function(x,Statistic="z",Weights,log_p.value=FALSE) {
#   #z is a vector of z-values
#   #Weights is typically sample size
#   #z=datatable[,RIF.MetaTau_tval]; Weights=datatable[,N];log_p.value=FALSE
#   #Weights=NA
#   if(missing(Weights)) {
#     w <- rep(1, length(x))/length(x)
#   } else if(length(Weights)==length(x)){
#     w <- sqrt(Weights)/sum(sqrt(Weights))
#   } else {
#     stop("Length of x and w must equal!")
#   }
#   if(Statistic=="z"){
#     Zi <- abs(x)
#     Z  <- sum(w*Zi)/sqrt(sum(w^2))
#   } else if(Statistic=="p"){
#     Zi <- qnorm(x/2,lower=FALSE) # note /2 is modification for two-sided pvalues
#     Z  <- sum(w*Zi)/sqrt(sum(w^2))
#   }
#   if(log_p.value){
#     p.val <- -1*pnorm(Z, lower.tail=FALSE, log=TRUE)
#   } else{
#     p.val <- pnorm(Z, lower.tail=FALSE)
#   }
#   return(c(Z,p.val))
#   #return(as.vector(p.val))
# }
# Fischer.test <- function(x,Statistic="z",log_p.value=FALSE) {
#   #z is a vector of z-values
#   #Weights is typically sample size
#   #z=datatable[,RIF.MetaTau_tval]
#   if(Statistic=="z"){
#     T <- sum(-2*pnorm(x,lower.tail=TRUE,log.p=TRUE))
#   } else if(Statistic=="p"){
#     T <- sum(-2*log(x))
#   }
#   if(log_p.value){
#     p.val <- pchisq(T, df=length(x)*2, lower.tail=FALSE, log=TRUE)
#   } else{
#     p.val <- pchisq(T, df=length(x)*2, lower.tail=FALSE)
#   }
#   return(c(T,p.val))
#   #return(as.vector(p.val))
# }
#
# MCQR.metafor.function <- function(datatable,fun1=Stouffer.test,fun2=Fischer.test,
#                                   log_p.value=FALSE){
#   #datatable=i.Results[SNP==i.SNPs[l]];fun1=Stouffer.test;fun2=Fischer.test;log_p.value=FALSE
#   #colnames(datatable)
#   # SUMMARIZE SNP INFO
#   Results <- c(SNP=unique(datatable[,SNP]),
#                Chr=unique(datatable[complete.cases(Chr)][,Chr]),
#                N=sum(datatable[,N],na.rm=TRUE),
#                avr.EAF=mean(datatable[,EAF],na.rm=TRUE),
#                min.EAF=min(datatable[,EAF],na.rm=TRUE),
#                max.EAF=max(datatable[,EAF],na.rm=TRUE))
#   Note <- NULL
#   # META-ANALYSIS OF MCQR.Median ESTIMATES
#   meta <- tryCatch({rma(yi=MCQR.Median_Beta, sei=MCQR.Median_SE,
#                         data=datatable[complete.cases(MCQR.Median_Beta,MCQR.Median_SE)],
#                         method="REML",control=list(maxiter=1000, stepadj=.5))},
#                    error=function(e) e)
#   if(!class(meta)[1]=="rma.uni"){
#     Results <- c(Results,MCQR.Median_Beta=NA,MCQR.Median_SE=NA,MCQR.Median_lower=NA,
#                  MCQR.Median_upper=NA,MCQR.Median_tval=NA,
#                  MCQR.Median_p.value=NA,MCQR.Median_I2=NA,MCQR.Median_QEp=NA)
#     Note <- c(Note,"MCQR.Median Meta-Analysis Failed")
#   } else {
#     Results <- c(Results,MCQR.Median_Beta=meta$b,MCQR.Median_SE=meta$se,
#                  MCQR.Median_lower=meta$ci.lb,MCQR.Median_upper=meta$ci.ub,
#                  MCQR.Median_tval=meta$zval,MCQR.Median_p.value=meta$pval,
#                  MCQR.Median_I2=meta$I2,MCQR.Median_QEp=meta$QEp)
#   }
#   # META-ANALYSIS OF abs_MCQR.Median ESTIMATES
#   meta <- tryCatch({rma(yi=abs_MCQR.Median_Beta, sei=MCQR.Median_SE,
#                         data=datatable[complete.cases(abs_MCQR.Median_Beta,MCQR.Median_SE)],
#                         method="REML",control=list(maxiter=1000, stepadj=.5))},
#                    error=function(e) e)
#   if(!class(meta)[1]=="rma.uni"){
#     Results <- c(Results,abs_MCQR.Median_Beta=NA,abs_MCQR.Median_SE=NA,abs_MCQR.Median_lower=NA,
#                  abs_MCQR.Median_upper=NA,abs_MCQR.Median_tval=NA,
#                  abs_MCQR.Median_p.value=NA,abs_MCQR.Median_I2=NA,abs_MCQR.Median_QEp=NA)
#     Note <- c(Note,"abs_MCQR.Median Meta-Analysis Failed")
#   } else {
#     Results <- c(Results,abs_MCQR.Median_Beta=meta$b,abs_MCQR.Median_SE=meta$se,
#                  abs_MCQR.Median_lower=meta$ci.lb,abs_MCQR.Median_upper=meta$ci.ub,
#                  abs_MCQR.Median_tval=meta$zval,abs_MCQR.Median_p.value=meta$pval,
#                  abs_MCQR.Median_I2=meta$I2,abs_MCQR.Median_QEp=meta$QEp)
#   }
#
#   # META-ANALYSIS OF MCQR.MetaTau ESTIMATES
#   meta <- tryCatch({rma(yi=MCQR.MetaTau_Beta, sei=MCQR.MetaTau_SE,
#                         data=datatable[complete.cases(MCQR.MetaTau_Beta,MCQR.MetaTau_SE)],
#                         method="REML",control=list(maxiter=1000, stepadj=.5))},
#                    error=function(e) e)
#   if(!class(meta)[1]=="rma.uni"){
#     Results <- c(Results,MCQR.MetaTau_Beta=NA,MCQR.MetaTau_SE=NA,MCQR.MetaTau_lower=NA,
#                  MCQR.MetaTau_upper=NA,MCQR.MetaTau_tval=NA,
#                  MCQR.MetaTau_p.value=NA,MCQR.MetaTau_I2=NA,MCQR.MetaTau_QEp=NA)
#     Note <- c(Note,"MCQR.MetaTau Meta-Analysis Failed")
#   } else {
#     Results <- c(Results,MCQR.MetaTau_Beta=meta$b,MCQR.MetaTau_SE=meta$se,
#                  MCQR.MetaTau_lower=meta$ci.lb,MCQR.MetaTau_upper=meta$ci.ub,
#                  MCQR.MetaTau_tval=meta$zval,MCQR.MetaTau_p.value=meta$pval,
#                  MCQR.MetaTau_I2=meta$I2,MCQR.MetaTau_QEp=meta$QEp)
#   }
#
#   # META-ANALYSIS OF abs_MCQR.MetaTau ESTIMATES
#   meta <- tryCatch({rma(yi=abs_MCQR.MetaTau_Beta, sei=MCQR.MetaTau_SE,
#                         data=datatable[complete.cases(abs_MCQR.MetaTau_Beta,MCQR.MetaTau_SE)],
#                         method="REML",control=list(maxiter=1000, stepadj=.5))},
#                    error=function(e) e)
#   if(!class(meta)[1]=="rma.uni"){
#     Results <- c(Results,abs_MCQR.MetaTau_Beta=NA,abs_MCQR.MetaTau_SE=NA,abs_MCQR.MetaTau_lower=NA,
#                  abs_MCQR.MetaTau_upper=NA,abs_MCQR.MetaTau_tval=NA,
#                  abs_MCQR.MetaTau_p.value=NA,abs_MCQR.MetaTau_I2=NA,abs_MCQR.MetaTau_QEp=NA)
#     Note <- c(Note,"abs_MCQR.MetaTau Meta-Analysis Failed")
#   } else {
#     Results <- c(Results,abs_MCQR.MetaTau_Beta=meta$b,abs_MCQR.MetaTau_SE=meta$se,
#                  abs_MCQR.MetaTau_lower=meta$ci.lb,abs_MCQR.MetaTau_upper=meta$ci.ub,
#                  abs_MCQR.MetaTau_tval=meta$zval,abs_MCQR.MetaTau_p.value=meta$pval,
#                  abs_MCQR.MetaTau_I2=meta$I2,abs_MCQR.MetaTau_QEp=meta$QEp)
#   }
#
#   if(is.null(Note)){
#     Results <- c(Results,Notes=NA)
#   } else {
#     Note <- paste(Note, sep="", collapse=", ")
#     Results <- c(Results,Notes=Note)
#   }
#   return(Results)
# }
# MUQR.metafor.function <- function(datatable,fun1=Stouffer.test,fun2=Fischer.test,
#                                   log_p.value=FALSE){
#   # datatable=i.Results[SNP==i.SNPs[l]];fun1=Stouffer.test;fun2=Fischer.test;log_p.value=FALSE
#   # colnames(datatable)
#
#   # SUMMARIZE SNP INFO
#   Results <- c(SNP=unique(datatable[,SNP]),
#                Chr=unique(datatable[complete.cases(Chr)][,Chr]),
#                N=sum(datatable[,N],na.rm=TRUE),
#                avr.EAF=mean(datatable[,EAF],na.rm=TRUE),
#                min.EAF=min(datatable[,EAF],na.rm=TRUE),
#                max.EAF=max(datatable[,EAF],na.rm=TRUE))
#   Note <- NULL
#   # META-ANALYSIS OF LN ESTIMATES
#   meta <- tryCatch({rma(yi=LN_Beta, sei=LN_SE,
#                         data=datatable[complete.cases(LN_Beta,LN_SE)],method="REML",
#                         control=list(maxiter=1000, stepadj=.5))},
#                    error=function(e) e)
#   if(!class(meta)[1]=="rma.uni"){
#     Results <- c(Results,LN_Beta=NA,LN_SE=NA,LN_lower=NA,LN_upper=NA,LN_tval=NA,
#                  LN_p.value=NA,LN_I2=NA,LN_QEp=NA)
#     Note <- c(Note,"LN Meta-Analysis Failed")
#   } else {
#     Results <- c(Results,LN_Beta=meta$b,LN_SE=meta$se,LN_lower=meta$ci.lb,
#                  LN_upper=meta$ci.ub,LN_tval=meta$zval,LN_p.value=meta$pval,
#                  LN_I2=meta$I2,LN_QEp=meta$QEp)
#   }
#
#   # META-ANALYSIS OF rt_LN ESTIMATES
#   meta <- tryCatch({rma(yi=rt_LN_Beta, sei=rt_LN_SE,
#                         data=datatable[complete.cases(rt_LN_Beta,rt_LN_SE)],method="REML",
#                         control=list(maxiter=1000, stepadj=.5))},
#                    error=function(e) e)
#   if(!class(meta)[1]=="rma.uni"){
#     Results <- c(Results,rt_LN_Beta=NA,rt_LN_SE=NA,rt_LN_lower=NA,rt_LN_upper=NA,rt_LN_tval=NA,
#                  rt_LN_p.value=NA,rt_LN_I2=NA,rt_LN_QEp=NA)
#     Note <- c(Note,"rt_LN Meta-Analysis Failed")
#   } else {
#     Results <- c(Results,rt_LN_Beta=meta$b,rt_LN_SE=meta$se,rt_LN_lower=meta$ci.lb,
#                  rt_LN_upper=meta$ci.ub,rt_LN_tval=meta$zval,rt_LN_p.value=meta$pval,
#                  rt_LN_I2=meta$I2,rt_LN_QEp=meta$QEp)
#   }
#
#   # META-ANALYSIS OF sc_LN ESTIMATES
#   meta <- tryCatch({rma(yi=sc_LN_Beta, sei=sc_LN_SE,
#                         data=datatable[complete.cases(sc_LN_Beta,sc_LN_SE)],
#                         method="REML",control=list(maxiter=1000, stepadj=.5))},
#                    error=function(e) e)
#   if(!class(meta)[1]=="rma.uni"){
#     Results <- c(Results,sc_LN_Beta=NA,sc_LN_SE=NA,
#                  sc_LN_lower=NA,sc_LN_upper=NA,
#                  sc_LN_tval=NA,sc_LN_p.value=NA,
#                  sc_LN_I2=NA,sc_LN_QEp=NA)
#     Note <- c(Note,"sc_LN Meta-Analysis Failed")
#   } else {
#     Results <- c(Results,sc_LN_Beta=meta$b,sc_LN_SE=meta$se,
#                  sc_LN_lower=meta$ci.lb,sc_LN_upper=meta$ci.ub,
#                  sc_LN_tval=meta$zval,sc_LN_p.value=meta$pval,
#                  sc_LN_I2=meta$I2,sc_LN_QEp=meta$QEp)
#   }
#
#   # META-ANALYSIS OF sc_LN ESTIMATES
#   meta <- tryCatch({rma(yi=log_LN_Beta, sei=log_LN_SE,
#                         data=datatable[complete.cases(log_LN_Beta,log_LN_SE)],
#                         method="REML",control=list(maxiter=1000, stepadj=.5))},
#                    error=function(e) e)
#   if(!class(meta)[1]=="rma.uni"){
#     Results <- c(Results,log_LN_Beta=NA,log_LN_SE=NA,
#                  log_LN_lower=NA,log_LN_upper=NA,
#                  log_LN_tval=NA,log_LN_p.value=NA,
#                  log_LN_I2=NA,log_LN_QEp=NA)
#     Note <- c(Note,"log_LN Meta-Analysis Failed")
#   } else {
#     Results <- c(Results,log_LN_Beta=meta$b,log_LN_SE=meta$se,
#                  log_LN_lower=meta$ci.lb,log_LN_upper=meta$ci.ub,
#                  log_LN_tval=meta$zval,log_LN_p.value=meta$pval,
#                  log_LN_I2=meta$I2,log_LN_QEp=meta$QEp)
#   }
#
#   # META-ANALYSIS OF MUQR.Median ESTIMATES
#   meta <- tryCatch({rma(yi=MUQR.Median_Beta, sei=MUQR.Median_SE,
#                         data=datatable[complete.cases(MUQR.Median_Beta,MUQR.Median_SE)],
#                         method="REML",control=list(maxiter=1000, stepadj=.5))},
#                    error=function(e) e)
#   if(!class(meta)[1]=="rma.uni"){
#     Results <- c(Results,MUQR.Median_Beta=NA,MUQR.Median_SE=NA,MUQR.Median_lower=NA,
#                  MUQR.Median_upper=NA,MUQR.Median_tval=NA,
#                  MUQR.Median_p.value=NA,MUQR.Median_I2=NA,MUQR.Median_QEp=NA)
#     Note <- c(Note,"MUQR.Median Meta-Analysis Failed")
#   } else {
#     Results <- c(Results,MUQR.Median_Beta=meta$b,MUQR.Median_SE=meta$se,
#                  MUQR.Median_lower=meta$ci.lb,MUQR.Median_upper=meta$ci.ub,
#                  MUQR.Median_tval=meta$zval,MUQR.Median_p.value=meta$pval,
#                  MUQR.Median_I2=meta$I2,MUQR.Median_QEp=meta$QEp)
#   }
#
#   # META-ANALYSIS OF rt_MUQR.Median ESTIMATES
#   meta <- tryCatch({rma(yi=rt_MUQR.Median_Beta, sei=rt_MUQR.Median_SE,
#                         data=datatable[complete.cases(rt_MUQR.Median_Beta,rt_MUQR.Median_SE)],
#                         method="REML",control=list(maxiter=1000, stepadj=.5))},
#                    error=function(e) e)
#   if(!class(meta)[1]=="rma.uni"){
#     Results <- c(Results,rt_MUQR.Median_Beta=NA,rt_MUQR.Median_SE=NA,rt_MUQR.Median_lower=NA,
#                  rt_MUQR.Median_upper=NA,rt_MUQR.Median_tval=NA,
#                  rt_MUQR.Median_p.value=NA,rt_MUQR.Median_I2=NA,rt_MUQR.Median_QEp=NA)
#     Note <- c(Note,"rt_MUQR.Median Meta-Analysis Failed")
#   } else {
#     Results <- c(Results,rt_MUQR.Median_Beta=meta$b,rt_MUQR.Median_SE=meta$se,
#                  rt_MUQR.Median_lower=meta$ci.lb,rt_MUQR.Median_upper=meta$ci.ub,
#                  rt_MUQR.Median_tval=meta$zval,rt_MUQR.Median_p.value=meta$pval,
#                  rt_MUQR.Median_I2=meta$I2,rt_MUQR.Median_QEp=meta$QEp)
#   }
#
#   # META-ANALYSIS OF sc_MUQR.Median ESTIMATES
#   meta <- tryCatch({rma(yi=sc_MUQR.Median_Beta, sei=sc_MUQR.Median_SE,
#                         data=datatable[complete.cases(sc_MUQR.Median_Beta,sc_MUQR.Median_SE)],
#                         method="REML",control=list(maxiter=1000, stepadj=.5))},
#                    error=function(e) e)
#   if(!class(meta)[1]=="rma.uni"){
#     Results <- c(Results,sc_MUQR.Median_Beta=NA,sc_MUQR.Median_SE=NA,sc_MUQR.Median_lower=NA,
#                  sc_MUQR.Median_upper=NA,sc_MUQR.Median_tval=NA,
#                  sc_MUQR.Median_p.value=NA,sc_MUQR.Median_I2=NA,sc_MUQR.Median_QEp=NA)
#     Note <- c(Note,"sc_MUQR.Median Meta-Analysis Failed")
#   } else {
#     Results <- c(Results,sc_MUQR.Median_Beta=meta$b,sc_MUQR.Median_SE=meta$se,
#                  sc_MUQR.Median_lower=meta$ci.lb,sc_MUQR.Median_upper=meta$ci.ub,
#                  sc_MUQR.Median_tval=meta$zval,sc_MUQR.Median_p.value=meta$pval,
#                  sc_MUQR.Median_I2=meta$I2,sc_MUQR.Median_QEp=meta$QEp)
#   }
#   # META-ANALYSIS OF log_MUQR.Median ESTIMATES
#   meta <- tryCatch({rma(yi=log_MUQR.Median_Beta, sei=log_MUQR.Median_SE,
#                         data=datatable[complete.cases(log_MUQR.Median_Beta,log_MUQR.Median_SE)],
#                         method="REML",control=list(maxiter=1000, stepadj=.5))},
#                    error=function(e) e)
#   if(!class(meta)[1]=="rma.uni"){
#     Results <- c(Results,log_MUQR.Median_Beta=NA,log_MUQR.Median_SE=NA,log_MUQR.Median_lower=NA,
#                  log_MUQR.Median_upper=NA,log_MUQR.Median_tval=NA,
#                  log_MUQR.Median_p.value=NA,log_MUQR.Median_I2=NA,log_MUQR.Median_QEp=NA)
#     Note <- c(Note,"log_MUQR.Median Meta-Analysis Failed")
#   } else {
#     Results <- c(Results,log_MUQR.Median_Beta=meta$b,log_MUQR.Median_SE=meta$se,
#                  log_MUQR.Median_lower=meta$ci.lb,log_MUQR.Median_upper=meta$ci.ub,
#                  log_MUQR.Median_tval=meta$zval,log_MUQR.Median_p.value=meta$pval,
#                  log_MUQR.Median_I2=meta$I2,log_MUQR.Median_QEp=meta$QEp)
#   }
#
#   # META-ANALYSIS OF abs_MUQR.Median ESTIMATES
#   meta <- tryCatch({rma(yi=abs_MUQR.Median_Beta, sei=MUQR.Median_SE,
#                         data=datatable[complete.cases(abs_MUQR.Median_Beta,MUQR.Median_SE)],
#                         method="REML",control=list(maxiter=1000, stepadj=.5))},
#                    error=function(e) e)
#   if(!class(meta)[1]=="rma.uni"){
#     Results <- c(Results,abs_MUQR.Median_Beta=NA,abs_MUQR.Median_SE=NA,abs_MUQR.Median_lower=NA,
#                  abs_MUQR.Median_upper=NA,abs_MUQR.Median_tval=NA,
#                  abs_MUQR.Median_p.value=NA,abs_MUQR.Median_I2=NA,abs_MUQR.Median_QEp=NA)
#     Note <- c(Note,"abs_MUQR.Median Meta-Analysis Failed")
#   } else {
#     Results <- c(Results,abs_MUQR.Median_Beta=meta$b,abs_MUQR.Median_SE=meta$se,
#                  abs_MUQR.Median_lower=meta$ci.lb,abs_MUQR.Median_upper=meta$ci.ub,
#                  abs_MUQR.Median_tval=meta$zval,abs_MUQR.Median_p.value=meta$pval,
#                  abs_MUQR.Median_I2=meta$I2,abs_MUQR.Median_QEp=meta$QEp)
#   }
#
#   # META-ANALYSIS OF MUQR.MetaTau ESTIMATES
#   meta <- tryCatch({rma(yi=MUQR.MetaTau_Beta, sei=MUQR.MetaTau_SE,
#                         data=datatable[complete.cases(MUQR.MetaTau_Beta,MUQR.MetaTau_SE)],
#                         method="REML",control=list(maxiter=1000, stepadj=.5))},
#                    error=function(e) e)
#   if(!class(meta)[1]=="rma.uni"){
#     Results <- c(Results,MUQR.MetaTau_Beta=NA,MUQR.MetaTau_SE=NA,MUQR.MetaTau_lower=NA,
#                  MUQR.MetaTau_upper=NA,MUQR.MetaTau_tval=NA,
#                  MUQR.MetaTau_p.value=NA,MUQR.MetaTau_I2=NA,MUQR.MetaTau_QEp=NA)
#     Note <- c(Note,"MUQR.MetaTau Meta-Analysis Failed")
#   } else {
#     Results <- c(Results,MUQR.MetaTau_Beta=meta$b,MUQR.MetaTau_SE=meta$se,
#                  MUQR.MetaTau_lower=meta$ci.lb,MUQR.MetaTau_upper=meta$ci.ub,
#                  MUQR.MetaTau_tval=meta$zval,MUQR.MetaTau_p.value=meta$pval,
#                  MUQR.MetaTau_I2=meta$I2,MUQR.MetaTau_QEp=meta$QEp)
#   }
#
#   # META-ANALYSIS OF rt_MUQR.MetaTau ESTIMATES
#   meta <- tryCatch({rma(yi=rt_MUQR.MetaTau_Beta, sei=rt_MUQR.MetaTau_SE,
#                         data=datatable[complete.cases(rt_MUQR.MetaTau_Beta,rt_MUQR.MetaTau_SE)],
#                         method="REML",control=list(maxiter=1000, stepadj=.5))},
#                    error=function(e) e)
#   if(!class(meta)[1]=="rma.uni"){
#     Results <- c(Results,rt_MUQR.MetaTau_Beta=NA,rt_MUQR.MetaTau_SE=NA,rt_MUQR.MetaTau_lower=NA,
#                  rt_MUQR.MetaTau_upper=NA,rt_MUQR.MetaTau_tval=NA,
#                  rt_MUQR.MetaTau_p.value=NA,rt_MUQR.MetaTau_I2=NA,rt_MUQR.MetaTau_QEp=NA)
#     Note <- c(Note,"rt_MUQR.MetaTau Meta-Analysis Failed")
#   } else {
#     Results <- c(Results,rt_MUQR.MetaTau_Beta=meta$b,rt_MUQR.MetaTau_SE=meta$se,
#                  rt_MUQR.MetaTau_lower=meta$ci.lb,rt_MUQR.MetaTau_upper=meta$ci.ub,
#                  rt_MUQR.MetaTau_tval=meta$zval,rt_MUQR.MetaTau_p.value=meta$pval,
#                  rt_MUQR.MetaTau_I2=meta$I2,rt_MUQR.MetaTau_QEp=meta$QEp)
#   }
#
#   # META-ANALYSIS OF sc_MUQR.MetaTau ESTIMATES
#   meta <- tryCatch({rma(yi=sc_MUQR.MetaTau_Beta, sei=sc_MUQR.MetaTau_SE,
#                         data=datatable[complete.cases(sc_MUQR.MetaTau_Beta,sc_MUQR.MetaTau_SE)],
#                         method="REML",control=list(maxiter=1000, stepadj=.5))},
#                    error=function(e) e)
#   if(!class(meta)[1]=="rma.uni"){
#     Results <- c(Results,sc_MUQR.MetaTau_Beta=NA,sc_MUQR.MetaTau_SE=NA,sc_MUQR.MetaTau_lower=NA,
#                  sc_MUQR.MetaTau_upper=NA,sc_MUQR.MetaTau_tval=NA,
#                  sc_MUQR.MetaTau_p.value=NA,sc_MUQR.MetaTau_I2=NA,sc_MUQR.MetaTau_QEp=NA)
#     Note <- c(Note,"sc_MUQR.MetaTau Meta-Analysis Failed")
#   } else {
#     Results <- c(Results,sc_MUQR.MetaTau_Beta=meta$b,sc_MUQR.MetaTau_SE=meta$se,
#                  sc_MUQR.MetaTau_lower=meta$ci.lb,sc_MUQR.MetaTau_upper=meta$ci.ub,
#                  sc_MUQR.MetaTau_tval=meta$zval,sc_MUQR.MetaTau_p.value=meta$pval,
#                  sc_MUQR.MetaTau_I2=meta$I2,sc_MUQR.MetaTau_QEp=meta$QEp)
#   }
#
#   # META-ANALYSIS OF log_MUQR.MetaTau ESTIMATES
#   meta <- tryCatch({rma(yi=log_MUQR.MetaTau_Beta, sei=log_MUQR.MetaTau_SE,
#                         data=datatable[complete.cases(log_MUQR.MetaTau_Beta,log_MUQR.MetaTau_SE)],
#                         method="REML",control=list(maxiter=1000, stepadj=.5))},
#                    error=function(e) e)
#   if(!class(meta)[1]=="rma.uni"){
#     Results <- c(Results,log_MUQR.MetaTau_Beta=NA,log_MUQR.MetaTau_SE=NA,log_MUQR.MetaTau_lower=NA,
#                  log_MUQR.MetaTau_upper=NA,log_MUQR.MetaTau_tval=NA,
#                  log_MUQR.MetaTau_p.value=NA,log_MUQR.MetaTau_I2=NA,log_MUQR.MetaTau_QEp=NA)
#     Note <- c(Note,"log_MUQR.MetaTau Meta-Analysis Failed")
#   } else {
#     Results <- c(Results,log_MUQR.MetaTau_Beta=meta$b,log_MUQR.MetaTau_SE=meta$se,
#                  log_MUQR.MetaTau_lower=meta$ci.lb,log_MUQR.MetaTau_upper=meta$ci.ub,
#                  log_MUQR.MetaTau_tval=meta$zval,log_MUQR.MetaTau_p.value=meta$pval,
#                  log_MUQR.MetaTau_I2=meta$I2,log_MUQR.MetaTau_QEp=meta$QEp)
#   }
#
#   # META-ANALYSIS OF abs_MUQR.MetaTau ESTIMATES
#   meta <- tryCatch({rma(yi=abs_MUQR.MetaTau_Beta, sei=MUQR.MetaTau_SE,
#                         data=datatable[complete.cases(abs_MUQR.MetaTau_Beta,MUQR.MetaTau_SE)],
#                         method="REML",control=list(maxiter=1000, stepadj=.5))},
#                    error=function(e) e)
#   if(!class(meta)[1]=="rma.uni"){
#     Results <- c(Results,abs_MUQR.MetaTau_Beta=NA,abs_MUQR.MetaTau_SE=NA,abs_MUQR.MetaTau_lower=NA,
#                  abs_MUQR.MetaTau_upper=NA,abs_MUQR.MetaTau_tval=NA,
#                  abs_MUQR.MetaTau_p.value=NA,abs_MUQR.MetaTau_I2=NA,abs_MUQR.MetaTau_QEp=NA)
#     Note <- c(Note,"abs_MUQR.MetaTau Meta-Analysis Failed")
#   } else {
#     Results <- c(Results,abs_MUQR.MetaTau_Beta=meta$b,abs_MUQR.MetaTau_SE=meta$se,
#                  abs_MUQR.MetaTau_lower=meta$ci.lb,abs_MUQR.MetaTau_upper=meta$ci.ub,
#                  abs_MUQR.MetaTau_tval=meta$zval,abs_MUQR.MetaTau_p.value=meta$pval,
#                  abs_MUQR.MetaTau_I2=meta$I2,abs_MUQR.MetaTau_QEp=meta$QEp)
#   }
#
#   # META-ANALYSIS OF MUQR.MetaTau p.values
#   # meta <- fun1(x=datatable[complete.cases(MUQR.MetaTau_tval)][,MUQR.MetaTau_tval],
#   #              Statistic="z",Weights=datatable[complete.cases(MUQR.MetaTau_tval)][,N],
#   #              log_p.value=log_p.value)
#   # Results <- c(Results,Stouf_MUQR.MetaTau_zval=meta[1],
#   #              Stouf_MUQR.MetaTau_p.value=meta[2])
#   # meta <- fun2(x=datatable[complete.cases(MUQR.MetaTau_tval)][,MUQR.MetaTau_tval],
#   #              Statistic="z",log_p.value=log_p.value)
#   # Results <- c(Results,Fisch_MUQR.MetaTau_Tval=meta[1],
#   #              Fisch_MUQR.MetaTau_p.value=meta[2])
#
#   # META-ANALYSIS OF LEVENE p.values
#   meta <- fun1(x=datatable[complete.cases(Levene_p.value)][,Levene_p.value],
#                Statistic="p",Weights=datatable[complete.cases(Levene_p.value)][,N],
#                log_p.value=log_p.value)
#   Results <- c(Results,Stouf_Levene_zval=meta[1],Stouf_Levene_p.value=meta[2])
#   meta <- fun2(x=datatable[complete.cases(Levene_p.value)][,Levene_p.value],
#                Statistic="p",log_p.value=log_p.value)
#   Results <- c(Results,Fisch_Levene_Tval=meta[1],Fisch_Levene_p.value=meta[2])
#
#   # META-ANALYSIS OF z2 ESTIMATES
#   meta <- tryCatch({rma(yi=z2_Beta, sei=z2_SE,
#                         data=datatable[complete.cases(z2_Beta,z2_SE)],method="REML",
#                         control=list(maxiter=1000, stepadj=.5))},
#                    error=function(e) e)
#   if(!class(meta)[1]=="rma.uni"){
#     Results <- c(Results,z2_Beta=NA,z2_SE=NA,z2_lower=NA,z2_upper=NA,z2_tval=NA,
#                  z2_p.value=NA,z2_I2=NA,z2_QEp=NA)
#     Note <- c(Note,"z2 Meta-Analysis Failed")
#   } else {
#     Results <- c(Results,z2_Beta=meta$b,z2_SE=meta$se,z2_lower=meta$ci.lb,
#                  z2_upper=meta$ci.ub,z2_tval=meta$zval,z2_p.value=meta$pval,
#                  z2_I2=meta$I2,z2_QEp=meta$QEp)
#   }
#   if(is.null(Note)){
#     Results <- c(Results,Notes=NA)
#   } else {
#     Note <- paste(Note, sep="", collapse=", ")
#     Results <- c(Results,Notes=Note)
#   }
#   return(Results)
# }
#
# # NEED A GENERIC FUNCTION
#
# # ============= UTILITY FUNCTIONS ====================
# rntransform <- function(x){qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
#
# comb <- function(x, ...) {
#   #rbind comb for data.table
#   mapply(function(...) rbind(...,fill=TRUE),x,...,SIMPLIFY=FALSE)
# }
#
