
mqr <- function(datatable,y,g,covariates=NULL,tau=seq(0.05, 0.95, by=0.05),
                mqr.method="UQR",boot.m="mcmb",boot.R=200,seed=31371,enable.dither=TRUE,
                Univariable=TRUE,fitOLS=TRUE){
  # datatable=DT;y="response";g="x";tau=seq(0.05,0.95,length.out=10);
  # mqr.method="UQR";boot.m="mcmb";boot.R=200;seed=31371;enable.dither=TRUE;
  # Univariable=TRUE;fitOLS=TRUE;covariates=NULL
  if(y=="y"){
    stop("The name of your respoonse variable is 'y' please rename it before calling MQR()")
  }
  if(is.null(covariates)){
    FML  <- as.formula(paste(y,g,sep=" ~ "))
  } else {
    covs <- paste(covariates,sep="",collapse=" + ")
    FML  <- as.formula(paste(paste(y,g,sep=" ~ "),covs,sep=" + "))
  }
  DT <- data.table(model.frame(FML, datatable))

  if(mqr.method=="CQR"){
    ifelse(nrow(DT)<2000,{rq_meth <- "br"},
           ifelse(nrow(DT)<100000,{rq_meth <- "fn"},
                  {rq_meth <- "pfn"}))

    # CALCULATE QUANTILE ESTIMATES
    boot.Bs <- matrix(NA,nrow=boot.R,ncol=length(tau),dimnames=list(1:boot.R,tau))
    Beta <- Notes <- rep(NA,length(tau))

    # BOOTSTRAP ERROR ESTIMATES UISNG boot.m
    g.ptm <- proc.time()
    for(i in 1:length(tau)){
      Q.mod <- tryCatch({quantreg::rq(FML, tau[i], DT, method=rq_meth)}, error=function(err){
        # print(paste0("rq model fitting at tau=",tau[i]," failed because : ",err))
        return(c("ERROR",paste0("rq(tau=",tau[i],") failed because : ",err[[1]])))})
      if(length(Q.mod)==1){
        Notes[i] <- Q.mod[[2]]
        next
      }
      Beta[i] <- coef(Q.mod)[g]
      i.Bs  <- tryCatch({
        bs.M <- boot.m
        set.seed(seed)
        quantreg::summary.rq(Q.mod, se="boot", covariance=TRUE, R=boot.R,bsmethod=bs.M)$B},
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
    Results <- tryCatch({mcqr(Beta,boot.Bs)}, error=function(err){
      # print(paste0("MR-CQR failed because : ",err))
      return(c("ERROR",paste0("MCQR failed because : ",err[[1]])))})
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
          Q.mod <- tryCatch({quantreg::rq(FML, tau[i], DT, method=rq_meth)}, error=function(err){
            # print(paste0("rq model fitting at tau=",tau[i]," failed because : ",err))
            return(c("ERROR",paste0("rq(tau=",tau[i],") dithered failed because : ",err[[1]])))})
          if(length(Q.mod)==1){
            Notes.2[i] <- Q.mod[[2]]
            next
          }
          Beta[i] <- coef(Q.mod)[g]
          i.Bs  <- tryCatch({
            bs.M <- boot.m
            set.seed(seed)
            quantreg::summary.rq(Q.mod, se="boot", covariance=TRUE, R=boot.R,bsmethod=bs.M)$B},
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
        Results <- tryCatch({mcqr(Beta,boot.Bs)}, error=function(err){
          # print(paste0("MR-CQR failed because : ",err))
          return(c("ERROR",paste0("MCQR dithered failed because : ",err[[1]])))})
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
    if(fitOLS){
      LN.sum <- coef(summary(lm(FML,DT)))
      LN.Results <- c(LN_Beta=LN.sum[g,1],LN_SE=LN.sum[g,2],
                      LN_tval=LN.sum[g,3],LN_p.value=LN.sum[g,4])
    } else {
      LN.Results <- NULL
    }

    # Add MEDIAN Regression Estimates
    Q.mod <- tryCatch({quantreg::rq(FML, tau=0.5, DT,method=rq_meth)},error=function(err){
      return(c("ERROR",paste0("rq(Median) using ",rq_meth, " , failed because : ",err[[1]])))})
    if(length(Q.mod)==1){
      Notes <- c(Notes,Q.mod[[2]])
      t.meth <- "br"
      Q.mod <- tryCatch({quantreg::rq(FML, tau=0.5, DT,method=t.meth)},error=function(err){
        return(c("ERROR",paste0("rq(Median) using ",t.meth, " , failed because : ",err[[1]])))})
    }
    if(length(Q.mod)==1){
      Notes <- c(Notes,Q.mod[[2]],"MCQR.Median fit failed")
      MEDIAN.Results <- c(MCQR.Median_Beta=NA,MCQR.Median_SE=NA,MCQR.Median_tval=NA,
                          MCQR.Median_p.value=NA)
    } else {
      Q.sum  <- tryCatch({
        bs.M <- boot.m
        set.seed(seed)
        quantreg::summary.rq(Q.mod, se="boot", covariance=TRUE, R=boot.R, bsmethod=bs.M)},
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

  } else if (mqr.method=="UQR"){
    DT[,y:=DT[[y]]]
    g.ptm <- proc.time()
    # Adjust for Median effects of g (i.e. residual scale)
    RIF <- rif(y=DT[,y], taus=0.5)
    Beta_0.5 <- coef(lm(as.formula(paste("RIF",g,sep=" ~ ")),data=DT))[g]
    DT[,y:= y - Beta_0.5*DT[[g]]]

    RIF <- rif(y=DT[,y], taus=tau) # compute RIF at taus of interest

    if(Univariable){
      FML <- as.formula(paste("RIF",g,sep=" ~ "))
    } else {
      FML <- as.formula(paste("RIF",as.character(FML)[1],as.character(FML)[3],
                              sep="",collapse=""))
    }

    Models <- lm(FML, DT)
    Beta <- Models$coefficients[2,]
    COV <- rif.cov(Models, pred=g)

    # re-centred tau
    taus <- tau-0.5
    Results <- tryCatch({muqr(Beta=Beta, COV=COV, taus=taus)}, error=function(err){
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
        RIF <- rif(y=DT[,Dithered_Response], taus=0.5)
        Beta_0.5 <- coef(lm(as.formula(paste("RIF",g,sep=" ~ ")),data=DT))[g]
        DT[,Residual_Dithered_Response:=Dithered_Response - Beta_0.5*DT[[g]]]

        RIF <- rif(y=DT[,Residual_Dithered_Response], taus=tau) # compute RIF at taus of interest
        Models <- lm(FML, DT)
        Beta <- Models$coefficients[2,]
        COV <- RIFmod.cov(Models, pred=g)
        Results <- tryCatch({muqr(Beta=Beta, COV=COV, taus=taus)}, error=function(err){
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
          RIF <- rif(y=DT[,Dithered_Response], taus=0.5)
          RIF.sum <- coef(summary(lm(FML,DT)))
          MEDIAN.Result <- c(MUQR.Median_Beta=RIF.sum[g,1],MUQR.Median_SE=RIF.sum[g,2],
                             MUQR.Median_tval=RIF.sum[g,3],MUQR.Median_p.value=RIF.sum[g,4])

          # Add MEAN Regression Estimates.
          if(fitOLS){
            FML <- as.formula(paste("Dithered_Response",FML[[1]],FML[[3]]))
            LN.sum <- coef(summary(lm(FML,DT)))
            LN.Results <- c(LN_Beta=LN.sum[g,1],LN_SE=LN.sum[g,2],
                            LN_tval=LN.sum[g,3],LN_p.value=LN.sum[g,4])
          } else {
            LN.Results <- NULL
          }

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
      RIF <- rif(y=DT[[y]], taus=0.5)
      RIF.sum <- coef(summary(lm(FML,DT)))
      MEDIAN.Result <- c(MUQR.Median_Beta=RIF.sum[g,1],MUQR.Median_SE=RIF.sum[g,2],
                         MUQR.Median_tval=RIF.sum[g,3],MUQR.Median_p.value=RIF.sum[g,4])

      # Add MEAN Regression Estimates.
      if(fitOLS){
        FML <- as.formula(paste(y,FML[[1]],FML[[3]]))
        LN.sum <- coef(summary(lm(FML,DT)))
        LN.Results <- c(LN_Beta=LN.sum[g,1],LN_SE=LN.sum[g,2],
                        LN_tval=LN.sum[g,3],LN_p.value=LN.sum[g,4])
      } else {
        LN.Results <- NULL
      }
    }
    g.etm <- proc.time()-g.ptm
    if(length(Notes)>1){
      Notes <- paste(Notes,sep="",collapse=". ")
    }
    Results <- c(N=nrow(DT),EAF=sum(DT[[g]])/(2*nrow(DT)),LN.Results,MEDIAN.Result,Results,
                 MUQR.TtC=g.etm[[3]],Notes=Notes)
    return(Results)
  }
}
