
# ============= FUNCTIONS AND PACKAGES ====================
library(data.table); library(quantreg); library(foreach); library(sn); library(deming);
library(lubridate); library(rmutil); library(doParallel);

#http://stackoverflow.com/questions/27279164/output-list-of-two-rbinded-data-frames-with-foreach-in-r
comb <- function(x, ...) {
  mapply(cbind,x,...,SIMPLIFY=FALSE)
}

# So right now i've settelled on 6 Different methods to skew the distribution still
# producing mean of zero and var of 1.
Error.Distribution <- function(N,Type="Normal",a,b,Skew="None"){
  #Type="Right-Outliers"; a=0.02;N=10000
  # METHOD 0 : Normal (not skewed)
  #            Type="Normal"
  # METHOD 1 : uses the Skew-T Distribution. This is pretty good but may not produce
  #            ideal shape of skewing. Increasing the term "alpha" term increases
  #            skewing. Dont bother with alpha < 15 because skewing is not
  #            substantial enough.
  #            Type="Skew-T"
  # METHOD 2 : uses the chi distribution which is really effective, increasing "df" term
  #            brings it closer and closer to normal distribution. 2 < df < 10 is ideal
  #            for testing.
  #            Type="Chi"
  # METHOD 3 : uses the beta distribution which allows independent control of right and
  #            left side of the distribution, "a" term specifies steepness on the left
  #            side of the distribution (increasing the term decreases the steepness)
  #            and "b" is the same except on the right side of the distribution.
  #            a=10,b=10 is pretty close-ish to normal (but not really). To skew,
  #            set a and b to lower or higher values. ie left skewed distribution a<<b
  #            try a=2,b=10 for only minor skewing ideally a=1 and 10<b<20
  #            SET a=2, b=2 FOR FAT DISTRIUBTION WITH SHORT TAILS
  #            Type="Beta"
  # METHOD 4 : log-normal distribution. Not much to this, no room to modify skewness
  #            throwing it in there because BMI is often log transformed
  #            Type="LogNormal"
  # METHOD 5 : Laplace distribution. Used to Specify SKINY DISTRIBUTION AND LONG TAILS
  #            Type="Laplace"
  # METHOD 6 : Uniform distribution. Used to Specify FAT DISTRIBUTION AND SHORT TAILS
  #            Type="Uniform"
  # METHOD 7 : Adding outliers to distribution. Specify left or right side or both. Specify
  #          : percent of samples to be made into outliers. a=0.02
  #            Type="Left-Outliers", "Right-Outliers", "Both-Outliers"
  ifelse(Skew=="Left",{dir <- -1}, {dir <- 1})
  if(Type=="Normal"){
    x <- rnorm(N)
    return(x)
  }
  if(Type=="Chi"){
    x <- dir*rchisq(N,df=a)
    x <- (x-mean(x))/sd(x)
    return(x)
  }
  if(Type=="Laplace"){
    x <- rlaplace(N)
    x <- (x-mean(x))/sd(x)
    return(x)
  }
  if(Type=="Uniform"){
    x <- runif(N)
    x <- (x-mean(x))/sd(x)
    return(x)
  }
  if(Type=="Skew-T"){
    x <- dir*rst(n=N,xi=0,omega=1,alpha=a)
    x <- (x-mean(x))/sd(x)
    return(x)
  }
  if(Type=="LogNormal"){
    x <- dir*rlnorm(N)
    x <- (x-mean(x))/sd(x)
    return(x)
  }
  if(Type=="Both-Outliers"){
    x <- rnorm(N)
    q <- quantile(x,c(0.25,0.75))
    x[sample(1:N,a*N)] <-c(runif(a*N/2,q[1]-3*(q[2]-q[1]),q[1]-1.6*(q[2]-q[1])),
                           runif(a*N/2,q[2]+1.6*(q[2]-q[1]),q[2]+3*(q[2]-q[1])))
    return(x)
  }
  if(Type=="Left-Outliers"){
    x <- rnorm(N)
    q <- quantile(x,c(0.25,0.75))
    x[sample(1:N,a*N)] <- runif(a*N,q[1]-3*(q[2]-q[1]),q[1]-1.6*(q[2]-q[1]))
    return(x)
  }
  if(Type=="Right-Outliers"){
    x <- rnorm(N)
    q <- quantile(x,c(0.25,0.75))
    x[sample(1:N,a*N)] <- runif(a*N,q[2]+1.6*(q[2]-q[1]),q[2]+3*(q[2]-q[1]))
    return(x)
  }
  if(Type=="Beta"){
    x <- dir*rbeta(N,a,b)
    x <- (x-mean(x))/sd(x)
    return(x)
  }
  #hist(x,breaks=100,prob=TRUE,main="",xlab="Outcome",col="grey",border=NA)
  #curve(dnorm(x, mean=mean(x), sd=sd(x)),add=TRUE, col="red")
  #qqnorm(x);qqline(x,col="red")
}

# So the overall idea here is to set genotypes based on hwe p-value. To get no
# hwe distortion hwe.p=1
Generate_Genotypes <- function(N=10000,MAF=0.05,hwe.p=1,min.grp.size=5,
                               Pare.Encoding=TRUE){
  #p=0.05;N=10000;p=0.05;hwe.p=10^-5;min.grp.size=5
  p <- MAF
  q <- 1-p
  ifelse(Pare.Encoding,{aa <- -2*p; ab <- 1-2*p; bb <- 2-2*p},{aa <- 0; ab <- 1; bb <- 2})

  # Genotype numbers based on perfect hwe are as follows
  pp <- N*p^2; pq <- 2*N*p*(1-p); qq <- N*(1-p)^2
  #p^2+2*p*q+q^2 # should be equal to 1

  if(hwe.p==1){
    # for no hwe distortions at all
    #obs.G <- c(rep(aa,qq), rep(ab,pq), rep(bb,pp)) # perfect hwe
    # Akram recommends doing sample instead
    obs.G <- sample(c(aa, ab, bb), N, prob=c(qq,pq,pp)/N, replace=TRUE)
  } else {
    # If hwe.p is specified then need to work backwards, from hwe.p to get chi value
    # at that hwe.p
    hwe.chi<- qchisq(hwe.p,1,lower.tail=FALSE)

    # obs.pq set to minimum (ie obs.pq=pq) and qq to max, ie pp will be minimized
    # Note this can mean that pp is bellow min.pp and even bellow 0 so need to take
    # precautions
    obs.pq <- pq
    obs.qq <- sqrt(hwe.chi*pp*qq/(pp+qq))+qq
    obs.pp <- N-obs.qq-obs.pq
    if(obs.pp<min.grp.size){
      # if obs.pp is too low for my liking ie. either less than 5 or less than 0
      # set obs.pp to my minimim group size threshold
      obs.pp <- min.grp.size
      obs.qq <- sqrt((hwe.chi-((obs.pp-pp)^2)/pp)*pq*qq/(pq+qq))+qq
      obs.pq <- N-obs.pp-obs.qq
    }
    #obs.G <- c(rep(aa,round(obs.qq)),rep(ab,round(obs.pq)),rep(bb,round(obs.pp)))
    # Akram recommends doing sample instead
    obs.G <- sample(c(aa, ab, bb), N, prob=c(obs.qq,obs.pq,obs.pp)/N, replace=TRUE)
    rm(hwe.chi,obs.pp,obs.qq,obs.pq)
  }
  if(length(obs.G)<N){
    obs.G <- c(obs.G,rep(NA, N-length(obs.G)))
  }
  #rm(p,q,aa,ab,bb,pp,pq,qq)
  return(obs.G)
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

MUQR.function <- function(datatable, tau=seq(0.05, 0.95, by=0.05), y, g,
                          covariates=NULL,Scale="Z_Score",Univariable=TRUE,
                          ufun1=RIF.Transformation,ufun2=cov.RIF,
                          ufun3=MetaReg.MUQR.function){
  #datatable=Batch.Data; tau=Taus; y=outcome; g=l.SNPs[j];
  #Univariable=TRUE;Scale="Z_Score"
  #covariates=c(covariates,"STUDY");
  #ufun1=RIF.Transformation; ufun2=cov.RIF;ufun3=MetaReg.MUQR.function

  g.ptm <- proc.time()
  if(is.null(covariates)){
    FML  <- as.formula(paste(y,g,sep=" ~ "))
  } else {
    covs <- paste(covariates,sep="",collapse=" + ")
    FML  <- as.formula(paste(paste(y,g,sep=" ~ "),covs,sep=" + "))
  }
  DT <- data.table(model.frame(FML, datatable))

  #if(Scale=="Median_Residuals"){
  #RIF <- ufun1(y=DT[[y]], taus=0.5)
  #Beta_0.5 <- coef(lm(as.formula(paste("RIF",g,sep=" ~ ")),data=DT))[g]
  #DT[,y:=DT[[y]] - Beta_0.5*DT[[g]]]
  #} else
  if(Scale=="Z_Score_ByStudy"){
    DT[,y:=DT[[y]]]
    DT[,y:=scale(y),by=STUDY]
  } else if(Scale=="Z_Score"){
    DT[,y:=DT[[y]]]
    DT[,y:=scale(y)]
  } else {
    DT[,y:=DT[[y]]]
  }

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
  meta.Results <- ufun3(Beta=Beta, COV=COV, taus=taus)
  Result <- c(meta.Results["(Intercept)",],meta.Results["taus",])
  names(Result) <- c("MUQR.MetaMarg_Beta","MUQR.MetaMarg_SE","MUQR.MetaMarg_tval",
                     "MUQR.MetaMarg_p.value","MUQR.MetaTau_Beta","MUQR.MetaTau_SE",
                     "MUQR.MetaTau_tval","MUQR.MetaTau_p.value")
  g.etm <- proc.time()-g.ptm

  # Add MEDIAN Regression Estimates
  RIF <- ufun1(y=DT[[y]], taus=0.5)
  RIF.sum <- coef(summary(lm(FML,DT)))
  MEDIAN.Result <- c(MUQR.Median_Beta=RIF.sum[g,1],MUQR.Median_SE=RIF.sum[g,2],
                     MUQR.Median_tval=RIF.sum[g,3],MUQR.Median_p.value=RIF.sum[g,4])

  # Add MEAN Regression Estimates on Raw Scale.
  FML <- as.formula(paste(y,as.character(FML)[1],as.character(FML)[3],
                          sep="",collapse=""))
  LN.sum <- coef(summary(lm(FML,DT)))
  LN.Results <- c(LN_Beta=LN.sum[g,1],LN_SE=LN.sum[g,2],
                  LN_tval=LN.sum[g,3],LN_p.value=LN.sum[g,4])

  Result <- c(N=nrow(DT),EAF=sum(DT[[g]])/(2*nrow(DT)),LN.Results,MEDIAN.Result,Result,
              MUQR.TtC=g.etm[[3]])
  return(Result)
}

#sfun1=RIF.Transformation,
#sfun2=cov.RIF,
#sfun3=MetaReg.MUQR.function
# Function takes 2 methods. "Classic-Levene" and "Brown-Forsythe". The latter is more
# robust to skewed distributions. The default is "classic.Levene". Note: Selecting
# "Brown-Forsythe" gives the exact same result as levene.test function in lawstat package
Levene.function <- function(datatable, y, g, covariates=NULL, method="Brown-Forsythe",
                            log.p=FALSE){
  #datatable=Data; method="Brown-Forsythe"
  if(is.null(covariates)){
    FML  <- as.formula(paste(y,g,sep=" ~ "))
  } else {
    covs <- paste(covariates,sep="",collapse=" + ")
    FML  <- as.formula(paste(paste(y,g,sep=" ~ "),covs,sep=" + "))
  }

  DT <- data.table(model.frame(FML, datatable))
  DT[,g:=ifelse(DT[[g]]<=0.5,0,
                ifelse(DT[[g]]<=1.5,1,
                       ifelse(DT[[g]]>1.5,2, NA)))]
  if(is.null(covariates)){
    DT[,y:=DT[[y]]]
  } else {
    FML.r <- paste(y, covs, sep=" ~ ")
    LMOD <- lm(FML,DT)
    DT[,y:=resid(LMOD)]
  }

  DT <- DT[,.(y,g)]

  # first calculate Zij (which is key) for each subject and add it to
  # data.table.
  # Use Median for "Brown–Forsythe" and Mean for "Classic-Levene". Otherwise tests
  # are the same.
  if(method=="Brown-Forsythe"){
    DT[,Zij:=abs(y-median(y)),by=g]
  } else if(method=="Classic-Levene"){
    DT[,Zij:=abs(y-mean(y)),by=g]
  }

  # calculate Z..
  Z.. <- mean(DT[,Zij])
  # find the group number
  k <- length(unique(DT[,g]))
  # find the sample size
  N <- nrow(DT)
  # calculate Zij_Zi, which is the difference between Zij (for each
  # subject) and the mean of Zij for each genotype, squared
  DT[,Zij_Zi:=(Zij-mean(Zij))^2,by=g]
  W <- (N-k)*sum(DT[,length(y)*(mean(Zij)-Z..)^2,by=g][,V1])/
    ((k-1)*sum(DT[,Zij_Zi]))
  if(log.p){
    p <- pf(W,df1=k-1,df2=N-k, lower.tail=FALSE, log.p=TRUE)
  } else {
    p <- pf(W,df1=k-1,df2=N-k, lower.tail=FALSE)
  }
  Result <- c(Levene_Fval=W,Levene_df1=k-1,Levene_df2=N-k,Levene_p.value=p)
  return(Result)
}

z2.function <- function(datatable, y, g,covariates=NULL){
  #datatable=Data; tau=Taus; y="BMI"; g=SNPs[j]; Adjusted.for=NA
  #covariates=c(covariates,"STUDY");
  #fun1=RIF.Transformation; fun2=cov.RIF;fun3=MetaReg.function

  if(is.null(covariates)){
    FML  <- as.formula(paste(y,g,sep=" ~ "))
  } else {
    covs <- paste(covariates,sep="",collapse=" + ")
    FML  <- as.formula(paste(paste(y,g,sep=" ~ "),covs,sep=" + "))
  }
  DT <- data.table(model.frame(FML, datatable))

  DT[,z2:=(scale(DT[[y]]))^2]
  FML <- as.formula(paste("z2",as.character(FML)[1],as.character(FML)[3],
                          sep="",collapse=""))

  # Add z2 Regression Estimates
  LN.sum <- coef(summary(lm(FML,DT)))
  Result <- c(z2_Beta=LN.sum[g,1],z2_SE=LN.sum[g,2],
              z2_tval=LN.sum[g,3],z2_p.value=LN.sum[g,4])
  return(Result)
}

predict.deming <- function(deming.object,x){
  #deming.object=dScale.Model; x=SNP.Results[,LN_Beta]
  y <- deming.object$coefficients[2]*x + deming.object$coefficients[1]
  Xp <- cbind(1,x)
  b <- deming.object$coefficients
  yh <- c(Xp %*% b)
  V <- deming.object$variance
  var.fit <- rowSums((Xp %*% V) * Xp)
  se <- sqrt(var.fit)
  return(c(Predicted_Beta=y,Predicted_SE=se))
}

# Extracting sparsity.rq from Qtools
sparsity.rq <- function (object, se="iid", hs=TRUE) {
  mt <- terms(object)
  m <- model.frame(object)
  y <- model.response(m)
  x <- model.matrix(mt, m, contrasts = object$contrasts)
  wt <- model.weights(object$model)
  taus <- object$tau
  nq <- length(taus)
  eps <- .Machine$double.eps^(2/3)
  vnames <- dimnames(x)[[2]]
  residm <- as.matrix(object$residuals)
  n <- length(y)
  p <- nrow(as.matrix(object$coef))
  rdf <- n - p
  if (!is.null(wt)) {
    residm <- residm * wt
    x <- x * wt
    y <- y * wt
  }
  if (is.null(se)) {
    se <- "nid"
  }
  spar <- dens <- matrix(NA, n, nq)
  for (i in 1:nq) {
    tau <- taus[i]
    if (se == "iid") {
      resid <- residm[, i]
      pz <- sum(abs(resid) < eps)
      h <- max(p + 1, ceiling(n * bandwidth.rq(tau, n,
                                               hs = hs)))
      ir <- (pz + 1):(h + pz + 1)
      ord.resid <- sort(resid[order(abs(resid))][ir])
      xt <- ir/(n - p)
      spar[, i] <- rq(ord.resid ~ xt)$coef[2]
      dens[, i] <- 1/spar[, i]
    }
    else if (se == "nid") {
      h <- bandwidth.rq(tau, n, hs = hs)
      if (tau + h > 1)
        stop("tau + h > 1:  error in summary.rq")
      if (tau - h < 0)
        stop("tau - h < 0:  error in summary.rq")
      bhi <- rq.fit.fnb(x, y, tau = tau + h)$coef
      blo <- rq.fit.fnb(x, y, tau = tau - h)$coef
      dyhat <- x %*% (bhi - blo)
      if (any(dyhat <= 0))
        warning(paste(sum(dyhat <= 0), "non-positive fis"))
      f <- pmax(0, (2 * h)/(dyhat - eps))
      dens[, i] <- f
      spar[, i] <- 1/f
    }
    else if (se == "ker") {
      h <- bandwidth.rq(tau, n, hs = hs)
      if (tau + h > 1)
        stop("tau + h > 1:  error in summary.rq")
      if (tau - h < 0)
        stop("tau - h < 0:  error in summary.rq")
      uhat <- c(residm[, i])
      h <- (qnorm(tau + h) - qnorm(tau - h)) * min(sqrt(var(uhat)),
                                                   (quantile(uhat, 0.75) - quantile(uhat, 0.25))/1.34)
      f <- dnorm(uhat/h)/h
      dens[, i] <- f
      spar[, i] <- 1/f
    }
  }
  colnames(dens) <- colnames(spar) <- taus
  return(list(density = dens, sparsity = spar, bandwidth = h))
}

MetaReg.MCQR.function <- function(Beta,boot.Bs,test="z"){
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

MCQR.function <- function(datatable, tau=seq(0.05, 0.95, by=0.05),y,g,
                          covariates=NULL, Scale="Z_Score",Univariable=TRUE,
                          boot.R=200,fun1=MetaReg.MCQR.function){
  #datatable=DT;tau=Taus;y=outcome;g=j.SNPs[l];covariates=MEGA.covs; Univariable=TRUE;
  #Scale="Z_Score";boot.R=N; fun1=MetaReg.MCQR.function

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
  bs.procedure <- "mcmb"

  #if(Scale=="Median_Residuals"){
  #Beta_0.5 <- coef(rq(as.formula(paste(y,g,sep=" ~ ")),tau=0.5,
  #                   data=DT,method=rq_meth))[g]
  #DT[,y:=DT[[y]] - Beta_0.5*DT[[g]]]
  #} else
  if(Scale=="Z_Score_ByStudy"){
    DT[,y:=DT[[y]]]
    DT[,y:=scale(y),by=STUDY]
  } else if(Scale=="Z_Score"){
    DT[,y:=DT[[y]]]
    DT[,y:=scale(y)]
  } else {
    DT[,y:=DT[[y]]]
  }

  if(Univariable){
    FML <- as.formula(paste("y",g,sep=" ~ "))
  } else {
    FML <- as.formula(paste("y",as.character(FML)[1],as.character(FML)[3],
                            sep="",collapse=""))
  }

  # CALCULATE QUANTILE ESTIMATES
  Q.mod <- rq(FML, tau, DT, method=rq_meth)
  Beta <- as.vector(coef(Q.mod)[g,])
  names(Beta) <- gsub("tau= ","", colnames(coef(Q.mod)))

  # BOOTSTRAP ERROR ESTIMATES
  Q.sum  <- tryCatch({
    bs.M <- bs.procedure
    summary.rqs(Q.mod, se="boot", covariance=TRUE, R=boot.R,
                bsmethod=bs.M)
  }, error=function(e) e)
  if(!class(Q.sum)[1]=="summary.rqs"){
    bs.M <- "pwy"
    Q.sum <- tryCatch({
      summary.rqs(Q.mod, se="boot", covariance=TRUE, R=boot.R,
                  bsmethod=bs.M)
    }, error=function(e) e)
  }
  if(!class(Q.sum)[1]=="summary.rqs"){
    Results <- "MCQR failed"
    return(Results)
  }
  boot.Bs <- matrix(NA,nrow=boot.R,ncol=length(Beta),
                    dimnames=list(1:boot.R,names(Beta)))
  for(h in 1:length(Beta)){
    #h <- 1
    colnames(Q.sum[[h]]$B) <- rownames(Q.sum[[h]]$coefficients)
    boot.Bs[,h] <- Q.sum[[h]]$B[,g]
  }

  # META-REGRESSION ANALYSIS
  Results <- fun1(Beta,boot.Bs)
  Results <- c(Results["(Intercept)",],Results["taus",])
  names(Results) <- c("MCQR.MetaMarg_Beta","MCQR.MetaMarg_SE",
                      "MCQR.MetaMarg_tval","MCQR.MetaMarg_p.value","MCQR.MetaTau_Beta",
                      "MCQR.MetaTau_SE","MCQR.MetaTau_tval","MCQR.MetaTau_p.value")
  g.etm <- proc.time()-g.ptm

  # Add MEAN Regression Estimates on Raw Scale.
  LN.sum <- coef(summary(lm(FML,DT)))
  LN.Results <- c(LN_Beta=LN.sum[g,1],LN_SE=LN.sum[g,2],
                  LN_tval=LN.sum[g,3],LN_p.value=LN.sum[g,4])

  # Add MEDIAN Regression Estimates
  Q.mod <- rq(FML, tau=0.5, DT, method=rq_meth)
  Q.sum  <- tryCatch({
    bs.M <- bs.procedure
    summary.rq(Q.mod, se="boot", covariance=TRUE, R=boot.R,
               bsmethod=bs.M)
  }, error=function(e) e)
  if(!class(Q.sum)[1]=="summary.rq"){
    bs.M <- "pwy"
    Q.sum <- tryCatch({
      summary.rq(Q.mod, se="boot", covariance=TRUE, R=boot.R,
                 bsmethod=bs.M)
    }, error=function(e) e)
  }
  if(!class(Q.sum)[1]=="summary.rq"){
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

  Results <- c(N=nrow(DT),EAF=sum(DT[[g]])/(2*nrow(DT)),LN.Results,MEDIAN.Results,
               Results,MCQR.TtC=g.etm[[3]])
  return(Results)
}


# Function bellow is same as above execpt all Qreg.Assessment options are implimented
# as is Levene.Test. Requires both qreg.function and Pare.Levene.function
Sim.function  <- function(Arguments, Special="None", N, MAF, Varying="v.GxE",
                          b0, v.G, v.E, v.GxE,Interaction.Dir="+ve", hwe.p=1,
                          min.grp.size=5, Pare.Encoding=TRUE,Type="Normal", a, b,
                          Skew="None",se, Adjusted.for=NULL,No.taus=19,
                          fun1=Generate_Genotypes,fun2=Error.Distribution,fun3=QREG.function,
                          fun4=Levene.Joint.function,fun5=RIF.function,
                          Levene.function=Levene.function,RIF.Transformation=RIF.Transformation,
                          sparsity.rq=sparsity.rq,cov.RIF=cov.RIF,cov.rq=cov.rq,
                          MetaReg.function=MetaReg.function){
  #Special="None";N=2000;MAF=0.05; Varying="v.GxE"; b0=0;
  #v.G=NA; v.E=NA; v.GxE=NA; v.GxE=0;
  #Interaction.Dir="+ve";hwe.p=1E-5;min.grp.size=5; Pare.Encoding=TRUE; Type="Normal";a=0.02;
  #b=c(1.5,3);Skew="None"; se="nid"; Adjusted.for=NULL
  #No.taus=19
  #fun1=Generate_Genotypes;fun2=Error.Distribution;fun3=MCQR.function;fun4=MUQR.function
  #Arguments <- arguments[l,]
  if(!missing(Arguments)){
    Special <- Arguments["Special"]; N <- as.numeric(Arguments["N"]);
    Varying <- Arguments["Varying"]; MAF <- as.numeric(Arguments["MAF"]);
    b0 <- as.numeric(Arguments["b0"]);v.G <- as.numeric(Arguments["v.G"]);
    v.E <- as.numeric(Arguments["v.E"]);v.GxE <- as.numeric(Arguments["v.GxE"]);
    Interaction.Dir <- Arguments["Interaction.Dir"];
    hwe.p <- as.numeric(Arguments["hwe.p"]); min.grp.size <- Arguments["min.grp.size"];
    Pare.Encoding <- Arguments["Pare.Encoding"]; Type <- Arguments["Type"];
    a <- as.numeric(Arguments["a"]); b <- as.numeric(Arguments["b"]);
    Skew <- Arguments["Skew"]; se <- Arguments["se"];
    Adjusted.for <- Arguments["Adjusted.for"]; No.taus <- as.numeric(Arguments["No.taus"])
  }
  p <- MAF
  L <- 2*p*(1-p)
  b1.G <- sqrt(v.G/L)
  b2.E <- sqrt(v.E)
  b3.GxE <- sqrt(v.GxE/L)
  C <- sqrt(1 - v.G - v.E - v.GxE)
  test.data <- data.table(data.frame(g=sample(fun1(N=N,MAF=p,hwe.p=hwe.p,
                                                   min.grp.size=min.grp.size,
                                                   Pare.Encoding=Pare.Encoding),
                                              replace=FALSE),
                                     E=rnorm(N),
                                     e=fun2(N,Type=Type,a,b,Skew)))
  if(Special=="None"){
    if(Interaction.Dir=="+ve"){
      test.data[,y:=b0 + b1.G*g + b2.E*E + b3.GxE*g*E + C*e]
    } else if (Interaction.Dir=="-ve"){
      test.data[,y:=b0 + b1.G*g + b2.E*E - b3.GxE*g*E + C*e]
    }
  }



  # For Joint Test of LN Model
  #lm1 <- lm(y ~ E, data=test.data)
  lm2 <- lm(y ~ g + E + g*E, data=test.data)
  #F <- as.data.frame(anova(lm2,lm1)) # the order or lm1 or lm2 does not matter
  #LN.Joint <- pf(F[2,"F"], df1=abs(F[2,"Df"]),
  #               df2=max(F[,"Res.Df"]), lower=FALSE)
  LN.Model <- coef(summary(lm2))
  #LN.Model
  Taus <- seq(0.05,0.95,length.out=No.taus)
  test.data[,pheno:=y]
  test.data[,geno:=g]
  MUQR.Time <- system.time(
    MUQR.Result <- fun4(test.data,tau=Taus,y="pheno", g="geno",Scale="Raw")
  )["elapsed"]

  Levene.Time <- system.time(
    Levene.Result <- fun4(test.data,tau=Taus,y="y", g="g",Scale="Raw")
  )["elapsed"]
  RIF.Time <- system.time(
    RIF.Result <- fun5(test.data,tau=Taus,y="y", g="g",Adjusted.for=Adjusted.for)
  )["elapsed"]
  Results <- cbind(b1.G=LN.Model["g",1],b2.E=LN.Model["E",1],b3.GxE=LN.Model["g:E",1],
                   GxE.LN_p.value=LN.Model["g:E",4],LN_Joint=LN.Joint,
                   QREG.Result, QREG.Time=QREG.Time,Levene.Result,
                   Levene.Time=Levene.Time,RIF.Result, RIF.Time=RIF.Time)
  #rm(test.data,LN.Model,lm1,lm2,F,LN.Joint)
  return(Results)
}

ScalingSim.function  <- function(Arguments,
                                 fun1=Generate_Genotypes,
                                 fun2=Error.Distribution,
                                 fun3=MUQR.function,
                                 fun4=predict.deming,
                                 sfun1=RIF.Transformation,
                                 sfun2=cov.RIF,
                                 sfun3=MetaReg.MUQR.function){
  #Special="None";N=2000;MAF=0.05; Varying="v.GxE"; b0=0;
  #v.G=NA; v.E=NA; v.GxE=NA; v.GxE=0;
  #Interaction.Dir="+ve";hwe.p=1E-5;min.grp.size=5; Pare.Encoding=TRUE; Type="Normal";a=0.02;
  #b=c(1.5,3);Skew="None"; se="nid"; Adjusted.for=NULL
  #No.taus=19; var.amp=0.5; Scaling.R=100
  #fun1=Generate_Genotypes; fun2=Error.Distribution;fun3=MUQR.function; fun4=predict.deming;
  #Arguments <- arguments[l,]
  #browser()
  if(!missing(Arguments)){
    Special <- Arguments["Special"]; N <- as.numeric(Arguments["N"]);
    Varying <- Arguments["Varying"]; MAF <- as.numeric(Arguments["MAF"]);
    b0 <- as.numeric(Arguments["b0"]);v.G <- as.numeric(Arguments["v.G"]);
    v.E <- as.numeric(Arguments["v.E"]);v.GxE <- as.numeric(Arguments["v.GxE"]);
    Interaction.Dir <- Arguments["Interaction.Dir"];
    hwe.p <- as.numeric(Arguments["hwe.p"]); min.grp.size <- as.numeric(Arguments["min.grp.size"]);
    Pare.Encoding <- Arguments["Pare.Encoding"]; Type <- Arguments["Type"];
    a <- as.numeric(Arguments["a"]); b <- as.numeric(Arguments["b"]);
    Skew <- Arguments["Skew"]; se <- Arguments["se"];
    Adjusted.for <- Arguments["Adjusted.for"]; No.taus <- as.numeric(Arguments["No.taus"]);
    var.amp <- as.numeric(Arguments["var.amp"]); Scaling.R <- as.numeric(Arguments["Scaling.R"]);
    Scaling.nodes <- as.numeric(Arguments["Scaling.nodes"])
  }

  p <- MAF
  L <- 2*p*(1-p)
  b1.G <- sqrt(v.G/L)
  b2.E <- sqrt(v.E)
  b3.GxE <- sqrt(v.GxE/L)

  #gamma_0 <- 0
  #gamma_1 <- 0.1
  #v.combo <- v.G*(1+gamma_1) + v.E*(1+gamma_1) + v.GxE*(1+gamma_1)
  #C <- sqrt(1-v.combo)
  C <- sqrt(1-v.G-v.E-v.GxE)
  #temp.1 <- NULL
  #for(i in 1:20){
  #print(i)
  test.data <- data.table(data.frame(g=sample(fun1(N,p,hwe.p,min.grp.size,Pare.Encoding),
                                              replace=FALSE),
                                     E=rnorm(N),
                                     e=fun2(N,Type,a,b,Skew)))
  if(Special=="None"){
    if(Interaction.Dir=="+ve"){
      test.data[,y:=b0 + b1.G*g + b2.E*E + b3.GxE*g*E + C*e]
    } else if (Interaction.Dir=="-ve"){
      test.data[,y:=b0 + b1.G*g + b2.E*E - b3.GxE*g*E + C*e]
    }
  } else if(Special=="Global_Scaling"){
    test.data[,pre_y:=b0 + b1.G*g + b2.E*E + b3.GxE*g*E + C*e]
    test.data <- test.data[order(pre_y)]
    amplifier <- runif(N,1,1+var.amp)
    amplifier <- amplifier[order(amplifier)]
    test.data[,Amp.Factor:=amplifier]
    test.data[,y:=pre_y*Amp.Factor]
    #if(Interaction.Dir=="+ve"){
    #test.data[,y:=b0 + b1.G*g + b2.E*E + b3.GxE*g*E + C*(gamma_0 + gamma_1*b1.G)*g*e +
    #            C*(gamma_0 + gamma_1*b2.E)*E*e + C*(gamma_0 + gamma_1*b3.GxE)*g*E*e]
    #} else if (Interaction.Dir=="-ve"){
    #test.data[,y:=b0 + b1.G*g + b2.E*E - b3.GxE*g*E + C*(gamma_0 + gamma_1*b1.G)*g*e +
    #            C*(gamma_0 + gamma_1*b2.E)*E*e + C*(gamma_0 + gamma_1*b3.GxE)*g*E*e]
    #}
    #test.data[,y:=b0 + b1.G*g + b2.E*E + b3.GxE*g*E + C*e]
    #test.data[,y:=scale(y)]
  }

  #test.data[,mean(pre_y)]
  #test.data[,var(pre_y)]
  #test.data[,mean(y)]
  #test.data[,var(y)]
  #hist(test.data[,pre_y],breaks=100,prob=TRUE,main="",xlab="Outcome",col="grey",border=NA)
  #curve(dnorm(x, mean=mean(test.data[,pre_y]), sd=sd(test.data[,pre_y])),add=TRUE, col="red")
  #qqnorm(test.data[,y]);qqline(test.data[,y],col="red")

  # hist(test.data[,y],breaks=100,prob=TRUE,main="",xlab="Outcome",col="grey",border=NA)
  # curve(dnorm(x, mean=mean(test.data[,y]), sd=sd(test.data[,y])),add=TRUE, col="red")
  # qqnorm(test.data[,y]);qqline(test.data[,y],col="red")

  # For Joint Test of LN Model
  #lm1 <- lm(y ~ E, data=test.data)
  lm2 <- lm(y ~ g + E + g*E, data=test.data)
  #F <- as.data.frame(anova(lm2,lm1)) # the order or lm1 or lm2 does not matter
  #LN.Joint <- pf(F[2,"F"], df1=abs(F[2,"Df"]),
  #               df2=max(F[,"Res.Df"]), lower=FALSE)
  LN.Model <- coef(summary(lm2))
  #LN.Model
  Taus <- seq(0.05,0.95,length.out=No.taus)
  test.data[,pheno:=y]
  test.data[,geno:=g]
  MUQR.Result <- fun3(test.data,tau=Taus,y="pheno", g="geno",Scale="Raw")
  #MUQR.Result
  #cor(temp.1[,as.numeric(MUQR.MetaTau_Beta)], temp.1[,as.numeric(LN_Beta)])
  #plot(x=temp.1[,as.numeric(LN_Beta)],y=temp.1[,as.numeric(MUQR.MetaTau_Beta)])
  #nrow(temp.1[as.numeric(LN_p.value)<0.05])/nrow(temp.1)
  #nrow(temp.1[as.numeric(MUQR.MetaTau_p.value)<0.05])/nrow(temp.1)

  #MUQR.Result

  #Levene.Time <- system.time(
  #  Levene.Result <- fun4(test.data,tau=Taus,y="y", g="g",Scale="Raw")
  #)["elapsed"]
  #RIF.Time <- system.time(
  #  RIF.Result <- fun5(test.data,tau=Taus,y="y", g="g",Adjusted.for=Adjusted.for)
  #)["elapsed"]
  Results <- c(b1.G=LN.Model["g",1],b2.E=LN.Model["E",1],b3.GxE=LN.Model["g:E",1],
               GxE.LN_p.value=LN.Model["g:E",4],MUQR.Result[3:length(MUQR.Result)])
  #rm(test.data,LN.Model,lm1,lm2,F,LN.Joint)
  #Scaling.R=100

  Scaling.batches <- Scaling.R/Scaling.nodes
  u.cols <- c("LN_Beta","LN_SE","LN_tval","LN_p.value","MUQR.Median_Beta","MUQR.Median_SE","MUQR.Median_tval",
              "MUQR.Median_p.value","MUQR.MetaMarg_Beta","MUQR.MetaMarg_SE","MUQR.MetaMarg_tval",
              "MUQR.MetaMarg_p.value","MUQR.MetaTau_Beta","MUQR.MetaTau_SE","MUQR.MetaTau_tval",
              "MUQR.MetaTau_p.value","MUQR.TtC")

  Scaling.Results <- foreach(u=1:Scaling.nodes, .combine="rbind",.multicombine=TRUE, .inorder=FALSE,
                             .export=c("RIF.Transformation","cov.RIF","MetaReg.MUQR.function",
                                       "MUQR.function")
  ) %dopar% {
    library(data.table)
    #library(quantreg); library(pracma);library(rmutil);library(sn)
    u.Result <- matrix(NA,Scaling.batches,length(u.cols),dimnames=list(1:Scaling.batches,u.cols))
    for(o in 1:Scaling.batches){
      #o <- 1
      o.DT <- copy(test.data)
      o.DT[,scale.geno:=sample(o.DT[["geno"]],N)]
      o.Result <- fun3(o.DT,tau=Taus,y="pheno",g="scale.geno",Scale="Raw")
      u.Result[o,] <- o.Result[3:length(o.Result)]
    }
    return(u.Result)
  }
  Scaling.Results <- data.table(Scaling.Results)

  Scale.Model <- lm(MUQR.MetaTau_tval~LN_tval,data=Scaling.Results)
  #Scale.Model <- lm(MUQR.MetaTau_Beta~LN_Beta,data=Scaling.Results)
  Scale.Summary <- coef(summary(Scale.Model))
  #Deming.Scale.Model <- deming(MUQR.MetaTau_Beta~LN_Beta,data=Scaling.Results,xstd=LN_SE,ystd=MUQR.MetaTau_SE)

  # cor(x=Scaling.Results[,LN_Beta], y=Scaling.Results[,MUQR.MetaTau_Beta])
  # Mean <- Scale.Summary[,"Estimate"]
  # Upper <- Mean+1.96*Scale.Summary[,"Std. Error"]
  # Lower <- Mean-1.96*Scale.Summary[,"Std. Error"]
  # dem.Mean <- Deming.Scale.Model$coefficients
  # dem.Upper <- Deming.Scale.Model$ci[,2]
  # dem.Lower <- Deming.Scale.Model$ci[,1]
  # png(file.path("/home/../media/StoreB/arkan","Scaling.Plot1.png"))
  # plot(x=Scaling.Results[,LN_Beta],y=Scaling.Results[,MUQR.MetaTau_Beta],pch=20,col="grey")
  # plot(x=Scaling.Results[,LN_tval],y=Scaling.Results[,MUQR.MetaTau_tval],pch=20,col="grey")
  # abline(a=Mean[[1]],b=Mean[[2]])
  # abline(a=Upper[[1]],b=Upper[[2]],lty=2)
  # abline(a=Lower[[1]],b=Lower[[2]],lty=2)
  #
  # abline(a=dem.Mean[[1]],b=dem.Mean[[2]],col="blue")
  # abline(a=dem.Upper[[1]],b=dem.Upper[[2]],lty=2,col="blue")
  # abline(a=dem.Lower[[1]],b=dem.Lower[[2]],lty=2,col="blue")
  # dev.off()


  # LM.Predicted <- predict(Scale.Model,data.frame(LN_Beta=Results["LN_Beta"]),se.fit=TRUE)
  # Results <- c(Results,Predicted_MUQR.MetaTau_Beta=LM.Predicted$fit[[1]],
  #              Predicted_MUQR.MetaTau_SE=LM.Predicted$se.fit)
  LM.Predicted <- predict(Scale.Model,data.frame(LN_tval=Results["LN_tval"]),se.fit=TRUE)
  Results <- c(Results,MUQR.MetaTau_tval_PredictedByScaling=LM.Predicted$fit[[1]])
  #tval <- diff(c(Results["MUQR.MetaTau_tval"],LM.Predicted$fit[[1]]))
  tval <- abs(Results["MUQR.MetaTau_tval"])-abs(LM.Predicted$fit[[1]])
  # tval <- (Results["MUQR.MetaTau_Beta"]-Results["Predicted_MUQR.MetaTau_Beta"])/
  #   (sqrt((Results["MUQR.MetaTau_SE"])^2+(Results["Predicted_MUQR.MetaTau_SE"])^2))
  #p.value <- 2*pt(abs(tval),df=Scaling.R-1, lower=FALSE)
  p.value <- 2*pnorm(abs(tval), lower=FALSE)
  # Results <- c(Results,MetaTau_vs_Predicted_tval=tval[[1]],MetaTau_vs_Predicted_p.value=p.value[[1]])
  Results <- c(Results,Scaling_Adjusted_tval=tval[[1]],Scaling_Adjusted_p.value=p.value[[1]])

  # DM.Predicted <- fun4(Deming.Scale.Model,Results["LN_Beta"])
  # names(DM.Predicted) <- c("dem.Predicted_MUQR.MetaTau_Beta","dem.Predicted_MUQR.MetaTau_SE")
  # Results <- c(Results,DM.Predicted)
  # tval <- (Results["MUQR.MetaTau_Beta"]-Results["dem.Predicted_MUQR.MetaTau_Beta"])/
  #   (sqrt((Results["MUQR.MetaTau_SE"])^2+(Results["dem.Predicted_MUQR.MetaTau_SE"])^2))
  # # p.value <- 2*pt(abs(tval),df=Scaling.R-1, lower=FALSE)
  # p.value <- 2*pnorm(abs(tval), lower=FALSE)
  # Results <- c(Results,MetaTau_vs_dem.Predicted_tval=tval[[1]],MetaTau_vs_dem.Predicted_p.value=p.value[[1]])
  return(Results)
}

# ============= Running Simulations: Varying number of Taus ============
# Conditions to test
Special <- "None"; N <- 10000; MAF <- NULL;
b0 <- 0; Interaction.Dir <- "+ve"; hwe.p <- 1; min.grp.size <- 3; Pare.Encoding <- TRUE;
Type <- "Normal"; a <- NULL ; b <- NULL; Skew <- "None"; tau.thresh <- NULL; bT.G <- NULL;
se <- "iid"; Adjusted.for <- NULL; No.taus <- 19

Varying <- "v.GxE"
v.G <- 0.002;
v.E <- 0.25
v.GxE <- c(0,0.006)
MAF <- 0.05
No.taus <- c(2,5,10,20,30,50)

# Setup Results organization
arguments<-matrix(NA, nrow=length(v.GxE)*length(MAF)*length(No.taus), ncol=19)
colnames(arguments) <- c("Special", "N", "MAF", "Varying","b0", "v.G", "v.E", "v.GxE",
                         "Interaction.Dir", "hwe.p", "min.grp.size","Pare.Encoding", "Type",
                         "a", "b", "Skew", "se","Adjusted.for","No.taus")
est <- list()
# Specify Simulations and Computations
R <- 1
nodes <- 1
batches <- R/nodes
Run.time <- NULL

#cl<-makeCluster(nodes)
#registerDoParallel(cl)
#cl <- startMPIcluster(count=nodes)
#registerDoMPI(cl)
#opt <- list(chunkSize=R/nodes, profile=TRUE, seed=159653)
for(i in 1:length(v.GxE)){
  for(j in 1:length(MAF)){
    print(paste("Analysis Started at: ", Sys.time(), sep=""))
    z.ptm <- proc.time()
    for(k in 1:length(No.taus)){
      #i <- 1; j <- 1; k <- 4
      l <- length(No.taus)*length(MAF)*(i-1)+(j-1)*length(No.taus)+k
      t <- cbind(Varying,v.E=v.E,MAF=MAF[j],v.GxE=v.GxE[i],N,b0,v.G,
                 Special,Interaction.Dir,hwe.p,min.grp.size,Pare.Encoding,Type,a,b,Skew,
                 se,Adjusted.for,No.taus=No.taus[k])
      arguments[l,match(colnames(t),colnames(arguments))] <- t

      est[[l]] <- foreach(z=1:nodes, .combine="rbind",.multicombine=TRUE, .inorder=FALSE
      ) %dopar% {
        library(data.table)
        #library(quantreg); library(pracma);library(sn); library(rmutil);
        z.Result <- matrix(NA,batches,16,
                           dimnames=list(
                             1:batches,
                             c("b1.G", "b2.E","b3.GxE","GxE.LN_p.value" ,"LN_Joint",
                               "G.LN_p.value","Levene_p.value","Levene_Joint",
                               "Levene.Time","G.RIF_p.value","RIF.Meta_Beta","RIF.Meta_SE",
                               "RIF.Meta_p.value","RIF.MetaMarg_p.value","RIF_Joint","RIF.Time")))
        for(h in 1:batches){
          #h <- 1
          z.Result[h,] <- Sim.function(Arguments=arguments[l,])
        }
        z.Result <- cbind(z.Result,matrix(rep(arguments[l,],batches),batches,
                                          length(arguments[l,]),byrow=TRUE,
                                          dimnames=list(1:batches,names(arguments[l,]))))
        return(z.Result)
      }
      print(paste(No.taus[k], " - finished at ", Sys.time(), sep=""))
    }
    z.etm <- proc.time()-z.ptm
    Run.time <- c(Run.time, z.etm["elapsed"])
    if(l/length(No.taus)<nrow(arguments)/length(No.taus)){
      TtC <- round(range(Run.time)*((nrow(arguments)-l)/length(No.taus)))
    } else {
      TtC <- round(range(Run.time))
    }
    print(paste("Analysis: ", l/length(No.taus), " of ", nrow(arguments)/length(No.taus),
                ". Completed in ", seconds_to_period(as.numeric(round(z.etm["elapsed"]))),
                ". Approx. Time to Completion: ",  seconds_to_period(as.numeric(TtC[1])),
                " - ",seconds_to_period(as.numeric(TtC[2])), sep=""))
  }
}
v.GxE_arguments <- arguments
v.GxE_est <- est
rm(arguments,est);gc()
#stopCluster(cl)
# mpi.quit()	# this will exit R as well
p <- "/home/alabada/pgpstore/R_Studio_Projects/Simulations_Quantile_Regression"

save.image(file=file.path(p,"Results","RData_Objects", "Number_of_Taus_Apr_1_2017.RData"))
list.files(file.path(p,"Results","RData_Objects"))


# ============= Running Simulations: Pare 2010 Fig3 : Varying v.GxE ============
# Conditions to test
Special <- "None"; N <- 10000; MAF <- NULL;
b0 <- 0; Interaction.Dir <- "+ve"; hwe.p <- 1; min.grp.size <- 3; Pare.Encoding <- TRUE;
Type <- "Normal"; a <- NULL ; b <- NULL; Skew <- "None"; tau.thresh <- NULL; bT.G <- NULL;
se <- "iid"; Adjusted.for <- NULL; No.taus <- 19

Varying <- "v.GxE"
v.G <- 0.002;
v.E <- c(0.15,0.25)
v.GxE <- seq(0,0.01, by=0.001)
MAF <- c(0.05,0.2,0.4)

# Setup Results organization
arguments<-matrix(NA, nrow=length(v.E)*length(MAF)*length(v.GxE), ncol=19)
colnames(arguments) <- c("Special", "N", "MAF", "Varying","b0", "v.G", "v.E", "v.GxE",
                         "Interaction.Dir", "hwe.p", "min.grp.size","Pare.Encoding", "Type",
                         "a", "b", "Skew", "se","Adjusted.for","No.taus")
est <- list()
# Specify Simulations and Computations
R <- 1
nodes <- 1
batches <- R/nodes
Run.time <- NULL

#cl<-makeCluster(nodes)
#registerDoParallel(cl)
#cl <- startMPIcluster(count=nodes)
#registerDoMPI(cl)
#opt <- list(chunkSize=R/nodes, profile=TRUE, seed=159653)
for(i in 1:length(v.E)){
  for(j in 1:length(MAF)){
    print(paste("Analysis Started at: ", Sys.time(), sep=""))
    z.ptm <- proc.time()
    for(k in 1:length(v.GxE)){
      #i <- 1; j <- 3; k <- 1
      l <- length(v.GxE)*length(MAF)*(i-1)+(j-1)*length(v.GxE)+k
      t <- cbind(Varying,v.E=v.E[i],MAF=MAF[j],v.GxE=v.GxE[k],N,b0,v.G,
                 Special,Interaction.Dir,hwe.p,min.grp.size,Pare.Encoding,Type,a,b,Skew,
                 se,Adjusted.for,No.taus)
      arguments[l,match(colnames(t),colnames(arguments))] <- t

      est[[l]] <- foreach(z=1:nodes, .combine="rbind",.multicombine=TRUE, .inorder=FALSE
      ) %dopar% {
        library(data.table)
        #library(quantreg); library(pracma);library(sn); library(rmutil);
        z.Result <- matrix(NA,batches,16,
                           dimnames=list(
                             1:batches,
                             c("b1.G", "b2.E","b3.GxE","GxE.LN_p.value" ,"LN_Joint",
                               "G.LN_p.value","Levene_p.value","Levene_Joint",
                               "Levene.Time","G.RIF_p.value","RIF.Meta_Beta","RIF.Meta_SE",
                               "RIF.Meta_p.value","RIF.MetaMarg_p.value","RIF_Joint","RIF.Time")))
        for(h in 1:batches){
          #h <- 1
          z.Result[h,] <- Sim.function(Arguments=arguments[l,])
        }
        z.Result <- cbind(z.Result,matrix(rep(arguments[l,],batches),batches,
                                          length(arguments[l,]),byrow=TRUE,
                                          dimnames=list(1:batches,names(arguments[l,]))))
        return(z.Result)
      }
      print(paste(v.GxE[k], " - finished at ", Sys.time(), sep=""))
    }
    z.etm <- proc.time()-z.ptm
    Run.time <- c(Run.time, z.etm["elapsed"])
    if(l/length(v.GxE)<nrow(arguments)/length(v.GxE)){
      TtC <- round(range(Run.time)*((nrow(arguments)-l)/length(v.GxE)))
    } else {
      TtC <- round(range(Run.time))
    }
    print(paste("Analysis: ", l/length(v.GxE), " of ", nrow(arguments)/length(v.GxE),
                ". Completed in ", seconds_to_period(as.numeric(round(z.etm["elapsed"]))),
                ". Approx. Time to Completion: ",  seconds_to_period(as.numeric(TtC[1])),
                " - ",seconds_to_period(as.numeric(TtC[2])), sep=""))
  }
}
v.GxE_arguments <- arguments
v.GxE_est <- est
rm(arguments,est);gc()
stopCluster(cl)
# mpi.quit()	# this will exit R as well
p <- "/home/alabada/pgpstore/R_Studio_Projects/Simulations_Quantile_Regression"

save.image(file=file.path(p,"Results","RData_Objects", "vGxE_Fig3Sim_Apr_1_2017.RData"))
list.files(file.path(p,"Results","RData_Objects"))
# ============= Running Simulations: Pare 2010 Fig3 : Varying v.G ============
# Conditions to test
Special <- "None"; N <- 10000; MAF <- NULL;
b0 <- 0; Interaction.Dir <- "+ve"; hwe.p <- 1; min.grp.size <- 3; Pare.Encoding <- TRUE;
Type <- "Normal"; a <- NULL ; b <- NULL; Skew <- "None"; tau.thresh <- NULL; bT.G <- NULL;
se <- "iid"; Adjusted.for <- NULL; No.taus <- 19

Varying <- "v.G"
v.G <- seq(0,0.006, by=0.0006)
v.E <- c(0.15,0.25)
v.GxE <- 0.006
MAF <- c(0.05,0.2,0.4)

# Setup Results organization
arguments<-matrix(NA, nrow=length(v.E)*length(MAF)*length(v.G), ncol=19)
colnames(arguments) <- c("Special", "N", "MAF", "Varying","b0", "v.G", "v.E", "v.GxE",
                         "Interaction.Dir", "hwe.p", "min.grp.size","Pare.Encoding", "Type",
                         "a", "b", "Skew", "se","Adjusted.for","No.taus")
est <- list()
# Specify Simulations and Computations
R <- 1
nodes <- 1
batches <- R/nodes
Run.time <- NULL

#cl<-makeCluster(nodes)
#registerDoParallel(cl)
#cl <- startMPIcluster(count=nodes)
#registerDoMPI(cl)
#opt <- list(chunkSize=R/nodes, profile=TRUE, seed=159653)
for(i in 1:length(v.E)){
  for(j in 1:length(MAF)){
    print(paste("Analysis Started at: ", Sys.time(), sep=""))
    z.ptm <- proc.time()
    for(k in 1:length(v.G)){
      #i <- 1; j <- 3; k <- 1
      l <- length(v.G)*length(MAF)*(i-1)+(j-1)*length(v.G)+k
      t <- cbind(Varying,v.E=v.E[i],MAF=MAF[j],v.GxE=v.GxE,N,b0,v.G=v.G[k],
                 Special,Interaction.Dir,hwe.p,min.grp.size,Pare.Encoding,Type,a,b,Skew,
                 se,Adjusted.for,No.taus)
      arguments[l,match(colnames(t),colnames(arguments))] <- t

      est[[l]] <- foreach(z=1:nodes, .combine="rbind",.multicombine=TRUE, .inorder=FALSE
      ) %dopar% {
        library(data.table)
        #library(quantreg); library(pracma);library(sn); library(rmutil);
        z.Result <- matrix(NA,batches,16,
                           dimnames=list(
                             1:batches,
                             c("b1.G", "b2.E","b3.GxE","GxE.LN_p.value" ,"LN_Joint",
                               "G.LN_p.value","Levene_p.value","Levene_Joint",
                               "Levene.Time","G.RIF_p.value","RIF.Meta_Beta","RIF.Meta_SE",
                               "RIF.Meta_p.value","RIF.MetaMarg_p.value","RIF_Joint","RIF.Time")))
        for(h in 1:batches){
          #h <- 1
          z.Result[h,] <- Sim.function(Arguments=arguments[l,])
        }
        z.Result <- cbind(z.Result,matrix(rep(arguments[l,],batches),batches,
                                          length(arguments[l,]),byrow=TRUE,
                                          dimnames=list(1:batches,names(arguments[l,]))))
        return(z.Result)
      }
      print(paste(v.G[k], " - finished at ", Sys.time(), sep=""))
    }
    z.etm <- proc.time()-z.ptm
    Run.time <- c(Run.time, z.etm["elapsed"])
    if(l/length(v.G)<nrow(arguments)/length(v.G)){
      TtC <- round(range(Run.time)*((nrow(arguments)-l)/length(v.G)))
    } else {
      TtC <- round(range(Run.time))
    }
    print(paste("Analysis: ", l/length(v.G), " of ", nrow(arguments)/length(v.G),
                ". Completed in ", seconds_to_period(as.numeric(round(z.etm["elapsed"]))),
                ". Approx. Time to Completion: ",  seconds_to_period(as.numeric(TtC[1])),
                " - ",seconds_to_period(as.numeric(TtC[2])), sep=""))
  }
}
v.G_arguments <- arguments
v.G_est <- est
rm(arguments,est);gc()
#stopCluster(cl)
# mpi.quit()	# this will exit R as well
#p <- "/home/alabada/pgpstore/R_Studio_Projects/Simulations_Quantile_Regression"
#save.image(file=file.path(p,"Results","RData_Objects", "vG_Fig3Sim_Apr_1_2017.RData"))
#list.files(file.path(p,"Results","RData_Objects"))
save.image(file="vG_Fig3Sim_Apr_1_2017.RData")


# ============= Running Simulations: Gauderman-Like FP analysis : Varying v.E ============
# Conditions to test
Special <- "None"; N <- 10000; MAF <- NULL;
b0 <- 0; Interaction.Dir <- "+ve"; hwe.p <- 1; min.grp.size <- 3; Pare.Encoding <- TRUE;
Type <- "Normal"; a <- NULL ; b <- NULL; Skew <- "None"; tau.thresh <- NULL; bT.G <- NULL;
se <- "iid"; Adjusted.for <- NULL; No.taus <- 19

Varying <- "v.E"
v.G <- 0
v.E <- seq(0,0.4, by=0.08)
v.GxE <- 0
Distributions <- c("Normal","Right-Skew")
MAF <- 0.1

# Setup Results organization
arguments<-matrix(NA, nrow=length(Distributions)*length(MAF)*length(v.E), ncol=19)
colnames(arguments) <- c("Special", "N", "MAF", "Varying","b0", "v.G", "v.E", "v.GxE",
                         "Interaction.Dir", "hwe.p", "min.grp.size","Pare.Encoding", "Type",
                         "a", "b", "Skew", "se","Adjusted.for","No.taus")
est <- list()
# Specify Simulations and Computations
R <- 5000
nodes <- 20
batches <- R/nodes
Run.time <- NULL

#cl<-makeCluster(nodes)
#registerDoParallel(cl)
#cl <- startMPIcluster(count=nodes)
#registerDoMPI(cl)
#opt <- list(chunkSize=R/nodes, profile=TRUE, seed=159653)
for(i in 1:length(Distributions)){
  #i <- 1; j <- 1; k <- 1
  if(Distributions[i]=="Normal"){
    Type <- "Normal"; Skew <- "None"; a <- b <- NULL
  } else if(Distributions[i]=="Right-Skew"){
    Type <- "Chi"; a <- 4; Skew <- "Right"; b <- NULL
  }
  for(j in 1:length(MAF)){
    print(paste("Analysis Started at: ", Sys.time(), sep=""))
    z.ptm <- proc.time()
    for(k in 1:length(v.E)){
      l <- length(v.E)*length(MAF)*(i-1)+(j-1)*length(v.E)+k
      t <- cbind(Varying,MAF=MAF[j],v.E=v.E[k],v.G,v.GxE,N,b0,
                 Special,Interaction.Dir,hwe.p,min.grp.size,Pare.Encoding,Type,a,b,Skew,
                 se,Adjusted.for,No.taus)

      arguments[l,match(colnames(t),colnames(arguments))] <- t
      est[[l]] <- foreach(z=1:nodes, .combine="rbind",.multicombine=TRUE, .inorder=FALSE
      ) %dopar% {
        library(data.table)
        #library(rmutil); library(quantreg); library(sn); library(pracma);
        z.Result <- matrix(NA,batches,16,
                           dimnames=list(
                             1:batches,
                             c("b1.G", "b2.E","b3.GxE","GxE.LN_p.value" ,"LN_Joint",
                               "G.LN_p.value","Levene_p.value","Levene_Joint",
                               "Levene.Time","G.RIF_p.value","RIF.Meta_Beta","RIF.Meta_SE",
                               "RIF.Meta_p.value","RIF.MetaMarg_p.value","RIF_Joint","RIF.Time")))
        for(h in 1:batches){
          #h <- 1
          z.Result[h,] <- Sim.function(Arguments=arguments[l,])
        }
        z.Result <- cbind(z.Result,matrix(rep(arguments[l,],batches),batches,
                                          length(arguments[l,]),byrow=TRUE,
                                          dimnames=list(1:batches,names(arguments[l,]))))
        return(z.Result)
      }
      print(paste(v.E[k], " - finished at ", Sys.time(), sep=""))
    }
    z.etm <- proc.time()-z.ptm
    Run.time <- c(Run.time, z.etm["elapsed"])
    if(l/length(v.E)<nrow(arguments)/length(v.E)){
      TtC <- round(range(Run.time)*((nrow(arguments)-l)/length(v.E)))
    } else {
      TtC <- round(range(Run.time))
    }
    print(paste("Analysis: ", l/length(v.E), " of ", nrow(arguments)/length(v.E),
                ". Completed in ", seconds_to_period(as.numeric(round(z.etm["elapsed"]))),
                ". Approx. Time to Completion: ",  seconds_to_period(as.numeric(TtC[1])),
                " - ",seconds_to_period(as.numeric(TtC[2])), sep=""))
  }
}
stopCluster(cl)

p <- "/home/alabada/pgpstore/R_Studio_Projects/Simulations_Quantile_Regression"
list.files(p)
save.image(file=file.path(p,"Results","RData_Objects", "unAdj_FalsePostives_Apr_1_2017.RData"))
list.files(file.path(p,"Results","RData_Objects"))


# ============= Running Simulations: Gauderman-Like FP analysis : Varying v.G ============
# Conditions to test
Special <- "None"; N <- 10000; MAF <- NULL;
b0 <- 0; Interaction.Dir <- "+ve"; hwe.p <- 1; min.grp.size <- 3; Pare.Encoding <- TRUE;
Type <- "Normal"; a <- NULL ; b <- NULL; Skew <- "None"; tau.thresh <- NULL; bT.G <- NULL;
se <- "iid"; Adjusted.for <- "E"; No.taus <- 19

Varying <- "v.G"
v.G <- seq(0,0.01, by=0.002)
v.E <- 0.2
v.GxE <- 0
Distributions <- c("Normal","Right-Skew")
MAF <- 0.1

# Setup Results organization
arguments<-matrix(NA, nrow=length(Distributions)*length(MAF)*length(v.G), ncol=19)
colnames(arguments) <- c("Special", "N", "MAF", "Varying","b0", "v.G", "v.E", "v.GxE",
                         "Interaction.Dir", "hwe.p", "min.grp.size","Pare.Encoding", "Type",
                         "a", "b", "Skew", "se","Adjusted.for","No.taus")
est <- list()
# Specify Simulations and Computations
R <- 1
nodes <- 1
batches <- R/nodes
Run.time <- NULL

#cl<-makeCluster(nodes)
#registerDoParallel(cl)
#cl <- startMPIcluster(count=nodes)
#registerDoMPI(cl)
#opt <- list(chunkSize=R/nodes, profile=TRUE, seed=159653)
for(i in 1:length(Distributions)){
  #i <- 1; j <- 1; k <- 1
  if(Distributions[i]=="Normal"){
    Type <- "Normal"; Skew <- "None"; a <- b <- NULL
  } else if(Distributions[i]=="Right-Skew"){
    Type <- "Chi"; a <- 4; Skew <- "Right"; b <- NULL
  }
  for(j in 1:length(MAF)){
    print(paste("Analysis Started at: ", Sys.time(), sep=""))
    z.ptm <- proc.time()
    for(k in 1:length(v.G)){
      l <- length(v.G)*length(MAF)*(i-1)+(j-1)*length(v.G)+k
      t <- cbind(Varying,MAF=MAF[j],v.E=v.E,v.G=v.G[k],v.GxE,N,b0,
                 Special,Interaction.Dir,hwe.p,min.grp.size,Pare.Encoding,Type,a,b,Skew,
                 se,Adjusted.for,No.taus)

      arguments[l,match(colnames(t),colnames(arguments))] <- t
      est[[l]] <- foreach(z=1:nodes, .combine="rbind",.multicombine=TRUE, .inorder=FALSE
      ) %dopar% {
        library(data.table)
        #library(rmutil); library(quantreg); library(sn); library(pracma);
        z.Result <- matrix(NA,batches,23,
                           dimnames=list(
                             1:batches,
                             c("b1.G", "b2.E","b3.GxE","GxE.LN_p.value" ,"LN_Joint",
                               "G.QREG_p.value","QREG.Meta_Beta","QREG.Meta_SE","QREG.Meta_p.value",
                               "QREG.MetaMarg_p.value","QREG_Joint","QREG.Time",
                               "G.LN_p.value","Levene_p.value","Levene_Joint",
                               "Levene.Time","G.RIF_p.value","RIF.Meta_Beta","RIF.Meta_SE",
                               "RIF.Meta_p.value","RIF.MetaMarg_p.value","RIF_Joint","RIF.Time")))
        for(h in 1:batches){
          #h <- 1
          z.Result[h,] <- Sim.function(Arguments=arguments[l,])
        }
        z.Result <- cbind(z.Result,matrix(rep(arguments[l,],batches),batches,
                                          length(arguments[l,]),byrow=TRUE,
                                          dimnames=list(1:batches,names(arguments[l,]))))
        return(z.Result)
      }
      print(paste(v.G[k], " - finished at ", Sys.time(), sep=""))
    }
    z.etm <- proc.time()-z.ptm
    Run.time <- c(Run.time, z.etm["elapsed"])
    if(l/length(v.G)<nrow(arguments)/length(v.G)){
      TtC <- round(range(Run.time)*((nrow(arguments)-l)/length(v.G)))
    } else {
      TtC <- round(range(Run.time))
    }
    print(paste("Analysis: ", l/length(v.G), " of ", nrow(arguments)/length(v.G),
                ". Completed in ", seconds_to_period(as.numeric(round(z.etm["elapsed"]))),
                ". Approx. Time to Completion: ",  seconds_to_period(as.numeric(TtC[1])),
                " - ",seconds_to_period(as.numeric(TtC[2])), sep=""))
  }
}
stopCluster(cl)

p <- "/home/alabada/pgpstore/R_Studio_Projects/Simulations_Quantile_Regression"
list.files(p)
save.image(file=file.path(p,"Results","RData_Objects", "unAdj_vG_FalsePostives_Apr_1_2017.RData"))
list.files(file.path(p,"Results","RData_Objects"))


# ============= Running Simulations: Skeweness ============
# Conditions to test
Special <- "None"; N <- 10000;
b0 <- 0; Interaction.Dir <- "+ve"; hwe.p <- 1; min.grp.size <- 3; Pare.Encoding <- TRUE;
Type <- "Normal"; a <- NULL ; b <- NULL; Skew <- "None"; tau.thresh <- NULL; bT.G <- NULL;
se <- "iid"; Adjusted.for <- NULL; No.taus <- 19

Varying <- "v.GxE"
v.G <- 0.002
v.E <- 0.25
v.GxE <- seq(0,0.01, by=0.001)
MAF <- 0.05
Distributions <- c("Normal","Left-Skew","Right-Skew","Short-Tail","Long-Tail",
                   "Right-Outliers","Both-Outliers")

# Setup Results organization
arguments<-matrix(NA, nrow=length(MAF)*length(Distributions)*length(v.GxE), ncol=19)
colnames(arguments) <- c("Special", "N", "MAF", "Varying","b0", "v.G", "v.E", "v.GxE",
                         "Interaction.Dir", "hwe.p", "min.grp.size","Pare.Encoding", "Type",
                         "a", "b", "Skew", "se","Adjusted.for","No.taus")
est <- list()

# Specity Simulations and Computations
R <- 1
nodes <- 1
batches <- R/nodes
Run.time <- NULL

#cl<-makeCluster(nodes)
#registerDoParallel(cl)
for(i in 1:length(MAF)){
  #i <- 1; j <- 1; k <- 1
  for(j in 1:length(Distributions)){
    print(paste("Analysis Started at: ", Sys.time(), sep=""))
    z.ptm <- proc.time()
    if(Distributions[j]=="Normal"){
      Type <- "Normal"; Skew <- "None"; a <- b <- NULL
    } else if(Distributions[j]=="Left-Skew"){
      Type <- "Chi"; a <- 4; Skew <- "Left"; b <- NULL
    } else if(Distributions[j]=="Right-Skew"){
      Type <- "Chi"; a <- 4; Skew <- "Right"; b <- NULL
    } else if(Distributions[j]=="Short-Tail"){
      Type <- "Uniform"; Skew <- "None"; a <- b <- NULL
    } else if(Distributions[j]=="Long-Tail"){
      Type <- "Laplace"; Skew <- "None"; a <- b <- NULL
    } else if(Distributions[j]=="Right-Outliers"){
      Type <- "Right-Outliers"; a <- 0.01; Skew <- "None"; b <- NULL
    } else if(Distributions[j]=="Both-Outliers"){
      Type <- "Both-Outliers"; a <- 0.02; Skew <- "None"; b <- NULL
    }
    for(k in 1:length(v.GxE)){
      l <- length(v.GxE)*length(Distributions)*(i-1)+(j-1)*length(v.GxE)+k
      t <- cbind(Varying,v.E=v.E,MAF=MAF[i],v.GxE=v.GxE[k],N,b0,v.G,
                 Special,Interaction.Dir,hwe.p,min.grp.size,Pare.Encoding,Type,a,b,Skew,
                 se,Adjusted.for,No.taus)
      arguments[l,match(colnames(t),colnames(arguments))] <- t

      est[[l]] <- foreach(z=1:nodes, .combine="rbind",.multicombine=TRUE, .inorder=FALSE
      ) %dopar% {
        library(data.table); library(sn); library(rmutil)
        #library(quantreg); library(pracma);
        z.Result <- matrix(NA,batches,16,
                           dimnames=list(
                             1:batches,
                             c("b1.G", "b2.E","b3.GxE","GxE.LN_p.value" ,"LN_Joint",
                               "G.LN_p.value","Levene_p.value","Levene_Joint",
                               "Levene.Time","G.RIF_p.value","RIF.Meta_Beta","RIF.Meta_SE",
                               "RIF.Meta_p.value","RIF.MetaMarg_p.value","RIF_Joint","RIF.Time")))
        for(h in 1:batches){
          #h <- 1
          z.Result[h,] <- Sim.function(Arguments=arguments[l,])
        }
        z.Result <- cbind(z.Result,matrix(rep(arguments[l,],batches),batches,
                                          length(arguments[l,]),byrow=TRUE,
                                          dimnames=list(1:batches,names(arguments[l,]))))
        return(z.Result)
      }
      print(paste(v.GxE[k], " - finished at ", Sys.time(), sep=""))
    }
    z.etm <- proc.time()-z.ptm
    Run.time <- c(Run.time, z.etm["elapsed"])
    if(l/length(v.GxE)<nrow(arguments)/length(v.GxE)){
      TtC <- round(range(Run.time)*((nrow(arguments)-l)/length(v.GxE)))
    } else {
      TtC <- round(range(Run.time))
    }
    print(paste("Analysis: ", l/length(v.GxE), " of ", nrow(arguments)/length(v.GxE),
                ". Completed in ", seconds_to_period(as.numeric(round(z.etm["elapsed"]))),
                ". Approx. Time to Completion: ",  seconds_to_period(as.numeric(TtC[1])),
                " - ",seconds_to_period(as.numeric(TtC[2])), sep=""))
  }
}
v.GxE_arguments <- arguments
v.GxE_est <- est
rm(arguments,est);gc()
stopCluster(cl)

p <- "/home/alabada/pgpstore/R_Studio_Projects/Simulations_Quantile_Regression"
save.image(file=file.path(p,"Results","RData_Objects", "Skew_SIM_Apr_1_2017.RData"))
list.files(file.path(p,"Results","RData_Objects"))
save.image(file="Skew_SIM_Apr_1_2017.RData")


