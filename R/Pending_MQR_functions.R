## ============= MQR FUNCTIONS ====================
MetaReg.CQR <- function(Beta,boot.Bs,test="z"){
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
RIF.Transform <- function(y,taus){
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
MetaReg.UQR <- function(Beta, COV, taus,test="z"){
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
