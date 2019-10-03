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
