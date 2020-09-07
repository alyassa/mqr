######### Compute FP Rate
rm(list=ls()); gc()
library(data.table); library(quantreg)
BootQR_fun<-function(DT, taus, method="pfn", B=10, Cores=4){
  # , seed=1
  # set.seed(seed)
  #Boot.Betas<-matrix(NA, ncol=length(taus), nrow=B)
  # Boot.Samp<-matrix(sample(1:n, size=n*B, replace=TRUE), ncol=B, nrow=n)

  require(foreach)
  require(doParallel)
  cl <- makeCluster(Cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl) )
  foreach(i = 1:B, .combine = "rbind", .inorder = FALSE, .multicombine = TRUE,
          .packages = "quantreg") %dopar% {
            QR.i<-quantreg::rq(y~x,tau=taus, data=DT[sample(1:nrow(DT), replace=TRUE),],method=method)
            return(coef(QR.i)[2,])
          }
}
MR_fun<-function(Beta, COV, taus){
  # COV=cov(Boot.Betas)
  # Beta=colMeans(Boot.Betas)
  m <- length(taus)
  A <- cbind(1, taus)
  k <- ncol(A)
  SigmaInv <- solve(COV)
  Beta.hat <- solve(t(A)%*%SigmaInv%*%A)%*%t(A)%*%SigmaInv%*%as.matrix(Beta)
  Beta.hat.var <- solve(t(A)%*%solve(COV)%*%A)
  Beta.hat.SE <- sqrt(diag(Beta.hat.var))
  tval <- Beta.hat/Beta.hat.SE
  pval<-2*pnorm(abs(tval), lower=FALSE)

  Result <- cbind(Beta.hat, Beta.hat.SE, tval, pval)
  dimnames(Result) <- list(c("(Intercept)","taus"),
                           c("MetaReg_Beta","MetaReg_SE","MetaReg_tval",
                             "MetaReg_p.value"))
  #rm(m,A,SigmaInv,Beta.hat.var,Beta.hat.SE,tval,pval)
  return(Result)
}

#MR_fun(Beta=colMeans(Boot.Betas), COV=cov(Boot.Betas), taus)
n<-1000; beta<-1
taus=seq(0.1,0.9,by=0.1)
Cores=6#; seed=1

R<-500; B<-1000
Results<-matrix(NA, ncol=4, nrow=R)
pb <- txtProgressBar(min = 0, max = R, style = 3)
for(i in 1:R){
  x<-rnorm(n)
  err<-rnorm(n)
  y<-beta*x + err
  DT<-data.table(y=y,x=x)

  Boot.Betas<-BootQR_fun(DT=DT, taus=taus, method="pfn", B=B, Cores=Cores) #, seed=seed)

  Results[i,]<-MR_fun(Beta=colMeans(Boot.Betas), COV=cov(Boot.Betas), taus)[2,]
  setTxtProgressBar(pb, i)
}

mean(Results[,4]<0.05)
saveRDS(Results, "results/BootCQR_FP.RDS")
# [1] 0.08 n<-1000; R<-100; B<-1000
