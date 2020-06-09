rif <- function(y,taus){
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
