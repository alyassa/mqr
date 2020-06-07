fischer.test <- function(x,Statistic="z",log_p.value=FALSE) {
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
