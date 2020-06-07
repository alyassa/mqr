
stouffer.test <- function(x,Statistic="z",Weights,log_p.value=FALSE) {
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
