# Function to extract covariance matrix from multivariate OLS
RIFmod.cov <- function(object, pred){
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
