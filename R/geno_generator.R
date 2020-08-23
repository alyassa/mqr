geno_generator <- function(N=10000,MAF=0.05,hwe.p=1,min.grp.size=5,
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
    obs.G <- sample(c(aa, ab, bb), N, prob=c(qq,pq,pp)/N, replace=TRUE)
  } else {
    # If hwe.p is specified then need to work backwards, from hwe.p to get chi value
    # at that hwe.p
    hwe.chi <- qchisq(hwe.p,1,lower.tail=FALSE)
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
    obs.G <- sample(c(aa, ab, bb), N, prob=c(obs.qq,obs.pq,obs.pp)/N, replace=TRUE)
    rm(hwe.chi,obs.pp,obs.qq,obs.pq)
  }
  if(length(obs.G)<N){
    obs.G <- c(obs.G,rep(NA, N-length(obs.G)))
  }
  return(obs.G)
}
