# ============= FUNCTIONS AND PACKAGES ====================
library(foreach); library(lubridate); library(doParallel);

# So right now i've settelled on 6 Different methods to skew the distribution still
# producing mean of zero and var of 1.
Error.Distribution <- function(N,Type="Normal",a,b,Skew.Dir="Right"){
  #Type="Right-Outliers"; a=0.02;N=10000
  # METHOD 0 : Normal (not skewed)
  #            Type="Normal"
  # METHOD 1 : uses the Skew-N Distribution. This is pretty good but may not produce
  #            ideal shape of skewing. Increasing the term "alpha" term increases
  #            skewing. Dont bother with alpha < 15 because skewing is not
  #            substantial enough.
  #            Type="Skew-N"
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
  # METHOD 6 : Exponential distribution, needed for more extreme skewing. a <=1
  #            Type="Exp"
  ifelse(Skew.Dir=="Left",{dir <- -1}, {dir <- 1})
  if(Type=="Normal"){
    x <- rnorm(N)
    return(x)
  }
  if(Type=="Chi"){
    x <- dir*rchisq(N,df=a)
    x <- scale(x)
    return(x)
  }
  if(Type=="Laplace"){
    x <- rlaplace(N)
    x <- scale(x)
    return(x)
  }
  if(Type=="Uniform"){
    x <- runif(N)
    x <- scale(x)
    return(x)
  }
  if(Type=="Skew-N"){
    # WAS RST!!! CHANGED TO RSN
    x <- dir*rsn(n=N,xi=0,omega=1,alpha=a)
    x <- scale(x)
    return(x)
  }
  if(Type=="LogNormal"){
    x <- dir*rlnorm(N)
    x <- scale(x)
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
    x <- scale(x)
    return(x)
  }
  if(Type=="Exp"){
    x <- dir*rexp(N,a)
    x <- scale(x)
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

Sim.function <- function(Arguments,add.GxE=FALSE,
                         S.fun1=Generate_Genotypes,
                         S.fun2=Error.Distribution,
                         S.fun3=rntransform,
                         S.fun4=MUQR.function,
                         S.fun5=RIF.Transformation,
                         S.fun6=cov.RIF,
                         S.fun7=MetaReg.MUQR.function,
                         S.fun8=MCQR.function,
                         S.fun9=MetaReg.MCQR.function){
  # S.fun1=Generate_Genotypes;S.fun2=Error.Distribution;S.fun3=rntransform;
  # S.fun4=MUQR.function;S.fun5=RIF.Transformation;S.fun6=cov.RIF;
  # S.fun7=MetaReg.MUQR.function;S.fun8=MCQR.function;
  # S.fun9=MetaReg.MCQR.function
  # Arguments=c(i.args[1,],Seed=l.s)
  if(!missing(Arguments)){
    N <- as.numeric(Arguments["N"]);
    MAF <- as.numeric(Arguments["MAF"]);
    b0 <- as.numeric(Arguments["b0"]);
    v.G <- as.numeric(Arguments["v.G"]);
    v.E <- as.numeric(Arguments["v.E"]);
    v.GxE <- as.numeric(Arguments["v.GxE"]);
    v.G_Buffer <- as.numeric(Arguments["v.G_Buffer"]);
    Interaction.Dir <- Arguments["Interaction.Dir"];
    hwe.p <- as.numeric(Arguments["hwe.p"]);
    min.grp.size <- as.numeric(Arguments["min.grp.size"]);
    Pare.Encoding <- Arguments["Pare.Encoding"];
    Type <- Arguments["Type"];
    a <- as.numeric(Arguments["a"]); b <- as.numeric(Arguments["b"]);
    Skew.Dir <- Arguments["Skew.Dir"];
    Seed <- as.numeric(Arguments["Seed"]);
    No.taus <- as.numeric(Arguments["No.taus"]);
    min.tau <- as.numeric(Arguments["min.tau"]);
    max.tau <- as.numeric(Arguments["max.tau"]);
  }
  if(!is.na(Seed)){
    set.seed(Seed)
  }
  # browser()
  b1.G <- sample(c(-1,1),1,prob=c(0.5,0.5))*sqrt(v.G/(2*MAF*(1-MAF)))
  b2.E <- sqrt(v.E)
  b3.GxE <- sign(b1.G)*sqrt(v.GxE/(2*MAF*(1-MAF)))
  b4.Buffer <- sqrt(v.G_Buffer-v.G-v.GxE)
  v.combo <- v.G_Buffer + v.E
  G <- S.fun1(N,MAF,hwe.p,min.grp.size,Pare.Encoding)
  E <- rnorm(N)
  err <- S.fun2(N,Type,a,b,Skew.Dir)
  C <- rnorm(N)
  if(Interaction.Dir=="+ve"){
    b3.Dir <- 1
  } else if (Interaction.Dir=="-ve"){
    b3.Dir <- -1
  }

  y <- b1.G*G + b2.E*E + b3.Dir*b3.GxE*G*E + sqrt(1-v.combo)*err + b4.Buffer*C
  DT <- data.table(y,G,E)
  setnames(DT,c("y","G","E"))
  DT[,pheno:=y]
  DT[,geno:=G]

  Taus <- seq(min.tau,max.tau,length.out=No.taus)
  MUQR.Result <- S.fun4(DT,tau=Taus,y="pheno", g="geno",ufun1=S.fun5,ufun2=S.fun6,ufun3=S.fun7)
  Notes <- MUQR.Result["Notes"][[1]]
  MCQR.Result <- S.fun8(DT,tau=Taus,y="pheno", g="geno",cfun1=S.fun9)
  Notes <- c(Notes,MCQR.Result["Notes"][[1]])
  Notes <- Notes[!is.na(Notes)]
  if(length(Notes)>0){
    Notes <- paste(Notes,sep="",collapse=". ")
  } else {
    Notes <- NA
  }
  MUQR.Result <- MUQR.Result[-length(MUQR.Result)]
  MCQR.Result <- MCQR.Result[-length(MCQR.Result)]
  Result <- c(SNP_v.G=v.G,SNP_MAF=MAF,MUQR.Result[3:length(MUQR.Result)],
              MCQR.Result[7:length(MCQR.Result)],Notes=Notes)

  # S.fun8(DT,tau=Taus,y="pheno", g="geno",seed=23472,fun1=S.fun9)
  # browser()
  if(add.GxE){
    lm1 <- coef(summary(lm(pheno~geno*E,data=DT)))
    Result <- c(Result,GxE_Beta=lm1["geno:E",1],GxE_SE=lm1["geno:E",2],
                GxE_tval=lm1["geno:E",3],GxE_p.value=lm1["geno:E",4])
  }
  return(Result)
}

# ============= MQR : RUNNING ANALYSIS ====================
# Conditions to test
Special <- "None"; N <- 10000; v.E <- 0.24; Interaction.Dir <- "+ve";
hwe.p <- 1; min.grp.size <- 3; Pare.Encoding <- TRUE; Adjusted.for <- NULL;
Type <- "Normal"; a <- NA  ; b <- NA; Skew.Dir <- "None"; b0 <- 25

N <- c(250,500,1000,2000,5000)
No.taus <- c(4,7,10,20)
Ranges <- rbind(c(0.02,0.98),c(0.05,0.95),c(0.1,0.9),c(0.2,0.8))
v.GxE <- seq(0,0.004,length.out=5)
v.G_Buffer <- 0.1
v.G <- 0.004
MAF <- 0.3
R <- 2000
nodes <- 40
batches <- split(1:R,cut(seq_along(1:R),nodes,labels=FALSE))
Seeding <- 3435 # if you dont want to seed, set to NA
Conditions <- data.table(expand.grid(N=N,No.taus=No.taus,Ranges=1:nrow(Ranges),v.GxE=v.GxE,
                                     stringsAsFactors=FALSE))
# Timer.Matrix <- data.table(expand.grid(N=N,No.taus=No.taus,stringsAsFactors=FALSE))
# Timer.List <- as.list(rep(NA,nrow(Timer.Matrix)))
# Timer.List <- as.list(rep(NA,length(N)*length(No.taus)))

cols <- c("SNP_v.G","SNP_MAF","LN_Beta","LN_SE","LN_tval","LN_p.value","MUQR.Median_Beta",
          "MUQR.Median_SE","MUQR.Median_tval","MUQR.Median_p.value","MUQR.MetaTau_Beta",
          "MUQR.MetaTau_SE","MUQR.MetaTau_tval","MUQR.MetaTau_p.value","MUQR.TtC",
          "MCQR.Median_Beta","MCQR.Median_SE","MCQR.Median_tval","MCQR.Median_p.value",
          "MCQR.MetaTau_Beta","MCQR.MetaTau_SE","MCQR.MetaTau_tval","MCQR.MetaTau_p.value",
          "Successful.Taus","No.Successful.Taus","rq.method","boot.method","MCQR.TtC","Notes")

Analysis.Timer <- NULL
Analysis.Results <- list()
cl <- makeCluster(nodes)
registerDoParallel(cl)
for(i in 1:nrow(Conditions)){
  # Timer.row <- Timer.Matrix[,which(Conditions[i,N]==N & Conditions[i,No.taus]==No.taus)]
  # if(length(Timer.List[[Timer.row]])==1 & is.na(Timer.List[[Timer.row]])){
  #   print(paste0("ANALYSIS OF CONDITIONS ",i, " OF ",nrow(Conditions),
  #                " STARTED AT: ", Sys.time()))
  # } else {
  #   print(paste0("ANALYSIS OF CONDITIONS ",i, " OF ",nrow(Conditions),
  #                " STARTED AT: ", Sys.time(),". ETC = ",
  #                seconds_to_period(as.numeric(round(mean(Timer.List[[Timer.row]]))))))
  # }
  if(is.null(Analysis.Results)){
    print(paste0("ANALYSIS OF CONDITIONS ",i, " OF ", nrow(Conditions),
                 " STARTED AT: ", Sys.time()))
  } else {
    print(paste0("ANALYSIS OF CONDITIONS ",i, " OF ", nrow(Conditions),
                 " STARTED AT: ", Sys.time(),". ETC = ",
                 seconds_to_period(as.numeric(round(mean(Analysis.Timer)*
                                                      (nrow(Conditions)-i+1))))," [",
                 seconds_to_period(as.numeric(round(min(Analysis.Timer)*
                                                      (nrow(Conditions)-i+1))))," - ",
                 seconds_to_period(as.numeric(round(max(Analysis.Timer)*
                                                      (nrow(Conditions)-i+1)))),"]"))
  }
  i.ptm <- proc.time()
  i.args <- cbind(N=Conditions[i,N],b0,v.G,v.E,v.GxE=Conditions[i,v.GxE],v.G_Buffer,
                  Interaction.Dir,MAF,hwe.p,min.grp.size,Pare.Encoding,Type,a,b,Skew.Dir,
                  Adjusted.for,
                  No.taus=Conditions[i,No.taus],
                  min.tau=Ranges[Conditions[i,Ranges],1],
                  max.tau=Ranges[Conditions[i,Ranges],2])
  # paste(names(Sim.function(Arguments=i.args[1,])),sep='',collapse='","')
  i.Results <- foreach(j=1:nodes, .combine=function(...) rbindlist(list(...)),
                     .multicombine=TRUE, .inorder=FALSE
  ) %dopar% {
    library(data.table);library(quantreg)
    Reps <- batches[[j]]
    j.Results <- matrix(NA,nrow=length(Reps),ncol=length(cols),
                        dimnames=list(1:length(Reps),cols))
    for(l in 1:length(Reps)){
      l.s <- nrow(Conditions)*R*(i-1) + (Reps[l]-1) + Seeding
      # print(l.s)
      j.Results[l,] <- Sim.function(Arguments=c(i.args[1,],Seed=l.s))
    }
    j.Results <- data.table(j.Results)
  }
  i.Results <- cbind(Conditions[i,],i.Results)
  Analysis.Results[[i]] <- i.Results
  i.etm <- proc.time()-i.ptm
  # if(length(Timer.List[[Timer.row]])==1 & is.na(Timer.List[[Timer.row]])){
  #   Timer.List[[Timer.row]] <- i.etm[[3]]
  # } else {
  #   Timer.List[[Timer.row]] <- c(Timer.List[[Timer.row]],i.etm[[3]])
  # }
  Analysis.Timer <- c(Analysis.Timer,i.etm[[3]])
  print(paste0("ANALYSIS OF CONDITIONS ",i, " OF ",nrow(Conditions),
               " COMPLETED IN ", seconds_to_period(round(i.etm[[3]]))))
}
stopCluster(cl)

p <- "/home/../media/StoreB/arkan/Simulations_Quantile_Regression/Results/RData_Objects"
save.image(file=file.path(p,"TestRangeNumber_Sims_July2_2019.RData"))

