# ============= FUNCTIONS AND PACKAGES ====================
library(foreach); library(doParallel); library(sn); library(lubridate); library(rmutil)

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

gamma_parameter_sampling <- function(y){
  gamma_mean <- ceiling(abs(min(y)))
  s <- (1-var(y))*exp(y)/(gamma_mean+y)
  y <- sapply(1:length(y),function(x) rgamma(1,shape=(gamma_mean+y[x])/s[x],scale=s[x]))
  return(y)
}

Sim.function <- function(Arguments,Model,Mode="Scaling.Model",add.GxE=FALSE,
                         S.fun1=Generate_Genotypes,
                         S.fun2=Error.Distribution,
                         S.fun3=rntransform,
                         S.fun4=MCQR.function,
                         S.fun7=MetaReg.MCQR.function,
                         S.fun8=Scaling.Adjustor,
                         S.fun9=gamma_parameter_sampling){
  # S.fun1=Generate_Genotypes;S.fun2=Error.Distribution;S.fun3=rntransform;
  # S.fun4=MCQR.function;
  # S.fun7=MetaReg.MCQR.function;S.fun8=Scaling.Adjustor
  # S.fun9=gamma_parameter_sampling
  # add.GxE=TRUE; Mode="Scaling.Adjustment"; #Mode="Scaling.Model.Fitting";
  # Arguments=c(j.args[1,],Seed=j.s)
  if(!missing(Arguments)){
    N <- as.numeric(Arguments["N"]);
    MAF <- as.numeric(Arguments["MAF"]);
    b0 <- as.numeric(Arguments["b0"]);
    v.G_Buffer <- as.numeric(Arguments["v.G_Buffer"]);
    v.G <- as.numeric(Arguments["v.G"]);
    v.E <- as.numeric(Arguments["v.E"]);
    v.GxE <- as.numeric(Arguments["v.GxE"]);
    Interaction.Dir <- Arguments["Interaction.Dir"];
    hwe.p <- as.numeric(Arguments["hwe.p"]);
    min.grp.size <- as.numeric(Arguments["min.grp.size"]);
    Pare.Encoding <- Arguments["Pare.Encoding"];
    Type <- Arguments["Type"];
    a <- as.numeric(Arguments["a"]); b <- as.numeric(Arguments["b"]);
    Skew.Dir <- Arguments["Skew.Dir"];
    gamma_0 <- as.numeric(Arguments["gamma_0"]);
    gamma_1 <- as.numeric(Arguments["gamma_1"]);
    Scaling.Method <- Arguments["Scaling.Method"]
    lambda <- as.numeric(Arguments["lambda"]);
    alpha <- as.numeric(Arguments["alpha"]);
    v.ParTot <- as.numeric(Arguments["v.ParTot"]);
    Rank.Multiplier <- as.numeric(Arguments["Rank.Multiplier"]);
    Seed <- as.numeric(Arguments["Seed"]);
    Adjusted.for <- Arguments["Adjusted.for"];
    No.taus <- as.numeric(Arguments["No.taus"]);
    Scale <- Arguments["Scale"];
    method <- Arguments["method"]
  }
  if(!is.na(Seed)){
    set.seed(Seed)
  }
  # browser()
  b1.G <- sample(c(-1,1),1,prob=c(0.5,0.5))*sqrt(v.G/(2*MAF*(1-MAF)))
  b2.E <- sqrt(v.E)
  b3.GxE <- sign(b1.G)*sqrt(v.GxE/(2*MAF*(1-MAF)))
  b4.Buffer <- sqrt(v.G_Buffer-v.G*(1+gamma_1)-v.GxE)
  v.combo <- v.G_Buffer + v.E
  if(is.na(v.ParTot)){
    v.ParTot <- 1
  }
  if(v.combo>1){
    stop("Variance components add up to more than 1 ")
  }
  if(v.ParTot-v.combo<0){
    stop("Total variance of parameters (v.ParTot) is less than sum of parameters variances (v.G, v.E, v.GxE). Increase v.ParTot so v.ParTot >= v.G + v.E + v.GxE")
  }
  G <- S.fun1(N,MAF,hwe.p,min.grp.size,Pare.Encoding)
  E <- rnorm(N)
  err <- S.fun2(N,Type,a,b,Skew.Dir)
  C <- rnorm(N)
  if(Interaction.Dir=="+ve"){
    b3.Dir <- 1
  } else if (Interaction.Dir=="-ve"){
    b3.Dir <- -1
  }

  pre_y <- b1.G*G + b2.E*E + b3.Dir*b3.GxE*G*E + (gamma_0 + gamma_1*b1.G)*G*err +
    sqrt(v.ParTot-v.combo)*err + b4.Buffer*C

  # mean(pre_y)
  # var(pre_y)
  # hist(pre_y,breaks=100,prob=TRUE,main="",xlab="Outcome",col="grey",border=NA)
  # curve(dnorm(x, mean=mean(pre_y), sd=sd(pre_y)),add=TRUE, col="red")
  # qqnorm(pre_y);qqline(pre_y,col="red")
  DT <- cbind(pre_y,G,E)
  DT <- data.table(DT)
  setnames(DT,c("pre_y","G","E"))
  if(Scaling.Method=="None"){
    DT[,y:=b0 + pre_y]
  } else if(Scaling.Method=="Rank.Multiplier"){
    DT <- DT[order(pre_y)]
    Multiplier <- seq(1,1+Rank.Multiplier,length.out=N)
    DT[,Rank_Multiplier:=Multiplier[order(Multiplier)]]
    DT[,y:=pre_y*Rank_Multiplier]
    DT[,y:=b0 + scale(y)]
  } else if (Scaling.Method=="EGM.Tranformation"){
    if(lambda!=0){
      remgtransform <- function(x,lambda){qemg((rank(x,na.last="keep")-0.5)/sum(!is.na(x)),
                                               mu=0,sigma=1,lambda=lambda)}
      DT[,y:=remgtransform(pre_y,lambda=lambda)]
    } else {
      DT[,y:=pre_y]
    }
    DT[,y:=b0 + scale(y)]
  } else if (Scaling.Method=="SN.Tranformation"){
    rskewtransform <- function(x,alpha){qsn((rank(x,na.last="keep")-0.5)/sum(!is.na(x)),
                                            xi=0,omega=1,alpha=alpha)}
    sn_y <- rskewtransform(y,alpha=4)
    DT[,y:=rskewnransform(pre_y,alpha=alpha)]
    DT[,y:=b0 + scale(y)]
  } else if(Scaling.Method=="gamma_parameter"){
    DT[,y:=S.fun9(pre_y)]
    DT[,y:=b0 + scale(y)]
  }

  # hist(DT[,y],breaks=100,prob=TRUE,main="",xlab="Outcome",col="grey",border=NA)
  # curve(dnorm(x, mean=mean(DT[,y]), sd=sd(DT[,y])),add=TRUE, col="red")
  # qqnorm(DT[,y]);qqline(DT[,y],col="red")
  if(Scale=="log"){
    DT[,sc_pre_y:=log(pre_y + b0)]
    DT[,sc_y:=log(y)]
  } else if(Scale=="RT"){
    DT[,sc_pre_y:=S.fun3(pre_y)]
    DT[,sc_y:=S.fun3(y)]
  } else {
    DT[,sc_pre_y:=pre_y]
    DT[,sc_y:=y]
  }

  DT[,pheno:=sc_y]
  DT[,geno:=G]
  if(Mode=="Data.Generator"){
    return(DT)
  }

  S.fun4=MCQR.function,
  S.fun7=MetaReg.MCQR.function,
  S.fun8=Scaling.Adjustor,
  Taus <- seq(0.05,0.95,length.out=No.taus)
  MCQR.Result <- muqr(DT,tau=Taus,y="pheno", g="geno",cfun1=S.fun7)
  MCQR.Result <- c(SNP_v.G=v.G,SNP_MAF=MAF,MCQR.Result)
  # browser()
  if(add.GxE){
    lm1 <- coef(summary(lm(pheno~geno*E,data=DT)))
    MCQR.Result <- c(MCQR.Result,GxE_Beta=lm1["geno:E",1],GxE_SE=lm1["geno:E",2],
                     GxE_tval=lm1["geno:E",3],GxE_p.value=lm1["geno:E",4])
  }

  if(Mode=="Scaling.Model.Fitting"){
    return(MCQR.Result)
  } else if(Mode=="Scaling.Adjustment"){
    return(c(MCQR.Result,
             as.matrix(S.fun8(data=data.table(t(MCQR.Result)),model=Model,
                              response="MCQR.MetaTau_Beta",
                              predictor="MCQR.Median_Beta",method=method))[1,]))
  }
}

Scaling.Estimator <- function(Arguments,return.full=TRUE,add.GxE=FALSE,
                                 D.fun1=Generate_Genotypes,
                                 D.fun2=Error.Distribution,
                                 D.fun3=rntransform,
                                 D.fun4=MCQR.function,
                                 D.fun7=MetaReg.MCQR.function,
                                 D.fun8=Sim.function,
                                 D.fun9=Scaling.Adjustor,
                                 D.fun10=gamma_parameter_sampling){
  # D.fun1=Generate_Genotypes;D.fun2=Error.Distribution;D.fun3=rntransform;
  # D.fun4=MCQR.function;
  # D.fun7=MetaReg.MCQR.function;D.fun8=Sim.function;D.fun9=Scaling.Adjustor
  # D.fun10=gamma_parameter_sampling
  # add.GxE=FALSE
  # Arguments=c(i.args[1,],Seed=Seeding1)
  if(!missing(Arguments)){
    Sc.R <- as.numeric(Arguments["Sc.R"]);
    v.G_Dist <- Arguments["v.G_Dist"];
    v.G_DistMean <- as.numeric(Arguments["v.G_DistMean"])
    v.G_DistQuant <- as.numeric(Arguments["v.G_DistQuant"])
    nodes <- as.numeric(Arguments["nodes"]);
    N <- as.numeric(Arguments["N"]); b0 <- as.numeric(Arguments["b0"]);
    v.G_Buffer <- as.numeric(Arguments["v.G_Buffer"]);
    v.G <- as.numeric(Arguments["v.G"]);
    v.E <- as.numeric(Arguments["v.E"]);
    v.GxE <- as.numeric(Arguments["v.GxE"]);
    Interaction.Dir <- Arguments["Interaction.Dir"];
    hwe.p <- as.numeric(Arguments["hwe.p"]); min.grp.size <- as.numeric(Arguments["min.grp.size"]);
    Pare.Encoding <- Arguments["Pare.Encoding"]; Type <- Arguments["Type"];
    a <- as.numeric(Arguments["a"]); b <- as.numeric(Arguments["b"]);
    Skew.Dir <- Arguments["Skew.Dir"];
    Scaling.Method <- Arguments["Scaling.Method"]
    lambda <- as.numeric(Arguments["lambda"]);
    alpha <- as.numeric(Arguments["alpha"]);
    Rank.Multiplier <- as.numeric(Arguments["Rank.Multiplier"]);
    gamma_0 <- as.numeric(Arguments["gamma_0"]);
    gamma_1 <- as.numeric(Arguments["gamma_1"]);
    v.ParTot <- as.numeric(Arguments["v.ParTot"]);
    Adjusted.for <- Arguments["Adjusted.for"]; No.taus <- as.numeric(Arguments["No.taus"]);
    Scale <- Arguments["Scale"]
    Seed <- as.numeric(Arguments["Seed"]);
    method <- Arguments["method"]
  }
  ptm <- proc.time()
  if(!is.na(Seed)){
    set.seed(Seed)
  }
  # browser()
  batches <- split(1:Sc.R,cut(seq_along(1:Sc.R),nodes,labels=FALSE))
  MAF <- runif(Sc.R,0.05,0.5)
  if(v.G_Dist=="chisq"){
    DEM_v.G <- rchisq(Sc.R,1)*v.G_DistMean
  } else if(v.G_Dist=="lomax"){
    DEM_v.G <- rparetoII(Sc.R,location=0,scale=1,shape=2)*v.G_DistMean
  } else if(v.G_Dist=="quantile-replace"){
    DEM_v.G <- rchisq(Sc.R,1)*v.G_DistMean
    DEM_v.G[sample(1:Sc.R,v.G_DistQuant*Sc.R)] <- runif(v.G_DistQuant*Sc.R,0.003,0.006)
  }

  # DEM_v.G[sample(1:Sc.R,1)] <- v.G
  # DEM_v.G <- rchisq(Sc.R,1)/5000
  # DEM_v.G <- DEM_v.G/sum(DEM_v.G)
  # DEM_v.G[DEM_v.G>0.001]
  # DEM_v.G[DEM_v.G>0.004]
  # summary(DEM_v.G)
  if(add.GxE){
    cols <- c("SNP_v.G","SNP_MAF","LN_Beta","LN_SE","LN_tval",
              "LN_p.value","MCQR.Median_Beta","MCQR.Median_SE","MCQR.Median_tval",
              "MCQR.Median_p.value","MCQR.MetaTau_Beta","MCQR.MetaTau_SE",
              "MCQR.MetaTau_tval","MCQR.MetaTau_p.value","MCQR.TtC","Notes","GxE_Beta",
              "GxE_SE","GxE_tval","GxE_p.value")
  } else {
    cols <- c("SNP_v.G","SNP_MAF","LN_Beta","LN_SE","LN_tval",
              "LN_p.value","MCQR.Median_Beta","MCQR.Median_SE","MCQR.Median_tval",
              "MCQR.Median_p.value","MCQR.MetaTau_Beta","MCQR.MetaTau_SE",
              "MCQR.MetaTau_tval","MCQR.MetaTau_p.value","MCQR.TtC","Notes")
  }
  # print("Generating Scaling Model")
  ptm <- proc.time()
  Results <- foreach(i=1:nodes, .combine=function(...) rbindlist(list(...)),
                     .multicombine=TRUE, .inorder=FALSE
  ) %dopar% {
    library(data.table);library(sn); library(emg); library(quantreg)
    i.Result <- matrix(NA,length(batches[[i]]),length(cols),
                       dimnames=list(1:length(batches[[i]]),cols))
    for(j in 1:length(batches[[i]])){
      if(!is.na(Seed)){
        se <- Sc.R/nodes*(i-1) + j + Seed
      } else {
        se <- NA
      }
      # print(se)
      j.args <- cbind(N,b0,v.G_Buffer,
                      v.G=DEM_v.G[batches[[i]][j]],v.E,v.GxE=0,
                      MAF=MAF[batches[[i]][j]],
                      Interaction.Dir,hwe.p,min.grp.size,Pare.Encoding,Type,a,b,
                      Skew.Dir,Adjusted.for,No.taus,Scaling.Method,lambda,alpha,Rank.Multiplier,
                      gamma_0,gamma_1,v.ParTot,Scale)

      i.Result[j,] <- D.fun8(Arguments=c(j.args[1,],Seed=se),Mode="Scaling.Model.Fitting",
                             add.GxE,
                             S.fun1=D.fun1,S.fun2=D.fun2,S.fun3=D.fun3,
                             S.fun4=D.fun4,S.fun7=D.fun7,S.fun8=D.fun9,S.fun9=D.fun10)
    }
    i.Result <- data.table(i.Result)
    return(i.Result)
  }

  Sc.Results <- D.fun9(MCQR.MetaTau_Beta ~ MCQR.Median_Beta,
                       data=Results,weights="MCQR.Median_SE",method=method)
  Results <- cbind(Results,Sc.Results[[2]])
  Model <- Sc.Results[[1]]
  # Model$model <- Model$residuals <- Model$weights <- NULL
  etm <- proc.time() - ptm
  # print(paste0("Completed in ", seconds_to_period(round(etm[[3]]))))
  if(return.full){
    Output <- list(Model,Results)
    names(Output) <- c("Model","Results")
    return(Output)
  } else {
    return(Results[SNP_v.G==v.G])
  }
}

Testing.Sim.function <- function(Arguments,Model,add.GxE=FALSE,
                              D.fun1=Generate_Genotypes,
                              D.fun2=Error.Distribution,
                              D.fun3=rntransform,
                              D.fun4=MCQR.function,
                              D.fun7=MetaReg.MCQR.function,
                              D.fun8=Sim.function,
                              D.fun9=Scaling.Adjustor,
                              D.fun10=gamma_parameter_sampling){
  # D.fun1=Generate_Genotypes;D.fun2=Error.Distribution;D.fun3=rntransform;
  # D.fun4=MCQR.function;D.fun7=MetaReg.MCQR.function;D.fun8=Sim.function;
  # D.fun9=Scaling.Adjustor;D.fun10=gamma_parameter_sampling
  # add.GxE=TRUE;Model=i.ScalingModel$Model
  # Arguments=c(j.args[1,],Seed=j.s)
  if(!missing(Arguments)){
    R <- as.numeric(Arguments["R"])
    nodes <- as.numeric(Arguments["nodes"]);
    Seed <- as.numeric(Arguments["Seed"]);
  }
  if(!is.na(Seed)){
    set.seed(Seed)
  }
  batches <- split(1:R,cut(seq_along(1:R),nodes,labels=FALSE))
  # browser()
  if(add.GxE){
    cols <- c("SNP_v.G","SNP_MAF","LN_Beta","LN_SE","LN_tval",
              "LN_p.value","MCQR.Median_Beta","MCQR.Median_SE","MCQR.Median_tval",
              "MCQR.Median_p.value","MCQR.MetaTau_Beta","MCQR.MetaTau_SE",
              "MCQR.MetaTau_tval","MCQR.MetaTau_p.value","MCQR.TtC","Notes","GxE_Beta",
              "GxE_SE","GxE_tval","GxE_p.value","Pred.MCQR.MetaTau_Beta",
              "Pred.MCQR.MetaTau_SE","Adj.MCQR.MetaTau_Beta",
              "Adj.MCQR.MetaTau_SE","Adj.MCQR.MetaTau_tval",
              "Adj.MCQR.MetaTau_p.value")

  } else {
    cols <- c("SNP_v.G","SNP_MAF","LN_Beta","LN_SE","LN_tval",
              "LN_p.value","MCQR.Median_Beta","MCQR.Median_SE","MCQR.Median_tval",
              "MCQR.Median_p.value","MCQR.MetaTau_Beta","MCQR.MetaTau_SE",
              "MCQR.MetaTau_tval","MCQR.MetaTau_p.value","MCQR.TtC","Notes",
              "Pred.MCQR.MetaTau_Beta","Pred.MCQR.MetaTau_SE",
              "Adj.MCQR.MetaTau_Beta","Adj.MCQR.MetaTau_SE","Adj.MCQR.MetaTau_tval",
              "Adj.MCQR.MetaTau_p.value")
  }
  Results <- foreach(i=1:nodes, .combine=function(...) rbindlist(list(...)),
                     .multicombine=TRUE, .inorder=FALSE
  ) %dopar% {
    library(data.table); library(MASS); library(quantreg); library(sn);
    library(emg); library(quantreg)
    i.Result <- matrix(NA,length(batches[[i]]),length(cols),
                       dimnames=list(1:length(batches[[i]]),cols))
    for(j in 1:length(batches[[i]])){
      if(!is.na(Seed)){
        se <- R/nodes*(i-1) + j + Seed
      } else {
        se <- NA
      }
      # print(se)
      j.args <- copy(Arguments)
      j.args["Seed"] <- as.character(se)
      i.Result[j,] <- D.fun8(Arguments=j.args,Model=Model,Mode="Scaling.Adjustment",
                             add.GxE,
                             S.fun1=D.fun1,S.fun2=D.fun2,S.fun3=D.fun3,
                             S.fun4=D.fun4,S.fun7=D.fun7,S.fun8=D.fun9,S.fun9=D.fun10)
    }
    i.Result <- data.table(i.Result)
    return(i.Result)
  }
  return(Results)
}
