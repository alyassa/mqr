
error_generator <- function(N=10000,Type="Normal",a,b,Skew.Dir="Right"){
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
