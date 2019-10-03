Levene <- function(datatable, y, g, covariates=NULL, method="Brown-Forsythe",
                   log.p=FALSE){
  #datatable=Data; method="Brown-Forsythe"
  if(is.null(covariates)){
    FML  <- as.formula(paste(y,g,sep=" ~ "))
  } else {
    covs <- paste(covariates,sep="",collapse=" + ")
    FML  <- as.formula(paste(paste(y,g,sep=" ~ "),covs,sep=" + "))
  }

  DT <- data.table(model.frame(FML, datatable))
  DT[,g:=ifelse(DT[[g]]<=0.5,0,
                ifelse(DT[[g]]<=1.5,1,
                       ifelse(DT[[g]]>1.5,2, NA)))]
  if(is.null(covariates)){
    DT[,y:=DT[[y]]]
  } else {
    FML.r <- paste(y, covs, sep=" ~ ")
    LMOD <- lm(FML,DT)
    DT[,y:=resid(LMOD)]
  }

  DT <- DT[,.(y,g)]

  # first calculate Zij (which is key) for each subject and add it to
  # data.table.
  # Use Median for "Brown–Forsythe" and Mean for "Classic-Levene". Otherwise tests
  # are the same.
  if(method=="Brown-Forsythe"){
    DT[,Zij:=abs(y-median(y)),by=g]
  } else if(method=="Classic-Levene"){
    DT[,Zij:=abs(y-mean(y)),by=g]
  }

  # calculate Z..
  Z.. <- mean(DT[,Zij])
  # find the group number
  k <- length(unique(DT[,g]))
  # find the sample size
  N <- nrow(DT)
  # calculate Zij_Zi, which is the difference between Zij (for each
  # subject) and the mean of Zij for each genotype, squared
  DT[,Zij_Zi:=(Zij-mean(Zij))^2,by=g]
  W <- (N-k)*sum(DT[,length(y)*(mean(Zij)-Z..)^2,by=g][,V1])/
    ((k-1)*sum(DT[,Zij_Zi]))
  if(log.p){
    p <- pf(W,df1=k-1,df2=N-k, lower.tail=FALSE, log.p=TRUE)
  } else {
    p <- pf(W,df1=k-1,df2=N-k, lower.tail=FALSE)
  }
  Result <- c(Levene_Fval=W,Levene_df1=k-1,Levene_df2=N-k,Levene_p.value=p)
  return(Result)
}
