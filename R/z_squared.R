z_squared <- function(datatable, y, g,covariates=NULL){
  #datatable=Data; tau=Taus; y="BMI"; g=SNPs[j]; Adjusted.for=NA
  #covariates=c(covariates,"STUDY");
  #fun1=RIF.Transformation; fun2=cov.RIF;fun3=MetaReg.function

  if(is.null(covariates)){
    FML  <- as.formula(paste(y,g,sep=" ~ "))
  } else {
    covs <- paste(covariates,sep="",collapse=" + ")
    FML  <- as.formula(paste(paste(y,g,sep=" ~ "),covs,sep=" + "))
  }
  DT <- data.table(model.frame(FML, datatable))

  DT[,z2:=(scale(DT[[y]]))^2]
  FML <- as.formula(paste("z2",as.character(FML)[1],as.character(FML)[3],
                          sep="",collapse=""))

  # Add z2 Regression Estimates
  LN.sum <- coef(summary(lm(FML,DT)))
  Result <- c(z2_Beta=LN.sum[g,1],z2_SE=LN.sum[g,2],
              z2_tval=LN.sum[g,3],z2_p.value=LN.sum[g,4])
  return(Result)
}
