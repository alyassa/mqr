scaling.adjustor <- function(FML,data,weights,model,response,predictor,
                             method="Median_rq"){
  # data=Results
  # browser()
  DT <- copy(data)
  if(missing(model)){
    DT[,w:=1/(DT[[weights]]^2)]
    if(method=="Iteratively_ReWeighted"){
      Model <- MASS::rlm(FML,data=DT,weights=w,method="M",wt.method="inv.var")
    } else if(method=="Median_rq"){
      Model <- quantreg::rq(FML,tau=0.5,data=DT,weights=w,method="fn")
    } else if(method=="ols"){
      Model <- lm(FML,data=DT)
    }
    y.name <- colnames(Model$model)[1]
    y.se <- gsub("Beta","SE",y.name)
    x.name <- colnames(Model$model)[2]
    x.se <- gsub("Beta","SE",x.name)

    Pred <- Model$model[,x.name]*Model$coefficients[[2]] + Model$coefficients[[1]]
    Pred.se <- DT[[x.se]]*Model$coefficients[[2]] + Model$coefficients[[1]]
    Beta <- Model$residuals
    se <- sqrt(DT[[y.se]]^2 + Pred.se^2)
    tval <- Beta/se
    p.value <- 2*pnorm(abs(tval),lower=FALSE)
    Results <- cbind(Pred,Pred.se,Beta,se,tval,p.value)
    r.name <- gsub("Beta","",y.name)
    colnames(Results) <- c(paste0("Pred.",r.name,"Beta"),paste0("Pred.",r.name,"SE"),
                           paste0("Adj.",r.name,"Beta"),paste0("Adj.",r.name,"SE"),
                           paste0("Adj.",r.name,"tval"),paste0("Adj.",r.name,"p.value"))
    Results <- data.table(Results)
    # Model$model <- Model$residuals <- Model$weights <- NULL
    return(list(Model,Results))
  } else {
    Model <- copy(model)
    x.se <- gsub("Beta","SE",predictor)
    y.se <- gsub("Beta","SE",response)
    Pred <- DT[[predictor]]*Model$coefficients[[2]] + Model$coefficients[[1]]
    Pred.se <- DT[[x.se]]*Model$coefficients[[2]] + Model$coefficients[[1]]
    Beta <- DT[[response]]-Pred
    se <- sqrt(DT[[y.se]]^2 + Pred.se^2)
    tval <- Beta/se
    p.value <- 2*pnorm(abs(tval),lower=FALSE)
    Results <- cbind(Pred,Pred.se,Beta,se,tval,p.value)
    r.name <- gsub("Beta","",response)
    colnames(Results) <- c(paste0("Pred.",r.name,"Beta"),paste0("Pred.",r.name,"SE"),
                           paste0("Adj.",r.name,"Beta"),paste0("Adj.",r.name,"SE"),
                           paste0("Adj.",r.name,"tval"),paste0("Adj.",r.name,"p.value"))
    Results <- data.table(Results)
    return(Results)
  }
}
