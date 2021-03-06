
DT <- data.table(outcome=rnorm(1000),predictor=rnorm(1000))

args(MQR)
t <- MQR(DT,y="outcome",g="predictor",mqr.method="CQR")
t
print(paste(names(t),sep="",collapse='","'))

# MCQR cols
c("N","EAF","LN_Beta","LN_SE","LN_tval","LN_p.value","MCQR.Median_Beta","MCQR.Median_SE",
  "MCQR.Median_tval","MCQR.Median_p.value","MCQR.MetaTau_Beta","MCQR.MetaTau_SE",
  "MCQR.MetaTau_tval","MCQR.MetaTau_p.value","Successful.Taus","No.Successful.Taus",
  "rq.method","boot.method","MCQR.TtC","Notes")

# MUQR cols
c("N","EAF","LN_Beta","LN_SE","LN_tval","LN_p.value","MUQR.Median_Beta","MUQR.Median_SE",
  "MUQR.Median_tval","MUQR.Median_p.value","MUQR.MetaTau_Beta","MUQR.MetaTau_SE",
  "MUQR.MetaTau_tval","MUQR.MetaTau_p.value","MUQR.TtC","Notes")



n <- 1000
taus <- seq(0.1, 0.9, 0.1)
R <- 1000

p <- NULL
for(i in 1:R){
  prs <- rnorm(n, mean=0, sd=1)
  noise <- rnorm(n, mean=0, sd=1)
  y <- prs + noise
  myDT <- data.table(score=prs,outcome=y)
  myDT
  # p <- c(p,MQR(DT,y="outcome",g="prs",tau=taus,mqr.method="CQR")["MUQR.MetaTau_p.value"])
  MQR(datatable=myDT,y="outcome",g="score",tau=taus,
      mqr.method="CQR")["MCQR.MetaTau_p.value"]
  myDT
  # p <- c(p,MQR(DT,y="outcome",g="prs",tau=taus,mqr.method="CQR")["MCQR.MetaTau_p.value"])
}
sum(p<0.05)/1000

summary(p)
