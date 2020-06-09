
DT <- data.table(outcome=rnorm(1000),predictor=rnorm(1000))

args(MQR)
MQR(DT,y="outcome",g="predictor",mqr.method="CQR")

