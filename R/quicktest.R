
DT <- data.table(outcome=rnorm(1000),predictor=rnorm(1000))

args(mqr)
mqr(DT,y="outcome",g="predictor")
