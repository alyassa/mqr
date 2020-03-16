
# ============= Plotting Simulations: Global Scaling : Varying v.G ============
library(data.table)
# p <- "/home/alabada/pgpstore/R_Studio_Projects/Simulations_Quantile_Regression"
p <- "/home/../media/StoreB/arkan/Simulations_Quantile_Regression"
#p <- "/Users/Arkan/Desktop"
list.files(file.path(p,"Results","RData_Objects"))
load(file.path(p,"Results","RData_Objects", "vG_Global_Scaling_Oct_2_2018.RData"),
     envir=.GlobalEnv)
# v.G_arguments
args <- data.table(v.G_arguments)
args[,Index:=1:nrow(args)]

threshold <- 0.05 #pvalue threshold for significance
head(v.G_est[[1]])
pdf(file.path(p,"Results","RData_Objects","vG_Global_Scaling_FP_Oct_2_2018.pdf"),
    width=6,height=4)
#par(oma=c(0,0,0,0)) # set outer margin to zero
par(mar=c(2,2,2,0.5)) # number of lines per margin (bottom,left,top,right)
# margins are 3 lines wide and 2 lines tall in total so plots are slightly wider not pefrect squares
mfrow=c(2,3)
layout(matrix(1:6, 2, 3, byrow=TRUE),
       widths=lcm(rep(2*2.54,2)), heights=lcm(rep(2*2.54,3)))
for(i in v.GxE){
  for(j in MAF){
    # i <- v.GxE[1]
    # j <- MAF[1]
    i.j.est <- data.table(do.call(rbind,v.G_est[args[v.GxE==i][MAF==j][,Index]]))
     i.j.Result <- i.j.est[,lapply(.SD, function(x)
      length(which(as.numeric(x)<threshold))/length(x)),by=v.G,
      .SDcols=colnames(i.j.est)[grep("p.value", colnames(i.j.est))]]

    plot(x=i.j.Result[,v.G], y=i.j.Result[,MUQR.MetaTau_p.value],
         type="n", xlab="", ylab="",xlim=c(0, 0.0042), ylim=c(0,0.42),
         xaxt="n", yaxt="n", bty="l")
    title(main=bquote(paste("MAF =",.(j), sep="")), line=1)
    axis(2, tick=TRUE, labels=FALSE)
    axis(2, tick=FALSE, labels=TRUE, line=-0.5)
    axis(1, tick=TRUE, labels=FALSE)
    axis(1, tick=FALSE, labels=TRUE, line=-0.5)
    lines(x=i.j.Result[,v.G], y=i.j.Result[,rt.GxE.LN_p.value],
          type="o", pch=20, col="black", lwd=1)
    lines(x=i.j.Result[,v.G], y=i.j.Result[,MUQR.MetaTau_p.value],
          type="o", pch=20, col="darkmagenta", lwd=1)
    lines(x=i.j.Result[,v.G], y=i.j.Result[,RT_MUQR.MetaTau_p.value],
          type="o", pch=20, col="mediumseagreen", lwd=1)
  }
}
dev.off()


# ============= Plotting Simulations: Global Scaling ============
library(data.table)
p <- "/home/../media/StoreB/arkan/Simulations_Quantile_Regression"
list.files(file.path(p,"Results","RData_Objects"))

attach(file.path(p,"Results","RData_Objects","Global_Scaling_Sept_26_2018a.RData"))
A <- v.GxE_est
head(A[[1]])

attach(file.path(p,"Results","RData_Objects","Global_Scaling_Sept_26_2018b.RData"))
B <- v.GxE_est
head(B[[1]])

attach(file.path(p,"Results","RData_Objects","Global_Scaling_Sept_26_2018c.RData"))
C <- v.GxE_est
head(C[[1]])

attach(file.path(p,"Results","RData_Objects","Global_Scaling_Sept_26_2018d.RData"))
D <- v.GxE_est
head(D[[1]])

length(A)
v.GxE_est <- NULL
v.GxE_est <- list()
for(i in 1:length(A)){
  v.GxE_est[[i]] <- rbind(A[[i]],B[[i]],C[[i]],D[[i]])
}
dim(v.GxE_est[[1]])

v.GxE_arguments
args <- data.table(v.GxE_arguments)
args[,Index:=1:nrow(args)]

threshold <- 0.05 #pvalue threshold for significance
head(v.GxE_est[[1]])
#pdf(file.path(p,"Results","Results_Plotted","Global_Scaling_Sept_20_2018.pdf"),
#    width=6,height=4)
#par(oma=c(0,0,0,0)) # set outer margin to zero
#par(mar=c(2,2,2,0.5)) # number of lines per margin (bottom,left,top,right)
# margins are 3 lines wide and 2 lines tall in total so plots are slightly wider not pefrect squares
#mfrow=c(2,3)
#layout(matrix(1:6, 2, 3, byrow=TRUE),
#       widths=lcm(rep(2*2.54,2)), heights=lcm(rep(2*2.54,3)))
#for (i in 1:length(MAF)){
#  for(j in 1:length(Distributions)){
    i <- MAF[1]
    j <- Distributions[1]
    #i <- 2; j <- 1
    i.j.est <- data.table(do.call(rbind,v.GxE_est[args[MAF==i][Distributions==j][,Index]]))
    i.j.Result <- i.j.est[,lapply(.SD, function(x)
      length(which(as.numeric(x)<threshold))/length(x)),by=v.GxE,
      .SDcols=colnames(i.j.est)[grep("p.value", colnames(i.j.est))]]
    i.j.Result[,.(v.GxE,GxE.LN_p.value,MUQR.MetaTau_p.value,Scaling.Adjusted.MUQR_p.value,
                  dem.Scaling.Adjusted.MUQR_p.value)]

    y.max <- max(as.numeric(i.j.Result[,RIF.Meta_p.value]))
    y.max <- y.max+0.1*y.max
    y.min <- min(as.numeric(i.j.Result[,RIF.Meta_p.value]))
    y.min <- y.min-0.1*y.min
    y.lim <- c(y.min,y.max)


    plot(x=i.j.Result[,No.taus], y=i.j.Result[,RIF.Meta_p.value],
         type="n", xlab="", ylab="", ylim=y.lim,xlim=c(0, 60),
         xaxt="n", yaxt="n", bty="l")
    title(main=paste0("v.GxE= ",100*i), line=1)
    axis(2, tick=TRUE, labels=FALSE)
    axis(2, tick=FALSE, labels=TRUE, line=-0.5)
    axis(1, tick=TRUE, labels=FALSE)
    axis(1, tick=FALSE, labels=TRUE, line=-0.5)
    lines(x=i.j.Result[,No.taus], y=i.j.Result[,RIF.Meta_p.value],
          type="o", pch=20, col="darkmagenta", lwd=1)
    if(i==0){
      lines(x=i.j.Result[,No.taus], y=i.j.Result[,Levene_p.value],
            type="o", pch=20, col="mediumseagreen", lwd=1)
    }
  }
}
dev.off()





# ============= Plotting Simulations: SIM_Number_of_Taus ============
library(data.table)
p <- "/home/alabada/pgpstore/R_Studio_Projects/Simulations_Quantile_Regression"
list.files(file.path(p,"Results","RData_Objects"))
load(file.path(p,"Results","RData_Objects","Number_of_Taus_Apr_12_2017.RData"))
v.GxE_arguments
args <- data.table(v.GxE_arguments)
args[,Index:=1:nrow(args)]

threshold <- 5E-8 #pvalue threshold for significance
head(v.GxE_est[[1]])
pdf(file.path(p,"Results","Results_Plotted","Number_of_Taus_Apr_12_2017.pdf"),
    width=6,height=4)
#par(oma=c(0,0,0,0)) # set outer margin to zero
par(mar=c(2,2,2,0.5)) # number of lines per margin (bottom,left,top,right)
# margins are 3 lines wide and 2 lines tall in total so plots are slightly wider not pefrect squares
mfrow=c(2,3)
layout(matrix(1:6, 2, 3, byrow=TRUE),
       widths=lcm(rep(2*2.54,2)), heights=lcm(rep(2*2.54,3)))
for (i in v.GxE){
  for(j in MAF){
    #i <- v.GxE[1]
    #j <- MAF[1]
    i.j.est <- data.table(do.call(rbind,v.GxE_est[args[v.GxE==i][MAF==j][,Index]]))
    if(i==0){
      i.j.Result <- i.j.est[,lapply(.SD, function(x)
        length(which(as.numeric(x)<0.05))/length(x)),by=No.taus,
        .SDcols=colnames(i.j.est)[grep("p.value", colnames(i.j.est))]]
      y.lim <- c(0,0.3)
    } else {
      i.j.Result <- i.j.est[,lapply(.SD, function(x)
        length(which(as.numeric(x)<threshold))/length(x)),by=No.taus,
        .SDcols=colnames(i.j.est)[grep("p.value", colnames(i.j.est))]]
      y.max <- max(as.numeric(i.j.Result[,RIF.Meta_p.value]))
      y.max <- y.max+0.1*y.max
      y.min <- min(as.numeric(i.j.Result[,RIF.Meta_p.value]))
      y.min <- y.min-0.1*y.min
      y.lim <- c(y.min,y.max)
    }

    plot(x=i.j.Result[,No.taus], y=i.j.Result[,RIF.Meta_p.value],
         type="n", xlab="", ylab="", ylim=y.lim,xlim=c(0, 60),
         xaxt="n", yaxt="n", bty="l")
    title(main=paste0("v.GxE= ",100*i), line=1)
    axis(2, tick=TRUE, labels=FALSE)
    axis(2, tick=FALSE, labels=TRUE, line=-0.5)
    axis(1, tick=TRUE, labels=FALSE)
    axis(1, tick=FALSE, labels=TRUE, line=-0.5)
    lines(x=i.j.Result[,No.taus], y=i.j.Result[,RIF.Meta_p.value],
          type="o", pch=20, col="darkmagenta", lwd=1)
    if(i==0){
      lines(x=i.j.Result[,No.taus], y=i.j.Result[,Levene_p.value],
            type="o", pch=20, col="mediumseagreen", lwd=1)
    }
  }
}
dev.off()





# ============= Plotting Simulations: SIM_Pare_Fig3 ============
library(data.table)
p <- "/home/alabada/pgpstore/R_Studio_Projects/Simulations_Quantile_Regression"
list.files(file.path(p,"Results","RData_Objects"))
load(file.path(p,"Results","RData_Objects","vGxE_Fig3Sim_Apr_12_2017.RData"))
v.GxE_arguments
args <- data.table(v.GxE_arguments)
args[,Index:=1:nrow(args)]

threshold <- 5E-8 #pvalue threshold for significance
head(v.GxE_est[[1]])
pdf(file.path(p,"Results","Results_Plotted","v.GxE_Fig3Sim_Apr_12_2017.pdf"),
   width=6,height=4)
#par(oma=c(0,0,0,0)) # set outer margin to zero
par(mar=c(2,2,2,0.5)) # number of lines per margin (bottom,left,top,right)
# margins are 3 lines wide and 2 lines tall in total so plots are slightly wider not pefrect squares
mfrow=c(2,3)
layout(matrix(1:6, 2, 3, byrow=TRUE),
       widths=lcm(rep(2*2.54,2)), heights=lcm(rep(2*2.54,3)))
for (i in v.E){
  for(j in MAF){
    #i <- v.E[1]
    #j <- MAF[1]
    i.j.est <- data.table(do.call(rbind,v.GxE_est[args[v.E==i][MAF==j][,Index]]))
    i.j.Result <- i.j.est[,lapply(.SD, function(x)
      length(which(as.numeric(x)<threshold))/length(x)),by=v.GxE,
      .SDcols=colnames(i.j.est)[grep("p.value", colnames(i.j.est))]]

    plot(x=i.j.Result[,v.GxE], y=i.j.Result[,GxE.LN_p.value],
         type="n", xlab="", ylab="", ylim=c(0,1),xlim=c(0, 0.01),
         xaxt="n", yaxt="n", bty="l")
    title(main=bquote(paste("MAF=", .(j), " and ", beta[2],"=",.(i), sep="")), line=1)
    axis(2, tick=TRUE, labels=FALSE)
    axis(2, tick=FALSE, labels=TRUE, line=-0.5)
    axis(1, tick=TRUE, labels=FALSE)
    axis(1, tick=FALSE, labels=TRUE, line=-0.5)
    lines(x=i.j.Result[,v.GxE], y=i.j.Result[,GxE.LN_p.value],
          type="o", pch=20, col="black", lwd=1)
    lines(x=i.j.Result[,v.GxE], y=i.j.Result[,RIF.Meta_p.value],
          type="o", pch=20, col="darkmagenta", lwd=1)
    lines(x=i.j.Result[,v.GxE], y=i.j.Result[,Levene_p.value],
          type="o", pch=20, col="mediumseagreen", lwd=1)
  }
}
dev.off()

threshold <- 5E-8 #pvalue threshold for significance
pdf(file.path(p,"Results","Results_Plotted","v.GxE_Joint_Fig3Sim_Apr_12_2017.pdf"),
    width=6,height=4)
#par(oma=c(0,0,0,0)) # set outer margin to zero
par(mar=c(2,2,2,0.5)) # number of lines per margin (bottom,left,top,right)
# margins are 3 lines wide and 2 lines tall in total so plots are slightly wider not pefrect squares
mfrow=c(2,3)
layout(matrix(1:6, 2, 3, byrow=TRUE),
       widths=lcm(rep(2*2.54,2)), heights=lcm(rep(2*2.54,3)))
for (i in v.E){
  for(j in MAF){
    #i <- v.E[1]
    #j <- MAF[1]
    i.j.est <- data.table(do.call(rbind,v.GxE_est[args[v.E==i][MAF==j][,Index]]))
    i.j.Result <- i.j.est[,lapply(.SD, function(x)
      length(which(as.numeric(x)<threshold))/length(x)),by=v.GxE,
      .SDcols=colnames(i.j.est)[grep("Joint", colnames(i.j.est))]]

    plot(x=i.j.Result[,v.GxE], y=i.j.Result[,LN_Joint],
         type="n", xlab="", ylab="", ylim=c(0,1),xlim=c(0, 0.01),
         xaxt="n", yaxt="n", bty="l")
    title(main=bquote(paste("MAF=", .(j), " and ", beta[2],"=",.(i), sep="")), line=1)
    axis(2, tick=TRUE, labels=FALSE)
    axis(2, tick=FALSE, labels=TRUE, line=-0.5)
    axis(1, tick=TRUE, labels=FALSE)
    axis(1, tick=FALSE, labels=TRUE, line=-0.5)
    lines(x=i.j.Result[,v.GxE], y=i.j.Result[,LN_Joint],
          type="o", pch=20, col="black", lwd=1)
    lines(x=i.j.Result[,v.GxE], y=i.j.Result[,RIF_Joint],
          type="o", pch=20, col="darkmagenta", lwd=1)
    lines(x=i.j.Result[,v.GxE], y=i.j.Result[,Levene_Joint],
          type="o", pch=20, col="mediumseagreen", lwd=1)
  }
}
dev.off()

threshold <- 5E-8 #pvalue threshold for significance
pdf(file.path(p,"Results","Results_Plotted","v.GxE_Marginal_Fig3Sim_Mar_26_2017.pdf"),
    width=6,height=4)
#par(oma=c(0,0,0,0)) # set outer margin to zero
par(mar=c(2,2,2,0.5)) # number of lines per margin (bottom,left,top,right)
# margins are 3 lines wide and 2 lines tall in total so plots are slightly wider not pefrect squares
mfrow=c(2,3)
layout(matrix(1:6, 2, 3, byrow=TRUE),
       widths=lcm(rep(2*2.54,2)), heights=lcm(rep(2*2.54,3)))
for (i in v.E){
  for(j in MAF){
    #i <- v.E[1]
    #j <- MAF[1]
    i.j.est <- data.table(do.call(rbind,v.GxE_est[args[v.E==i][MAF==j][,Index]]))
    i.j.Result <- i.j.est[,lapply(.SD, function(x)
      length(which(as.numeric(x)<threshold))/length(x)),by=v.GxE,
      .SDcols=c("G.LN_p.value","G.RIF_p.value","RIF.MetaMarg_p.value")]

    plot(x=i.j.Result[,v.GxE], y=i.j.Result[,G.LN_p.value],
         type="n", xlab="", ylab="", ylim=c(0,1),xlim=c(0, 0.01),
         xaxt="n", yaxt="n", bty="l")
    title(main=bquote(paste("MAF=", .(j), " and ", beta[2],"=",.(i), sep="")), line=1)
    axis(2, tick=TRUE, labels=FALSE)
    axis(2, tick=FALSE, labels=TRUE, line=-0.5)
    axis(1, tick=TRUE, labels=FALSE)
    axis(1, tick=FALSE, labels=TRUE, line=-0.5)
    lines(x=i.j.Result[,v.GxE], y=i.j.Result[,G.LN_p.value],
          type="o", pch=20, col="black", lwd=1)
    lines(x=i.j.Result[,v.GxE], y=i.j.Result[,RIF.MetaMarg_p.value],
          type="o", pch=20, col="darkmagenta", lwd=1)
    lines(x=i.j.Result[,v.GxE], y=i.j.Result[,G.RIF_p.value],
          type="o", pch=20, col="darkmagenta", lwd=1,lty=3)
  }
}
dev.off()

load(file.path(p,"Results","RData_Objects","vG_Fig3Sim_Apr_12_2017.RData"))
v.G_arguments
args <- data.table(v.G_arguments)
args[,Index:=1:nrow(args)]
threshold <- 5E-8 #pvalue threshold for significance
pdf(file.path(p,"Results","Results_Plotted","v.G_Joint_Fig3Sim_Apr_12_2017.pdf"),
    width=6,height=4)
#par(oma=c(0,0,0,0)) # set outer margin to zero
par(mar=c(2,2,2,0.5)) # number of lines per margin (bottom,left,top,right)
# margins are 3 lines wide and 2 lines tall in total so plots are slightly wider not pefrect squares
mfrow=c(2,3)
layout(matrix(1:6, 2, 3, byrow=TRUE),
       widths=lcm(rep(2*2.54,2)), heights=lcm(rep(2*2.54,3)))
for (i in v.E){
  for(j in MAF){
    #i <- v.E[1]
    #j <- MAF[1]
    i.j.est <- data.table(do.call(rbind,v.G_est[args[v.E==i][MAF==j][,Index]]))
    i.j.Result <- i.j.est[,lapply(.SD, function(x)
      length(which(as.numeric(x)<threshold))/length(x)),by=v.G,
      .SDcols=colnames(i.j.est)[grep("Joint", colnames(i.j.est))]]

    plot(x=i.j.Result[,v.G], y=i.j.Result[,LN_Joint],
         type="n", xlab="", ylab="", ylim=c(0,1),xlim=c(0, 0.006),
         xaxt="n", yaxt="n", bty="l")
    title(main=bquote(paste("MAF=", .(j), " and ", beta[2],"=",.(i), sep="")), line=1)
    axis(2, tick=TRUE, labels=FALSE)
    axis(2, tick=FALSE, labels=TRUE, line=-0.5)
    axis(1, tick=TRUE, labels=FALSE)
    axis(1, tick=FALSE, labels=TRUE, line=-0.5)
    lines(x=i.j.Result[,v.G], y=i.j.Result[,LN_Joint],
          type="o", pch=20, col="black", lwd=1)
    lines(x=i.j.Result[,v.G], y=i.j.Result[,RIF_Joint],
          type="o", pch=20, col="darkmagenta", lwd=1)
    lines(x=i.j.Result[,v.G], y=i.j.Result[,Levene_Joint],
          type="o", pch=20, col="mediumseagreen", lwd=1)
  }
}
dev.off()


# ============= Plotting Simulations: SIM_Skew ============
library(data.table)
path <- "/home/alabada/pgpstore/R_Studio_Projects/Simulations_Quantile_Regression"
#path <- "/Users/Arkan/Desktop"
list.files(file.path(path,"Results","RData_Objects"))

load(file.path(path,"Results","RData_Objects","Skew_SIM_Apr_12_2017.RData"))
#load(file.path(path,"Results","RData_Objects","N_Skew_HWE_Sim_Dec_5_2016.RData"))
path <- "/home/alabada/pgpstore/R_Studio_Projects/Simulations_Quantile_Regression"
#path <- "/Users/Arkan/Desktop"

v.GxE_arguments
args <- data.table(v.GxE_arguments)
args[,Index:=1:nrow(args)]
#setnames(args,c("Directions","Distributions"),c("Interaction_Direction","error.setting"))
args[,Index:=1:nrow(args)]
args[,Distributions:=ifelse(Type=="Normal", "Normal",
                            ifelse(Type=="Chi" & Skew=="Left", "Left-Skew",
                                   ifelse(Type=="Chi" & Skew=="Right", "Right-Skew",
                                          ifelse(Type=="Uniform", "Short-Tail",
                                                 ifelse(Type=="Laplace", "Long-Tail",
                                                        ifelse(Type=="Right-Outliers","Right-Outliers",
                                                               ifelse(Type=="Both-Outliers","Both-Outliers",NA)))))))]

table(args[,Distributions])
threshold <- 5E-8 #pvalue threshold for significance
pdf(file.path(path,"Results","Results_Plotted","Skew_Sim_Apr_12_2017.pdf"),
    width=7.2,height=9)
#par(oma=c(0,0,0,0)) # set outer margin to zero
par(mar=c(2,1.6,2,0.2)) # number of lines per margin (bottom,left,top,right)
# margins are 3 lines wide and 2 lines tall in total so plots are slightly wider not pefrect squares
mfrow=c(4,3)
layout(matrix(1:12, 4, 3, byrow = TRUE),
       widths=lcm(rep(2*2.54,3)), heights=lcm(rep(2*2.54,4)))
for(j in Distributions){
  #i <- Directions[2]
  #j <- Distributions[1]
  t <- args[Distributions==j][v.GxE==0][
    ,.(Special,N,MAF,Varying,b0,v.G,v.E,v.GxE,Interaction.Dir,hwe.p,min.grp.size,
      Pare.Encoding,Type,a,b,Skew,se,Adjusted.for,No.taus)]
  t <- as.matrix(data.frame(t))[1,]
  test.data <- Test.Data.function(Arguments=t)
  hist(test.data[,y],breaks=100,prob=TRUE,main=j,xlab="",col="grey",border=NA)
  curve(dnorm(x, mean=mean(test.data[,y]), sd=sd(test.data[,y])),add=TRUE, col="red")
  qqnorm(test.data[,y],bty="l",xlab="", ylab="");
  qqline(test.data[,y],col="red")

  i.j.est <- data.table(do.call(rbind,v.GxE_est[args[Distributions==j][,Index]]))
  i.j.Result <- i.j.est[,lapply(.SD, function(x)
    length(which(as.numeric(x)<threshold))/length(x)),by=v.GxE,
    .SDcols=colnames(i.j.est)[grep("p.value",colnames(i.j.est))]]

  plot(x=i.j.Result[,v.GxE], y=i.j.Result[,GxE.LN_p.value],
       type="n", xlab="", ylab="", ylim=c(0,1),xlim=c(0, 0.01),
       xaxt="n", yaxt="n", bty="l")
  title(main=j,line=1)
  axis(2, tick=TRUE, labels=FALSE)
  axis(2, tick=FALSE, labels=TRUE, line=-0.5)
  axis(1, tick=TRUE, labels=FALSE)
  axis(1, tick=FALSE, labels=TRUE, line=-0.5)
  lines(x=i.j.Result[,v.GxE], y=i.j.Result[,GxE.LN_p.value],
        type="o", pch=20, col="black", lwd=1)
  lines(x=i.j.Result[,v.GxE], y=i.j.Result[,RIF.Meta_p.value],
        type="o", pch=20, col="darkmagenta", lwd=1)
  lines(x=i.j.Result[,v.GxE], y=i.j.Result[,Levene_p.value],
        type="o", pch=20, col="mediumseagreen", lwd=1)
}
dev.off()

threshold <- 5E-8 #pvalue threshold for significance
pdf(file.path(path,"Results","Results_Plotted","Joint_Skew_Sim_Apr_12_2017.pdf"),
    width=7.2,height=9)
#par(oma=c(0,0,0,0)) # set outer margin to zero
par(mar=c(2,1.6,2,0.2)) # number of lines per margin (bottom,left,top,right)
# margins are 3 lines wide and 2 lines tall in total so plots are slightly wider not pefrect squares
mfrow=c(4,3)
layout(matrix(1:12, 4, 3, byrow = TRUE),
       widths=lcm(rep(2*2.54,3)), heights=lcm(rep(2*2.54,4)))
for(j in Distributions){
  #i <- Directions[2]
  #j <- Distributions[1]
  t <- args[Distributions==j][v.GxE==0][
    ,.(Special,N,MAF,Varying,b0,v.G,v.E,v.GxE,Interaction.Dir,hwe.p,min.grp.size,
       Pare.Encoding,Type,a,b,Skew,se,Adjusted.for,No.taus)]
  t <- as.matrix(data.frame(t))[1,]
  test.data <- Test.Data.function(Arguments=t)
  hist(test.data[,y],breaks=100,prob=TRUE,main=j,xlab="",col="grey",border=NA)
  curve(dnorm(x, mean=mean(test.data[,y]), sd=sd(test.data[,y])),add=TRUE, col="red")
  qqnorm(test.data[,y],bty="l",xlab="", ylab="");
  qqline(test.data[,y],col="red")

  i.j.est <- data.table(do.call(rbind,v.GxE_est[args[Distributions==j][,Index]]))
  i.j.Result <- i.j.est[,lapply(.SD, function(x)
    length(which(as.numeric(x)<threshold))/length(x)),by=v.GxE,
    .SDcols=colnames(i.j.est)[grep("Joint",colnames(i.j.est))]]

  plot(x=i.j.Result[,v.GxE], y=i.j.Result[,LN_Joint],
       type="n", xlab="", ylab="", ylim=c(0,1),xlim=c(0, 0.01),
       xaxt="n", yaxt="n", bty="l")
  title(main=j,line=1)
  axis(2, tick=TRUE, labels=FALSE)
  axis(2, tick=FALSE, labels=TRUE, line=-0.5)
  axis(1, tick=TRUE, labels=FALSE)
  axis(1, tick=FALSE, labels=TRUE, line=-0.5)
  lines(x=i.j.Result[,v.GxE], y=i.j.Result[,LN_Joint],
        type="o", pch=20, col="black", lwd=1)
  lines(x=i.j.Result[,v.GxE], y=i.j.Result[,RIF_Joint],
        type="o", pch=20, col="darkmagenta", lwd=1)
  lines(x=i.j.Result[,v.GxE], y=i.j.Result[,Levene_Joint],
        type="o", pch=20, col="mediumseagreen", lwd=1)
}
dev.off()




# ============= Plotting Simulations: SIM_FalsePositives : Varying v.E ============
library(data.table)
p <- "/home/alabada/pgpstore/R_Studio_Projects/Simulations_Quantile_Regression"
#p <- "/Users/Arkan/Desktop"
list.files(file.path(p,"Results","RData_Objects"))
load(file.path(p,"Results","RData_Objects","Adj_E_vE_FalsePostives_Apr_12_2017.RData"),
     envir=.GlobalEnv)
Adj.arguments <- arguments
Adj.est <- est
rm(arguments,est);gc()

p <- "/home/alabada/pgpstore/R_Studio_Projects/Simulations_Quantile_Regression"
list.files(file.path(p,"Results","RData_Objects"))
load(file.path(p,"Results","RData_Objects","unAdj_vE_FalsePostives_Apr_12_2017.RData"),
     envir=.GlobalEnv)
unAdj.arguments <- arguments
unAdj.est <- est
rm(arguments,est);gc()

args <- data.table(Adj.arguments)
args[,Index:=1:nrow(args)]
head(Adj.est[[1]])

threshold <- 0.05 #pvalue threshold for significance
pdf(file.path(p,"Results","Results_Plotted","vE_FalsePostives_Apr_12_2017.pdf"),
   width=6,height=4)
#par(oma=c(0,0,0,0)) # set outer margin to zero
par(mar=c(2,2,2,0.5)) # number of lines per margin (bottom,left,top,right)
# margins are 3 lines wide and 2 lines tall in total so plots are slightly wider not pefrect squares
mfrow=c(2,3)
layout(matrix(1:6, 2, 3, byrow=TRUE),
       widths=lcm(rep(2*2.54,2)), heights=lcm(rep(2*2.54,3)))
for (i in Distributions){
  for(j in MAF){
    #i <- Distributions[1]
    #j <- MAF[1]
    if(i=="Normal"){T <- "Normal"}else{T="Chi"}
    i.j.Adj.est <- data.table(do.call(rbind,Adj.est[args[Type==T][MAF==j][,Index]]))
    i.j.Adj.Result <- i.j.Adj.est[,lapply(.SD, function(x)
      length(which(as.numeric(x)<threshold))/length(x)),by=v.E,
      .SDcols=colnames(i.j.Adj.est)[grep("p.value", colnames(i.j.Adj.est))]]

    i.j.unAdj.est <- data.table(do.call(rbind,unAdj.est[args[Type==T][MAF==j][,Index]]))
    i.j.unAdj.Result <- i.j.unAdj.est[,lapply(.SD, function(x)
      length(which(as.numeric(x)<threshold))/length(x)),by=v.E,
      .SDcols=colnames(i.j.unAdj.est)[grep("p.value", colnames(i.j.unAdj.est))]]

    plot(x=i.j.unAdj.Result[,v.E], y=i.j.unAdj.Result[,GxE.LN_p.value],
         type="n", xlab="", ylab="", ylim=c(0,0.2),xlim=c(0, 0.4),
         xaxt="n", yaxt="n", bty="l")
    title(main=bquote(paste(.(i), ", MAF =",.(j), sep="")), line=1)
    axis(2, tick=TRUE, labels=FALSE)
    axis(2, tick=FALSE, labels=TRUE, line=-0.5)
    axis(1, tick=TRUE, labels=FALSE)
    axis(1, tick=FALSE, labels=TRUE, line=-0.5)
    lines(x=i.j.unAdj.Result[,v.E], y=i.j.unAdj.Result[,GxE.LN_p.value],
          type="o", pch=20, col="black", lwd=1)
    lines(x=i.j.unAdj.Result[,v.E], y=i.j.unAdj.Result[,RIF.Meta_p.value],
          type="o", pch=20, col="darkmagenta", lwd=1)
    lines(x=i.j.unAdj.Result[,v.E], y=i.j.unAdj.Result[,Levene_p.value],
          type="o", pch=20, col="mediumseagreen", lwd=1)

    #lines(x=i.j.Adj.Result[,v.E], y=i.j.Adj.Result[,GxE.LN_p.value],
     #     type="o", pch=20, col="black", lwd=1,lty=2)
    #lines(x=i.j.Adj.Result[,v.E], y=i.j.Adj.Result[,RIF.Meta_p.value],
     #     type="o", pch=20, col="darkmagenta", lwd=1,lty=2)
    #lines(x=i.j.Adj.Result[,v.E], y=i.j.Adj.Result[,Levene_p.value],
     #     type="o", pch=20, col="mediumseagreen", lwd=1,lty=2)
  }
}
dev.off()

threshold <- 0.05
pdf(file.path(p,"Results","Results_Plotted","vE_2STEP_FalsePostives_Apr_12_2017.pdf"),
    width=6,height=4)
#par(oma=c(0,0,0,0)) # set outer margin to zero
par(mar=c(2,2,2,0.5)) # number of lines per margin (bottom,left,top,right)
# margins are 3 lines wide and 2 lines tall in total so plots are slightly wider not pefrect squares
mfrow=c(2,3)
layout(matrix(1:6, 2, 3, byrow=TRUE),
       widths=lcm(rep(2*2.54,2)), heights=lcm(rep(2*2.54,3)))
for (i in Distributions){
  for(j in MAF){
    #i <- Distributions[1]
    #j <- MAF[1]
    if(i=="Normal"){T <- "Normal"}else{T="Chi"}
    i.j.Adj.est <- data.table(do.call(rbind,Adj.est[args[Type==T][MAF==j][,Index]]))
    i.j.Adj.ORs <- i.j.Adj.est[,lapply(.SD, function(x)
      log(fisher.test(table(as.numeric(GxE.LN_p.value)<threshold,
                        as.numeric(x)<threshold))$estimate)),by=v.E,
      .SDcols=colnames(i.j.Adj.est)[grep("p.value", colnames(i.j.Adj.est))]]
    i.j.Adj.lower <- i.j.Adj.est[,lapply(.SD, function(x)
      log(fisher.test(table(as.numeric(GxE.LN_p.value)<threshold,
                        as.numeric(x)<threshold))$conf.int[[1]])),by=v.E,
      .SDcols=colnames(i.j.Adj.est)[grep("p.value", colnames(i.j.Adj.est))]]
    i.j.Adj.upper <- i.j.Adj.est[,lapply(.SD, function(x)
      log(fisher.test(table(as.numeric(GxE.LN_p.value)<threshold,
                        as.numeric(x)<threshold))$conf.int[[2]])),by=v.E,
      .SDcols=colnames(i.j.Adj.est)[grep("p.value", colnames(i.j.Adj.est))]]

    i.j.unAdj.est <- data.table(do.call(rbind,unAdj.est[args[Type==T][MAF==j][,Index]]))
    i.j.unAdj.ORs <- i.j.unAdj.est[,lapply(.SD, function(x)
      log(fisher.test(table(as.numeric(GxE.LN_p.value)<threshold,
                        as.numeric(x)<threshold))$estimate)),by=v.E,
      .SDcols=colnames(i.j.unAdj.est)[grep("p.value", colnames(i.j.unAdj.est))]]
    i.j.unAdj.lower <- i.j.unAdj.est[,lapply(.SD, function(x)
      log(fisher.test(table(as.numeric(GxE.LN_p.value)<threshold,
                        as.numeric(x)<threshold))$conf.int[[1]])),by=v.E,
      .SDcols=colnames(i.j.unAdj.est)[grep("p.value", colnames(i.j.unAdj.est))]]
    i.j.unAdj.upper <- i.j.unAdj.est[,lapply(.SD, function(x)
      log(fisher.test(table(as.numeric(GxE.LN_p.value)<threshold,
                        as.numeric(x)<threshold))$conf.int[[2]])),by=v.E,
      .SDcols=colnames(i.j.unAdj.est)[grep("p.value", colnames(i.j.unAdj.est))]]

    y.MAX <- rbind(data.frame(i.j.Adj.upper[,lapply(.SD, function(x) max(as.numeric(x)))]),
                   data.frame(i.j.unAdj.upper[,lapply(.SD, function(x) max(as.numeric(x)))]))
    y.MAX <- y.MAX[,c("RIF.Meta_p.value","Levene_p.value")]
    y.MAX <- max(as.numeric(as.matrix(y.MAX)))
    y.MAX <- y.MAX+0.1*y.MAX

    y.MIN <- rbind(data.frame(i.j.Adj.lower[,lapply(.SD, function(x) min(as.numeric(x)))]),
                   data.frame(i.j.unAdj.lower[,lapply(.SD, function(x) min(as.numeric(x)))]))
    y.MIN <- y.MIN[,c("RIF.Meta_p.value","Levene_p.value")]
    y.MIN <- min(as.numeric(as.matrix(y.MIN)))
    y.MIN <- y.MIN+0.1*y.MIN

    plot(x=i.j.unAdj.ORs[,v.E], y=log(i.j.unAdj.ORs[,RIF.Meta_p.value]),
         type="n", xlab="", ylab="", ylim=c(y.MIN,y.MAX),xlim=c(0, 0.4),
         xaxt="n", yaxt="n", bty="l")
    title(main=bquote(paste(.(i), ", MAF =",.(j), sep="")), line=1)
    axis(2, tick=TRUE, labels=FALSE)
    axis(2, tick=FALSE, labels=TRUE, line=-0.5)
    axis(1, tick=TRUE, labels=FALSE)
    axis(1, tick=FALSE, labels=TRUE, line=-0.5)
    abline(h=0,col="black")
    lines(x=i.j.Adj.ORs[,v.E], y=i.j.Adj.ORs[,RIF.Meta_p.value],
          type="o", pch=20, col="darkmagenta", lwd=1,lty=2)
    arrows(as.numeric(i.j.Adj.ORs[,v.E]), i.j.Adj.lower[,RIF.Meta_p.value],
           as.numeric(i.j.Adj.ORs[,v.E]), i.j.Adj.upper[,RIF.Meta_p.value], length=0.02,
           angle=90, code=3,lwd=0.5)
    lines(x=i.j.Adj.ORs[,v.E], y=i.j.Adj.ORs[,Levene_p.value],
          type="o", pch=20, col="mediumseagreen", lwd=1,lty=2)
    arrows(as.numeric(i.j.Adj.ORs[,v.E]), i.j.Adj.lower[,Levene_p.value],
           as.numeric(i.j.Adj.ORs[,v.E]), i.j.Adj.upper[,Levene_p.value], length=0.02,
           angle=90, code=3,lwd=0.5)

    lines(x=i.j.unAdj.ORs[,v.E], y=i.j.unAdj.ORs[,RIF.Meta_p.value],
          type="o", pch=20, col="darkmagenta", lwd=1)
    arrows(as.numeric(i.j.unAdj.ORs[,v.E]), i.j.unAdj.lower[,RIF.Meta_p.value],
           as.numeric(i.j.unAdj.ORs[,v.E]), i.j.unAdj.upper[,RIF.Meta_p.value], length=0.02,
           angle=90, code=3,lwd=0.5)
    lines(x=i.j.unAdj.ORs[,v.E], y=i.j.unAdj.ORs[,Levene_p.value],
          type="o", pch=20, col="mediumseagreen", lwd=1)
    arrows(as.numeric(i.j.unAdj.ORs[,v.E]), i.j.unAdj.lower[,Levene_p.value],
           as.numeric(i.j.unAdj.ORs[,v.E]), i.j.unAdj.upper[,Levene_p.value], length=0.02,
           angle=90, code=3,lwd=0.5)
  }
}
dev.off()

# ============= Plotting Simulations: SIM_FalsePositives : Varying v.G ============
library(data.table)
p <- "/home/alabada/pgpstore/R_Studio_Projects/Simulations_Quantile_Regression"
#p <- "/Users/Arkan/Desktop"
list.files(file.path(p,"Results","RData_Objects"))
load(file.path(p,"Results","RData_Objects","Adj_E_vG_FalsePostives_Apr_12_2017.RData"),
     envir=.GlobalEnv)
Adj.arguments <- arguments
Adj.est <- est
rm(arguments,est);gc()

p <- "/home/alabada/pgpstore/R_Studio_Projects/Simulations_Quantile_Regression"
list.files(file.path(p,"Results","RData_Objects"))
load(file.path(p,"Results","RData_Objects","unAdj_vG_FalsePostives_Apr_12_2017.RData"),
     envir=.GlobalEnv)
unAdj.arguments <- arguments
unAdj.est <- est
rm(arguments,est);gc()

args <- data.table(Adj.arguments)
args[,Index:=1:nrow(args)]
head(Adj.est[[1]])

threshold <- 0.05 #pvalue threshold for significance
pdf(file.path(p,"Results","Results_Plotted","vG_FalsePostives_Apr_12_2017.pdf"),
    width=6,height=4)
#par(oma=c(0,0,0,0)) # set outer margin to zero
par(mar=c(2,2,2,0.5)) # number of lines per margin (bottom,left,top,right)
# margins are 3 lines wide and 2 lines tall in total so plots are slightly wider not pefrect squares
mfrow=c(2,3)
layout(matrix(1:6, 2, 3, byrow=TRUE),
       widths=lcm(rep(2*2.54,2)), heights=lcm(rep(2*2.54,3)))
for (i in Distributions){
  for(j in MAF){
    #i <- Distributions[1]
    #j <- MAF[1]
    if(i=="Normal"){T <- "Normal"}else{T="Chi"}
    i.j.Adj.est <- data.table(do.call(rbind,Adj.est[args[Type==T][MAF==j][,Index]]))
    i.j.Adj.Result <- i.j.Adj.est[,lapply(.SD, function(x)
      length(which(as.numeric(x)<threshold))/length(x)),by=v.G,
      .SDcols=colnames(i.j.Adj.est)[grep("p.value", colnames(i.j.Adj.est))]]

    i.j.unAdj.est <- data.table(do.call(rbind,unAdj.est[args[Type==T][MAF==j][,Index]]))
    i.j.unAdj.Result <- i.j.unAdj.est[,lapply(.SD, function(x)
      length(which(as.numeric(x)<threshold))/length(x)),by=v.G,
      .SDcols=colnames(i.j.unAdj.est)[grep("p.value", colnames(i.j.unAdj.est))]]

    plot(x=i.j.unAdj.Result[,v.G], y=i.j.unAdj.Result[,RIF.Meta_p.value],
         type="n", xlab="", ylab="",xlim=c(0, 0.01), ylim=c(0,0.2),
         xaxt="n", yaxt="n", bty="l")
    title(main=bquote(paste(.(i), ", MAF =",.(j), sep="")), line=1)
    axis(2, tick=TRUE, labels=FALSE)
    axis(2, tick=FALSE, labels=TRUE, line=-0.5)
    axis(1, tick=TRUE, labels=FALSE)
    axis(1, tick=FALSE, labels=TRUE, line=-0.5)
    lines(x=i.j.unAdj.Result[,v.G], y=i.j.unAdj.Result[,GxE.LN_p.value],
          type="o", pch=20, col="black", lwd=1)
    lines(x=i.j.unAdj.Result[,v.G], y=i.j.unAdj.Result[,RIF.Meta_p.value],
          type="o", pch=20, col="darkmagenta", lwd=1)
    lines(x=i.j.unAdj.Result[,v.G], y=i.j.unAdj.Result[,Levene_p.value],
          type="o", pch=20, col="mediumseagreen", lwd=1)

    #lines(x=i.j.Adj.Result[,v.G], y=i.j.Adj.Result[,GxE.LN_p.value],
     #     type="o", pch=20, col="black", lwd=1,lty=2)
    #lines(x=i.j.Adj.Result[,v.G], y=i.j.Adj.Result[,RIF.Meta_p.value],
     #     type="o", pch=20, col="darkmagenta", lwd=1,lty=2)
    #lines(x=i.j.Adj.Result[,v.G], y=i.j.Adj.Result[,Levene_p.value],
    #      type="o", pch=20, col="mediumseagreen", lwd=1,lty=2)
  }
}
dev.off()

threshold <- 0.05
pdf(file.path(p,"Results","Results_Plotted","vG_2STEP_FalsePostives_Apr_12_2017.pdf"),
    width=6,height=4)
#par(oma=c(0,0,0,0)) # set outer margin to zero
par(mar=c(2,2,2,0.5)) # number of lines per margin (bottom,left,top,right)
# margins are 3 lines wide and 2 lines tall in total so plots are slightly wider not pefrect squares
mfrow=c(2,3)
layout(matrix(1:6, 2, 3, byrow=TRUE),
       widths=lcm(rep(2*2.54,2)), heights=lcm(rep(2*2.54,3)))
for (i in Distributions){
  for(j in MAF){
    #i <- Distributions[1]
    #j <- MAF[1]
    if(i=="Normal"){T <- "Normal"}else{T="Chi"}
    i.j.Adj.est <- data.table(do.call(rbind,Adj.est[args[Type==T][MAF==j][,Index]]))
    #STOPPED HERE
    i.j.Adj.ORs <- i.j.Adj.est[,lapply(.SD, function(x)
      log(fisher.test(table(as.numeric(GxE.LN_p.value)<threshold,
                        as.numeric(x)<threshold))$estimate)),by=v.G,
      .SDcols=c("QREG.Meta_p.value","RIF.Meta_p.value","Levene_p.value")]
    i.j.Adj.lower <- i.j.Adj.est[,lapply(.SD, function(x)
      log(fisher.test(table(as.numeric(GxE.LN_p.value)<threshold,
                        as.numeric(x)<threshold))$conf.int[[1]])),by=v.G,
      .SDcols=c("QREG.Meta_p.value","RIF.Meta_p.value","Levene_p.value")]
    i.j.Adj.upper <- i.j.Adj.est[,lapply(.SD, function(x)
      log(fisher.test(table(as.numeric(GxE.LN_p.value)<threshold,
                        as.numeric(x)<threshold))$conf.int[[2]])),by=v.G,
      .SDcols=c("QREG.Meta_p.value","RIF.Meta_p.value","Levene_p.value")]

    i.j.unAdj.est <- data.table(do.call(rbind,unAdj.est[args[Type==T][MAF==j][,Index]]))
    i.j.unAdj.ORs <- i.j.unAdj.est[,lapply(.SD, function(x)
      log(fisher.test(table(as.numeric(GxE.LN_p.value)<threshold,
                        as.numeric(x)<threshold))$estimate)),by=v.G,
      .SDcols=c("QREG.Meta_p.value","RIF.Meta_p.value","Levene_p.value")]
    i.j.unAdj.lower <- i.j.unAdj.est[,lapply(.SD, function(x)
      log(fisher.test(table(as.numeric(GxE.LN_p.value)<threshold,
                        as.numeric(x)<threshold))$conf.int[[1]])),by=v.G,
      .SDcols=c("QREG.Meta_p.value","RIF.Meta_p.value","Levene_p.value")]
    i.j.unAdj.upper <- i.j.unAdj.est[,lapply(.SD, function(x)
      log(fisher.test(table(as.numeric(GxE.LN_p.value)<threshold,
                        as.numeric(x)<threshold))$conf.int[[2]])),by=v.G,
      .SDcols=c("QREG.Meta_p.value","RIF.Meta_p.value","Levene_p.value")]

    y.MAX <- rbind(data.frame(i.j.Adj.upper[,lapply(.SD, function(x) max(as.numeric(x)))]),
                   data.frame(i.j.unAdj.upper[,lapply(.SD, function(x) max(as.numeric(x)))]))
    y.MAX <- y.MAX[,c("RIF.Meta_p.value","Levene_p.value")]
    y.MAX <- max(as.numeric(as.matrix(y.MAX)))
    y.MAX <- y.MAX+0.1*y.MAX

    y.MIN <- rbind(data.frame(i.j.Adj.lower[,lapply(.SD, function(x) min(as.numeric(x)))]),
                   data.frame(i.j.unAdj.lower[,lapply(.SD, function(x) min(as.numeric(x)))]))
    y.MIN <- y.MIN[,c("RIF.Meta_p.value","Levene_p.value")]
    y.MIN <- min(as.numeric(as.matrix(y.MIN)))
    y.MIN <- y.MIN+0.1*y.MIN

    plot(x=i.j.unAdj.ORs[,v.G], y=log(i.j.unAdj.ORs[,RIF.Meta_p.value]),
         type="n", xlab="", ylab="", ylim=c(y.MIN,y.MAX),xlim=c(0, 0.01),
         xaxt="n", yaxt="n", bty="l")
    title(main=bquote(paste(.(i), ", MAF =",.(j), sep="")), line=1)
    axis(2, tick=TRUE, labels=FALSE)
    axis(2, tick=FALSE, labels=TRUE, line=-0.5)
    axis(1, tick=TRUE, labels=FALSE)
    axis(1, tick=FALSE, labels=TRUE, line=-0.5)
    abline(h=0,col="black")
    lines(x=i.j.Adj.ORs[,v.G], y=i.j.Adj.ORs[,RIF.Meta_p.value],
          type="o", pch=20, col="darkmagenta", lwd=1,lty=2)
    arrows(as.numeric(i.j.Adj.ORs[,v.G]), i.j.Adj.lower[,RIF.Meta_p.value],
           as.numeric(i.j.Adj.ORs[,v.G]), i.j.Adj.upper[,RIF.Meta_p.value], length=0.02,
           angle=90, code=3,lwd=0.5)
    lines(x=i.j.Adj.ORs[,v.G], y=i.j.Adj.ORs[,Levene_p.value],
          type="o", pch=20, col="mediumseagreen", lwd=1,lty=2)
    arrows(as.numeric(i.j.Adj.ORs[,v.G]), i.j.Adj.lower[,Levene_p.value],
           as.numeric(i.j.Adj.ORs[,v.G]), i.j.Adj.upper[,Levene_p.value], length=0.02,
           angle=90, code=3,lwd=0.5)

    lines(x=i.j.unAdj.ORs[,v.G], y=i.j.unAdj.ORs[,RIF.Meta_p.value],
          type="o", pch=20, col="darkmagenta", lwd=1)
    arrows(as.numeric(i.j.unAdj.ORs[,v.G]), i.j.unAdj.lower[,RIF.Meta_p.value],
           as.numeric(i.j.unAdj.ORs[,v.G]), i.j.unAdj.upper[,RIF.Meta_p.value], length=0.02,
           angle=90, code=3,lwd=0.5)
    lines(x=i.j.unAdj.ORs[,v.G], y=i.j.unAdj.ORs[,Levene_p.value],
          type="o", pch=20, col="mediumseagreen", lwd=1)
    arrows(as.numeric(i.j.unAdj.ORs[,v.G]), i.j.unAdj.lower[,Levene_p.value],
           as.numeric(i.j.unAdj.ORs[,v.G]), i.j.unAdj.upper[,Levene_p.value], length=0.02,
           angle=90, code=3,lwd=0.5)
  }
}
dev.off()

threshold <- 0.05
pdf(file.path(p,"Results","Results_Plotted","vG_2STEP_Joint_FalsePostives_Apr_12_2017.pdf"),
    width=6,height=4)
#par(oma=c(0,0,0,0)) # set outer margin to zero
par(mar=c(2,2,2,0.5)) # number of lines per margin (bottom,left,top,right)
# margins are 3 lines wide and 2 lines tall in total so plots are slightly wider not pefrect squares
mfrow=c(2,3)
layout(matrix(1:6, 2, 3, byrow=TRUE),
       widths=lcm(rep(2*2.54,2)), heights=lcm(rep(2*2.54,3)))
for (i in Distributions){
  for(j in MAF){
    #i <- Distributions[1]
    #j <- MAF[1]
    if(i=="Normal"){T <- "Normal"}else{T="Chi"}
    i.j.Adj.est <- data.table(do.call(rbind,Adj.est[args[Type==T][MAF==j][,Index]]))
    i.j.Adj.ORs <- i.j.Adj.est[,lapply(.SD, function(x)
      fisher.test(table(as.numeric(GxE.LN_p.value)<threshold,
                        as.numeric(x)<threshold))$estimate),by=v.G,
      .SDcols=colnames(i.j.Adj.est)[grep("Joint", colnames(i.j.Adj.est))]]
    i.j.Adj.lower <- i.j.Adj.est[,lapply(.SD, function(x)
      fisher.test(table(as.numeric(GxE.LN_p.value)<threshold,
                        as.numeric(x)<threshold))$conf.int[[1]]),by=v.G,
      .SDcols=colnames(i.j.Adj.est)[grep("Joint", colnames(i.j.Adj.est))]]
    i.j.Adj.upper <- i.j.Adj.est[,lapply(.SD, function(x)
      fisher.test(table(as.numeric(GxE.LN_p.value)<threshold,
                        as.numeric(x)<threshold))$conf.int[[2]]),by=v.G,
      .SDcols=colnames(i.j.Adj.est)[grep("Joint", colnames(i.j.Adj.est))]]

    i.j.unAdj.est <- data.table(do.call(rbind,unAdj.est[args[Type==T][MAF==j][,Index]]))
    i.j.unAdj.ORs <- i.j.unAdj.est[,lapply(.SD, function(x)
      fisher.test(table(as.numeric(GxE.LN_p.value)<threshold,
                        as.numeric(x)<threshold))$estimate),by=v.G,
      .SDcols=colnames(i.j.unAdj.est)[grep("Joint", colnames(i.j.unAdj.est))]]
    i.j.unAdj.lower <- i.j.unAdj.est[,lapply(.SD, function(x)
      fisher.test(table(as.numeric(GxE.LN_p.value)<threshold,
                        as.numeric(x)<threshold))$conf.int[[1]]),by=v.G,
      .SDcols=colnames(i.j.unAdj.est)[grep("Joint", colnames(i.j.unAdj.est))]]
    i.j.unAdj.upper <- i.j.unAdj.est[,lapply(.SD, function(x)
      fisher.test(table(as.numeric(GxE.LN_p.value)<threshold,
                        as.numeric(x)<threshold))$conf.int[[2]]),by=v.G,
      .SDcols=colnames(i.j.unAdj.est)[grep("Joint", colnames(i.j.unAdj.est))]]

    y.MAX <- rbind(data.frame(i.j.Adj.upper[,lapply(.SD, function(x) max(as.numeric(x)))]),
                   data.frame(i.j.unAdj.upper[,lapply(.SD, function(x) max(as.numeric(x)))]))
    y.MAX <- y.MAX[,c("RIF_Joint","Levene_Joint")]
    y.MAX <- max(as.numeric(as.matrix(y.MAX)))
    y.MAX <- y.MAX+0.1*y.MAX

    y.MIN <- rbind(data.frame(i.j.Adj.lower[,lapply(.SD, function(x) min(as.numeric(x)))]),
                   data.frame(i.j.unAdj.lower[,lapply(.SD, function(x) min(as.numeric(x)))]))
    y.MIN <- y.MIN[,c("RIF_Joint","Levene_Joint")]
    y.MIN <- min(as.numeric(as.matrix(y.MIN)))
    y.MIN <- y.MIN-0.1*y.MIN

    plot(x=i.j.unAdj.ORs[,v.G], y=log(i.j.unAdj.ORs[,RIF_Joint]),
         type="n", xlab="", ylab="", ylim=c(0,y.MAX),xlim=c(0, 0.01),
         xaxt="n", yaxt="n", bty="l")
    title(main=bquote(paste(.(i), ", MAF =",.(j), sep="")), line=1)
    axis(2, tick=TRUE, labels=FALSE)
    axis(2, tick=FALSE, labels=TRUE, line=-0.5)
    axis(1, tick=TRUE, labels=FALSE)
    axis(1, tick=FALSE, labels=TRUE, line=-0.5)
    abline(h=1,col="black")
    lines(x=i.j.Adj.ORs[,v.G], y=i.j.Adj.ORs[,RIF_Joint],
          type="o", pch=20, col="darkmagenta", lwd=1,lty=2)
    arrows(as.numeric(i.j.Adj.ORs[,v.G]), i.j.Adj.lower[,RIF_Joint],
           as.numeric(i.j.Adj.ORs[,v.G]), i.j.Adj.upper[,RIF_Joint], length=0.02,
           angle=90, code=3,lwd=0.5)
    lines(x=i.j.Adj.ORs[,v.G], y=i.j.Adj.ORs[,Levene_Joint],
          type="o", pch=20, col="mediumseagreen", lwd=1,lty=2)
    arrows(as.numeric(i.j.Adj.ORs[,v.G]), i.j.Adj.lower[,Levene_Joint],
           as.numeric(i.j.Adj.ORs[,v.G]), i.j.Adj.upper[,Levene_Joint], length=0.02,
           angle=90, code=3,lwd=0.5)

    lines(x=i.j.unAdj.ORs[,v.G], y=i.j.unAdj.ORs[,RIF_Joint],
          type="o", pch=20, col="darkmagenta", lwd=1)
    arrows(as.numeric(i.j.unAdj.ORs[,v.G]), i.j.unAdj.lower[,RIF_Joint],
           as.numeric(i.j.unAdj.ORs[,v.G]), i.j.unAdj.upper[,RIF_Joint], length=0.02,
           angle=90, code=3,lwd=0.5)
    lines(x=i.j.unAdj.ORs[,v.G], y=i.j.unAdj.ORs[,Levene_Joint],
          type="o", pch=20, col="mediumseagreen", lwd=1)
    arrows(as.numeric(i.j.unAdj.ORs[,v.G]), i.j.unAdj.lower[,Levene_Joint],
           as.numeric(i.j.unAdj.ORs[,v.G]), i.j.unAdj.upper[,Levene_Joint], length=0.02,
           angle=90, code=3,lwd=0.5)
  }
}
dev.off()

