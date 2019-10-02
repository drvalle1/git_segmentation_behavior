rm(list=ls(all=TRUE))

library(ggplot2)
library(dplyr)
library(coda)


source('gibbs functions2.R')
source('helper functions.R')
source('gibbs sampler2.R')
dat<- read.csv("Snail Kite Gridded Data.csv", header = T, sep = ",")
names(dat)[7]<- "dist"  #change to generic form

behav.list<- behav.prep(dat)


#################################
#### Run Gibbs Sampler by ID ####
#################################

### ID 1


length(breakpt)
#write.csv(breakpt, "ID1 Breakpoints (Behavior).csv", row.names = F)

#PLOT
png("ID 1 behav seg plot.png", width = 8, height = 6, units = "in", res = 400)
par(mfrow=c(3,1))
plot(dat1$SL, xlab = "", ylab = "SL"); abline(v=breakpt,col='red'); title(main="ID 1 (n=7)")
plot(dat1$TA, xlab = "", ylab = "TA"); abline(v=breakpt,col='red')
plot(dat1$TAA, xlab = "Time", ylab = "TAA"); abline(v=breakpt,col='red')
dev.off()


#IMAGE
behav.list<- behav.seg.image(dat1)

png("ID 1 behav seg image.png", width = 8, height = 6, units = "in", res = 400)
par(mfrow=c(3,1))
image(behav.list$SL, xlab = "", ylab = "SL"); abline(v=breakpt/nrow(dat1),col='blue'); title(main="ID 1 (n=7)")
image(behav.list$TA, xlab = "", ylab = "TA"); abline(v=breakpt/nrow(dat1),col='blue')
image(behav.list$TAA, xlab = "Time", ylab = "TAA"); abline(v=breakpt/nrow(dat1),col='blue')
dev.off()



## MCMC traceplots
store.param.mcmc<- as.mcmc(store.param)

par(mfrow=c(1,1))
#breakpts
traceplot(store.param.mcmc[,1]); title(y="# of Breakpoints", main = "ID 1")

#LML
traceplot(store.param.mcmc[,2]); title(y="Log Marginal Likelihood", main = "ID 1")
