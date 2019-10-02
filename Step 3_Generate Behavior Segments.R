rm(list=ls(all=TRUE))

library(ggplot2)
library(dplyr)
library(coda)


source('gibbs functions2.R')
source('helper functions.R')
source('gibbs sampler2.R')

dat<- read.csv("Snail Kite Gridded Data.csv", header = T, sep = ",")
names(dat)[7]<- "dist"  #change to generic form
behav.list<- behav.prep(dat=dat, tstep = 3600)  #add move params and filter by 3600 s interval


#################################
#### Run Gibbs Sampler by ID ####
#################################

### ID 1
dat1.res<- behav.gibbs.sampler(behav.list$`1`,50000)

length(dat1.res$breakpt)
#write.csv(dat1.res$breakpt, "ID1 Breakpoints (Behavior).csv", row.names = F)

#PLOT
png("ID 1 behav seg plot.png", width = 8, height = 6, units = "in", res = 400)
par(mfrow=c(3,1))
plot(behav.list$`1`$SL, xlab = "", ylab = "SL"); abline(v=dat1.res$breakpt,col='red'); title(main="ID 1 (n=26)")
plot(behav.list$`1`$TA, xlab = "", ylab = "TA"); abline(v=dat1.res$breakpt,col='red')
plot(behav.list$`1`$TAA, xlab = "Time", ylab = "TAA"); abline(v=dat1.res$breakpt,col='red')
dev.off()


#IMAGE
behav.heat<- behav.seg.image(behav.list$`1`)

png("ID 1 behav seg image.png", width = 8, height = 6, units = "in", res = 400)
par(mfrow=c(3,1))
image(behav.heat$SL, xlab = "", ylab = "SL"); abline(v=dat1.res$breakpt/nrow(behav.list$`1`),col='blue'); title(main="ID 1 (n=26)")
image(behav.heat$TA, xlab = "", ylab = "TA"); abline(v=dat1.res$breakpt/nrow(behav.list$`1`),col='blue')
image(behav.heat$TAA, xlab = "Time", ylab = "TAA"); abline(v=dat1.res$breakpt/nrow(behav.list$`1`),col='blue')
dev.off()


## MCMC traceplots
par(mfrow=c(1,1))
#breakpts
plot(dat1.res$store.param[,1], type = "l", xlab = "Iterations", ylab = "# of Breakpoints",
     main = "ID 1")
#LML
plot(dat1.res$store.param[,2], type = "l", xlab = "Iterations", ylab = "Log Marginal Likelihood",
     main = "ID 1")




### ID 12
dat12.res<- behav.gibbs.sampler(behav.list$`12`,50000)

length(dat12.res$breakpt)
#write.csv(dat12.res$breakpt, "ID12 Breakpoints (Behavior).csv", row.names = F)

#PLOT
png("ID 12 behav seg plot.png", width = 8, height = 6, units = "in", res = 400)
par(mfrow=c(3,1))
plot(behav.list$`12`$SL, xlab = "", ylab = "SL"); abline(v=dat12.res$breakpt,col='red'); title(main="ID 12 (n=31)")
plot(behav.list$`12`$TA, xlab = "", ylab = "TA"); abline(v=dat12.res$breakpt,col='red')
plot(behav.list$`12`$TAA, xlab = "Time", ylab = "TAA"); abline(v=dat12.res$breakpt,col='red')
dev.off()


#IMAGE
behav.heat<- behav.seg.image(behav.list$`12`)

png("ID 12 behav seg image.png", width = 8, height = 6, units = "in", res = 400)
par(mfrow=c(3,1))
image(behav.heat$SL, xlab = "", ylab = "SL"); abline(v=dat12.res$breakpt/nrow(behav.list$`12`),col='blue'); title(main="ID 12 (n=31)")
image(behav.heat$TA, xlab = "", ylab = "TA"); abline(v=dat12.res$breakpt/nrow(behav.list$`12`),col='blue')
image(behav.heat$TAA, xlab = "Time", ylab = "TAA"); abline(v=dat12.res$breakpt/nrow(behav.list$`12`),col='blue')
dev.off()


## MCMC traceplots
par(mfrow=c(1,1))
#breakpts
plot(dat12.res$store.param[,1], type = "l", xlab = "Iterations", ylab = "# of Breakpoints",
     main = "ID 12")
#LML
plot(dat12.res$store.param[,2], type = "l", xlab = "Iterations", ylab = "Log Marginal Likelihood",
     main = "ID 12")




### ID 19
dat19.res<- behav.gibbs.sampler(behav.list$`19`,50000)

length(dat19.res$breakpt)
#write.csv(dat19.res$breakpt, "ID19 Breakpoints (Behavior).csv", row.names = F)

#PLOT
png("ID 19 behav seg plot.png", width = 8, height = 6, units = "in", res = 400)
par(mfrow=c(3,1))
plot(behav.list$`19`$SL, xlab = "", ylab = "SL"); abline(v=dat19.res$breakpt,col='red'); title(main="ID 19 (n=5)")
plot(behav.list$`19`$TA, xlab = "", ylab = "TA"); abline(v=dat19.res$breakpt,col='red')
plot(behav.list$`19`$TAA, xlab = "Time", ylab = "TAA"); abline(v=dat19.res$breakpt,col='red')
dev.off()


#IMAGE
behav.heat<- behav.seg.image(behav.list$`19`)

png("ID 19 behav seg image.png", width = 8, height = 6, units = "in", res = 400)
par(mfrow=c(3,1))
image(behav.heat$SL, xlab = "", ylab = "SL"); abline(v=dat19.res$breakpt/nrow(behav.list$`19`),col='blue'); title(main="ID 19 (n=5)")
image(behav.heat$TA, xlab = "", ylab = "TA"); abline(v=dat19.res$breakpt/nrow(behav.list$`19`),col='blue')
image(behav.heat$TAA, xlab = "Time", ylab = "TAA"); abline(v=dat19.res$breakpt/nrow(behav.list$`19`),col='blue')
dev.off()


## MCMC traceplots
par(mfrow=c(1,1))
#breakpts
plot(dat19.res$store.param[,1], type = "l", xlab = "Iterations", ylab = "# of Breakpoints",
     main = "ID 19")
#LML
plot(dat19.res$store.param[,2], type = "l", xlab = "Iterations", ylab = "Log Marginal Likelihood",
     main = "ID 19")




### ID 27
dat27.res<- behav.gibbs.sampler(behav.list$`27`,50000)

length(dat27.res$breakpt)
#write.csv(dat27.res$breakpt, "ID27 Breakpoints (Behavior).csv", row.names = F)

#PLOT
png("ID 27 behav seg plot.png", width = 8, height = 6, units = "in", res = 400)
par(mfrow=c(3,1))
plot(behav.list$`27`$SL, xlab = "", ylab = "SL"); abline(v=dat27.res$breakpt,col='red'); title(main="ID 27 (n=2)")
plot(behav.list$`27`$TA, xlab = "", ylab = "TA"); abline(v=dat27.res$breakpt,col='red')
plot(behav.list$`27`$TAA, xlab = "Time", ylab = "TAA"); abline(v=dat27.res$breakpt,col='red')
dev.off()


#IMAGE
behav.heat<- behav.seg.image(behav.list$`27`)

png("ID 27 behav seg image.png", width = 8, height = 6, units = "in", res = 400)
par(mfrow=c(3,1))
image(behav.heat$SL, xlab = "", ylab = "SL"); abline(v=dat27.res$breakpt/nrow(behav.list$`27`),col='blue'); title(main="ID 27 (n=2)")
image(behav.heat$TA, xlab = "", ylab = "TA"); abline(v=dat27.res$breakpt/nrow(behav.list$`27`),col='blue')
image(behav.heat$TAA, xlab = "Time", ylab = "TAA"); abline(v=dat27.res$breakpt/nrow(behav.list$`27`),col='blue')
dev.off()


## MCMC traceplots
par(mfrow=c(1,1))
#breakpts
plot(dat27.res$store.param[,1], type = "l", xlab = "Iterations", ylab = "# of Breakpoints",
     main = "ID 27")
#LML
plot(dat27.res$store.param[,2], type = "l", xlab = "Iterations", ylab = "Log Marginal Likelihood",
     main = "ID 27")










#############################################
#### Run Gibbs Sampler w Velocity Params ####
#############################################

behav.list2<- behav.prep_veloc(dat=dat, tstep = 3600)  #add move params and filter by 3600 s interval



### ID 19
dat19.res2<- behav.gibbs.sampler_veloc(behav.list2$`19`,50000)

length(dat19.res2$breakpt)
#write.csv(dat19.res2$breakpt, "ID19 Breakpoints (Behavior) Vp_Vt.csv", row.names = F)

#PLOT
png("ID 19 behav seg plot Vp_Vt.png", width = 8, height = 6, units = "in", res = 400)
par(mfrow=c(2,1))
plot(behav.list2$`19`$Vp, xlab = "", ylab = "Vp"); abline(v=dat19.res2$breakpt,col='red'); title(main="ID 19 (n=3)")
plot(behav.list2$`19`$Vt, xlab = "", ylab = "Vt"); abline(v=dat19.res2$breakpt,col='red')
dev.off()


#IMAGE
behav.heat<- behav.seg.image_veloc(behav.list2$`19`)

png("ID 19 behav seg image Vp_Vt.png", width = 8, height = 6, units = "in", res = 400)
par(mfrow=c(2,1))
image(behav.heat$Vp, xlab = "", ylab = "Vp"); abline(v=dat19.res2$breakpt/nrow(behav.list2$`19`),
                                                     col='blue'); title(main="ID 19 (n=3)")
image(behav.heat$Vt, xlab = "", ylab = "Vt"); abline(v=dat19.res2$breakpt/nrow(behav.list2$`19`),
                                                     col='blue')
dev.off()


## MCMC traceplots
par(mfrow=c(1,1))
#breakpts
plot(dat19.res2$store.param[,1], type = "l", xlab = "Iterations", ylab = "# of Breakpoints",
     main = "ID 19")
#LML
plot(dat19.res2$store.param[,2], type = "l", xlab = "Iterations", ylab = "Log Marginal Likelihood",
     main = "ID 19")



### ID 27
dat27.res2<- behav.gibbs.sampler_veloc(behav.list2$`27`,50000)

length(dat27.res2$breakpt)
#write.csv(dat27.res2$breakpt, "ID27 Breakpoints (Behavior) Vp_Vt.csv", row.names = F)

#PLOT
png("ID 27 behav seg plot Vp_Vt.png", width = 8, height = 6, units = "in", res = 400)
par(mfrow=c(2,1))
plot(behav.list2$`27`$Vp, xlab = "", ylab = "Vp"); abline(v=dat27.res2$breakpt,col='red'); title(main="ID 27 (n=3)")
plot(behav.list2$`27`$Vt, xlab = "", ylab = "Vt"); abline(v=dat27.res2$breakpt,col='red')
dev.off()


#IMAGE
behav.heat<- behav.seg.image_veloc(behav.list2$`27`)

png("ID 27 behav seg image Vp_Vt.png", width = 8, height = 6, units = "in", res = 400)
par(mfrow=c(2,1))
image(behav.heat$Vp, xlab = "", ylab = "Vp"); abline(v=dat27.res2$breakpt/nrow(behav.list2$`27`),
                                                     col='blue'); title(main="ID 27 (n=3)")
image(behav.heat$Vt, xlab = "", ylab = "Vt"); abline(v=dat27.res2$breakpt/nrow(behav.list2$`27`),
                                                     col='blue')
dev.off()


## MCMC traceplots
par(mfrow=c(1,1))
#breakpts
plot(dat27.res2$store.param[,1], type = "l", xlab = "Iterations", ylab = "# of Breakpoints",
     main = "ID 27")
#LML
plot(dat27.res2$store.param[,2], type = "l", xlab = "Iterations", ylab = "Log Marginal Likelihood",
     main = "ID 27")
