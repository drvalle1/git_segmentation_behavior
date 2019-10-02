rm(list=ls(all=TRUE))

library(ggplot2)
library(dplyr)
library(coda)

set.seed(1)

source('gibbs functions2.R')
source('helper functions.R')
dat<- read.csv("Snail Kite Gridded Data.csv", header = T, sep = ",")
names(dat)[7]<- "dist"  #change to generic form

#discretize step length, turning angle, and turning angle autocorr by ID
#only subset obs w/ 1 hr time step; removal of 5815 total obs
dat1<- dat %>% filter(id==1) %>% mutate(obs = 1:nrow(.)) %>% filter(dt==3600) %>% mutate(time1 = 1:nrow(.)) 
dat12<- dat %>% filter(id==12) %>% mutate(obs = 1:nrow(.)) %>% filter(dt==3600) %>% mutate(time1 = 1:nrow(.)) 
dat19<- dat %>% filter(id==19) %>% mutate(obs = 1:nrow(.)) %>% filter(dt==3600) %>% mutate(time1 = 1:nrow(.)) 
dat27<- dat %>% filter(id==27) %>% mutate(obs = 1:nrow(.)) %>% filter(dt==3600) %>% mutate(time1 = 1:nrow(.))


dat1<- dat1 %>% assign.rel_angle.bin() %>% assign.dist.bin() %>% chng.rel_angle.sign()
dat12<- dat12 %>% assign.rel_angle.bin() %>% assign.dist.bin() %>% chng.rel_angle.sign()
dat19<- dat19 %>% assign.rel_angle.bin() %>% assign.dist.bin() %>% chng.rel_angle.sign()
dat27<- dat27 %>% assign.rel_angle.bin() %>% assign.dist.bin() %>% chng.rel_angle.sign()


dat.list<- list(dat1=dat1,dat12=dat12,dat19=dat19,dat27=dat27)

#for (i in 1:length(dat.list)) {   WORK ON LOOPING OVER IDs

#################################
#### Run Gibbs Sampler by ID ####
#################################

### ID 1

#priors
alpha=0.01
bern.a=bern.b=1

#useful stuff
max.time=max(dat1$time1)
max.SL=max(dat1$SL, na.rm = T)
max.TA=max(dat1$TA, na.rm = T)

#starting values
breakpt=mean(dat1$time1)

#number of iterations
ngibbs=20000

#matrix to store results
store.param=matrix(NA,ngibbs,2)

for (i in 1:ngibbs){
  print(i)
  vals=samp.move(breakpt=breakpt,max.time=max.time,dat=dat1,
                    alpha=alpha,bern.a=bern.a,bern.b=bern.b,
                    max.SL=max.SL,max.TA=max.TA)   
  
  breakpt=vals[[1]]
  
  #store results
  store.param[i,]=c(length(vals[[1]]), vals[[2]])
}
###Takes ~ 6 min to run 20000 iterations; identified 7 breakpoints


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








### ID 12

#priors
alpha=0.01
bern.a=bern.b=1

#useful stuff
max.time=max(dat12$time1)
max.SL=max(dat12$SL, na.rm = T)
max.TA=max(dat12$TA, na.rm = T)

#starting values
breakpt=mean(dat12$time1)

#number of iterations
ngibbs=40000

#matrix to store results
store.param=matrix(NA,ngibbs,2)

for (i in 1:ngibbs){
  print(i)
  vals=samp.move(breakpt=breakpt,max.time=max.time,dat=dat12,
                 alpha=alpha,bern.a=bern.a,bern.b=bern.b,
                 max.SL=max.SL,max.TA=max.TA)   
  
  breakpt=vals[[1]]
  
  #store results
  store.param[i,]=c(length(vals[[1]]), vals[[2]])
}
###Takes ~ 6 min to run 20000 iterations; identified 1 breakpoint


length(breakpt)
#write.csv(breakpt, "ID12 Breakpoints (Behavior).csv", row.names = F)

#PLOT
png("ID 12 behav seg plot.png", width = 8, height = 6, units = "in", res = 400)
par(mfrow=c(3,1))
plot(dat12$SL, xlab = "", ylab = "SL"); abline(v=breakpt,col='red')
plot(dat12$TA, xlab = "", ylab = "TA"); abline(v=breakpt,col='red')
plot(dat12$TAA, xlab = "Time", ylab = "TAA"); abline(v=breakpt,col='red')
dev.off()


#IMAGE
behav.list<- behav.seg.image(dat12)

png("ID 12 behav seg image.png", width = 8, height = 6, units = "in", res = 400)
par(mfrow=c(3,1))
image(behav.list$SL, xlab = "", ylab = "SL"); abline(v=breakpt/nrow(dat12),col='blue'); title(main="ID 12 (n=1)")
image(behav.list$TA, xlab = "", ylab = "TA"); abline(v=breakpt/nrow(dat12),col='blue')
image(behav.list$TAA, xlab = "Time", ylab = "TAA"); abline(v=breakpt/nrow(dat12),col='blue')
dev.off()



## MCMC traceplots
store.param.mcmc<- as.mcmc(store.param)

par(mfrow=c(1,1))
#breakpts
traceplot(store.param.mcmc[,1]); title(y="# of Breakpoints", main = "ID 12")

#LML
traceplot(store.param.mcmc[,2]); title(y="Log Marginal Likelihood", main = "ID 12")








### ID 19

#priors
alpha=0.01
bern.a=bern.b=1

#useful stuff
max.time=max(dat19$time1)
max.SL=max(dat19$SL, na.rm = T)
max.TA=max(dat19$TA, na.rm = T)

#starting values
breakpt=mean(dat19$time1)

#number of iterations
ngibbs=20000

#matrix to store results
store.param=matrix(NA,ngibbs,2)

for (i in 1:ngibbs){
  print(i)
  vals=samp.move(breakpt=breakpt,max.time=max.time,dat=dat19,
                 alpha=alpha,bern.a=bern.a,bern.b=bern.b,
                 max.SL=max.SL,max.TA=max.TA)   
  
  breakpt=vals[[1]]
  
  #store results
  store.param[i,]=c(length(vals[[1]]), vals[[2]])
}
###Takes ~ 2 min to run 20000 iterations; identified 1 breakpoint


length(breakpt)
#write.csv(breakpt, "ID19 Breakpoints (Behavior).csv", row.names = F)

png("ID 19 behav seg plot.png", width = 8, height = 6, units = "in", res = 400)
par(mfrow=c(3,1))
plot(dat19$SL, xlab = "", ylab = "SL"); abline(v=breakpt,col='red')
plot(dat19$TA, xlab = "", ylab = "TA"); abline(v=breakpt,col='red')
plot(dat19$TAA, xlab = "Time", ylab = "TAA"); abline(v=breakpt,col='red')
dev.off()



## MCMC traceplots
store.param.mcmc<- as.mcmc(store.param)

par(mfrow=c(1,1))
#breakpts
traceplot(store.param.mcmc[,1]); title(y="# of Breakpoints", main = "ID 19")

#LML
traceplot(store.param.mcmc[,2]); title(y="Log Marginal Likelihood", main = "ID 19")








### ID 27

#priors
alpha=0.01
bern.a=bern.b=1

#useful stuff
max.time=max(dat27$time1)
max.SL=max(dat27$SL, na.rm = T)
max.TA=max(dat27$TA, na.rm = T)

#starting values
breakpt=mean(dat27$time1)

#number of iterations
ngibbs=20000

#matrix to store results
store.param=matrix(NA,ngibbs,2)

for (i in 1:ngibbs){
  print(i)
  vals=samp.move(breakpt=breakpt,max.time=max.time,dat=dat27,
                 alpha=alpha,bern.a=bern.a,bern.b=bern.b,
                 max.SL=max.SL,max.TA=max.TA)   
  
  breakpt=vals[[1]]
  
  #store results
  store.param[i,]=c(length(vals[[1]]), vals[[2]])
}
###Takes ~ 1 min to run 20000 iterations; identified 1 breakpoint


length(breakpt)
#write.csv(breakpt, "ID27 Breakpoints (Behavior).csv", row.names = F)

png("ID 27 behav seg plot.png", width = 8, height = 6, units = "in", res = 400)
par(mfrow=c(3,1))
plot(dat27$SL, xlab = "", ylab = "SL"); abline(v=breakpt,col='red')
plot(dat27$TA, xlab = "", ylab = "TA"); abline(v=breakpt,col='red')
plot(dat27$TAA, xlab = "Time", ylab = "TAA"); abline(v=breakpt,col='red')
dev.off()


## MCMC traceplots
store.param.mcmc<- as.mcmc(store.param)

par(mfrow=c(1,1))
#breakpts
traceplot(store.param.mcmc[,1]); title(y="# of Breakpoints", main = "ID 27")

#LML
traceplot(store.param.mcmc[,2]); title(y="Log Marginal Likelihood", main = "ID 27")

