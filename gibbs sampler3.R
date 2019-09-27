rm(list=ls(all=TRUE))

library(ggplot2)
library(dplyr)
library(coda)

set.seed(1)

source('gibbs functions3.R')
dat<- read.csv("Snail Kite Gridded Data.csv", header = T, sep = ",")
names(dat)[7]<- "dist" #change to generic form

#discretize step length, turning angle, and turning angle autocorr by ID
dat1<- dat %>% filter(id==1) %>% veloc.persist() %>% veloc.turn
dat12<- dat %>% filter(id==12) %>% veloc.persist() %>% veloc.turn()
dat19<- dat %>% filter(id==19) %>% veloc.persist() %>% veloc.turn()
dat27<- dat %>% filter(id==27) %>% veloc.persist() %>% veloc.turn()

dat1$time1<- 1:nrow(dat1)
dat12$time1<- 1:nrow(dat12)
dat19$time1<- 1:nrow(dat19)
dat27$time1<- 1:nrow(dat27)


#################################
#### Run Gibbs Sampler by ID ####
#################################

### ID 1

#priors
alpha=0.01

#useful stuff
max.time=max(dat1$time1)
max.Vp=max(dat1$Vp, na.rm = T)
max.Vt=max(dat1$Vt, na.rm = T)

#starting values
breakpt=mean(dat1$time1)

#number of iterations
ngibbs=10000

#matrix to store results
store.param=matrix(NA,ngibbs,2)

for (i in 1:ngibbs){
  print(i)
  vals=samp.move(breakpt=breakpt,max.time=max.time,dat=dat1,
                  alpha=alpha,max.Vp=max.Vp,max.Vt=max.Vt)   
  
  breakpt=vals[[1]]
  
  #store results
  store.param[i,]=c(length(vals[[1]]), vals[[2]])
}
###Takes ~ 4 min to run 10000 iterations; identified 7 breakpoints


length(breakpt)
#write.csv(breakpt, "ID1 Breakpoints (Behavior).csv", row.names = F)

par(mfrow=c(2,1))
plot(dat1$Vp); abline(v=breakpt,col='red')
plot(dat1$Vt); abline(v=breakpt,col='red')


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

#useful stuff
max.time=max(dat12$time1)
max.SL=max(dat12$SL, na.rm = T)
max.TA=max(dat12$TA, na.rm = T)
max.Vp=max(dat12$Vp, na.rm = T)

#starting values
breakpt=mean(dat12$time1)

#number of iterations
ngibbs=10000

#matrix to store results
store.param=matrix(NA,ngibbs,2)

for (i in 1:ngibbs){
  print(i)
  vals=samp.move2(breakpt=breakpt,max.time=max.time,dat=dat12,
                  alpha=alpha,max.Vp=max.Vp,max.Vt=max.Vt)   
  
  breakpt=vals[[1]]
  
  #store results
  store.param[i,]=c(length(vals[[1]]), vals[[2]])
}
###Takes ~ 2 min to run 10000 iterations; identified 5 breakpoints


length(breakpt)
#write.csv(breakpt, "ID12 Breakpoints (Behavior).csv", row.names = F)

par(mfrow=c(3,1))
plot(dat12$SL); abline(v=breakpt,col='red')
plot(dat12$TA); abline(v=breakpt,col='red')
plot(dat12$Vp); abline(v=breakpt,col='red')



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

#useful stuff
max.time=max(dat19$time1)
max.SL=max(dat19$SL, na.rm = T)
max.TA=max(dat19$TA, na.rm = T)
max.Vp=max(dat19$Vp, na.rm = T)

#starting values
breakpt=mean(dat19$time1)

#number of iterations
ngibbs=10000

#matrix to store results
store.param=matrix(NA,ngibbs,2)

for (i in 1:ngibbs){
  print(i)
  vals=samp.move2(breakpt=breakpt,max.time=max.time,dat=dat19,
                  alpha=alpha,max.Vp=max.Vp,max.Vt=max.Vt)   
  
  breakpt=vals[[1]]
  
  #store results
  store.param[i,]=c(length(vals[[1]]), vals[[2]])
}
###Takes ~ 2 min to run 10000 iterations; identified 2 breakpoints


length(breakpt)
#write.csv(breakpt, "ID19 Breakpoints (Behavior).csv", row.names = F)

par(mfrow=c(3,1))
plot(dat19$SL); abline(v=breakpt,col='red')
plot(dat19$TA); abline(v=breakpt,col='red')
plot(dat19$Vp); abline(v=breakpt,col='red')



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

#useful stuff
max.time=max(dat27$time1)
max.SL=max(dat27$SL, na.rm = T)
max.TA=max(dat27$TA, na.rm = T)
max.Vp=max(dat27$Vp, na.rm = T)

#starting values
breakpt=mean(dat27$time1)

#number of iterations
ngibbs=10000

#matrix to store results
store.param=matrix(NA,ngibbs,2)

for (i in 1:ngibbs){
  print(i)
  vals=samp.move2(breakpt=breakpt,max.time=max.time,dat=dat27,
                  alpha=alpha,max.Vp=max.Vp,max.Vt=max.Vt)   
  
  breakpt=vals[[1]]
  
  #store results
  store.param[i,]=c(length(vals[[1]]), vals[[2]])
}
###Takes ~ 1 min to run 10000 iterations; identified 1 breakpoint


length(breakpt)
#write.csv(breakpt, "ID27 Breakpoints (Behavior).csv", row.names = F)

par(mfrow=c(3,1))
plot(dat27$SL); abline(v=breakpt,col='red')
plot(dat27$TA); abline(v=breakpt,col='red')
plot(dat27$Vp); abline(v=breakpt,col='red')



## MCMC traceplots
store.param.mcmc<- as.mcmc(store.param)

par(mfrow=c(1,1))
#breakpts
traceplot(store.param.mcmc[,1]); title(y="# of Breakpoints", main = "ID 27")

#LML
traceplot(store.param.mcmc[,2]); title(y="Log Marginal Likelihood", main = "ID 27")

