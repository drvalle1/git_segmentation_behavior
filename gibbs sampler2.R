rm(list=ls(all=TRUE))
set.seed(1)

setwd('U:\\GIT_models\\git_segmentation_behavior')
source('gibbs functions2.R')
dat=read.csv('fake data.csv',as.is=T)

#discretize step length, turning angle, and turning angle autocorr
dat<- dat %>% assign.rel_angle.bin() %>% assign.dist.bin() %>% chng.rel_angle.sign()

#priors
alpha=0.01
bern.a=bern.b=1

#useful stuff
max.time=max(dat$time1)
max.SL=max(dat$SL, na.rm = T)
max.TA=max(dat$TA, na.rm = T)

#starting values
breakpt=mean(dat$time1)

ngibbs=10000
for (i in 1:ngibbs){
  print(i)
  vals=samp.move(breakpt=breakpt,max.time=max.time,dat=dat,
                    alpha=alpha,bern.a=bern.a,bern.b=bern.b,
                    max.SL=max.SL,max.TA=max.TA)   
  
  breakpt=vals[[1]]
  
  #store results
  store.param[i,]=c(length(vals[[1]]), vals[[2]])
}
length(breakpt)
abline(v=breakpt,lty=3,col='green')
