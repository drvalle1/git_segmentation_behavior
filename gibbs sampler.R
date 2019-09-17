rm(list=ls(all=TRUE))
set.seed(1)

setwd('U:\\GIT_models\\git_segmentation_behavior')
source('gibbs functions.R')
dat=read.csv('fake data.csv',as.is=T)

#priors
alpha=0.01
bern.a=bern.b=1

#useful stuff
max.time=max(dat$time1)
max.SL=max(dat$SL)
max.TA=max(dat$TA)

#starting values
breakpt=mean(dat$time1)

ngibbs=10000
for (i in 1:ngibbs){
  print(i)
  breakpt=samp.move(breakpt=breakpt,max.time=max.time,dat=dat,
                    alpha=alpha,bern.a=bern.a,bern.b=bern.b,
                    max.SL=max.SL,max.TA=max.TA)   
}
length(breakpt)
abline(v=breakpt,lty=3,col='green')