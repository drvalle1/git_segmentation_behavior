rm(list=ls(all=TRUE))
library('MCMCpack')
set.seed(2)

nobs=1000
nseg=10
tmp=runif(nseg)
prob=tmp/sum(tmp); prob
partition=rmultinom(1,size=nobs,prob=prob)
seg.index=rep(1:nseg,times=partition)
nloc=100

#step length (multinomial with 4 categories)
prob.SL=matrix(c(0.6,0.2,0.2,0.0,
                 0.2,0.1,0.5,0.2,
                 0.0,0.1,0.2,0.7),3,4,byrow=T) 
#apply(prob.SL,1,sum) #basic checking
ind=sample(3,size=nseg,replace=T)
prob.SL1=prob.SL[ind,]

#turning angle  (multinomial with 4 categories)
prob.TA=matrix(c(0.6,0.3,0.1,0.0,
                 0.3,0.3,0.2,0.2,
                 0.0,0.1,0.7,0.2,
                 0.0,0.1,0.2,0.7),4,4,byrow=T)
#apply(prob.TA,1,sum)  #basic checking
ind=sample(4,size=nseg,replace=T)
prob.TA1=prob.TA[ind,]

#turning angle autocorrelation (bernoulli)
prob.TAA=c(0.3,0.9)
ind=sample(2,size=nseg,replace=T)
prob.TAA1=prob.TAA[ind]

#generate observations
obs=matrix(NA,nobs,3)
for (i in 1:nobs){
  tmp.SL=rmultinom(1,size=1,prob=prob.SL1[seg.index[i],])
  tmp.TA=rmultinom(1,size=1,prob=prob.TA1[seg.index[i],])
  tmp.TAA=rbinom(1,size=1,prob=prob.TAA1[seg.index[i]])
  obs[i,]=c(which(tmp.SL==1),which(tmp.TA==1),tmp.TAA)
}

#plot these data
par(mfrow=c(3,1))
for (i in 1:3){
  plot(obs[,i])
  abline(v=cumsum(partition))
}

obs1=data.frame(SL=obs[,1],TA=obs[,2],TAA=obs[,3])
obs1$time1=1:nobs

setwd('U:\\GIT_models\\git_segmentation_behavior')
write.csv(obs1,'fake data.csv',row.names=F)
