behav.gibbs.sampler=function(dat,ngibbs) {
set.seed(1)

#priors
alpha=1
bern.a=bern.b=1

#useful stuff
max.time=max(dat$time1)
max.SL=max(dat$SL, na.rm = T)
max.TA=max(dat$TA, na.rm = T)

#starting values
breakpt=mean(dat$time1)

#matrix to store results
store.param=matrix(NA,ngibbs,2)

for (i in 1:ngibbs){
  print(i)
  vals=samp.move(breakpt=breakpt,max.time=max.time,dat=dat,
                    alpha=alpha,bern.a=bern.a,bern.b=bern.b,
                    max.SL=max.SL,max.TA=max.TA)   
  
  breakpt=vals[[1]]
  
  #store results
  store.param[i,]=c(length(vals[[1]]), vals[[2]])

  }
  list(breakpt=breakpt,store.param=store.param)
}