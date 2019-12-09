behav.gibbs.sampler=function(dat,ngibbs) {
  set.seed(1)
  
  uni.id=unique(dat$id)
  
  #priors
  alpha=1
  
  #useful stuff
  max.time=max(dat$time1)
  max.SL=max(dat$SL, na.rm = T)
  max.TA=max(dat$TA, na.rm = T)
  
  #starting values
  breakpt=mean(dat$time1)
  
  #to store results
  res.brks=vector("list", ngibbs)
  res.LML=matrix(NA,1,(ngibbs+1))
  res.nbrks=matrix(NA,1,(ngibbs+1))
  store.param=matrix(NA,ngibbs,2)
  
  for (i in 1:ngibbs){
    vals=samp.move(breakpt=breakpt,max.time=max.time,dat=dat,
                   alpha=alpha,max.SL=max.SL,max.TA=max.TA)   
    
    breakpt=vals[[1]]
    
    #store results
    res.brks[[i]]<- breakpt
    store.param[i,]=c(length(breakpt), vals[[2]])
    
  }
  
  tmp=store.param[,1]
  res.nbrks[1,]=c(uni.id,tmp)
  colnames(res.nbrks)<- c('id', paste0("Iter_",1:ngibbs))
  
  tmp=store.param[,2]
  res.LML[1,]=c(uni.id,tmp)
  colnames(res.LML)<- c('id', paste0("Iter_",1:ngibbs))
  
  list(breakpt=res.brks, nbrks=res.nbrks, LML=res.LML)
}
#----------------------------------------------------
behavior_segment=function(data, ngibbs) {
  
  tic()  #start timer
  mod<- future_map(data, function(x) behav.gibbs.sampler(dat = x, ngibbs = ngibbs), .progress = TRUE)
  toc()  #provide elapsed time
  
  
  brkpts<- map(mod, 1)  #create list of all sets breakpoints by ID
  
  nbrks<- map_dfr(mod, 2) %>% t() %>% data.frame()  #create DF of number of breakpoints by ID
  names(nbrks)<- c('id', paste0("Iter_",1:ngibbs))

  LML<- map_dfr(mod, 3) %>% t() %>% data.frame()  #create DF of LML by ID
  names(LML)<- c('id', paste0("Iter_",1:ngibbs))


  list(brkpts = brkpts, nbrks = nbrks, LML = LML)
}
