#discretize step length, turning angle, and turning angle autocorr by ID
#only subset obs w/ 1 hr (3600 s) time step; removal of 5815 total obs
behav.prep=function(dat,tstep) {
  id<- unique(dat$id)
  n=length(id)
  dat.behav<- vector("list", n)
  names(dat.behav)<- id
  
  for (i in 1:length(id)) {
    dat.behav[[i]]<- dat %>% filter(id==id[i]) %>% mutate(obs = 1:nrow(.)) %>% filter(dt==tstep) %>%
                  mutate(time1 = 1:nrow(.)) %>% assign.rel_angle.bin() %>% assign.dist.bin() %>%
                  chng.rel_angle.sign()
  }
  dat.behav
}
#---------------------------------------------
behav.seg.image=function(dat) {  #Transform single var vectors into pres/abs matrices for heatmap
  SL<- matrix(0, nrow(dat), length(unique(dat$SL)))
  for (i in 1:nrow(dat)){
    SL[i,dat$SL[i]]=1
  }
  
  TA<- matrix(0, nrow(dat), length(unique(dat[!is.na(dat$TA),"TA"])))
  for (i in 1:nrow(dat)){
    TA[i,dat$TA[i]]=1
  }
  
  TAA<- matrix(0, nrow(dat), length(unique(dat[!is.na(dat$TAA),"TAA"])))
  for (i in 1:nrow(dat)){
    TAA[i,dat$TAA[i]+1]=1
  }
  
  behav.list<- list(SL=SL,TA=TA,TAA=TAA)
  
  behav.list
}
#---------------------------------------------

