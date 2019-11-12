assign.rel_angle.bin=function(dat){  #creates 8 bins
  angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)
  dat$TA<- NA
  
  for(i in 1:length(angle.bin.lims)) {
    tmp=which(dat$rel.angle >= angle.bin.lims[i] & dat$rel.angle < angle.bin.lims[i+1])
    dat[tmp,"TA"]=i
  }
  dat
}
#---------------------------------------------
assign.dist.bin=function(dat){  #creates 6 bins
  max.dist=max(dat$dist, na.rm = T) #using value from entire dataset, not specific time segment
  upper90.thresh=as.numeric(quantile(dat$dist, 0.90, na.rm=T)) #using value from entire dataset
  
  dist.bin.lims=seq(from=0, to=upper90.thresh, length.out = 6)
  dist.bin.lims=c(dist.bin.lims, max.dist)
  dat$SL<- NA
  #names(dat)[7]<- "dist"
  
  for(i in 1:length(dist.bin.lims)) {
    tmp=which(dat$dist >= dist.bin.lims[i] & dat$dist < dist.bin.lims[i+1])
    dat[tmp,"SL"]=i
  }
  tmp=which(dat$dist == max.dist)
  dat[tmp,"SL"]=6
  
  dat
}
#---------------------------------------------
chng.rel_angle.sign=function(dat){
  dat$TAA<- NA
  
  for(i in 1:(nrow(dat)-1)) {
    dat$TAA[i+1]<- ifelse(sign(dat$rel.angle[i]) == sign(dat$rel.angle[i+1]), 1, 0)
  }
  dat
}
#---------------------------------------------
#discretize step length, turning angle, and turning angle autocorr by ID
#only subset obs w/ 1 hr (3600 s) time step; removal of 5815 total obs
behav.prep=function(dat,tstep) {
  id<- unique(dat$id)
  n=length(id)
  behav.list<- vector("list", n)
  names(behav.list)<- id
  
  for (i in 1:length(id)) {
    behav.list[[i]]<- dat[dat$id==id[i],] %>% mutate(obs = 1:nrow(.)) %>% filter(dt==tstep) %>%
                  mutate(time1 = 1:nrow(.)) %>% assign.rel_angle.bin() %>% assign.dist.bin() %>%
                  chng.rel_angle.sign()
  }
  behav.list
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
#---------------------------------------
assign.behav.seg=function(breakpt,dat){
  breakpt1=c(0,breakpt,Inf)
  n=length(breakpt1)
  dat$behav.seg<- NA
  
  for (i in 2:n){
    ind=which(breakpt1[i-1]<dat$time1 & dat$time1<breakpt1[i])
    dat[ind,"behav.seg"]=i-1
  }
  dat
}
#------------------------------------------------
df.to.list=function(dat) {  #only for id as col in dat
  id<- unique(dat$id)
  n=length(id)
  dat.list<- vector("list", n)
  names(dat.list)<- id
  
  for (i in 1:length(id)) {
    dat.list[[i]]<- dat[dat$id==id[i],]
  }
  dat.list
}
#------------------------------------------------
get.summary.stats_behav=function(dat){  #dat must have time.seg assigned; for all IDs
  
  #create list of input and to store output
  dat.list<- df.to.list(dat = dat)
  id<- unique(dat$id)
  n<- length(id)
  obs.list<- vector("list", n)
  names(obs.list)<- id
  
  #calculate # of obs in each bin (per move param) by behav.seg
  for (i in 1:length(dat.list)) {
    dat.ind=dat.list[[i]]
    ntseg=max(dat.ind$behav.seg)
    
    
    #SL
    SL<- matrix(0, ntseg, max(dat.ind$SL, na.rm = T))
    colnames(SL)<- paste0("SL.",1:max(dat.ind$SL, na.rm = T))
    for (j in 1:ntseg){
      tmp<- dat.ind %>% filter(behav.seg == j) %>% dplyr::select(SL) %>% table()
      SL[j,as.numeric(names(tmp))]<- tmp
    }
    
    #TA
    TA<- matrix(0, ntseg, max(dat.list[[i]]$TA, na.rm = T))
    colnames(TA)<- paste0("TA.",1:max(dat.list[[i]]$TA, na.rm = T))
    for (j in 1:ntseg){
      tmp<- dat.list[[i]] %>% filter(behav.seg == j) %>% dplyr::select(TA) %>% table()
      TA[j,as.numeric(names(tmp))]<- tmp
    }
    
    #TAA
    TAA<- matrix(0, ntseg, 1)
    colnames(TAA)<- 'TAA.1'
    for (j in 1:ntseg){
      tmp<- dat.list[[i]] %>% filter(behav.seg == j) %>% dplyr::select(TAA) %>% table()
      tmp1<- tmp %>% data.frame() %>% filter(names(tmp) == 1) %>% dplyr::select(Freq) %>%
                     as.numeric()
      TAA[j,]<- tmp1
    }
    
    
    id<- rep(unique(dat.ind$id), ntseg)
    behav.res<- cbind(id, TA, SL, TAA) %>% data.frame()
    obs.list[[i]]<- behav.res
  }
  #obs<- do.call(rbind.data.frame, obs.list)
  obs<- map_dfr(obs.list, `[`)
  obs
}