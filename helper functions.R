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
veloc.persist=function(dat) {  #persistence velocity via Gurarie et al 2009; create 8 bins
  dist=dat$dist*1000  #from km to m 
  t.step=dat$dt
  theta=dat$abs.angle
  
  Vp<- (dist/t.step)*cos(theta)
  veloc.persist.lims=seq(from=min(Vp, na.rm = T), to=max(Vp, na.rm = T), length.out = 9)
  dat$Vp<- NA
  
  for(i in 1:length(veloc.persist.lims)) {
    tmp=which(Vp >= veloc.persist.lims[i] & Vp < veloc.persist.lims[i+1])
    dat[tmp,"Vp"]=i
  }
  tmp=which(Vp == veloc.persist.lims[9])
  dat[tmp,"Vp"]=8
  dat
}
#---------------------------------------------
veloc.turn=function(dat) {  #turning velocity via Gurarie et al 2009; create 8 bins
  dist=dat$dist*1000  #from km to m 
  t.step=dat$dt
  theta=dat$abs.angle
  
  Vt<- (dist/t.step)*sin(theta)
  veloc.turn.lims=seq(from=min(Vt, na.rm = T), to=max(Vt, na.rm = T), length.out = 9)
  dat$Vt<- NA
  
  for(i in 1:length(veloc.turn.lims)) {
    tmp=which(Vt >= veloc.turn.lims[i] & Vt < veloc.turn.lims[i+1])
    dat[tmp,"Vt"]=i
  }
  tmp=which(Vt == veloc.turn.lims[9])
  dat[tmp,"Vt"]=8
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
#discretize persistence (Vp) and turning velocity (Vt) by ID
#only subset obs w/ 1 hr (3600 s) time step; removal of 5815 total obs
behav.prep_veloc=function(dat,tstep) {
  id<- unique(dat$id)
  n=length(id)
  behav.list<- vector("list", n)
  names(behav.list)<- id
  
  for (i in 1:length(id)) {
    behav.list[[i]]<- dat[dat$id==id[i],] %>% mutate(obs = 1:nrow(.)) %>% filter(dt==tstep) %>%
      mutate(time1 = 1:nrow(.)) %>% veloc.persist() %>% veloc.turn
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
#---------------------------------------------
behav.seg.image_veloc=function(dat) {  #Transform single var vectors into pres/abs matrices for heatmap
  Vp<- matrix(0, nrow(dat), length(unique(dat$Vp)))
  for (i in 1:nrow(dat)){
    Vp[i,dat$Vp[i]]=1
  }
  
  Vt<- matrix(0, nrow(dat), length(unique(dat$Vt)))
  for (i in 1:nrow(dat)){
    Vt[i,dat$Vt[i]]=1
  }
  
  behav.list<- list(Vp=Vp,Vt=Vt)
  
  behav.list
}