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
  
  dat$SL<- NA
  
  for(i in 1:length(dist.bin.lims)) {
    tmp=which(dat$dist >= dist.bin.lims[i] & dat$dist < dist.bin.lims[i+1])
    dat[tmp,"SL"]=i
  }
  tmp=which(dat$dist == max.dist)
  dat[tmp,"SL"]=6
  
  dat
}
#---------------------------------------------
chng.rel_angle.sign=function(dat){  #turning angle autocorr
  dat$TAA<- NA
  
  for(i in 1:(nrow(dat)-1)) {
    dat$TAA[i+1]<- ifelse(sign(dat$rel.angle[i]) == sign(dat$rel.angle[i+1]), 1, 0)
  }
  dat
}
#---------------------------------------------
round_track_time=function(dat, int, tol) {  #replacement for sett0() when wanting to only round some of the times
  dat<- df.to.list(dat)
  
  for (i in 1:length(dat)) {
    tmp=matrix(NA,nrow(dat[[i]]),2)
    
    for (j in 1:nrow(dat[[i]])) {
      if (is.na(dat[[i]]$dt[j])) {
        tmp[j, 1:2]<- NA
      } else if (dat[[i]]$dt[j] >= (int - tol) & dat[[i]]$dt[j] <= (int + tol)) {
        tmp[j, 1:2]<- c(int, as.numeric(round(dat[[i]]$date[j], units = "hours")))
      } else {
        tmp[j, 1:2]<- c(dat[[i]]$dt[j], dat[[i]]$date[j])
      }
    }
    dat[[i]]$dt<- tmp[,1]
    dat[[i]]$date<- tmp[,2] %>% as.POSIXct(origin = '1970-01-01', tz = "UTC")
  }
  dat<- map_dfr(dat, `[`)
  dat
}
#---------------------------------------------
#discretize step length and turning angle by ID
#only subset obs w/ 1 hr (3600 s) time step; results in loss of irregular obs
behav.prep=function(dat,tstep) {
  id<- unique(dat$id)
  
  #set bin limits for SL
  max.dist=max(dat[dat$dt == 3600,]$dist, na.rm = T) #using value from entire dataset, not specific time segment
  upper90.thresh=as.numeric(quantile(dat$dist, 0.90, na.rm=T)) #using value from entire dataset
  
  dist.bin.lims=seq(from=0, to=upper90.thresh, length.out = 6)
  dist.bin.lims=c(dist.bin.lims, max.dist)
  
  cond<- matrix(0, length(dat.list), 1)
  for (i in 1:length(dat.list)) {  #don't include IDs w < 1 obs of dt == 3600
    cond[i,]<- if(nrow(dat.list[[i]][dat.list[[i]]$dt == 3600,]) > 1) {
      1
    } else {
      0
    }
  }
  tmp<- cbind(id,cond)
  ind<- which(tmp[,2] == 1)
  
  n=length(ind)
  behav.list<- vector("list", n)
  names(behav.list)<- id[ind]
  
  for (i in 1:n) {
    behav.list[[i]]<- dat[dat$id==id[ind[i]],] %>% mutate(obs = 1:nrow(.)) %>% filter(dt==tstep) %>%
                  mutate(time1 = 1:nrow(.)) %>% assign.rel_angle.bin() %>% assign.dist.bin()
  }
  behav.list
}
#---------------------------------------------
behav.seg.image=function(dat) {  #Transform single var vectors into pres/abs matrices for heatmap
  SL<- matrix(0, nrow(dat), max(length(unique(dat$SL)), max(dat$SL), na.rm = T))
  for (i in 1:nrow(dat)){
    SL[i,dat$SL[i]]=1
  }
  
  TA<- matrix(0, nrow(dat), max(length(unique(dat[!is.na(dat$TA),]$TA)), max(dat$TA), na.rm = T))
  for (i in 1:nrow(dat)){
    TA[i,dat$TA[i]]=1
  }
  
  behav.list<- list(SL=SL,TA=TA)
  behav.list
}
#---------------------------------------
assign.time.seg=function(dat){  #add tseg assignment to each obs
  tmp=which(unique(dat$id) == brkpts$id)
  breakpt<- brkpts[tmp,-1] %>% discard(is.na) %>% as.numeric(.[1,])
  breakpt1=c(0,breakpt,Inf)
  tmp1<- which(diff(breakpt1) < 1)  #check for impossible time units
  breakpt1[tmp1+1]<- breakpt1[(tmp1+1)] + 1  #fix impossible time units
  n=length(breakpt1)
  res=matrix(NA,nrow(dat),1)
  for (i in 2:n){
    ind=which(breakpt1[i-1]<=dat$time1 & dat$time1<breakpt1[i])
    res[ind,]=i-1
  }
  dat$tseg<- as.vector(res)
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
traceplot=function(data, type, identity) {  #create traceplots for nbrks or LML for all IDs
  for (i in 1:length(identity)) {
    par(ask=TRUE)
    plot(x=1:ngibbs, y=data[i,-1], type = "l", xlab = "Iteration",
         ylab = ifelse(type == "nbrks", "# of Breakpoints",
                       ifelse(type == "LML","Log Marginal Likelihood",
                              stop("Need to select one of 'nbrks' or 'LML' for plotting"))),
         main = paste("ID",data[i,1]))
  }
  on.exit(par(ask = FALSE))
}
#---------------------------------------------
getMAP=function(dat,nburn) {  #select MAP value that is beyond burn-in phase
  if (which.max(dat[-1]) < nburn) {
    MAP<- dat[-1] %>% order(decreasing = T) %>% subset(. > nburn) %>% first()
  } else {
    MAP<- which.max(dat[-1])
  }
  return(MAP)
}
#---------------------------------------------
getBreakpts=function(dat,MAP,brk.cols) {  #extract breakpoints of MAP per ID
  tmp<- matrix(NA, length(dat), brk.cols)
  
  for(i in 1:length(dat)) {
    ind<- MAP[i]
    tmp[i,1:length(dat[[i]][[ind]])]<- round(dat[[i]][[ind]], 2)
  }
  
  tmp<- cbind(id = identity, tmp) %>% data.frame()
  names(tmp)<- c('id', paste0("Brk_",1:brk.cols))
  tmp<- mutate_all(tmp, function(x) as.numeric(as.character(x)))
  tmp
}
#------------------------------------------------
plot.heatmap.behav=function(data, brkpts, dat.res) {
  
  behav.heat<- behav.seg.image(data)
  
  SL<- data.frame(behav.heat$SL)
  names(SL)<- 1:6
  SL<- SL %>% gather(key, value) %>% mutate(time=rep(data$time1, times=6),
                                            behav=rep("SL", nrow(data)*6))
  
  TA<- data.frame(behav.heat$TA)
  names(TA)<- 1:8
  TA<- TA %>% gather(key, value) %>% mutate(time=rep(data$time1, times=8),
                                            behav=rep("TA", nrow(data)*8))
  
  behav.heat_long<- rbind(SL,TA)
  behav.heat_long$value<- factor(behav.heat_long$value)
  levels(behav.heat_long$value)<- c("Unoccupied","Occupied")
  
  ind=which(unique(data$id) == brkpts$id)
  breakpt<- brkpts[ind,-1] %>% discard(is.na) %>% t() %>% data.frame()
  names(breakpt)<- "breaks"
  
  print(
  ggplot(behav.heat_long, aes(x=time, y=key, fill=value)) +
    geom_tile() +
    facet_wrap(~behav, scales = 'free', nrow = 2) +
    scale_fill_viridis_d('') +
    scale_y_discrete(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    geom_vline(data = breakpt, aes(xintercept = breaks), color = viridis(n=9)[7],
               size = 0.4) +
    labs(x = "Observations", y = "Bin", title = paste("ID", unique(data$id))) +
    theme_bw() +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12),
          strip.text = element_text(size = 12, face = 'bold'),
          title = element_text(size = 20))
  )
}
#------------------------------------------------
heatmap=function(data, brkpts, dat.res, type) {  #type can either be 'loc' or 'behav'
  
  if (type == "loc") {
    par(ask = TRUE)
    map(data, ~plot.heatmap.loc(., brkpts = brkpts, dat.res = dat.res))
    par(ask = FALSE)
  } else if (type == "behav") {
    par(ask = TRUE)
    map(data, ~plot.heatmap.behav(., brkpts = brkpts, dat.res = dat.res))
    par(ask = FALSE)
  } else {
    warning("Need to select type as either 'loc' or 'behav'")
  }
}
