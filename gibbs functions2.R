get.summary.stats=function(breakpt,dat,max.SL,max.TA){
  breakpt1=c(0,breakpt,Inf)
  n=length(breakpt1)
  
  #get results for SL
  res.SL=matrix(0,n-1,max.SL)
  res.TA=matrix(0,n-1,max.TA)
  res.TAA=matrix(0,n-1,2)
  for (i in 2:n){
    ind=which(breakpt1[i-1]<dat$time1 & dat$time1<breakpt1[i])
    tmp=dat[ind,]
    
    #get SL results
    tmp1=table(tmp[,'SL'])
    ind=as.numeric(names(tmp1))
    res.SL[i-1,ind]=tmp1
    
    #get TA results
    tmp1=table(tmp[,'TA'])
    ind=as.numeric(names(tmp1))
    res.TA[i-1,ind]=tmp1
    
    #get TAA results
    tmp1=tmp[,'TAA']
    tmp1=if (length(tmp1[is.na(tmp1)])>0) {
      tmp1[!is.na(tmp1)]
      } else { tmp1
    }
    soma=sum(tmp1)
    res.TAA[i-1,]=c(soma,length(tmp1)-soma)
  }
  colnames(res.TAA)=c('n1','n0')
  list(res.TA=res.TA,res.SL=res.SL,res.TAA=res.TAA)
}
#---------------------------------------------
log.marg.likel=function(alpha,summary.stats,bern.a,bern.b,max.SL,max.TA){
  #get ratio for SL
  lnum=rowSums(lgamma(alpha+summary.stats$res.SL))
  lden=lgamma(max.SL*alpha+rowSums(summary.stats$res.SL))
  p2=sum(lnum)-sum(lden)
  p1=nrow(summary.stats$res.SL)*(lgamma(max.SL*alpha)-max.SL*lgamma(alpha))
  p.SL=p1+p2
  
  #get ratio for TA
  lnum=rowSums(lgamma(alpha+summary.stats$res.TA))
  lden=lgamma(max.TA*alpha+rowSums(summary.stats$res.TA))
  p2=sum(lnum)-sum(lden)
  p1=nrow(summary.stats$res.TA)*(lgamma(max.TA*alpha)-max.TA*lgamma(alpha))
  p.TA=p1+p2
  
  #get ratio for TAA
  lnum=lgamma(bern.a+summary.stats$res.TAA[,'n1'])+lgamma(bern.b+summary.stats$res.TAA[,'n0'])
  lden=lgamma(rowSums(summary.stats$res.TAA)+bern.a+bern.b)
  p2=sum(lnum)-sum(lden)
  p1=nrow(summary.stats$res.TAA)*(lgamma(bern.a+bern.b)-lgamma(bern.a)-lgamma(bern.b))
  p.TAA=p1+p2
  
  p.SL+p.TA+p.TAA
}
#---------------------------------------------
samp.move=function(breakpt,max.time,dat,alpha,bern.a,bern.b,max.SL,max.TA){
  breakpt.old=breakpt
  p=length(breakpt)
  rand1=runif(1)	
  p0=1
  new.brk=runif(1,min=0,max=max.time)
  
  if (p == 1) {
    #birth
    if (rand1 < 1/2){
      breakpt.new=sort(c(breakpt.old,new.brk))
      p0=2/3 #death prob 2 -> 1 is (1/3) and birth prob 1 -> 2 is 1/2. 
    }
    #swap
    if (rand1 > 1/2) breakpt.new=new.brk
  }
  if (p > 1) {
    #birth
    if (rand1 < 1/3) {
      breakpt.new=sort(c(breakpt.old,new.brk))
    }
    #death
    if (rand1 > 1/3 & rand1 < 2/3) {
      ind=sample(1:length(breakpt.old),size=1)
      breakpt.new=breakpt.old[-ind]
      if (p==2) p0=3/2 #birth prob from 1 -> 2 is 1/2 and death prob from 2 -> 1 is 1/3
    }
    #swap
    if (rand1 > 2/3) {
      ind=sample(1:length(breakpt.old),size=1)
      breakpt.new=sort(c(breakpt.old[-ind],new.brk))
    }
  }
  
  #get sufficient statistics
  stats.old=get.summary.stats(breakpt=breakpt.old,dat=dat,max.SL=max.SL,max.TA=max.TA)
  stats.new=get.summary.stats(breakpt=breakpt.new,dat=dat,max.SL=max.SL,max.TA=max.TA)
  
  #get marginal loglikel
  pold=log.marg.likel(alpha=alpha,summary.stats=stats.old,bern.a=bern.a,bern.b=bern.b,max.SL=max.SL,max.TA=max.TA)
  pnew=log.marg.likel(alpha=alpha,summary.stats=stats.new,bern.a=bern.a,bern.b=bern.b,max.SL=max.SL,max.TA=max.TA)+log(p0)
  prob=exp(pnew-pold)
  rand2=runif(1)
  
  if (rand2<prob) return(list(breakpt.new, pnew))
  return(list(breakpt.old, pold))
}
#---------------------------------------------
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
chng.SL=function(dat){ #Step Length Autocorr w/ ref to median
  dat$SLA<- NA
  
  for(i in 1:(nrow(dat)-1)) {
    tmp=dat$dist-median(dat$dist, na.rm = T)
    dat$SLA[i+1]<- ifelse(sign(tmp[i]) == sign(tmp[i+1]), 1, 0)
  }
  dat
}
#---------------------------------------------
assign.move.persist=function(dat) {
  dat$MP<- NA
  SL=dat$dist-median(dat$dist, na.rm = T) #"centered" SL data
  SL[which(is.na(SL))]=0
  dat[which(is.na(dat$rel.angle)),"rel.angle"]=0

  
  for (i in 1:(nrow(dat)-1)) {
    if (sign(SL[i]) == sign(SL[i+1]) & sign(dat$rel.angle[i]) == sign(dat$rel.angle[i+1])) {
      dat[i,"MP"]=1
    } else if (sign(SL[i]) != sign(SL[i+1]) & sign(dat$rel.angle[i]) == sign(dat$rel.angle[i+1])) {
      dat[i,"MP"]=2
    } else if (sign(SL[i]) == sign(SL[i+1]) & sign(dat$rel.angle[i]) != sign(dat$rel.angle[i+1])) {
      dat[i,"MP"]=3
    } else {
      dat[i,"MP"]=4
    }
  }
  dat
}
#---------------------------------------------
direct.persist=function(dat){  #creates 4 bins
  diff.angle=diff(dat$abs.angle)
  tmp<- matrix(NA, nrow(dat), 1)
  
  for(i in 1:length(diff.angle)) {
    tmp[i]=cos(diff.angle[i])
  }
  direct.persist.lims=seq(from=-1, to=1, by=0.5)
  dat$DP<- NA
  
  for(i in 1:length(direct.persist.lims)) {
    tmp1=which(tmp >= direct.persist.lims[i] & tmp < direct.persist.lims[i+1])
    dat[tmp1,"DP"]=i
  }
  dat
}
#---------------------------------------------
veloc.persist=function(dat) {  # create 8 bins
  dist=dat$dist*1000  #from km to m 
  t.step=dat$dt
  theta=dat$rel.angle
  
  Vp<- (dist/t.step)*cos(theta)
  veloc.persist.lims=seq(from=min(Vp, na.rm = T), to=max(Vp, na.rm = T), length.out = 9)
  dat$Vp<- NA
  
  for(i in 1:length(veloc.persist.lims)) {
    tmp=which(Vp >= veloc.persist.lims[i] & Vp < veloc.persist.lims[i+1])
    dat[tmp,"Vp"]=i
  }
  dat
}
#---------------------------------------------









get.summary.stats2=function(breakpt,dat,max.SL,max.TA,max.Vp){ #to change Vp to multinomial
  breakpt1=c(0,breakpt,Inf)
  n=length(breakpt1)
  
  #get results for SL
  res.SL=matrix(0,n-1,max.SL)
  res.TA=matrix(0,n-1,max.TA)
  res.Vp=matrix(0,n-1,max.Vp)
  for (i in 2:n){
    ind=which(breakpt1[i-1]<dat$time1 & dat$time1<breakpt1[i])
    tmp=dat[ind,]
    
    #get SL results
    tmp1=table(tmp[,'SL'])
    ind=as.numeric(names(tmp1))
    res.SL[i-1,ind]=tmp1
    
    #get TA results
    tmp1=table(tmp[,'TA'])
    ind=as.numeric(names(tmp1))
    res.TA[i-1,ind]=tmp1
    
    #get Vp results
    tmp1=table(tmp[,'Vp'])
    ind=as.numeric(names(tmp1))
    res.Vp[i-1,ind]=tmp1
  }
  list(res.TA=res.TA,res.SL=res.SL,res.Vp=res.Vp)
}
#---------------------------------------------
log.marg.likel2=function(alpha,summary.stats,max.SL,max.TA,max.Vp){ #to change Vp to multinomial
  #get ratio for SL
  lnum=rowSums(lgamma(alpha+summary.stats$res.SL))
  lden=lgamma(max.SL*alpha+rowSums(summary.stats$res.SL))
  p2=sum(lnum)-sum(lden)
  p1=nrow(summary.stats$res.SL)*(lgamma(max.SL*alpha)-max.SL*lgamma(alpha))
  p.SL=p1+p2
  
  #get ratio for TA
  lnum=rowSums(lgamma(alpha+summary.stats$res.TA))
  lden=lgamma(max.TA*alpha+rowSums(summary.stats$res.TA))
  p2=sum(lnum)-sum(lden)
  p1=nrow(summary.stats$res.TA)*(lgamma(max.TA*alpha)-max.TA*lgamma(alpha))
  p.TA=p1+p2
  
  #get ratio for Vp
  lnum=rowSums(lgamma(alpha+summary.stats$res.Vp))
  lden=lgamma(max.Vp*alpha+rowSums(summary.stats$res.Vp))
  p2=sum(lnum)-sum(lden)
  p1=nrow(summary.stats$res.Vp)*(lgamma(max.Vp*alpha)-max.Vp*lgamma(alpha))
  p.Vp=p1+p2
  
  p.SL+p.TA+p.Vp
}
#---------------------------------------------
samp.move2=function(breakpt,max.time,dat,alpha,max.SL,max.TA,max.Vp){ #to change Vp to multinomial
  breakpt.old=breakpt
  p=length(breakpt)
  rand1=runif(1)	
  p0=1
  new.brk=runif(1,min=0,max=max.time)
  
  if (p == 1) {
    #birth
    if (rand1 < 1/2){
      breakpt.new=sort(c(breakpt.old,new.brk))
      p0=2/3 #death prob 2 -> 1 is (1/3) and birth prob 1 -> 2 is 1/2. 
    }
    #swap
    if (rand1 > 1/2) breakpt.new=new.brk
  }
  if (p > 1) {
    #birth
    if (rand1 < 1/3) {
      breakpt.new=sort(c(breakpt.old,new.brk))
    }
    #death
    if (rand1 > 1/3 & rand1 < 2/3) {
      ind=sample(1:length(breakpt.old),size=1)
      breakpt.new=breakpt.old[-ind]
      if (p==2) p0=3/2 #birth prob from 1 -> 2 is 1/2 and death prob from 2 -> 1 is 1/3
    }
    #swap
    if (rand1 > 2/3) {
      ind=sample(1:length(breakpt.old),size=1)
      breakpt.new=sort(c(breakpt.old[-ind],new.brk))
    }
  }
  
  #get sufficient statistics
  stats.old=get.summary.stats2(breakpt=breakpt.old,dat=dat,max.SL=max.SL,max.TA=max.TA,max.Vp=max.Vp)
  stats.new=get.summary.stats2(breakpt=breakpt.new,dat=dat,max.SL=max.SL,max.TA=max.TA,max.Vp=max.Vp)
  
  #get marginal loglikel
  pold=log.marg.likel2(alpha=alpha,summary.stats=stats.old,max.SL=max.SL,max.TA=max.TA,max.Vp=max.Vp)
  pnew=log.marg.likel2(alpha=alpha,summary.stats=stats.new,max.SL=max.SL,max.TA=max.TA,max.Vp=max.Vp)+log(p0)
  prob=exp(pnew-pold)
  rand2=runif(1)
  
  if (rand2<prob) return(list(breakpt.new, pnew))
  return(list(breakpt.old, pold))
}
#---------------------------------------------










get.summary.stats3=function(breakpt,dat,max.SL,max.TA){  #include Step Length Autocorr
  breakpt1=c(0,breakpt,Inf)
  n=length(breakpt1)
  
  #get results for SL
  res.SL=matrix(0,n-1,max.SL)
  res.TA=matrix(0,n-1,max.TA)
  res.TAA=matrix(0,n-1,2)
  res.SLA=matrix(0,n-1,2)
  for (i in 2:n){
    ind=which(breakpt1[i-1]<dat$time1 & dat$time1<breakpt1[i])
    tmp=dat[ind,]
    
    #get SL results
    tmp1=table(tmp[,'SL'])
    ind=as.numeric(names(tmp1))
    res.SL[i-1,ind]=tmp1
    
    #get TA results
    tmp1=table(tmp[,'TA'])
    ind=as.numeric(names(tmp1))
    res.TA[i-1,ind]=tmp1
    
    #get TAA results
    tmp1=tmp[,'TAA']
    tmp1=if (length(tmp1[is.na(tmp1)])>0) {
      tmp1[!is.na(tmp1)]
    } else { tmp1
    }
    soma=sum(tmp1)
    res.TAA[i-1,]=c(soma,length(tmp1)-soma)
  
    #get SLA results
    tmp1=tmp[,'SLA']
    tmp1=if (length(tmp1[is.na(tmp1)])>0) {
      tmp1[!is.na(tmp1)]
    } else { tmp1
  }
  soma=sum(tmp1)
  res.SLA[i-1,]=c(soma,length(tmp1)-soma)
  }
colnames(res.TAA)=c('n1','n0')
colnames(res.SLA)=c('n1','n0')

  list(res.TA=res.TA,res.SL=res.SL,res.TAA=res.TAA,res.SLA=res.SLA)
}
#---------------------------------------------
log.marg.likel3=function(alpha,summary.stats,bern.a,bern.b,max.SL,max.TA){
  #get ratio for SL
  lnum=rowSums(lgamma(alpha+summary.stats$res.SL))
  lden=lgamma(max.SL*alpha+rowSums(summary.stats$res.SL))
  p2=sum(lnum)-sum(lden)
  p1=nrow(summary.stats$res.SL)*(lgamma(max.SL*alpha)-max.SL*lgamma(alpha))
  p.SL=p1+p2
  
  #get ratio for TA
  lnum=rowSums(lgamma(alpha+summary.stats$res.TA))
  lden=lgamma(max.TA*alpha+rowSums(summary.stats$res.TA))
  p2=sum(lnum)-sum(lden)
  p1=nrow(summary.stats$res.TA)*(lgamma(max.TA*alpha)-max.TA*lgamma(alpha))
  p.TA=p1+p2
  
  #get ratio for TAA
  lnum=lgamma(bern.a+summary.stats$res.TAA[,'n1'])+lgamma(bern.b+summary.stats$res.TAA[,'n0'])
  lden=lgamma(rowSums(summary.stats$res.TAA)+bern.a+bern.b)
  p2=sum(lnum)-sum(lden)
  p1=nrow(summary.stats$res.TAA)*(lgamma(bern.a+bern.b)-lgamma(bern.a)-lgamma(bern.b))
  p.TAA=p1+p2
  
  #get ratio for SLA
  lnum=lgamma(bern.a+summary.stats$res.SLA[,'n1'])+lgamma(bern.b+summary.stats$res.SLA[,'n0'])
  lden=lgamma(rowSums(summary.stats$res.SLA)+bern.a+bern.b)
  p2=sum(lnum)-sum(lden)
  p1=nrow(summary.stats$res.SLA)*(lgamma(bern.a+bern.b)-lgamma(bern.a)-lgamma(bern.b))
  p.SLA=p1+p2
  
  p.SL+p.TA+p.TAA+p.SLA
}
#---------------------------------------------
samp.move3=function(breakpt,max.time,dat,alpha,bern.a,bern.b,max.SL,max.TA){
  breakpt.old=breakpt
  p=length(breakpt)
  rand1=runif(1)	
  p0=1
  new.brk=runif(1,min=0,max=max.time)
  
  if (p == 1) {
    #birth
    if (rand1 < 1/2){
      breakpt.new=sort(c(breakpt.old,new.brk))
      p0=2/3 #death prob 2 -> 1 is (1/3) and birth prob 1 -> 2 is 1/2. 
    }
    #swap
    if (rand1 > 1/2) breakpt.new=new.brk
  }
  if (p > 1) {
    #birth
    if (rand1 < 1/3) {
      breakpt.new=sort(c(breakpt.old,new.brk))
    }
    #death
    if (rand1 > 1/3 & rand1 < 2/3) {
      ind=sample(1:length(breakpt.old),size=1)
      breakpt.new=breakpt.old[-ind]
      if (p==2) p0=3/2 #birth prob from 1 -> 2 is 1/2 and death prob from 2 -> 1 is 1/3
    }
    #swap
    if (rand1 > 2/3) {
      ind=sample(1:length(breakpt.old),size=1)
      breakpt.new=sort(c(breakpt.old[-ind],new.brk))
    }
  }
  
  #get sufficient statistics
  stats.old=get.summary.stats3(breakpt=breakpt.old,dat=dat,max.SL=max.SL,max.TA=max.TA)
  stats.new=get.summary.stats3(breakpt=breakpt.new,dat=dat,max.SL=max.SL,max.TA=max.TA)
  
  #get marginal loglikel
  pold=log.marg.likel3(alpha=alpha,summary.stats=stats.old,bern.a=bern.a,bern.b=bern.b,max.SL=max.SL,max.TA=max.TA)
  pnew=log.marg.likel3(alpha=alpha,summary.stats=stats.new,bern.a=bern.a,bern.b=bern.b,max.SL=max.SL,max.TA=max.TA)+log(p0)
  prob=exp(pnew-pold)
  rand2=runif(1)
  
  if (rand2<prob) return(list(breakpt.new, pnew))
  return(list(breakpt.old, pold))
}