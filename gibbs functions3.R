get.summary.stats=function(breakpt,dat,max.Vp,max.Vt){ #to change Vp to multinomial
  breakpt1=c(0,breakpt,Inf)
  n=length(breakpt1)
  
  #get results for SL
  res.Vp=matrix(0,n-1,max.Vp)
  res.Vt=matrix(0,n-1,max.Vt)
  for (i in 2:n){
    ind=which(breakpt1[i-1]<dat$time1 & dat$time1<breakpt1[i])
    tmp=dat[ind,]
    
    #get Vp results
    tmp1=table(tmp[,'Vp'])
    ind=as.numeric(names(tmp1))
    res.Vp[i-1,ind]=tmp1
    
    #get Vt results
    tmp1=table(tmp[,'Vt'])
    ind=as.numeric(names(tmp1))
    res.Vt[i-1,ind]=tmp1
  }
  list(res.Vp=res.Vp,res.Vt=res.Vt)
}
#---------------------------------------------
log.marg.likel=function(alpha,summary.stats,max.Vp,max.Vt){ #to change Vp to multinomial
  #get ratio for Vp
  lnum=rowSums(lgamma(alpha+summary.stats$res.Vp))
  lden=lgamma(max.Vp*alpha+rowSums(summary.stats$res.Vp))
  p2=sum(lnum)-sum(lden)
  p1=nrow(summary.stats$res.Vp)*(lgamma(max.Vp*alpha)-max.Vp*lgamma(alpha))
  p.Vp=p1+p2
  
  #get ratio for Vt
  lnum=rowSums(lgamma(alpha+summary.stats$res.Vt))
  lden=lgamma(max.Vt*alpha+rowSums(summary.stats$res.Vt))
  p2=sum(lnum)-sum(lden)
  p1=nrow(summary.stats$res.Vt)*(lgamma(max.Vt*alpha)-max.Vt*lgamma(alpha))
  p.Vt=p1+p2
  
  p.Vp+p.Vt
}
#---------------------------------------------
samp.move=function(breakpt,max.time,dat,alpha,max.Vp,max.Vt){ #to change Vp to multinomial
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
  stats.old=get.summary.stats(breakpt=breakpt.old,dat=dat,max.Vp=max.Vp,max.Vt=max.Vt)
  stats.new=get.summary.stats(breakpt=breakpt.new,dat=dat,max.Vp=max.Vp,max.Vt=max.Vt)
  
  #get marginal loglikel
  pold=log.marg.likel(alpha=alpha,summary.stats=stats.old,max.Vp=max.Vp,max.Vt=max.Vt)
  pnew=log.marg.likel(alpha=alpha,summary.stats=stats.new,max.Vp=max.Vp,max.Vt=max.Vt)+log(p0)
  prob=exp(pnew-pold)
  rand2=runif(1)
  
  if (rand2<prob) return(list(breakpt.new, pnew))
  return(list(breakpt.old, pold))
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