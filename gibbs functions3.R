get.summary.stats=function(breakpt,dat,max.SL,max.TA,max.Vp){ #to change Vp to multinomial
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
log.marg.likel=function(alpha,summary.stats,max.SL,max.TA,max.Vp){ #to change Vp to multinomial
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
samp.move=function(breakpt,max.time,dat,alpha,max.SL,max.TA,max.Vp){ #to change Vp to multinomial
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
  stats.old=get.summary.stats(breakpt=breakpt.old,dat=dat,max.SL=max.SL,max.TA=max.TA,max.Vp=max.Vp)
  stats.new=get.summary.stats(breakpt=breakpt.new,dat=dat,max.SL=max.SL,max.TA=max.TA,max.Vp=max.Vp)
  
  #get marginal loglikel
  pold=log.marg.likel(alpha=alpha,summary.stats=stats.old,max.SL=max.SL,max.TA=max.TA,max.Vp=max.Vp)
  pnew=log.marg.likel(alpha=alpha,summary.stats=stats.new,max.SL=max.SL,max.TA=max.TA,max.Vp=max.Vp)+log(p0)
  prob=exp(pnew-pold)
  rand2=runif(1)
  
  if (rand2<prob) return(list(breakpt.new, pnew))
  return(list(breakpt.old, pold))
}
#---------------------------------------------
veloc.persist=function(dat) {  #persistence velocity via Gurarie et al 2009; create 8 bins
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
