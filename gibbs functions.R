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
  
  if (rand2<prob) return(breakpt.new)
  return(breakpt.old)
}