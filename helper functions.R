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
