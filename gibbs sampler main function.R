time.segm.behavior=function(dat,ngibbs,breakpt=NULL){
  #useful stuff
  max.time=nrow(dat)
  ncat.data=apply(dat,2,max)
  
  #starting values
  if (is.null(breakpt)) breakpt=floor(max.time/2)
  
  for (i in 1:ngibbs){
    print(i)
    breakpt=samp.move(breakpt=breakpt,max.time=max.time,dat=dat,
                      alpha=alpha,ncat.data=ncat.data,ndata.types=ndata.types)
  }
  breakpt
}
