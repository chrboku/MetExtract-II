suppressWarnings(suppressMessages(library(ptw)))
suppressWarnings(suppressMessages(library(signal)))
#source("./MassSpecWaveletIdentification.r")

filtSmooth=5

getAlignXICs<-function(eics, ncol, align=TRUE, npoly=2, ra=c(1:length(eics[0]))){
  mdat<-matrix(eics, ncol=ncol)
  len<-dim(mdat)[1]

  eicCum=mdat[,1]

  ret=matrix(eics, ncol=ncol)

  if(align){
    for(i in 1:ncol){
      poly<-rep(0, npoly)
      poly[2]=1

      w<-ptw(eicCum, mdat[,i], init.coef=poly)

      #warp.samp <- w$warped.sample
      #warp.samp[is.na(warp.samp)] <- 0
      #w <- ptw(eicCum, warp.samp, init.coef=c(0,1,0,0,0,0))

      ret[,i]=w$warped.sample
      ret[is.na(ret[,i]),i]<-0
    }
  }

  return(ret)
}

alignPeaks<-function(eics, refTimes, peaks, ncol, scanDuration=1, align=TRUE, npoly=2, ra=c(1000:1200), maxGroupDist=22){

  mpeaks<-matrix(peaks, ncol=ncol)
  allpeaks=rep(0, length(peaks))
  k=1

  if(align){
      mdat<-matrix(eics, ncol=ncol)
      len<-dim(mdat)[1]

      for(i in 1:ncol){
        #mdat[,i]=filtfilt(rep(1, filtSmooth)/filtSmooth,1,mdat[,i])
      }

      eicCum=mdat[,1]
    for(i in 1:ncol){
      poly<-rep(0, npoly)
      poly[2]=1

      w<-ptw(eicCum, mdat[,i], init.coef=poly)

      #warp.samp <- w$warped.sample
      #warp.samp[is.na(warp.samp)] <- 0
      #w <- ptw(eicCum, warp.samp, init.coef=c(0,1,0,0,0,0))

      #x=rep(0, len)
      #x=w$warped.sample
      #x[is.na(x)]<-0

      for(peak in mpeaks[,i]){
        allpeaks[k]=refTimes[which.min(abs(w$warp.fun[1,]-peak[1]))]
        k=k+1
      }
    }
  }else{
    for(i in 1:ncol){
      for(peak in mpeaks[,i]){
      	allpeaks[k]=peak[1]
	    k=k+1
      }
    }
  }

  hcl=hclust(dist(allpeaks, method = "euclidean"), method="average")
  cuted<-cutree(hcl, h=maxGroupDist)

  ret=c()
  for(i in 1:length(cuted)){
    ret=c(ret, allpeaks[i], cuted[i])
  }

  return((array(ret, c(2,length(ret)/2))))
}
