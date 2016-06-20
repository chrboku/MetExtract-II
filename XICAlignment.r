suppressWarnings(suppressMessages(library(ptw)))
suppressWarnings(suppressMessages(library(signal)))
#source("./MassSpecWaveletIdentification.r")

filtSmooth=5

getAlignXICs<-function(eics, ncol, align=TRUE, npoly=2, ra=c(1:length(eics[0]))){
  mdat<-matrix(eics, ncol=ncol)
  len<-dim(mdat)[1]
  #ra=c(1100:1300)
  ra=1:len
  plotRes<-FALSE

  
  eicCum=rep(0, len)
  for(i in 1:len){
    eicCum[i]=0
    for(j in 1:1){
      eicCum[i]=eicCum[i]+mdat[i,j]
    }
  }

  if(plotRes){
    plot(eicCum, type="l")
    for(i in 1:(dim(mdat)[2])){
      lines(mdat[,i], col="green")
    }
    lines(eicCum, col="black")
    eicCum=filtfilt(rep(1, filtSmooth)/filtSmooth,1,eicCum)
    lines(eicCum, col="red")
  }


  if(plotRes){
    plot(eicCum[ra], type="l")
    k=1
    for(i in 1:ncol){
      if(i==1){
      x11()
        plot(mdat[,i][ra], col="grey", type="l", main="not aligned")
      }else{
        lines(mdat[,i][ra], col="grey")
      }
    }
  }
  
  ret=matrix(eics, ncol=ncol)
  if(align){
    for(i in 1:ncol){
      poly<-rep(0, npoly)
      poly[2]=1
      w<-ptw(eicCum, mdat[,i], init.coef=poly)
      ret[,i]=w$warped.sample
      ret[is.na(ret[,i]),i]<-0
      for(j in 1:len){
	    eicCum[j]=eicCum[j]+ret[j,i]
      }
      if(plotRes){
        if(i==1){
          x11()
          plot(w$warped.sample[ra], col="grey", type="l", main="aligned")
        }else{
          lines(w$warped.sample[ra], col="grey")
        }
      }
    }
  }
  return(ret)
}

alignPeaks<-function(eics, peaks, ncol, align=TRUE, npoly=2, ra=c(1000:1200), maxGroupDist=22){
  mdat<-matrix(eics, ncol=ncol)
  len<-dim(mdat)[1]
  mpeaks<-matrix(peaks, ncol=ncol)
  #ra=c(350:800)
  
  for(i in 1:ncol){
    #mdat[,i]=filtfilt(rep(1, filtSmooth)/filtSmooth,1,mdat[,i])  
  }
  
  eicCum=rep(0, len)
  for(i in 1:len){
    eicCum[i]=0
    for(j in 1:1){ 
      eicCum[i]=eicCum[i]+mdat[i,j]
    }
  }
  
  
  #eicCum=filtfilt(rep(1, filtSmooth)/filtSmooth,1,eicCum)
  #plot(eicCum[ra], type="l")
  allpeaks=rep(0, length(peaks))
  k=1
  for(i in 1:ncol){
    #if(i==1){
    #  plot( mdat[,i][ra], col="grey", type="l", main="not aligned")
    #}else{
    #  lines( mdat[,i][ra], col="grey")
    #}
    for(peak in mpeaks[,i]){
      #cat("adding peak", peak, peak[1], "\n")
      #abline(v=peak[1]-ra[1], col="orange")
    }
  }
  #x11()
  #eicd=rep(0, len)
  if(align){
    for(i in 1:ncol){
      poly<-rep(0, npoly)
      poly[2]=1
      w<-ptw(eicCum, mdat[,i], init.coef=poly)
      #for(j in 1:len){
      #  eicd[j]=eicd[j]+w$warped.sample[j]
      #}
      x=rep(0, len)
      x=w$warped.sample
      x[is.na(x)]<-0
      for(j in 1:min(len, length(x))){
        eicCum[j]=eicCum[j]+x[j]
      }
      #if(i==1){
      #  plot(w$warped.sample[ra], col="grey", type="l", main="aligned")
      #}else{
      #  lines(w$warped.sample[ra], col="grey")
      #}
    
      for(peak in mpeaks[,i]){
        #abline(v=which.min(abs(w$warp.fun-peak[1]))-ra[1], col="orange")
        #cat(peak, " translated to ", which.min(abs(w$warp.fun-peak[1])), "\n")
        allpeaks[k]=which.min(abs(w$warp.fun-peak[1]))
        #allpeaks[k]=peak[1]
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
  #x11()
  #plot(eicd, type="l")
  #eicd[is.na(eicd)]<-0 
  #peaks=getMajorPeaks(eicd)
  #for(peak in peaks){
  #  cat(peak, "\n")
  #  abline(v=peak[1])
  #}
  
  #x11()
  hcl=hclust(dist(allpeaks))
  #plot(hcl)
  #abline(h=maxGroupDist, col="red")
  #Sys.sleep(3)
  cuted<-cutree(hcl, h=maxGroupDist)
  ret=c()

  #for(j in unique(cuted)){
  #  centers=allpeaks[cuted==j]
  #
  #  n=length(centers)
  #  print(paste(c("group", j, "Count", n, ":", paste(centers, collapse=",")), collapse=" "))
  #}
  #print("---------")


  for(i in 1:length(cuted)){
    ret=c(ret, allpeaks[i], cuted[i])
  }
  return((array(ret, c(2,length(ret)/2))))
}
