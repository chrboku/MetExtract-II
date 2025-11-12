if(FALSE){
    source("https://bioconductor.org/biocLite.R")
    biocLite("MassSpecWavelet")
}

cat(paste0("\n\nLoading packages from", paste0(.libPaths(), collapse = ", "), "\n\n"))

suppressWarnings(suppressMessages(library(waveslim)))
suppressWarnings(suppressMessages(library(MassSpecWavelet)))
suppressWarnings(suppressMessages(library(signal)))
suppressWarnings(suppressMessages(library(baseline)))
#suppressWarnings(suppressMessages(library(rgl)))

pack<-installed.packages()
#cat("MassSpecWavelet version", pack[which(pack[,1]=="MassSpecWavelet"),3], "\n")

tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
   type <- match.arg(type)
   assign(".type", type, envir=baseenv())
   if(gcFirst) gc(FALSE)
   tic <- proc.time()[type]
   assign(".tic", tic, envir=baseenv())
   invisible(tic)
}

toc <- function()
{
   type <- get(".type", envir=baseenv())
   toc <- proc.time()[type]
   tic <- get(".tic", envir=baseenv())
   print(toc - tic)
   invisible(toc)
}

getMajorPeaksd<-function(eic, scales=c(5,33), snrTh=1, eicSmoothing="None"){

#	if(eicSmoothing=="None"){
#		#do not perform any eic smoothing
#	}else if(eicSmoothing=="Moving Average (3)"){
#		eic=filtfilt(rep(1, 3)/3,1,eic)
#	}else if(eicSmoothing=="Moving Average (5)"){
#		eic=filtfilt(rep(1, 5)/5,1,eic)
#	}else if(eicSmoothing=="Savitzky Golay Filter (3)"){
#		eic=sgolayfilt(eic, p=3)
#	}else if(eicSmoothing=="Savitzky Golay Filter (5)"){
#		eic=sgolayfilt(eic, p=5)
#	}else{
#		print("Unsupported smoothing")
#	}

    #I don't know why, but the scale factor 1 is required
    if(scales[1]==1){
        scales <- seq(scales[1], scales[2], 1)
    }else{
        scales <- c(1, seq(scales[1], scales[2], 1))
    }

    plot<-T

    majorPeakInfo<-1
    peakIndex<-1

    if(T){
        peakInfo<-peakDetectionCWT(eic, scales=scales, SNR.Th=snrTh, peakScaleRange = 2)

        majorPeakInfo = peakInfo$majorPeakInfo
        peakIndex <- majorPeakInfo$peakIndex

    }else{
        #cat("using own method\n")
        wCoefs <- cwt(eic, scales = scales, wavelet = "mexh")
        wCoefs <- cbind(as.vector(eic), wCoefs)
        #save(wCoefs, file="c:/temp/wavelet.txt")
        colnames(wCoefs) <- c(0, scales)
        localMax <- getLocalMaximumCWT(wCoefs)

        plotRange <- c(length(eic)*2/16., length(eic)*4/16.)
        plotRange <- c(1, length(eic))
        #cat("damnit", plotRange[1], " ", plotRange[2])
        #cat(plotRange[1]:plotRange[2], "\n", c(scales,max(scales)+1), "\n")
        #cat("Wavelet Max", max(wCoefs[plotRange[1]:plotRange[2],]), " Wavelet Min", min(wCoefs[plotRange[1]:plotRange[2],]), "\n")
        #plot3d(plotRange[1]:plotRange[2], c(scales,max(scales)+1), wCoefs[plotRange[1]:plotRange[2],])
        if(plot){
            xTickInterval <- 10
            image(plotRange[1]:plotRange[2], c(scales,max(scales)+1), wCoefs[plotRange[1]:plotRange[2],
                  ], col = terrain.colors(256), axes = FALSE, xlab = "m/z index",
                  ylab = "CWT coefficient scale", main = "CWT coefficients")
            axis(1, at = seq(plotRange[1], plotRange[2], by = xTickInterval))
            axis(2, at = c(1, seq(10, 64, by = 10)))
            box()
            image(plotRange[1]:plotRange[2], c(scales,max(scales)+1), localMax[plotRange[1]:plotRange[2],
                  ], col = terrain.colors(256), axes = FALSE, xlab = "m/z index",
                  ylab = "CWT coefficient scale", main = "CWT coefficients")

        }
        ridgeList <- getRidge(localMax)

        nearbyPeak<-TRUE
        majorPeakInfo<-identifyMajorPeaks(eic, ridgeList, wCoefs, SNR.Th=snrTh, nearbyPeak=nearbyPeak, nearbyWinSize=300)
        majorPeakInfo<-tuneInPeakInfo(eic, majorPeakInfo)
        peakIndex<-majorPeakInfo$peakIndex
    }

    ret=c()

    for(p in peakIndex){
      ind=paste("1_", toString(p), sep="")
      if(!is.na(majorPeakInfo$peakValue[ind])){
        ret=c(ret, majorPeakInfo$peakIndex[ind], majorPeakInfo$peakScale[ind], majorPeakInfo$peakSNR[ind], majorPeakInfo$peakValue[ind])
      }
    }

    if(length(peakIndex)==0){
        return(c(-1))
    }

    return(array(ret, c(4,length(ret)/4)))
}








getMajorPeaks<-function(eic, times, scales=c(5,33), snrTh=1, eicSmoothing="None", minCorr=0.85, maxSquaredError=20.5, testForGaussPeak=FALSE, removeBaseline=FALSE){

    if(removeBaseline){
      bl=baseline(t(eic), method='medianWindow', hwm=40)
      eic=vapply(bl@corrected, function(x){max(0, x)}, c(1))
    }

    x<-try(getMajorPeaksd(eic, scales, snrTh, eicSmoothing), silent=TRUE)

    if(testForGaussPeak){
      ret=c()

      peakInd=1
      while(peakInd<length(x)){
        peakCenter=x[peakInd]
        peakScale=x[peakInd+1]
        peakSNR=x[peakInd+2]
        peakArea=x[peakInd+3]

        bestFit=testPeakForGauss(eic, times, peakCenter, peakScale)
        if(!is.na(bestFit[1]) && bestFit[1]>=minCorr && bestFit[2]<=maxSquaredError){
          ret=c(ret, peakCenter, peakScale, peakSNR, peakArea)
        }

        peakInd=peakInd+4
      }

      if(length(ret)==0){
        return(c(-1))
      }
      return(array(ret, c(4, length(ret)/4)))
    }else{
      return(x)
    }
}










genName <- function(n = 1) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

testPeakForGauss<-function(eic, times, peakCenter, peakScale, shiftallowed=2, widths=seq(0.1, 0.7, 0.05),
                           debug=FALSE){

  if(debug){
    png(paste0("C:/tmp/", genName(), ".png"))
  }

  psm=2
  bestFit=c(NA,NA,NA,NA)
  for(i in (-shiftallowed):shiftallowed){
    for(j in widths){
      model=dnorm(times/60., mean=times[peakCenter+i]/60., sd=j)

      co=cor(eic[(peakCenter-peakScale*psm):(peakCenter+peakScale*psm)],
             model[(peakCenter-peakScale*psm):(peakCenter+peakScale*psm)],
             method="pearson")

      v1=eic[(peakCenter-peakScale*psm):(peakCenter+peakScale*psm)]
      v2=model[(peakCenter-peakScale*psm):(peakCenter+peakScale*psm)]

      v1=v1/max(v1)
      v2=v2/max(v2)

      sumOfSquaredDiff=sum((v1-v2)^2)

      if(is.na(bestFit[1]) || co>bestFit[1]){
        bestFit=c(co, sumOfSquaredDiff, i, j)
      }
    }
  }

  if(debug){
    model=dnorm(times/60., mean=times[peakCenter+bestFit[3]]/60., sd=bestFit[4])

    v1=eic[(peakCenter-peakScale*psm):(peakCenter+peakScale*psm)]
    v2=model[(peakCenter-peakScale*psm):(peakCenter+peakScale*psm)]
    vt=times[(peakCenter-peakScale*psm):(peakCenter+peakScale*psm)]/60.

    v1=v1/max(v1)
    v2=v2/max(v2)

    plot(vt, v1, type="l", lwd=5)
    lines(vt, v2, col="red", lwd=5)
    title(paste0("Cor: ", bestFit[1], " SSE: ", bestFit[2]))
  }

  if(debug){
    dev.off()
  }

  return(bestFit)
}









