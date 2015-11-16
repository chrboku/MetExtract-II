suppressWarnings(suppressMessages(library(waveslim)))
suppressWarnings(suppressMessages(library(MassSpecWavelet)))
suppressWarnings(suppressMessages(library(signal)))
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
        peakInfo<-peakDetectionCWT(eic, scales=scales, SNR.Th=snrTh)
        
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

getMajorPeaks<-function(eic, scales=c(5,33), snrTh=1, eicSmoothing="None"){
    #cat("scales ", scales[1], " ", scales[2], "  snrTh ", snrTh, "\n")
    x<-try(getMajorPeaksd(eic, scales, snrTh, eicSmoothing), silent=TRUE)
    return(x)
}