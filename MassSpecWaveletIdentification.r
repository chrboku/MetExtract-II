
suppressWarnings(suppressMessages(library(waveslim)))
suppressWarnings(suppressMessages(library(MassSpecWavelet)))
suppressWarnings(suppressMessages(library(signal)))

pack<-installed.packages()



getMajorPeaksimpl<-function(eic, scales=c(5,33), snrTh=1){

    ## Prepare scales
    #I don't know why, but the scale factor 1 is required
    if(scales[1]==1){
        scales <- seq(scales[1], scales[2], 1)
    }else{
        scales <- c(1, seq(scales[1], scales[2], 1))
    }

    ## Detect peaks
    peakInfo <- peakDetectionCWT(eic, scales=scales, SNR.Th=snrTh)#, tuneIn=FALSE, peakScaleRange = scales[1], amp.Th = 0.00000000000001, winSize.noise=20, SNR.method="mad")
    majorPeakInfo <- peakInfo$majorPeakInfo
    peakIndex <- majorPeakInfo$peakIndex


    ## No peaks detected
    if(length(peakIndex)==0){
        return(c(-1))
    }


    ## Peaks detected, format return array
    ret=c()
    for(p in peakIndex){
      ind=paste("1_", toString(p), sep="")
      if(!is.na(majorPeakInfo$peakValue[ind])){
        ret=c(ret, majorPeakInfo$peakIndex[ind], majorPeakInfo$peakScale[ind], majorPeakInfo$peakSNR[ind], majorPeakInfo$peakValue[ind])
      }
    }
    return(array(ret, c(4,length(ret)/4)))
}

getMajorPeaks<-function(eic, scales=c(5,33), snrTh=1, eicSmoothing="None"){
    x<-try(getMajorPeaksimpl(eic, scales, snrTh), silent=TRUE)
    return(x)
}