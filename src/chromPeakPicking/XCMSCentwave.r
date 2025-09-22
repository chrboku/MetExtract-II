if(FALSE){
    source("https://bioconductor.org/biocLite.R")
    biocLite("MassSpecWavelet")
}
suppressWarnings(suppressMessages(library(xcms)))


## Centwave algorithm developed and published by
## Tautenhahn, R., Böttcher, C. & Neumann, S. Highly sensitive feature detection for high resolution LC/MS. BMC Bioinformatics 9, 504 (2008). https://doi.org/10.1186/1471-2105-9-504

## Code for integration adapted from MzMine 2
## T. Pluskal, S. Castillo, A. Villar-Briones, M. Orešič, MZmine 2: Modular framework for processing, visualizing, and analyzing mass spectrometry-based molecular profile data, BMC Bioinformatics 11:395 (2010). PMID: 20650010
## https://github.com/mzmine/mzmine2/blob/master/src/main/java/net/sf/mzmine/modules/peaklistmethods/peakpicking/deconvolution/centwave/CentWaveDetector.java
getCentwavePeaks<-function(scanTime, intensity, peakwidth=c(5, 30)){

    mz = 300
    numPoints = length(scanTime)

    xRaw = new("xcmsRaw")
    xRaw@tic = intensity
    xRaw@scantime = scanTime
    xRaw@scanindex = 0:(numPoints-1)
    xRaw@env$mz = rep(mz, numPoints)
    xRaw@env$intensity = intensity

    ROIs = list()
    roi=1

    start=1
    while(start <= length(scanTime)){
      if(intensity[start]>0){
        end=start+1
        while(end <= length(scanTime) && intensity[end] > 0){
          end=end+1
        }

        ROIs[[roi]]=list(scmin=start, scmax=end-1, mzmin=mz, mzmax=mz)
        start=end
        roi=roi+1
      }else{
        start=start+1
      }
    }

    mtx <- findPeaks.centWave(xRaw, ppm=1, mzdiff=-100, verbose=TRUE, peakwidth=c(5,50), snthresh=3, integrate=2, ROI.list=ROIs, firstBaselineCheck=TRUE)

    m=matrix(mtx, nrow=length(mtx)/22, byrow=F)
    m=m[,c(4:13, 17, 18, 19, 20), drop=F]
    colnames(m)=c("rt", "rtmin", "rtmax", "into", "intb", "maxo", "sn", "egauss", "mu", "sigma", "scale", "scpos", "scmin", "scmax")
    for(i in 1:nrow(m)){
        m[i, "scmin"]=min(length(scanTime), max(1, which.min(abs(scanTime-m[i, "rtmin"]))))
        m[i, "scmax"]=min(length(scanTime), max(1, which.min(abs(scanTime-m[i, "rtmax"]))))
        m[i, "scpos"]=min(length(scanTime), max(1, m[i, "scmin"]+which.max(intensity[m[i, "scmin"]:m[i, "scmax"]])))-1
        cat(m[i, "scmin"], m[i, "scpos"], m[i, "scmax"], length(scanTime), "\n")

        m[i, "rt"]=scanTime[m[i, "scpos"]]
        m[i, "rtmin"]=scanTime[m[i, "scmin"]]
        m[i, "rtmax"]=scanTime[m[i, "scmax"]]
    }

    return(m)
}











