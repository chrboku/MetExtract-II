is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])

missingWaveslim=FALSE
missingSignal=FALSE
missingPTW=FALSE
missingMassSpecWavelet=FALSE
missingMASS=FALSE
#missingRGL=FALSE

if(!is.installed('waveslim')){
    cat("Required Package 'waveslim' missing\n")
    missingWaveslim=TRUE
}
if(!is.installed('signal')){
    cat("Required Package 'signal' missing\n")
    missingSignal=TRUE
}
if(!is.installed('ptw')){
    cat("Required Package 'ptw' missing\n")
    missingPTW=TRUE
}
if(!is.installed('MASS')){
    cat("Required Package 'MASS' missing\n")
    missingMASS=TRUE
}
#if(!is.installed('rgl')){
#    cat("Required Package 'rgl' missing\n")
#    missingRGL=TRUE
#}
if(!is.installed('MassSpecWavelet')){
    cat("Required Package 'MassSpecWavelet' missing\n")
    missingMassSpecWavelet=TRUE
}


if(missingWaveslim){
    install.packages("waveslim")
    cat("Package 'waveslim' installed")
}
if(missingSignal){
    install.packages("signal")
    cat("Package 'signal' installed")
}
if(missingPTW){
    install.packages("ptw")
    cat("Package 'ptw' installed")
}
if(missingMASS){
    install.packages("MASS")
    cat("Package 'MASS' installed")
}
#if(missingRGL){
#    install.packages("rgl")
#    cat("Package 'rgl' installed")
#}
if(missingMassSpecWavelet){
    source("http://bioconductor.org/biocLite.R")
    biocLite("MassSpecWavelet")
    cat("Package 'MassSpecWavelet' installed")
}
