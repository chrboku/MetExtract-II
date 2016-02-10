# example
# Calculate isotopic enrichment using M+1 and M (or M'-1 and M')
# Note: Use determined number of carbon atoms (delta M' and M) only in AllExtract data processing or 
#       to determine the isotopic enrichment of the 13C-labelled tracer
Cn=11
ratioOfMpOneToM=.19
substitutions=1


################################################################################
################################################################################
######################## Function definitions

# calculate ratio of M+s relative to M (or M'-s relative to M')
getNormalisedIsotopePatternRatio<-function(a,s,p){
  return( (p^(a-s)) * ((1.-p)^s) * choose(a,s) / (p^a) )
}
# caluclate the absolute, relative ratio of M+s (or M'-s)
# Note: similar to getNormalisedIsotopePatternRatio but no normalisation to M (or M') is performed
getAbsoluteIsotopePatternRatio<-function(a,s,p){
  return( (p^(a-s)) * ((1.-p)^s) * choose(a,s))
}
# calculates the isotopic enrichment using M+s and M (or M'-s and M')
getIsotopicEnrichment<-function(a,s,r){
  return( choose(a,s)^(1./s) / ( choose(a,s)^(1./s)+r^(1./s) ) )
}



################################################################################
################################################################################
######################## Calculate isotopic enrichment

enrichment=getIsotopicEnrichment(Cn, substitutions, ratioOfMpOneToM)




################################################################################
################################################################################
######################## Display results

print(sprintf("Enrichment is %.2f%% at a ratio of %.4f with Cn %d at M+%d/M (or M'-%d/M')", enrichment*100, ratioOfMpOneToM, Cn, substitutions, substitutions))









################################################################################
################################################################################
######################## Show results as a theoretical mass spectra
if(FALSE){
  CnPos=0:Cn
  natRatios=getAbsoluteIsotopePatternRatio(Cn, CnPos, 0.9893)     ## 98.93% is the natural enrichment with 13C
  labRatios=getAbsoluteIsotopePatternRatio(Cn, CnPos, enrichment)
  
  op <- par(mfrow = c(2,1),
            oma = c(1,2,0,0) + 0.1,
            mar = c(0,0,2,0) + 0.1)
  
  ## theoretical abundances
  m=rbind(natRatios, rev(labRatios))
  rownames(m)=c("Natural", "Labelled")
  
  barplot(m, beside=T, col=c("Olivedrab", "Firebrick"), main="Theoretical ratios")
  legend("top", legend=rownames(m), col=c("Olivedrab", "Firebrick"), lwd=5, box.lwd = 0, box.col="white")
  
  ## normalised theoretical abundances
  m=rbind(natRatios/natRatios[1], rev(labRatios)/labRatios[1])
  rownames(m)=c("Natural", "Labelled")
  
  barplot(m, beside=T, col=c("Olivedrab", "Firebrick"), main="Normalised, theoretical ratios")
  par(op)
}