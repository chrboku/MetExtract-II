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


# example
# Calculate isotopic enrichment using M+1 and M (or M'-1 and M')
# Note: Use determined number of carbon atoms (delta M' and M) only in AllExtract data processing or 
#       to determine the isotopic enrichment of the 13C-labelled tracer
substitutions=1

Cn=10
ratioOfMpOneToM=.15
print(sprintf("Enrichment is %.2f%% at a ratio of %.4f with Cn %d at M+%d/M (or M'-%d/M')", getIsotopicEnrichment(Cn, substitutions, ratioOfMpOneToM)*100, ratioOfMpOneToM, Cn, substitutions, substitutions))

CnPos=0:Cn
natRatios=getAbsoluteIsotopePatternRatio(Cn, CnPos, 0.9893)
labRatios=getAbsoluteIsotopePatternRatio(Cn, CnPos, 0.96)

labRatios2=c(0,0,getAbsoluteIsotopePatternRatio(Cn-2, CnPos[0:(length(CnPos)-2)], 0.96))*.5+labRatios
labRatios2=labRatios2/max(labRatios2)

m=rbind(natRatios, rev(labRatios), rev(labRatios2))
rownames(m)=c("Natural", "Labelled", "Labelled-Mixed")

barplot(m, beside=T, col=c("olivedrab", "firebrick", "slategrey"))
legend("top", legend=rownames(m), col=c("olivedrab", "firebrick", "slategrey"), lwd=5)


