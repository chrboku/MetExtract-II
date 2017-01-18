# MetExtract-II

MetExtract II is a software tool for stable isotope labeling (SIL) assisted untargeted metabolomics research using liquid chromatography high resolution mass spectrometry (LC-HRMS). It comprises of 3 modules:

* ### AllExtract 
This module requires LC-HRMS data obtained from samples containing both native and uniformly as well as highly (≥98%) isotopically enriched samples material. AllExtract detects metabolite ions using the mirror-symmetric isotopolog patterns. Any substances that do not show these charaterisitica are discarded as non-sample derived compounds. Thus, this module efficiently detects only metabolites of the studied biological system. Additionally, using experiment- and metabolome wide internal standardisation, the analytical precision of the measurements and consequently the subsequently performed statistical analysis of the experiment are improved.

* ### TracExtract 
This module is designed to investigate the metabolic fate of precursor substances (tracers) and requires the tracer to be applied as both a native and a highly isotopically enriched form. It mainly supports biotransformation experiments using secondary metabolites or such compounds, which are not heavily fragmented in the biological system and incorporated as larger parts in downstream secondary metabolites. Other than the AllExtract module, it is not designed to detect all metabolites of the biological system but will only report the biotransformation products of the studied tracer.

* ### FragExtract 
The last module of MetExtract II is designed for processing LC-HRMS/MS datasets derived from MS/MS experiments of native and highly isotopically enriched metabolite ions. Using the distinct characteristics of SIL, FragExtract is able to efficiently clean and annotate the fragmentation spectra of unknown metabolites therefor supporting metabolite annotation and identification.

For more information and pre-compiled binaries please visit https://metabolomics-ifa.boku.ac.at/metextractII/
