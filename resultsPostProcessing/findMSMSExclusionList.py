version=0.1

#################################################################################
################################################################################
### Import dependencies
##
#
from utils import Bunch
from Chromatogram import Chromatogram
import csv

from matplotlib import pyplot as plt





#################################################################################
################################################################################
### Parameters
##
#
inF="E:/170308_posneg_326_337_335_MSMSmeasurements/12C13C_leaf_forexclusionlist.mzXML"
ppm=25
minWallScans=25









#################################################################################
################################################################################
### Load data
##
#
chrom=Chromatogram()
chrom.parse_file(inF)

data=[]
for scan in chrom.MS2_list:
    data.append([scan.retention_time, scan.precursor_mz, None])




#################################################################################
################################################################################
### Process
##
#
data=sorted(data, key=lambda x:x[1])

bins=[]
curbin=[]
lastS=None
for si, s in enumerate(data):
    rt, mz, groupID=s

    if len(curbin)==0:
        curbin.append([rt, mz, si])
    else:
        if abs(mz-lastS[1])*1000000./mz<=ppm:
            curbin.append([rt, mz, si])
        else:
            bins.append(curbin)
            curbin=[[rt, mz, si]]

    data[si][2] = len(bins)
    lastS=[rt, mz]


exclusionList=[]
for bini, bin in enumerate(bins):
    for si, s in enumerate(bin):
        data[s[2]][2]=1 if len(bin)>=minWallScans else 0
    if len(bin)>=minWallScans:
        start=min([s[0] for s in bin])
        stop=max([s[0] for s in bin])
        avmz=sum([s[1] for s in bin])*1./len(bin)
        exclusionList.append([start, stop, avmz])






fileLineArray = []

fileLineArray.append('<?xml version="1.0" encoding="ISO-8859-1"?>')
from time import gmtime, strftime
fileLineArray.append(
    '<featureMap version="1.4" id="fm_16311276685788915066" xsi:noNamespaceSchemaLocation="http://open-ms.sourceforge.net/schemas/FeatureXML_1_4.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">')
fileLineArray.append('	<dataProcessing completion_time="%s">' % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
fileLineArray.append('		<software name="findMSMSExclusionList" version="%s" />' % (version))
fileLineArray.append('		<processingAction name="Feature pair detection" />')
# fileLineArray.append('		<UserParam type="string" name="parameter: key" value="value"/>')
fileLineArray.append('	</dataProcessing>')
fileLineArray.append('	<featureList count="%d">")' % (len(exclusionList)))

for featurei, feature in enumerate(exclusionList):

    fileLineArray.append('		<feature id="%s">' % (featurei))
    fileLineArray.append('			<position dim="0">%f</position>' % ((feature[0]+feature[1])/2.))
    fileLineArray.append('			<position dim="1">%f</position>' % (feature[2]))
    fileLineArray.append('			<intensity>1</intensity>')
    fileLineArray.append('			<quality dim="0">-1</quality>')
    fileLineArray.append('			<quality dim="1">-1</quality>')
    fileLineArray.append('			<overallquality>-1</overallquality>')
    fileLineArray.append('			<convexhull nr="1">')
    fileLineArray.append('				<pt x="%f" y="%f" />' % (feature[0], feature[2]))
    fileLineArray.append('				<pt x="%f" y="%f" />' % (feature[1], feature[2]))
    fileLineArray.append('				<pt x="%f" y="%f" />' % (feature[0], feature[2]))
    fileLineArray.append('				<pt x="%f" y="%f" />' % (feature[1], feature[2]))
    fileLineArray.append('			</convexhull>')
    fileLineArray.append('		</feature>')

fileLineArray.append('	</featureList>')
fileLineArray.append('</featureMap>')





with open(inF+".featureML", "wb") as fOut:
    for line in fileLineArray:
        fOut.write(line)
        fOut.write("\r\n")




plt.scatter([s[0] for s in data], [s[1] for s in data], c=[s[2]%8 for s in data], cmap=plt.get_cmap("Accent"))
plt.show()