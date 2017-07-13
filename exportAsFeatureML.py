from utils import Bunch

from time import gmtime, strftime
import csv

from MetExtractII_Main import MetExtractVersion


def writeFeatureListToFeatureML(features, toFile, ppmPM=5, rtPM=0.25*60):

    fileLineArray=[]

    fileLineArray.append('<?xml version="1.0" encoding="ISO-8859-1"?>')
    fileLineArray.append('<featureMap version="1.4" id="fm_16311276685788915066" xsi:noNamespaceSchemaLocation="http://open-ms.sourceforge.net/schemas/FeatureXML_1_4.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">')
    fileLineArray.append('	<dataProcessing completion_time="%s">'%(strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
    fileLineArray.append('		<software name="MetExtract II" version="%s" />'%(MetExtractVersion))
    fileLineArray.append('		<processingAction name="Feature pair detection" />')
    #fileLineArray.append('		<UserParam type="string" name="parameter: key" value="value"/>')
    fileLineArray.append('	</dataProcessing>')
    fileLineArray.append('	<featureList count="%d">")'%(len(features)))

    for feature in features:
        num    =feature.id
        grpNum =feature.ogroup
        mz     =feature.mz
        lmz    =feature.lmz
        rt     =feature.rt
        z      =feature.charge
        name   =feature.name
        xn     =feature.Xn

        fileLineArray.append('		<feature id="%s">'%(num))
        fileLineArray.append('			<UserParam type="string" name="label" value="%s (Num: %s, OGroup: %s, MZ: %.5f, RT (min): %.2f, Charge: %d, Xn: %s)"/>'%(name, num, grpNum, mz, rt, z, xn))
        fileLineArray.append('			<position dim="0">%f</position>'%(rt))
        fileLineArray.append('			<position dim="1">%f</position>'%(mz))
        fileLineArray.append('			<charge>%d</charge>'%(z))
        fileLineArray.append('			<intensity>1</intensity>')
        fileLineArray.append('			<quality dim="0">-1</quality>')
        fileLineArray.append('			<quality dim="1">-1</quality>')
        fileLineArray.append('			<overallquality>-1</overallquality>')
        fileLineArray.append('			<convexhull nr="1">')
        fileLineArray.append('				<pt x="%f" y="%f" />'%(rt-rtPM, mz*(1-ppmPM/1000000.)))
        fileLineArray.append('				<pt x="%f" y="%f" />'%(rt+rtPM, mz*(1-ppmPM/1000000.)))
        fileLineArray.append('				<pt x="%f" y="%f" />'%(rt-rtPM, mz*(1+ppmPM/1000000.)))
        fileLineArray.append('				<pt x="%f" y="%f" />'%(rt+rtPM, mz*(1+ppmPM/1000000.)))
        fileLineArray.append('			</convexhull>')
        fileLineArray.append('			<convexhull nr="2_lab">')
        fileLineArray.append('				<pt x="%f" y="%f" />'%(rt-rtPM, lmz*(1-ppmPM/1000000.)))
        fileLineArray.append('				<pt x="%f" y="%f" />'%(rt+rtPM, lmz*(1-ppmPM/1000000.)))
        fileLineArray.append('				<pt x="%f" y="%f" />'%(rt-rtPM, lmz*(1+ppmPM/1000000.)))
        fileLineArray.append('				<pt x="%f" y="%f" />'%(rt+rtPM, lmz*(1+ppmPM/1000000.)))
        fileLineArray.append('			</convexhull>')
        fileLineArray.append('		</feature>')

    fileLineArray.append('	</featureList>')
    fileLineArray.append('</featureMap>')


    with open(toFile, "wb") as fOut:
        for line in fileLineArray:
            fOut.write(line)
            fOut.write("\r\n")


def convertMEMatrixToFeatureML(meMatrixFile, featureMLFile=None):
    if featureMLFile is None:
        featureMLFile=meMatrixFile.replace(".tsv", ".txt").replace(".txt", "")+".featureML"

    with open(meMatrixFile) as fIn:
        csvReader=csv.reader(fIn, delimiter="\t", quotechar="\"")

        headers={}
        features=[]
        for linei, row in enumerate(csvReader):
            if linei==0:
                for colInd, header in enumerate(row):
                    headers[header]=colInd
            elif row[0].startswith("#"):
                pass
            else:
                b=Bunch(id=row[headers["Num"]], ogroup=row[headers["OGroup"]], mz=float(row[headers["MZ"]]), rt=float(row[headers["RT"]])*60, Xn=int(row[headers["Xn"]]), lmz=float(row[headers["L_MZ"]]),
                        charge=int(row[headers["Charge"]]), name=row[headers["Num"]], ionMode=row[headers["Ionisation_Mode"]])
                features.append(b)

    writeFeatureListToFeatureML(features, featureMLFile, ppmPM=5., rtPM=0.25*60)


def convertMEMatrixToFeatureMLSepPolarities(meMatrixFile, featureMLFile=None, postfix=""):
    if featureMLFile is None:
        featureMLFile=meMatrixFile.replace(".tsv", ".txt").replace(".txt", "")+".featureML"

    features = {'+': [], '-': []}
    with open(meMatrixFile) as fIn:
        csvReader=csv.reader(fIn, delimiter="\t", quotechar="\"")

        headers={}
        for linei, row in enumerate(csvReader):
            if linei==0:
                for colInd, header in enumerate(row):
                    headers[header]=colInd
            elif row[0].startswith("#"):
                pass
            else:
                b=Bunch(id=row[headers["Num"]], ogroup=row[headers["OGroup"]], mz=float(row[headers["MZ"]]), rt=float(row[headers["RT"]])*60, Xn=int(row[headers["Xn"]]), lmz=float(row[headers["L_MZ"]]),
                        charge=int(row[headers["Charge"]]), name=row[headers["Num"]], ionMode=row[headers["Ionisation_Mode"]])
                features[b.ionMode].append(b)

    if len(features['-'])>0:
        writeFeatureListToFeatureML(features['-'], featureMLFile.replace(".featureML", "_negMode%s.featureML"%(postfix)), ppmPM=5., rtPM=0.25*60)
    if len(features['+'])>0:
        writeFeatureListToFeatureML(features['+'], featureMLFile.replace(".featureML", "_posMode%s.featureML"%(postfix)), ppmPM=5., rtPM=0.25 * 60)


