import sys
import base64
import zlib
import struct
import xml.parsers.expat
from xml.dom.minidom import parse

from MSScan import MS1Scan, MS2Scan

from utils import Bunch

# reads and holds the information of one MZXML file
class Chromatogram():
    def __init__(self):
        self.msLevel = 0
        self.current_tag = ''
        self.tag_level = 0
        self.MS1_list = []
        self.MS2_list = []
        self.msInstruments = {}

    def getMS1ScanCount(self):
        return len(self.MS1_list)

    def getMS1SignalCount(self):
        c = 0
        for scan in self.MS1_list:
            c += len(scan.mz_list)
        return c

    def getAllMS1Intensities(self):
        h = []
        for scan in self.MS1_list:
            h.extend(scan.intensity_list)
        return h

    def getIthMS1Scan(self, index, filterLine):
        i = 0
        for scan in self.MS1_list:
            if scan.filter_line == filterLine:
                if i == index:
                    return scan
                i += 1
        return None

    def getMS1ScanByNum(self, scanNum):
        return self.MS1_list[scanNum]

    def getClosestMS1Scan(self, scanTimeMin, filterLine=""):
        fscan = self.MS1_list[0]
        for scan in self.MS1_list:
            if scan.filter_line==filterLine and abs(scan.retention_time / 60. - scanTimeMin) < abs(fscan.retention_time / 60. - scanTimeMin):
                    fscan = scan
        return fscan

    def getScanByID(self, id):
        for scan in self.MS1_list:
            if scan.id == id:
                return scan
        for scan in self.MS2_list:
            if scan.id == id:
                return scan

    def getFilterLines(self, includeMS1=True, includeMS2=False, includePosPolarity=True, includeNegPolarity=True):
        filterLines = set()
        if includeMS1:
            for scan in self.MS1_list:
                if (includePosPolarity and scan.polarity=="+") or (includeNegPolarity and scan.polarity=="-"):
                    filterLines.add(scan.filter_line)
        if includeMS2:
            for scan in self.MS2_list:
                if (includePosPolarity and scan.polarity=="+") or (includeNegPolarity and scan.polarity=="-"):
                    filterLines.add(scan.filter_line)

        return filterLines

    def getFilterLinesExtended(self, includeMS1=True, includeMS2=False, includePosPolarity=True, includeNegPolarity=True):
        filterLines = {}
        if includeMS1:
            for scan in self.MS1_list:
                if (includePosPolarity and scan.polarity=="+") or (includeNegPolarity and scan.polarity=="-"):
                    if scan.filter_line not in filterLines.keys():
                        filterLines[scan.filter_line]=Bunch(scanType="MS1", polarity=scan.polarity, targetStartTime=10000000, targetEndTime=0)
                    filterLines[scan.filter_line].targetStartTime=min(filterLines[scan.filter_line].targetStartTime, scan.retention_time)
                    filterLines[scan.filter_line].targetEndTime=min(filterLines[scan.filter_line].targetEndTime, scan.retention_time)

        if includeMS2:
            for scan in self.MS2_list:
                if (includePosPolarity and scan.polarity=="+") or (includeNegPolarity and scan.polarity=="-"):
                    if scan.filter_line not in filterLines.keys():
                        filterLines[scan.filter_line]=Bunch(scanType="MS2", polarity=scan.polarity, targetStartTime=10000000, targetEndTime=0, preCursorMz=[], colisionEnergy=0, scanTimes=[])
                    filterLines[scan.filter_line].targetStartTime=min(filterLines[scan.filter_line].targetStartTime, scan.retention_time)
                    filterLines[scan.filter_line].targetEndTime=max(filterLines[scan.filter_line].targetEndTime, scan.retention_time)
                    filterLines[scan.filter_line].preCursorMz.append(scan.precursor_mz)
                    filterLines[scan.filter_line].colisionEnergy=scan.colisionEnergy
                    filterLines[scan.filter_line].scanTimes.append(scan.retention_time)

            for k, v in filterLines.items():
                if v.scanType=="MS2":
                    v.preCursorMz=sum(v.preCursorMz)/len(v.preCursorMz)

        return filterLines

    def getFilterLinesPerPolarity(self, includeMS1=True, includeMS2=False):
        filterLines = {'+':set(), '-':set()}
        if includeMS1:
            for scan in self.MS1_list:
                filterLines[scan.polarity].add(scan.filter_line)
        if includeMS2:
            for scan in self.MS2_list:
                filterLines[scan.polarity].add(scan.filter_line)

        return filterLines

    def getPolarities(self):
        polarities = set()
        for scan in self.MS1_list:
            polarities.add(scan.polarity)
        for scan in self.MS2_list:
            polarities.add(scan.polarity)

        return polarities

    def getTIC(self, filterLine="", useMS2=False):
        TIC = []
        times = []
        scanIds = []

        msLevelList=self.MS1_list
        if useMS2:
            msLevelList=self.MS2_list

        for scan in msLevelList:
            if filterLine == "" or scan.filter_line == filterLine:
                TIC.append(scan.total_ion_current)
                #TIC.append(sum(scan.intensity_list))
                times.append(scan.retention_time)
                scanIds.append(scan.id)

        return TIC, times, scanIds

    # returns a specific area (2-dimensionally bound in rt and mz direction) of the LC-HRMS data
    def getArea(self, startScan, endScan, mz, ppm, filterLine="", intThreshold=0):
        scans = []
        times = []
        scanIDs = []

        scannum = 0
        for scan in self.MS1_list:
            if filterLine == "" or scan.filter_line == filterLine:
                if endScan >= scannum >= startScan:
                    bounds = scan.findMZ(mz, ppm)
                    curScan = []
                    if bounds[0] != -1:
                        for cur in range(bounds[0], bounds[1] + 1):
                            intensity = scan.intensity_list[cur]
                            if intensity >= intThreshold:
                                curScan.append((scan.mz_list[cur], intensity))
                    scans.append(curScan)
                    times.append(scan.retention_time)
                    scanIDs.append(scan.id)

                scannum += 1

        return scans, times, scanIDs

    # returns a specific area (2-dimensionally bound in rt and mz direction) of the LC-HRMS data
    def getSpecificArea(self, startTime, endTime, mzmin, mzmax, filterLine="", intThreshold=0):

        scans = []
        times = []
        scanIDs = []

        scannum = 0
        for scan in self.MS1_list:
            if filterLine == "" or scan.filter_line == filterLine:
                if startTime <= scan.retention_time <= endTime:
                    bounds = scan._findMZGeneric(mzmin, mzmax)
                    curScan = []
                    if bounds[0] != -1:
                        for cur in range(bounds[0], bounds[1] + 1):
                            intensity = scan.intensity_list[cur]
                            if intensity >= intThreshold:
                                curScan.append((scan.mz_list[cur], intensity))
                    scans.append(curScan)
                    times.append(scan.retention_time)
                    scanIDs.append(scan.id)

                scannum += 1

        return scans, times, scanIDs

    # returns an eic. Single MS peaks and such below a certain threshold may be removed
    def getEIC(self, mz, ppm, filterLine="", removeSingles=True, intThreshold=0, useMS1=True, useMS2=False, startTime=0, endTime=1000000):
        eic = []
        times = []
        scanIds = []
        mzs = []

        eicAppend = eic.append
        timesAppend = times.append
        scanIdsAppend = scanIds.append

        if useMS1:
            for scan in self.MS1_list:
                if (scan.filter_line == filterLine or filterLine == "") and startTime <= scan.retention_time <= endTime:
                    bounds = scan.findMZ(mz, ppm)
                    if bounds[0] != -1:
                        df=scan.intensity_list[bounds[0]:(bounds[1] + 1)]
                        uIndex=df.index(max(scan.intensity_list[bounds[0]:(bounds[1] + 1)]))
                        eicAppend(scan.intensity_list[bounds[0]+uIndex])
                        mzs.append(scan.mz_list[bounds[0]+uIndex])
                    else:
                        eicAppend(0)
                        mzs.append(-1)
                    timesAppend(scan.retention_time)
                    scanIdsAppend(scan.id)

        if useMS2:
            for scan in self.MS2_list:
                if (scan.filter_line == filterLine or filterLine == "") and startTime <= scan.retention_time <= endTime:
                    bounds = scan.findMZ(mz, ppm)
                    if bounds[0] != -1:
                        df=scan.intensity_list[bounds[0]:(bounds[1] + 1)]
                        uIndex=df.index(max(scan.intensity_list[bounds[0]:(bounds[1] + 1)]))
                        eicAppend(scan.intensity_list[bounds[0]+uIndex])
                        mzs.append(scan.mz_list[bounds[0]+uIndex])
                    else:
                        eicAppend(0)
                        mzs.append(-1)
                    timesAppend(scan.retention_time)
                    scanIdsAppend(scan.id)

        for i in range(1, len(eic)):
            if eic[i] < intThreshold:
                eic[i] = 0
                mzs[i] = -1

        if removeSingles:
            for i in range(1, len(eic) - 1):
                if eic[i - 1] == 0 and eic[i + 1] == 0:
                    eic[i] = 0
                    mzs[i] = -1

        return eic, times, scanIds, mzs


    def getSignalCount(self, filterLine="", removeSingles=True, intThreshold=0, useMS1=True, useMS2=False, startTime=0, endTime=1000000):
        totSignals=0
        for scan in self.MS1_list+self.MS2_list:
            if (scan.filter_line == filterLine or filterLine == "") and startTime <= scan.retention_time <= endTime:
                totSignals+=len([i for i in scan.intensity_list if i>=intThreshold])

        return totSignals


    def getMinMaxAvgSignalIntensities(self, filterLine="", removeSingles=True, intThreshold=0, useMS1=True, useMS2=False, startTime=0, endTime=1000000):
        minInt=10000000000
        maxInt=0
        avgsum=0
        avgcount=0
        for scan in self.MS1_list+self.MS2_list:
            if (scan.filter_line == filterLine or filterLine == "") and startTime <= scan.retention_time <= endTime:

                uVals=[i for i in scan.intensity_list if i>=intThreshold]
                if len(uVals)==0:
                    uVals=[0]

                minInt=min(minInt, min(uVals))
                maxInt=max(minInt, max(uVals))
                avgsum+=sum(uVals)
                avgcount+=len(uVals)

        if avgcount==0:
            return 0, 0, 0
        else:
            return minInt, maxInt, avgsum/avgcount

    # converts base64 coded spectrum in an mz and intensity array
    def decode_spectrum(self, line, compression=None, precision=32):
        idx = 0
        mz_list = []
        intensity_list = []

        if line is not None and len(line)>0:
            decoded = base64.decodestring(line)
            if compression=="zlib":
                decoded = zlib.decompress(decoded)

            if precision == 64:
                tmp_size = len(decoded) / 8
                unpack_format1 = ">%dd" % tmp_size
            else:
                tmp_size = len(decoded) / 4
                unpack_format1 = ">%df" % tmp_size


            for tmp in struct.unpack(unpack_format1, decoded):
                if idx % 2 == 0:
                    mz_list.append(float(tmp))
                else:
                    intensity_list.append(float(tmp))
                idx += 1

        return mz_list, intensity_list

    # xml parser - start element handler
    def _start_element(self, name, attrs):
        self.tag_level += 1
        self.current_tag = name
        #print "start tag", name

        if name == 'parentFile' and attrs.has_key('fileName'):
            self.parentFile = attrs['fileName']

        if name == 'msInstrument' and attrs.has_key('msInstrumentID'):
            self.msInstruments[attrs['msInstrumentID']] = {}
            self.lastMSInstrument = self.msInstruments[attrs['msInstrumentID']]

        if name == 'msModel' and len(self.msInstruments) > 0:
            self.lastMSInstrument['msModel'] = attrs['value']

        if name == 'msMassAnalyzer' and len(self.msInstruments) > 0:
            self.lastMSInstrument['msMassAnalyzer'] = attrs['value']

        if name == 'precursorMz':
            self.MS2_list[-1].precursor_intensity = float(attrs['precursorIntensity'])
            self.MS2_list[-1].precursor_charge = 0
            if attrs.has_key('precursorCharge'):
                self.MS2_list[-1].precursor_charge = int(attrs['precursorCharge'])
            if attrs.has_key("activationMethod"):
                self.MS2_list[-1].activationMethod = str(attrs['activationMethod'])

        if name == 'scan':
            self.curScan = self.curScan + 1

            self.msLevel = int(attrs['msLevel'])
            if self.msLevel == 1:
                tmp_ms = MS1Scan()
            elif self.msLevel == 2:
                tmp_ms = MS2Scan()
            else:
                print "What is it?", attrs
                sys.exit(1)

            tmp_ms.id = int(attrs['num'])

            tmp_ms.peak_count = int(attrs['peaksCount'])
            tmp_ms.peak_count_tag = int(attrs['peaksCount'])
            tmp_ms.retention_time = float(attrs['retentionTime'].strip('PTS'))
            if tmp_ms.peak_count > 0:
                if attrs.has_key('lowMz'):
                    tmp_ms.low_mz = float(attrs['lowMz'])
                if attrs.has_key('highMz'):
                    tmp_ms.high_mz = float(attrs['highMz'])
                if attrs.has_key('basePeakMz'):
                    tmp_ms.base_peak_mz = float(attrs['basePeakMz'])
                if attrs.has_key('basePeakIntensity'):
                    tmp_ms.base_peak_intensity = float(attrs['basePeakIntensity'])
            if attrs.has_key('totIonCurrent'):
                tmp_ms.total_ion_current = float(attrs['totIonCurrent'])
            else:
                tmp_ms.total_ion_current = 0
            tmp_ms.list_size = 0
            tmp_ms.polarity = str(attrs['polarity'])
            tmp_ms.encoded_mz = ''
            tmp_ms.encoded_intensity = ''
            tmp_ms.encodedData = ''
            tmp_ms.mz_list = []
            tmp_ms.intensity_list = []
            tmp_ms.msInstrumentID = ""
            if attrs.has_key('msInstrumentID'):
                tmp_ms.msInstrumentID = attrs['msInstrumentID']
            tmp_ms.filter_line = "N/A"
            if attrs.has_key('filterLine'):
                tmp_ms.filter_line = attrs['filterLine'] + " (pol: %s)" % tmp_ms.polarity
            elif len(self.msInstruments) == 1:
                tmp_ms.filter_line = "%s (MS lvl: %d, pol: %s)" % (
                    self.msInstruments.values()[0]["msModel"], self.msLevel, tmp_ms.polarity)
            else:
                tmp_ms.filter_line = "%s (MS lvl: %d, pol: %s)" % (
                    self.msInstruments[tmp_ms.msInstrumentID]["msModel"], self.msLevel, tmp_ms.polarity)

            if self.msLevel == 1:
                self.MS1_list.append(tmp_ms)
            elif self.msLevel == 2:
                if attrs.has_key("collisionEnergy"):
                    tmp_ms.collisionEnergy=float(attrs["collisionEnergy"])
                else:
                    tmp_ms.collisionEnergy=-1
                tmp_ms.ms1_id = self.MS1_list[-1].id
                self.MS2_list.append(tmp_ms)

        if name == "peaks":
            if self.msLevel==1:
                curScan=self.MS1_list[-1]
            elif self.msLevel==2:
                curScan=self.MS2_list[-1]

            curScan.compression=None
            if attrs.has_key("compressionType"):
                curScan.compression=str(attrs["compressionType"])
            curScan.precision = int(attrs['precision'])

    def _findMZGeneric(self, mzleft, mzright, mzlist):
        if len(mzlist) == 0:
            return -1, -1

        min = 0
        max = len(mzlist) - 1
        peakCount = len(mzlist)

        while min <= max:
            cur = int((max + min) // 2)

            if mzleft <= mzlist[cur] <= mzright:
                leftBound = cur
                while leftBound > 0 and mzlist[leftBound - 1] >= mzleft:
                    leftBound -= 1

                rightBound = cur
                while (rightBound + 1) < peakCount and mzlist[rightBound + 1] <= mzright:
                    rightBound += 1

                return leftBound, rightBound

            if mzlist[cur] > mzright:
                max = cur - 1
            else:
                min = cur + 1

        return -1, -1


    # xml parser - end element handler
    def _end_element(self, name):
        #print "end tag", name

        if name == 'scan':
            if self.msLevel==1:
                curScan=self.MS1_list[-1]
            elif self.msLevel==2:
                curScan=self.MS2_list[-1]

            mz_list, intensity_list = self.decode_spectrum(curScan.encodedData, compression=curScan.compression, precision=curScan.precision)

            ## Remove peaks below threshold
            if self.intensityCutoff>0:
                mz_list = [mz_list[i] for i in range(len(intensity_list)) if intensity_list[i] > self.intensityCutoff]
                intensity_list = [intensity_list[i] for i in range(len(intensity_list)) if intensity_list[i] > self.intensityCutoff]

            ## Remove peaks with unwanted mz values
            if self.mzFilter != None:
                keep = []
                for ps, pe in self.mzFilter:
                    a,b=self._findMZGeneric(ps, pe, mz_list)
                    if a!=-1:
                        keep.extend(range(a,b+1))
                mz_list = [mz_list[i] for i in sorted(keep)]
                intensity_list = [intensity_list[i] for i in sorted(keep)]

            assert len(mz_list) == len(intensity_list)

            curScan.mz_list=mz_list
            curScan.intensity_list=intensity_list
            curScan.peak_count=len(mz_list)
            curScan.encodedData=None


        if name == "precursorMz":
            self.MS2_list[-1].precursor_mz = float(self.MS2_list[-1].precursor_mz_data)
            #self.MS2_list[-1].filter_line="%s (PreMZ: %.2f ColEn: %.1f ActMet: %s)"%(self.MS2_list[-1].filter_line, self.MS2_list[-1].precursor_mz, self.MS2_list[-1].collisionEnergy, self.MS2_list[-1].activationMethod if hasattr(self.MS2_list[-1], "activationMethod") else "-")

        self.tag_level -= 1
        self.current_tag = ''
        self.msLevel == 0

    # xml parser - CDATA handler
    def _char_data(self, data):

        if self.current_tag == 'precursorMz':
            self.MS2_list[-1].precursor_mz_data+=data

        if self.current_tag == 'peaks' and not self.ignorePeaksData:
            if self.msLevel == 1:
                self.MS1_list[-1].encodedData+=data

            elif self.msLevel == 2:
                self.MS2_list[-1].encodedData+=data

    # parses an mzxml file and stores its information in the respective object.
    # If intensityCutoff is set, all MS peaks below this threshold will be discarded.
    # If ignoreCharacterData is set, only the metainformation will be parsed but not the actual MS data
    def parseMZXMLFile(self, filename_xml, intensityCutoff=-1, ignoreCharacterData=False):
        self.intensityCutoff = intensityCutoff
        self.curScan = 1

        expat = xml.parsers.expat.ParserCreate()
        expat.StartElementHandler = self._start_element
        expat.EndElementHandler = self._end_element
        self.ignorePeaksData=ignoreCharacterData
        expat.CharacterDataHandler = self._char_data

        ##expat.Parse(content_listd)
        expat.ParseFile(open(filename_xml, 'r'))

        if not ignoreCharacterData:
            for scan in self.MS1_list:

                assert len(scan.mz_list) == len(scan.intensity_list) == scan.peak_count
                #if intensityCutoff<0:
                #    assert scan.peak_count == scan.peak_count_tag

        signals=0
        for msscan in self.MS1_list:
            signals+=len(msscan.mz_list)


        for ionMode in ["+", "-"]:

            precursors=[]
            for msmsscan in self.MS2_list:
                if msmsscan.polarity==ionMode:
                    precursors.append(msmsscan.precursor_mz)

            precursors=sorted(precursors)

            packets=[]
            curpack=[]

            for precursorMZ in precursors:
                if len(curpack)==0:
                    curpack.append(precursorMZ)
                else:
                    mzmean=sum(curpack)/len(curpack)
                    if abs(mzmean-precursorMZ)/precursorMZ*1E6>10:
                        packets.append(curpack)
                        curpack=[]
                    curpack.append(precursorMZ)
            if len(curpack)>0:
                packets.append(curpack)

            for curpack in packets:
                meanmz = sum(curpack) / len(curpack)
                #print ionMode, meanmz, (meanmz-min(curpack))/meanmz*1E6, (max(curpack)-meanmz)/meanmz*1E6, curpack

            for msmsscan in self.MS2_list:
                if msmsscan.polarity == ionMode:
                    for curpack in packets:
                        meanmz=sum(curpack)/len(curpack)
                        if abs(msmsscan.precursor_mz-meanmz)/msmsscan.precursor_mz*1E6<=10:
                            #msmsscan.precursor_mz=meanmz
                            msmsscan.filter_line = "%s (PreMZ: %.4f [%.6f-%6f; %.1f ppm] ColEn: %.1f ActMet: %s)" % (msmsscan.filter_line, meanmz, min(curpack), max(curpack), (max(curpack)-min(curpack))/meanmz*1E6, msmsscan.collisionEnergy, msmsscan.activationMethod if hasattr(msmsscan, "activationMethod") else "-")




    def parseMzMLFile(self, filename_xml, intensityCutoff, ignoreCharacterData):

        import pymzml

        run=pymzml.run.Reader(filename_xml)

        for specturm in run:

            if specturm.has_key("time array") or specturm["id"] is None:
                continue

            try:
                msLevel=int(specturm["ms level"])
            except :
                print "Error: What is it?", specturm["id"], type(specturm), specturm
                continue


            if msLevel == 1:
                tmp_ms = MS1Scan()
            elif msLevel == 2:
                tmp_ms = MS2Scan()
            else:
                print "What is it?", msLevel, specturm["id"]
                sys.exit(1)

            tmp_ms.id = int(specturm["id"])
            tmp_ms.filter_line = specturm["filter string"]

            tmp_ms.peak_count = len(specturm.peaks)
            if "scan time" in specturm.keys():
                tmp_ms.retention_time = specturm["scan time"]*60
            elif "scan start time" in specturm.keys():
                tmp_ms.retention_time = specturm["scan start time"]*60
            else:
                raise Exception("no scan retention time found")
            #if tmp_ms.peak_count > 0:
            tmp_ms.total_ion_current = specturm["total ion current"]
            tmp_ms.list_size = 0

            if specturm.has_key("positive scan"):
                tmp_ms.polarity = "+"
            elif specturm.has_key("negative scan"):
                tmp_ms.polarity = "-"
            else:
                raise RuntimeError("No polarity for scan available")
            tmp_ms.mz_list = [p[0] for p in specturm.peaks]
            tmp_ms.intensity_list = [p[1] for p in specturm.peaks]
            tmp_ms.msInstrumentID = ""


            if msLevel == 1:
                self.MS1_list.append(tmp_ms)
            elif msLevel == 2:
                tmp_ms.precursor_mz = specturm["precursors"][0]["mz"]
                tmp_ms.precursor_intensity = 0
                tmp_ms.precursor_charge = specturm["precursors"][0]["charge"]
                self.MS2_list.append(tmp_ms)

        precursors=[]
        i=0
        polarities=set()
        for ms2Scan in self.MS2_list:
            precursors.append(Bunch(precursormz=ms2Scan.precursor_mz, id=i, polarity=ms2Scan.polarity))
            polarities.add(ms2Scan.polarity)
            i=i+1

        from utils import mean

        for polarity in polarities:
            precursorsTemp=[pc for pc in precursors if pc.polarity==polarity]
            precursorsTemp=sorted(precursorsTemp, key=lambda x: x.precursormz)

            lastmz=-1000
            lastmzs=[]

            for i in range(len(precursorsTemp)):
                if (precursorsTemp[i].precursormz-lastmz)>=25*lastmz/1000000.:
                    if len(lastmzs)>0:
                        minMZ=lastmzs[0].precursormz
                        maxMZ=lastmzs[-1].precursormz
                        #print polarity, minMZ, maxMZ, (maxMZ-minMZ)*1000000./minMZ, "  \n --> ", (precursorsTemp[i].precursormz-lastmz)*1000000./lastmz, ": ", precursorsTemp[i].precursormz

                        for j in range(len(lastmzs)):
                            self.MS2_list[lastmzs[j].id].filter_line="MSn %s mean mz precursor: %.5f"%(polarity, mean([pc.precursormz for pc in lastmzs]))
                            self.MS2_list[lastmzs[j].id].precursor_mz=float(mean([pc.precursormz for pc in lastmzs]))
                    lastmz=precursorsTemp[i].precursormz
                    lastmzs=[precursorsTemp[i]]
                else:
                    lastmzs.append(precursorsTemp[i])
            if len(lastmzs)>0:
                for j in range(len(lastmzs)):
                    self.MS2_list[lastmzs[j].id].filter_line="MSn %s mean mz precursor: %.5f"%(polarity, mean([pc.precursormz for pc in lastmzs]))
                    self.MS2_list[lastmzs[j].id].precursor_mz=float(mean([pc.precursormz for pc in lastmzs]))






    def parse_file(self, filename_xml, intensityCutoff=-1, ignoreCharacterData=False, mzFilter=None):
        self.mzFilter=mzFilter
        if filename_xml.lower().endswith(".mzxml"):
            return self.parseMZXMLFile(filename_xml, intensityCutoff, ignoreCharacterData)
        elif filename_xml.lower().endswith(".mzml"):
            return self.parseMzMLFile(filename_xml, intensityCutoff, ignoreCharacterData)
        else:
            raise RuntimeError("Invalid file type")


    # updates the data of this MS scan
    # WARNING: unstable, horridly implemented method. SUBJECT OF CHANGE
    def updateScan(self, scan, data):
        scanID = int(scan.getAttribute("num"))

        #  tmp_ms.retention_time = float(attrs['retentionTime'].strip('PTS'))
        if data.has_key(scanID) and hasattr(data[scanID], "rt"):
            scan.setAttribute("retentionTime", "PT%.3fS"%(data[scanID].rt))

        peaks = None
        for kid in scan.childNodes:
            if kid.nodeName == "peaks":
                peaks = kid
        assert peaks is not None and len(peaks.childNodes) == 1

        if data.has_key(scanID) and len(data[scanID].mzs) > 0:
            h = data[scanID]
            assert len(h.mzs) == len(h.ints)
            peaks.childNodes[0].nodeValue = base64.encodestring("".join([struct.pack(">f", h.mzs[i]) + struct.pack(">f", h.ints[i]) for i in range(len(h.mzs))])).strip().replace("\n", "")
            scan.setAttribute("peaksCount", str(len(h.mzs)))
            scan.setAttribute("lowMz", str(min(h.mzs)))
            scan.setAttribute("highMz", str(max(h.mzs)))
            scan.setAttribute("basePeakMz", "0")
            scan.setAttribute("basePeakIntensity", "0")
            scan.setAttribute("totIonCurrent", str(sum(h.ints)))
        else:
            peaks.childNodes[0].nodeValue = ""
            scan.setAttribute("peaksCount", "0")
            scan.setAttribute("lowMz", "0")
            scan.setAttribute("highMz", "0")
            scan.setAttribute("basePeakMz", "0")
            scan.setAttribute("basePeakIntensity", "0")
            scan.setAttribute("totIonCurrent", "0")

    # updates the data of this mzxml file
    # WARNING: unstable, horridly implemented method. SUBJECT OF CHANGE
    def updateData(self, node, data):
        for kid in node.childNodes:
            if kid.nodeName != '#text':
                if kid.nodeName == "scan":
                    self.updateScan(kid, data)
                else:
                    self.updateData(kid, data)

    # updates the data of this mzxml file
    # WARNING: unstable, horridly implemented method. SUBJECT OF CHANGE
    def resetMZData(self, mzxmlfile, toFile, data):
        dom = parse(mzxmlfile)
        self.updateData(dom, data)
        dom.writexml(open(toFile, "wb"))

    def freeMe(self):
        for scan in self.MS2_list:
            scan.freeMe()
        for scan in self.MS1_list:
            scan.freeMe()


if __name__=="__main__":
    f="F:/MaxPlanck_FrederikDethloff/exp4/POS/pos-MF-1_25_01_4001.mzXML"
    f="E:/___Backup/Publications/Manuscripts/MetExtractII/_assets/geoRge/MTBLS213_20170613_073311/CELL_Glc13_05mM_Normo_01.mzXML"
    f="H:/20181016_472_Alternaria_ExperimentForDFGProposal/FTICRMS/F_negNS300147444mL1zu109_000001.mzXML"

    t=Chromatogram()
    t.parse_file(f)

    for i in t.getFilterLinesPerPolarity(includeMS1=False, includeMS2=True):
        print i
        for j in i:
            print "   ", j

    print t.getPolarities()

    x=Chromatogram()
    x.parse_file(f)




    scan=x.getMS1ScanByNum(0)
    print scan.retention_time/60


    for mz, inte in zip(scan.mz_list, scan.intensity_list):
        print mz, inte

    import matplotlib.pyplot as plt

    plt.vlines(x=scan.mz_list, ymin=0, ymax=scan.intensity_list)

    plt.show()