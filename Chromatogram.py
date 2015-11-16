import sys
import base64
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

    def getClosestMS1Scan(self, scanTimeMin):
        fscan = self.MS1_list[0]
        for scan in self.MS1_list:
            if abs(scan.retention_time / 60. - scanTimeMin) < abs(fscan.retention_time / 60. - scanTimeMin):
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
                        filterLines[scan.filter_line]=Bunch(scanType="MS2", polarity=scan.polarity, targetStartTime=10000000, targetEndTime=0, preCursorMz=[], colisionEnergy=0)
                    filterLines[scan.filter_line].targetStartTime=min(filterLines[scan.filter_line].targetStartTime, scan.retention_time)
                    filterLines[scan.filter_line].targetEndTime=max(filterLines[scan.filter_line].targetEndTime, scan.retention_time)
                    filterLines[scan.filter_line].preCursorMz.append(scan.precursor_mz)
                    filterLines[scan.filter_line].colisionEnergy=scan.colisionEnergy

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
    def getEIC(self, mz, ppm, filterLine="", removeSingles=True, intThreshold=0):
        ret = []
        times = []
        scanIds = []
        for scan in self.MS1_list:
            if filterLine == "" or scan.filter_line == filterLine:
                bounds = scan.findMZ(mz, ppm)
                if bounds[0] != -1:
                    ret.append(max(scan.intensity_list[bounds[0]:(bounds[1] + 1)]))
                else:
                    ret.append(0)
                times.append(scan.retention_time)
                scanIds.append(scan.id)

        for i in range(1, len(ret)):
            if ret[i] < intThreshold:
                ret[i] = 0

        if removeSingles:
            for i in range(1, len(ret) - 1):
                if ret[i - 1] == 0 and ret[i + 1] == 0:
                    ret[i] = 0

        return ret, times, scanIds

    # converts base64 coded spectrum in an mz and intensity array
    def decode_spectrum(self, line):
        #print line
        decoded = base64.decodestring(line)
        tmp_size = len(decoded) / 4
        unpack_format1 = ">%dL" % tmp_size

        idx = 0
        mz_list = []
        intensity_list = []

        for tmp in struct.unpack(unpack_format1, decoded):
            tmp_i = struct.pack("I", tmp)
            tmp_f = struct.unpack("f", tmp_i)[0]
            if idx % 2 == 0:
                mz_list.append(float(tmp_f))
            else:
                intensity_list.append(float(tmp_f))
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
                tmp_ms.low_mz = float(attrs['lowMz'])
                tmp_ms.high_mz = float(attrs['highMz'])
                tmp_ms.base_peak_mz = float(attrs['basePeakMz'])
                tmp_ms.base_peak_intensity = float(attrs['basePeakIntensity'])
            tmp_ms.total_ion_current = float(attrs['totIonCurrent'])
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
                tmp_ms.ms1_id = self.MS1_list[-1].id
                self.MS2_list.append(tmp_ms)

    # xml parser - end element handler
    def _end_element(self, name):
        #print "end tag", name

        if name == 'scan':
            if self.msLevel==1:
                curScan=self.MS1_list[-1]
            elif self.msLevel==2:
                curScan=self.MS2_list[-1]

            mz_list, intensity_list = self.decode_spectrum(curScan.encodedData)
            assert len(mz_list) == len(intensity_list)
            mz_list = [mz_list[i] for i in range(len(intensity_list)) if intensity_list[i] > self.intensityCutoff]
            intensity_list = [intensity_list[i] for i in range(len(intensity_list)) if intensity_list[i] > self.intensityCutoff]
            assert len(mz_list) == len(intensity_list)
            curScan.mz_list=mz_list
            curScan.intensity_list=intensity_list
            curScan.peak_count=len(mz_list)
            curScan.encodedData=None


        if name == "precursorMz":
            self.MS2_list[-1].precursor_mz = float(self.MS2_list[-1].precursor_mz_data)
            self.MS2_list[-1].filter_line="%s (PreMZ: %.2f ColEn: %.1f ActMet: %s)"%(self.MS2_list[-1].filter_line, self.MS2_list[-1].precursor_mz, self.MS2_list[-1].collisionEnergy, self.MS2_list[-1].activationMethod if hasattr(self.MS2_list[-1], "activationMethod") else "-")

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
                if intensityCutoff<0:
                    assert scan.peak_count == scan.peak_count_tag

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
            tmp_ms.retention_time = specturm["scan time"]*60
            #if tmp_ms.peak_count > 0:
            #    tmp_ms.low_mz = float(specturm['lowest observed m/z'])
            #    tmp_ms.high_mz = float(specturm['highest observed m/z'])
            #    tmp_ms.base_peak_mz = float(specturm['base peak m/z'])
            #    tmp_ms.base_peak_intensity = float(specturm['base peak intensity'])
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


    def parse_file(self, filename_xml, intensityCutoff=-1, ignoreCharacterData=False):
        if filename_xml.lower().endswith(".mzxml"):
            return self.parseMZXMLFile(filename_xml, intensityCutoff, ignoreCharacterData)
        elif filename_xml.lower().endswith(".mzml"):
            return self.parseMzMLFile(filename_xml, intensityCutoff, ignoreCharacterData)
        else:
            return RuntimeError("Invalid file type")


    # updates the data of this MS scan
    # WARNING: unstable, horridly implemented method. SUBJECT OF CHANGE
    def updateScan(self, scan, data):
        scanID = int(scan.getAttribute("num"))

        # Do not use. Only if you want to simulate a positive and negative recorded full scan.
        if False:
            if scanID % 2:
                scan.setAttribute("polarity", "+")
            else:
                scan.setAttribute("polarity", "-")

        peaks = None
        for kid in scan.childNodes:
            if kid.nodeName == "peaks":
                peaks = kid
        assert peaks is not None and len(peaks.childNodes) == 1

        if data.has_key(scanID) and len(data[scanID][0]) > 0:
            h = data[scanID]
            assert len(h) == 2 and len(h[0]) == len(h[1])
            peaks.childNodes[0].nodeValue = base64.encodestring("".join([struct.pack(">f", h[0][i]) + struct.pack(">f", h[1][i]) for i in range(len(h[0]))])).strip().replace("\n", "")
            scan.setAttribute("peaksCount", str(len(h[0])))
            scan.setAttribute("lowMz", str(min(h[0])))
            scan.setAttribute("highMz", str(max(h[0])))
            scan.setAttribute("basePeakMz", "0")
            scan.setAttribute("basePeakIntensity", "0")
            scan.setAttribute("totIonCurrent", str(sum(h[1])))
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
    f="E:/150807_pos_213_MSMS_12C13C_PPAs_Exp36/12C13C_Wheat_Fg_conc_MSMS_CID20_Set1.mzML"
    f="C:/Users/cbueschl/Desktop/implTest/Exactive_plus/WheatEar_DON_posneg.mzML"

    t=Chromatogram()
    t.parse_file(f)

    for i in t.getFilterLines(includeMS2=False):
        print i

    print t.getPolarities()

    x=Chromatogram()
    x.parse_file(f.replace(".mzML", ".mzXML"))


    import matplotlib.pyplot as plt


    TIC, times, scanIds=t.getTIC()
    plt.plot(times, TIC)


    TIC, times, scanIds=x.getTIC()
    plt.plot(times, [t*1 for t in TIC])

    plt.show()