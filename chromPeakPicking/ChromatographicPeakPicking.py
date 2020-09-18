from MassSpecWavelet import MassSpecWavelet
from GradientPeaks import GradientPeaks
from XCMSCentwave import XCMSCentwave

class ChromPeakPicking:

    method=None

    def __init__(self, algorithm="MassSpecWavelet", *argv):
        self.algorithm=algorithm

        if self.algorithm.lower()=="MassSpecWavelet".loweR():
            return MassSpecWavelet(*argv)
        if self.algorithm.lower()=="GradientPeaks".lower():
            return GradientPeaks(*argv)
        if self.algorithm.lower()=="XCMSCentwave".lower():
            return XCMSCentwave(*argv)


    @staticmethod
    def getCurrentPeakPicking(self, *argv):
        if ChromPeakPicking.method==None:
            method=ChromPeakPicking(*argv)
        return method