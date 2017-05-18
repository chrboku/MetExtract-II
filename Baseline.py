
if __name__=="__main__":
    import MExtract

import rpy2
import rpy2.robjects as ro

from utils import Bunch


r = ro.r

class Baseline:

    def __init__(self):
        r("suppressWarnings(suppressMessages(library(baseline))); getBaseline<-function(eic){  return(baseline(eic, method='medianWindow', hwm=40)@baseline[1,])  };")

    def getBaseline(self, eic, times):

        eicR="c(" + ",".join(str(e) for e in eic) + ")"

        try:
            ret = r("getBaseline(t("+eicR+"));")

            ret=[float(i) for i in ret]
            return ret

        except Exception as e:
            raise Exception("Baseline could not be calculated successfully")

