USEGRADIENTDESCENDPEAKPICKING=False

from math import pow, sqrt
from copy import copy as copycopy


import numpy
from numpy import *
import numpy as np

from scipy.ndimage import gaussian_filter1d



from operator import itemgetter



import numpy as np
from math import factorial
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    y=np.array(y)
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')



def __smooth(x, window_len=11, window='bartlett'):
    """
    (c) http://www.scipy.org/Cookbook/SignalSmooth (last accessed: 30. May 2012)
    
    smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is non of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s = numpy.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat':  #moving average
        w = numpy.ones(window_len, 'd')
    else:
        w = eval('numpy.' + window + '(window_len)')

    y = numpy.convolve(w / w.sum(), s, mode='valid')
    return y

def _smooth(x, window_len=11, window='hanning'):
    x = list(__smooth(array(x), window_len, window))
    x = x[int(window_len / 2):int((len(x) - window_len / 2 + (1 if not (window_len % 2) else 0)))]
    return x

def _smoothTriangle(data, degree, dropVals=False):
    """performs moving triangle smoothing with a variable degree."""
    """note that if dropVals is False, output length will be identical
    to input length, but with copies of data at the flanking regions"""
    triangle = numpy.array(range(degree) + [degree] + range(degree)[::-1]) + 1
    smoothed = []
    for i in range(degree, len(data) - degree * 2):
        point = data[i:i + len(triangle)] * triangle
        smoothed.append(sum(point) / sum(triangle))
    if dropVals: return smoothed
    smoothed = [smoothed[0]] * (degree + degree / 2) + smoothed
    while len(smoothed) < len(data): smoothed.append(smoothed[-1])
    return smoothed

# smooth data series
def smoothDataSeries(x, y, windowLen=2, polynom=3, window="gaussian", removeNegIntensities=True):
    if window == None or window.lower() == "none":
        return y

    window = window.lower()

    addD = int(round(windowLen / 2.))

    yi = copycopy(y)
    for i in range(addD):
        yi.insert(0, yi[0])
        yi.append(yi[-1])

    if window == "triangle":
        ys = _smoothTriangle(yi, windowLen)
        ys = ys[0:len(y)]

    elif window == "gaussian":
        ys = gaussian_filter1d(yi, windowLen)
        ys = ys[addD:(len(y) + addD)]

    elif window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        ys = _smooth(yi, window_len=windowLen, window=window)
        ys = ys[addD:(len(y) + addD)]
    elif window == "savitzkygolay":
        ys = savitzky_golay(yi, windowLen, order=polynom)
        ys = ys[addD:(len(y) + addD)]
    else:
        raise Exception("Unsupported smoothing window")

    if removeNegIntensities:
        ys=[max(0, i) for i in ys]

    return ys






# HELPER METHOD
# taken from http://www.py2exe.org/index.cgi/HowToDetermineIfRunningFromExe
import imp, os, sys
def main_is_frozen():
    return (hasattr(sys, "frozen") or # new py2exe
            hasattr(sys, "importers") # old py2exe
            or imp.is_frozen("__main__")) # tools/freeze

# method to determine location of python script / executable (not working directory)
def get_main_dir():
    if main_is_frozen():
        return os.path.dirname(sys.executable).replace("\\","/")
    #return os.path.dirname(sys.argv[0]).replace("\\","/")
    return os.path.abspath(__file__).replace("\\", "/")[:os.path.abspath(__file__).replace("\\", "/").rfind("/")]

import threading
class FuncThread(threading.Thread):
    def __init__(self, _target, *args, **kwds):
        threading.Thread.__init__(self)

        self._target = _target
        self._args = args
        self._kwds = kwds

        self._stop = threading.Event()

    def addArg(self, val):
        self._args.append(val)

    def addKwd(self, param, val):
        self._kwds[param]=val

    def stop(self):
        self._stop.set()

    def stopped(self):
        return self._stop.isSet()

    def run(self):
        self._target(*self._args, **self._kwds)

if __name__=="__main__" and False:
    # Example usage
    def someOtherFunc(data, key):
        print "someOtherFunc was called : data=%s; key=%s" % (str(data), str(key))

    t1 = FuncThread(_target=someOtherFunc, data=[1,2], key=6)
    t1.start()
    t1.join()





from multiprocessing import Process, Queue

def process_func(q, _target, *args, **kwds):
    try:
        ret = _target(*args, **kwds)
    except Exception:
        from sys import exc_info
        from traceback import format_tb
        ex_type, ex_value, tb = exc_info()
        error = ex_type, ex_value, ''.join(format_tb(tb))
        ret = None
    else:
        error = None

    q.put((ret, error))

class FuncProcess():
    def __init__(self, _target, *args, **kwds):
        self._target = _target
        self._args = args
        self._kwds = kwds


        self.q=None

    def getQueue(self):
        if self.q is None:
            self.q=Queue()
        return self.q

    def addArg(self, val):
        self._args.append(val)

    def addKwd(self, param, val):
        self._kwds[param]=val

    def start(self):
        self.p=Process(target=process_func, args=[self.q, self._target] + list(self._args), kwargs=self._kwds)
        self.p.start()

    def isAlive(self):
        if self.p is None:
            return False
        return self.p.is_alive()

    def stopped(self):
        return not self.isAlive()

    def join(self):
        if self.p is None:
            return
        self.p.join()

    def terminate(self):
        self.p.terminate()

    def getReturn(self):
        self.p.join()
        ret, error = self.q.get()

        if error:
            ex_type, ex_value, tb_str = error
            message = '%s (in subprocess)\n%s' % (ex_value.message, tb_str)
            raise ex_type(message)

        return ret


def someOtherFunc(data, key):
    print "someOtherFunc was called in separate process (%s): data=%s; key=%s" % (os.getpid(), str(data), str(key))
    return "someReturn"

if __name__=="__main__" and False:
    import os
    # Example usage
    t1 = FuncProcess(_target=someOtherFunc, data=[1,2], key=7)
    t1.start()
    print "main process", os.getpid(), "t1 is alive:", t1.isAlive()

    print "t1 is alive:", t1.isAlive()
    t1.join()
    print "t1 is alive:", t1.isAlive()
    print t1.getReturn()
    t1.terminate()









class CallBackMethod():
    def __init__(self, _target, *args, **kwds):
        self._target = _target
        self._args = args
        self._kwds = kwds

    def run(self, *args, **kwds):
        ar=[]
        ar.extend(self._args)
        ar.extend(args)
        kw={}
        kw.update(self._kwds)
        kw.update(kwds)
        return self._target(*ar, **kw)

    def getRunMethod(self):
        return self.run

if False:
    def callbackMethod(a="hello", b="world"):
        print a, b

    cm=CallBackMethod(_target=callbackMethod, b="Joe").getRunMethod()
    cm()






# GENERIC CLASS that is composed of member variables supplied during object construction.
# It implements a generic string representation
class Bunch:
    def __init__(self, _addFollowing=None, _addFollowingDefaultValue=None, **kwds):
        if _addFollowing is not None:
            if isinstance(_addFollowing, list):
                for att in _addFollowing:
                    setattr(self, att, _addFollowingDefaultValue)
            elif isinstance(_addFollowing, dict):
                for att, val in _addFollowing.items():
                    setattr(self, att, val)
            else:
                raise Exception("Unknown Parameter Type: Only list and dict are supported!")
        self.__dict__.update(kwds)

    def __str__(self):
        return "(%s)"%(",".join(["%s:%s"%(key, str(val)) for key, val in self.__dict__.items()]))

    def __repr__(self):
        return self.__str__()

    def hasMember(self, member):
        return self.__dict__.has_key(member)
    def hasVar(self, var):
        return self.hasMember(var)

    def dumpAsJSon(self):
        import json
        return json.dumps(self.__dict__)
    @staticmethod
    def loadFromJSon(jsonString):
        import json
        d=json.loads(jsonString)
        b=Bunch()
        for k, v in d.items():
            b.__dict__[k]=v
        return b













# class for storing a detected chromatographic peak pair.
# Thus, one ChromPeakPair instance stores two chromatographic peaks of a native and a labelled metabolite ion.
# Member variables are only set, if provides as parameters
class ChromPeakPair:
    def __init__(self, id=-1, fGroupID=-1, eicID=-1, massSpectrumID=-1, assignedName=-1, tracer=-1, tracerName="",
                 mz=-1, lmz=-1, xCount=-1, loading=-1, ionMode="", NPeakCenter=-1, NPeakCenterMin=-1, NPeakScale=-1,
                 NSNR=-1, NPeakArea=-1, NPeakAbundance=-1, LPeakCenter=-1, LPeakCenterMin=-1, LPeakScale=-1, LSNR=-1,
                 LPeakArea=-1, LPeakAbundance=-1, heteroIsotoplogues={}, assignedMZs=[], comments=[], artificialEICLShift=0, **args):
        argsUsed = 0

        self.id = id
        self.fGroupID = fGroupID
        self.eicID = eicID
        self.massSpectrumID = massSpectrumID
        self.assignedName = assignedName
        self.tracer = tracer
        self.tracerName = tracerName

        self.mz = mz
        self.lmz = lmz
        self.xCount = xCount
        self.loading = loading
        self.ionMode = ionMode

        self.NPeakCenter = NPeakCenter
        self.NPeakCenterMin = NPeakCenterMin
        self.NPeakScale = NPeakScale
        self.NSNR = NSNR
        self.NPeakArea = NPeakArea
        self.NPeakAbundance = NPeakAbundance

        self.LPeakCenter = LPeakCenter
        self.LPeakCenterMin = LPeakCenterMin
        self.LPeakScale = LPeakScale
        self.LSNR = LSNR
        self.LPeakArea = LPeakArea
        self.LPeakAbundance = LPeakAbundance

        self.heteroIsotopologues = heteroIsotoplogues or {}
        self.assignedMZs = assignedMZs or []

        self.comments = comments

        self.artificialEICLShift = artificialEICLShift

        if args.has_key("tmz"):
            self.tmz = args["tmz"]
            argsUsed += 1

        if args.has_key("peaksCorr"):
            self.peaksCorr = args["peaksCorr"]
            argsUsed += 1
        if args.has_key("silRatios"):
            self.silRatios = args["silRatios"]
            argsUsed += 1
        if args.has_key("peaksRatio"):
            self.peaksRatio = args["peaksRatio"]
            argsUsed += 1

        if args.has_key("NXIC"):
            self.NXIC = args["NXIC"]
            argsUsed += 1
        if args.has_key("LXIC"):
            self.LXIC = args["LXIC"]
            argsUsed += 1
        if args.has_key("NXICSmoothed"):
            self.NXICSmoothed = args["NXICSmoothed"]
            argsUsed += 1
        if args.has_key("LXICSmoothed"):
            self.LXICSmoothed = args["LXICSmoothed"]
            argsUsed += 1
        if args.has_key("times"):
            self.times = args["times"]
            argsUsed += 1

        if args.has_key("fDesc"):
            self.fDesc = args["fDesc"]
            argsUsed += 1
        if args.has_key("adducts"):
            self.adducts = args["adducts"]
            argsUsed += 1
        if args.has_key("heteroAtomsFeaturePairs"):
            self.heteroAtomsFeaturePairs = args["heteroAtomsFeaturePairs"]
            argsUsed += 1
        if args.has_key("Ms"):
            self.Ms = args["Ms"]
            argsUsed +=1
        if args.has_key("heteroAtoms"):
            self.heteroAtoms = args["heteroAtoms"]
            argsUsed += 1

        if args.has_key("NBorderLeft"):
            self.NBorderLeft = args["NBorderLeft"]
            argsUsed += 1
        if args.has_key("NBorderRight"):
            self.NBorderRight = args["NBorderRight"]
            argsUsed += 1
        if args.has_key("LBorderLeft"):
            self.LBorderLeft = args["LBorderLeft"]
            argsUsed += 1
        if args.has_key("LBorderRight"):
            self.LBorderRight = args["LBorderRight"]
            argsUsed += 1

        if args.has_key("isotopeRatios"):
            self.isotopeRatios = args["isotopeRatios"]
            argsUsed += 1
        if args.has_key("mzDiffErrors"):
            self.mzDiffErrors = args["mzDiffErrors"]
            argsUsed += 1

        if args.has_key("correlationsToOthers"):
            self.correlationsToOthers = args["correlationsToOthers"]
            argsUsed += 1

        assert argsUsed == len(args), "Not all agruments used %d %d" % (argsUsed, len(args))


    def __str__(self):
        return "chrompeak: %f, %d (%f min), %d" % (self.mz, self.NPeakCenter, self.NPeakCenterMin / 60., self.id)


#Faktorial of x
def fak(x):
    if x < 2:
        return 1.
    return x * fak(x - 1)


#binomial coefficient n over k
def binom(n, k):
    return fak(n) / (fak(k) * fak(n - k))


#relative ratio of intensity for peak with a total of x 12-carbon atoms having s 13-carbon substitutions
#with p being the purity for 12-carbon and 1-p being the purity for 13-carbon
def getRatio(p, x, s):
    return (p ** (x - s)) * ((1. - p) ** s) * binom(x, s)


#normalized ratio to s=0:=1
def getNormRatio(p, x, s):
    return getRatio(p, x, s) / (p ** x)


#Get most intense isotopologue for purity p and Cx
def getAtomAdd(p, x):
    m = 5
    minIndex, minValue = max(enumerate([getNormRatio(p, x, j) for j in range(m)]), key=itemgetter(1))
    return minIndex


# maps and eic to a certain reference time array.
# Used to reduce the number of scans in EICs before chromatographic alignment
def mapArrayToRefTimes(data, times, refTimes):
    assert len(data) == len(times), "Data arrays have different length"

    ref = [0]
    refTimes.insert(0, 0)
    data.insert(0, 0)
    times.insert(0, 0)
    dlen = len(data)

    j = 1
    for i in range(1, len(refTimes)):
        refTime = refTimes[i]

        while j < dlen and times[j] < refTime:
            j = j + 1

        if j == dlen:
            d = data[j - 1]
        elif refTime == times[j - 1]:
            d = data[j - 1]
        elif refTime == times[j]:
            d = data[j]
        elif j < dlen:
            d = data[j - 1] + (data[j] - data[j - 1]) * (refTime - times[j - 1]) / (times[j] - times[j - 1])

        ref.append(d)

    assert len(ref) == len(refTimes), "Error: Ref arrays have different length"

    data.pop(0)
    times.pop(0)
    refTimes.pop(0)
    ref.pop(0)
    return ref


def mean(x):
    if len(x) == 0:
        return None

    s = sum(x)
    return s / (len(x) * 1.)

def weightedMean(x, w):
    assert len(x) == len(w)

    if len(x)==0:
        return None

    s = 1. * sum(w)
    w = [i / s for i in w]

    xw = [x[i] * w[i] for i in range(len(x))]

    return sum(xw)


def sd(x):
    if len(x) <= 1:
        return None

    m = mean(x)
    return sqrt(sum([pow(xi - m, 2) for xi in x]) / len(x))


def weightedSd(x, w):
    assert len(x) == len(w)

    if len(x) <= 1:
        return None

    s = 1. * max(w)
    w = [i / s for i in w]

    m=weightedMean(x, w)
    return sqrt(sum([w[i]*pow(x[i] - m, 2) for i in range(len(x))]) / len(x))




def corr(a, b):
    assert (len(a) == len(b))

    n = len(a)

    if n < 2:
        return -1

    a = [float(i) for i in a]
    b = [float(i) for i in b]

    meana = sum(a) / n
    meanb = sum(b) / n
    #if meanb == 1 or meana == 1:
    #    return -1

    try:
        rxy = (sum([(a[i] - meana) * (b[i] - meanb) for i in range(n)])) / (
            sqrt(sum([pow(a[i] - meana, 2) for i in range(n)])) * sqrt(sum([pow(b[i] - meanb, 2) for i in range(n)])))
        return rxy
    except:
        return -1

def getLastTimeBefore(times, refTime):
    return min([(i, abs(s-refTime)) for i, s in zip(range(len(times)), times)], key=lambda x:x[1])[0]

# returns non-connected sub-graphs in a graph
def getSubGraphs(nodes):
    used = []
    groups = []

    for k in nodes.keys():
        if k not in used:
            current = []
            group = []
            current.append(k)
            while len(current) > 0:
                cur = current.pop()
                if cur not in used:
                    group.append(cur)
                    used.append(cur)
                    current.extend(nodes[cur])
            groups.append(group)
    return groups

def getSubGraphsFromDictDict(nodes):
    used = []
    groups = []

    for k in nodes.keys():
        if k not in used:
            current = []
            group = []
            current.append(k)
            while len(current) > 0:
                cur = current.pop()
                if cur not in used:
                    group.append(cur)
                    used.append(cur)
                    current.extend(nodes[cur].keys())
            groups.append(group)
    return groups


#taken from http://www.codinghorror.com/blog/2007/12/sorting-for-humans-natural-sort-order.html
#use for a=["1", "2", "10", "11", "3"]
#natSort(a)
#for a=[("1", 1), ("2", 2), ("10", 3), ("11", 4), ("3", 5)]
#natSort(a, key=itemgetter(0))
import re
def natSort(l, key=lambda ent: ent):
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda ent, key=key: [convert(c) for c in re.split('([0-9]+)', str(key(ent)))]
    l.sort(key=alphanum_key)
    return l


# HELPER METHOD used to calculate all possible (tupleSize-)combinations of elements
def _getCombinations(elements, tupleSize=2, start=0):
    if tupleSize > 1:
        j = []

        s = start
        while s < (len(elements) - tupleSize + 1):
            h = _getCombinations(elements, tupleSize - 1, start=s + 1)
            for i in h:
                a = [elements[s]]
                a.extend(i)
                j.append(a)
            s += 1

        return j

    elif tupleSize == 1:
        return [[x] for x in elements[start:(len(elements))]]

# method for calculating all possible (startTupleSize-)combinations of elements
def getXCombinations(elements, startTupleSize=2, endTupleSize=-1):
    if endTupleSize == -1:
        endTupleSize = len(elements)

    j = []
    for r in range(startTupleSize, endTupleSize + 1):
        j.extend(_getCombinations(elements, r))
    return j














# HELPER CLASS for storing a defined group
class SampleGroup:
    def __init__(self, name, files, minFound, omitFeatures, useForMetaboliteGrouping, removeAsFalsePositive, color):
        self.name = name
        self.files = files
        self.minFound = minFound
        self.omitFeatures = omitFeatures
        self.useForMetaboliteGrouping = useForMetaboliteGrouping
        self.removeAsFalsePositive = removeAsFalsePositive
        self.color = color

















def is_int(s):
    try:
        int(s)
    except ValueError:
        return False

    return True
def is_float(s):
    try:
        float(s)
    except ValueError:
        return False

    return True
def merge_dicts(*dict_args):
    '''
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    '''
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result







########################################################################################################################
########################################################################################################################
########################################################################################################################
##############                                                                                            ##############
##############                      HELPER METHODS for FILE IO, SQLITE and PRINTING                       ##############
##############                                                                                            ##############
########################################################################################################################
########################################################################################################################
########################################################################################################################



def readTSVFileAsBunch(file, delim="\t", parseToNumbers=True, useColumns=None, omitFirstNRows=0, renameRows=None):
    if renameRows is None:
        renameRows={}

    with open(file, "rb") as fin:
        curRowi=0
        headers={}

        bunchs=[]

        for line in fin:
            if omitFirstNRows>0:
                continue

            hs=line.split(delim)
            hs=[h.strip() for h in hs]
            if curRowi==0:
                for j, h in enumerate(hs):
                    if h in renameRows.keys():
                        headers[j]=renameRows[h]
                    else:
                        headers[j]=h
            else:
                b=Bunch(_addFollowing=headers.values(), _addFollowingDefaultValue="")
                for j, h in enumerate(hs):

                    if useColumns is None or headers[j] in useColumns:
                        setattr(b, headers[j], h)

                bunchs.append(b)
            curRowi+=1

    headTypes={}

    if parseToNumbers:
        for header in headers.values():
            headTypes[header]="str"

            minValTypeFound="int"
            for b in bunchs:
                attr=getattr(b, header)
                if minValTypeFound is "int" and is_int(attr):
                    minValTypeFound="int"

                elif (minValTypeFound is "float" or (minValTypeFound is "int" and not is_int(attr))) and is_float(attr):
                    minValTypeFound="float"

                else:
                    minValTypeFound="str"

            if minValTypeFound is not "str":
                for b in bunchs:
                    val=eval("%s('%s')"%(minValTypeFound, getattr(b, header)))
                    setattr(b, header, val)
                    headTypes[header]=minValTypeFound
    else:
        for header in headers.values():
            headTypes[header]="str"

    return headTypes, bunchs

from xlrd import open_workbook
def readXLSXFileAsBunch(file, sheetName="", useColumns=None, omitFirstNRows=0):
    if sheetName == "":
        bn = "Sheet1"
        rf = bn.rfind(".")
        if rf > -1:
            sheetName = bn[:rf]
        else:
            sheetName = bn
    assert omitFirstNRows >= 0

    ret = None
    rb = open_workbook(file)

    ind = 0
    sheetInd = -1
    for s in rb.sheet_names():
        if s == sheetName:
            sheetInd = ind
        ind = ind + 1
    assert sheetInd != -1, "Sheet not found"
    sh = rb.sheet_by_name(sheetName)

    i = 0
    remainingCols = True
    headers = []
    colIDs = []
    while remainingCols:
        col = i
        row = omitFirstNRows
        try:
            if sh.cell(row, col).value is None or str(sh.cell(row, col).value) == "":
                remainingCols = False
                continue
        except:
            remainingCols = False
            continue

        header = str(sh.cell(row, col).value)
        if useColumns is None or header in useColumns:
            headers.append(header)
            colIDs.append(i)
        i = i + 1

    remainingRows = True
    rows = []
    r = 1 + omitFirstNRows
    allEmpty = False
    while not allEmpty:
        rowValues = []
        allEmpty = True
        for j in colIDs:
            col = j
            row = r
            try:
                if sh.cell(row, col).value is None or str(sh.cell(row, col).value) == "":
                    rowValues.append("")
                    j = j + 1
                    continue
            except:
                rowValues.append("")
                j = j + 1
                continue
            rowValues.append(str(sh.cell(row, col).value))
            allEmpty = False
        if not allEmpty:
            rows.append(rowValues)

        r = r + 1

    dat=[]
    for row in rows:
        b=Bunch()
        for j, head in enumerate(headers):
            setattr(b, head, row[j])
        dat.append(b)

    return headers, dat

def writeBunchsAsTSVFile(bunchs, outFile, delim="\t", writeAttributes=None, useAllAttributes=False, renameRows={}):
    if writeAttributes is None:
        writeAttributes=[]

        attributes={}
        for b in bunchs:
            for att in b.__dict__.keys():
                if att not in attributes.keys():
                    attributes[att]=0
                attributes[att]+=1

        for att in sorted(attributes.keys()):
            if attributes[att]==len(bunchs) or useAllAttributes:
                writeAttributes.append(att)

    with open(outFile, "wb") as fo:
        fo.write(delim.join([att if att not in renameRows.keys() else renameRows[att] for att in writeAttributes]))
        fo.write("\n")
        for b in bunchs:
            fo.write(delim.join([str(getattr(b, att)) if hasattr(b, att) else "" for att in writeAttributes]))
            fo.write("\n")




def SQLSelectAsObject(curs, selectStatement, returnColumns=False, newObject=Bunch):
    curs.execute(selectStatement)
    names = list(map(lambda x: x[0], curs.description))
    for row in curs:
        b=newObject()
        for i in range(len(names)):
            setattr(b, names[i], row[i])
        if returnColumns:
            yield names, b
        else:
            yield b

def SQLInsert(curs, tableName, **kwargs):
    keys=kwargs.keys()
    vals=[]

    for key in keys:
        val=kwargs[key]

        if isinstance(val, int) or isinstance(val, float):
            vals.append(val)
        else:
            vals.append("%s"%val)

    assert len(keys)==len(vals)

    parts=["INSERT INTO"]
    parts.append(tableName)
    parts.append("(")

    needsComa=False
    for key in keys:
        if needsComa:
            parts.append(",")
        needsComa=True
        parts.append(key)

    parts.append(")")

    parts.append("VALUES (")

    needsComa=False
    for val in vals:
        if needsComa:
            parts.append(",")
        needsComa=True
        parts.append("?")

    parts.append(")")

    sqlCommand=" ".join(parts)
    curs.execute(sqlCommand, vals)

def writeObjectsAsSQLInsert(curs, toTable, objs, writeFields, fieldsMappingToColumns=None):
    if fieldsMappingToColumns is None:
        fieldsMappingToColumns={}
    rowIDs=[]
    for obj in objs:
        kwargs={}
        for field in writeFields:
            param=field
            if field in fieldsMappingToColumns.keys():
                param=fieldsMappingToColumns[field]
            kwargs[param]=getattr(obj, field)
        SQLInsert(curs, toTable, **kwargs)
        rowIDs.append(curs.lastrowid)
    return rowIDs
def writeObjectAsSQLInsert(curs, toTable, obj, writeFields, fieldsMappingToColumns=None):
    if fieldsMappingToColumns is None:
        fieldsMappingToColumns
    return writeObjectsAsSQLInsert(curs, toTable, [obj], writeFields, fieldsMappingToColumns)[0]


def createTableFromBunch(tableName, bunchObject, cursor, useAttrs=None, primaryKeys=None, autoIncrements=None, ifNotExists=False):
    if useAttrs is None:
        useAttrs=bunchObject.__dict__.keys()
    if primaryKeys is None:
        primaryKeys=[]
    if autoIncrements is None:
        autoIncrements=[]

    sqlCommand=["CREATE TABLE "]
    if ifNotExists:
        sqlCommand.append("IF NOT EXISTS ")
    sqlCommand.extend([tableName, " ("])

    needsComma=False
    for attr in useAttrs:
        typ=getattr(bunchObject, attr)
        if (type(typ) is type and typ is int) or (type(typ) is int):
            if needsComma:
                sqlCommand.append(", ")
            sqlCommand.extend([attr, " ", "INTEGER"])
            if attr in primaryKeys:
                sqlCommand.extend([" ", "PRIMARY KEY"])
            if attr in autoIncrements:
                sqlCommand.extend([" ", "AUTOINCREMENT"])
            needsComma=True
        elif (type(typ) is type and typ is float) or (type(typ) is float):
            if needsComma:
                sqlCommand.append(", ")
            sqlCommand.extend([attr, " ", "FLOAT"])
            if attr in primaryKeys:
                sqlCommand.extend([" ", "PRIMARY KEY"])
            if attr in autoIncrements:
                sqlCommand.extend([" ", "AUTOINCREMENT"])
            needsComma=True
        elif (type(typ) is type and typ is str) or (type(typ) is str) or True:
            if needsComma:
                sqlCommand.append(", ")
            sqlCommand.extend([attr, " ", "TEXT"])
            if attr in primaryKeys:
                sqlCommand.extend([" ", "PRIMARY KEY"])
            if attr in autoIncrements:
                sqlCommand.extend([" ", "AUTOINCREMENT"])
            needsComma=True

    sqlCommand.append(")")

    cursor.execute("".join(sqlCommand))



def printAsTable(heads, data, printInExcelFormat=False, excelSep="\t"):
    widths = [0 for i in range(len(heads))]
    for i in range(len(heads)):
        widths[i] = len(heads[i])

    for i in range(len(data)):
        if len(data[i]) != len(heads):
            raise Exception("data not consistent")
        for j in range(len(data[i])):
            widths[j] = max(widths[j], len(data[i][j]))

    if printInExcelFormat:
        row_format = excelSep.join("{%d}" % (i) for i in range(len(heads)))
    else:
        row_format = "  ".join("{%d:>%d}" % (i, widths[i]) for i in range(len(heads)))
    print row_format.format(*heads)
    if not printInExcelFormat:
        hrow = row_format.replace(" ", "-").format(*["".join(["-" for j in range(widths[i])]) for i in range(len(widths))])
        print hrow

    j = 0
    for row in data:
        print row_format.format(*row)
        j += 1
        if j == 15 and not printInExcelFormat:
            print hrow
            j = 0

def printObjectsAsTable(objs, attrs, printInExcelFormat=False, excelSep="\t"):
    widths = [0 for i in range(len(attrs))]
    for i in range(len(attrs)):
        widths[i] = len(attrs[i])

    for obj in objs:
        for j, attr in enumerate(attrs):
            if not hasattr(obj, attr):
                raise Exception("data not consistent")
            widths[j] = max(widths[j], len(str(getattr(obj, attr))))

    if printInExcelFormat:
        row_format = excelSep.join("{%d}" % (i) for i in range(len(attrs)))
    else:
        row_format = "  ".join("{%d:>%d}" % (i, widths[i]) for i in range(len(attrs)))
    print row_format.format(*attrs)
    if not printInExcelFormat:
        hrow = row_format.replace(" ", "-").format(*["".join(["-" for j in range(widths[i])]) for i in range(len(widths))])
        print hrow

    j = 0
    for obj in objs:
        row=[]
        for attr in attrs:
            row.append(getattr(obj, attr))
        print row_format.format(*row)
        j += 1
        if j == 15 and not printInExcelFormat:
            print hrow
            j = 0

def printArraysAsTable(objs, positions, printInExcelFormat=False, excelSep="\t"):
    widths = [0 for i in range(len(positions))]
    for i in range(len(positions)):
        widths[i] = len(str(positions[i]))

    for obj in objs:
        for j, pos in enumerate(positions):
            if len(obj)<=pos:
                raise Exception("data not consistent")
            widths[j] = max(widths[j], len(str(obj[pos])))

    if printInExcelFormat:
        row_format = excelSep.join("{%d}" % (i) for i in range(len(positions)))
    else:
        row_format = "  ".join("{%d:>%d}" % (i, widths[i]) for i in range(len(positions)))
    print row_format.format(*positions)
    if not printInExcelFormat:
        hrow = row_format.replace(" ", "-").format(*["".join(["-" for j in range(widths[i])]) for i in range(len(widths))])
        print hrow

    j = 0
    for obj in objs:
        row=[]
        for pos in positions:
            row.append(obj[pos])
        print row_format.format(*row)
        j += 1
        if j == 15 and not printInExcelFormat:
            print hrow
            j = 0

def printDictDictAsTable(rowDict, colHeads=[], rowsOrd=None, printInExcelFormat=False, excelSep="\t"):
    if rowsOrd is None:
        rowsOrd=rowDict.keys()

    colHeads.insert(0, "")
    for row in rowDict.keys():
        if len(colHeads[0])<len(row):
            colHeads[0]=row

    widths=[len(str(colHead)) for colHead in colHeads]

    for row in rowsOrd:
        cells=rowDict[row]
        for i in range(1, len(colHeads)):
            widths[i]=max(widths[i], len(str(cells[colHeads[i]])))

    colHeads[0]=""

    if printInExcelFormat:
        row_format = excelSep.join("{%d}" % (i) for i in range(len(colHeads)))
    else:
        row_format = "  ".join("{%d:>%d}" % (i, widths[i]) for i in range(len(colHeads)))
    print row_format.format(*colHeads)
    if not printInExcelFormat:
        hrow = row_format.replace(" ", "-").format(*["".join(["-" for j in range(widths[i])]) for i in range(len(widths))])
        print hrow


    j = 0
    for rowName in rowsOrd:
        row=[]
        row.append(rowName)
        for i in range(1,len(colHeads)):
            row.append(rowDict[rowName][colHeads[i]])
        print row_format.format(*row)
        j += 1
        if j == 4 and not printInExcelFormat:
            print hrow
            j = 0



if __name__=="__main__" and True:

    import sqlite3
    conn=sqlite3.connect(":memory:")
    curs=conn.cursor()

    createTableFromBunch(tableName="foo", bunchObject=Bunch(id=int, name=str, address=str, sex=int), cursor=curs,
                         primaryKeys=["id"], autoIncrements=["id"])

    b=[Bunch(name="Peter", address="New York", sex=1), Bunch(name="Anna", address="Sydney", sex=2), Bunch(name="Norbert", address="Linz", sex=1)]
    writeObjectsAsSQLInsert(curs, "foo", b, writeFields=["name", "address", "sex"])
    conn.commit()

    objs=[obj for obj in SQLSelectAsObject(curs, "select id, name, address, sex from foo")]
    printObjectsAsTable(objs, attrs=["id", "name", "address", "sex"])

    curs.close()
    conn.close()























