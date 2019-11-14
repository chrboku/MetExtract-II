from utils import Bunch, CallBackMethod

from multiprocessing import Pool, Manager, Queue
from PyQt4 import QtGui
import time


from mePyGuis import ProgressWrapper

def runImapUnordered(functionToCall, parameters, processes=1, pw=None):

    # initialise multiprocessing queue
    p = Pool(processes=processes, maxtasksperchild=1) # only in python >=2.7; experimental

    if pw == "createNewPW":
        pw=ProgressWrapper.ProgressWrapper(1, showProgressBars=True, showLog=False, showIndProgress=False)
        pw.show()

        pw.getCallingFunction()("max")(len(parameters))
        pw.getCallingFunction()("value")(0)
        pw.getCallingFunction()("log")("")

    startProc=time.time()
    if pw!=None: pw.getCallingFunction()("log")("Starting evaluation..")

    # start the multiprocessing pool
    res = p.imap_unordered(functionToCall, parameters)

    # wait until all subprocesses have finished re-integrating their respective LC-HRMS data file
    loop = True
    while loop:
        completed = res._index
        if completed == len(parameters):
            loop = False
        else:
            if pw!=None: pw.getCallingFunction()("text")("%d of %d processes completed\n(%d in parallel, %.2f minutes running)"%(completed, len(parameters), processes, (time.time()-startProc)/60.))
            if pw!=None: pw.getCallingFunction()("value")(completed)

            QtGui.QApplication.processEvents();
            time.sleep(5)

    if pw!=None: pw.getCallingFunction()("log")("Finished evalution (%.2f minutes)"%((time.time()-startProc)/60.))

    ret=[]
    for re in res:
        ret.append(re)

    return ret