import sys

from PyQt4 import QtGui
from PyQt4.QtCore import Qt

from ColoredProgressBar import ColoredProgressBar
from utils import CallBackMethod, natSort


def setLabelBackground(qlabel, colorName=None, r=255, g=0, b=0, alpha=255):
    qlabel.setAutoFillBackground(True)

    if colorName is not None:
        color = QtGui.QColor(colorName)
    else:
        color = QtGui.QColor(r, g, b)

    alpha = alpha
    values = "{r},{g},{b},{a}".format(r=color.red(), g=color.green(), b=color.blue(), a=alpha)
    qlabel.setStyleSheet(
        "QLabel { background-color: rgba(" + values + "); border-width: 1px;border-style: solid;border-color: rgb(170, 170, 170);}")


# a dialog for showing the progress of individual calculations each having individual subtasks
class ProgressWrapper(QtGui.QDialog):

    ## Disable closing of dialog with the ESC key
    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Escape:
            pass

    def __init__(self, pwCount=1, showProgressBars=True, showLog=False, showIndProgress=False, indGroups={},
                 indProgColumns=15, parent=None, closeCallback=None, skipCallBack=False):
        super(ProgressWrapper, self).__init__(parent)
        self.setModal(True)
        self.setWindowTitle("Progress Wrapper")


        l = QtGui.QGridLayout(self)

        self.hasProgressBars = False
        if showProgressBars:
            self.hasProgressBars = True
            self.texts = []
            self.ps = []

            for pw in range(pwCount):
                text = QtGui.QLabel("#%d" % pw)
                l.addWidget(text)
                self.texts.append(text)

                p = ColoredProgressBar()
                l.addWidget(p)
                self.ps.append(p)

        self.hasIndFileGroups = False
        if showIndProgress:
            self.hasIndFileGroups = True
            self.statuss = {}

            if self.hasProgressBars:
                line = QtGui.QFrame()
                line.setFrameShape(QtGui.QFrame.HLine)
                l.addWidget(line)

            p = QtGui.QScrollArea()
            p.setStyleSheet("QScrollArea { border-width: 0px;border-style: solid;border-color: rgb(170, 170, 170);}")
            p.setWidgetResizable(True)
            p.setVerticalScrollBarPolicy(2)
            p.setMaximumHeight(min(900, 35*len(indGroups))-(10 if len(indGroups)>1 else 0))
            p.setMinimumHeight(min(600, 35*len(indGroups))-(10 if len(indGroups)>1 else 0))
            p.setContentsMargins(0,0,0,0)
            l.addWidget(p)

            k = QtGui.QWidget()
            k.setContentsMargins(0,0,0,0)
            p.setWidget(k)

            o = QtGui.QGridLayout(k)
            o.setContentsMargins(0,0,0,0)
            for i, indGroup in enumerate(natSort(indGroups.keys())):
                text = QtGui.QLabel(indGroup)
                o.addWidget(text,i,0)

                w = QtGui.QWidget()
                o.addWidget(w,i,1)
                t = QtGui.QGridLayout(w)
                t.setContentsMargins(0, 0, 0, 0)

                crow = 0
                ccol = 0
                for f in indGroups[indGroup]:
                    text = QtGui.QLabel("")
                    if f not in self.statuss:
                        self.statuss[f] = [text]
                    else:
                        self.statuss[f].append(text)
                    text.setFixedSize(15, 15)
                    text.setToolTip("File: " + str(f) + "\nStatus: pending")
                    setLabelBackground(text, colorName="slategrey")
                    t.addWidget(text, crow, ccol)

                    ccol = ccol + 1

                    if ccol == indProgColumns:
                        spacerItem = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
                        t.addItem(spacerItem, crow, ccol)
                        crow = crow + 1
                        ccol = 0

                spacerItem = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
                t.addItem(spacerItem, crow, indProgColumns)

                if ccol < indProgColumns:
                    spacerItem = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
                    t.addItem(spacerItem, crow, ccol + 1, 1, indProgColumns - ccol - 1)

                if i < (len(indGroups) - 1):
                    line = QtGui.QFrame()
                    line.setFrameShape(QtGui.QFrame.HLine)
                    line.setFrameShadow(QtGui.QFrame.Sunken)
                    t.addWidget(line, crow + 1, 0, 1, indProgColumns + 1)

        self.hasLog = False
        if showLog:
            self.hasLog = True

            if self.hasProgressBars or self.hasIndFileGroups:
                line = QtGui.QFrame()
                line.setFrameShape(QtGui.QFrame.HLine)
                l.addWidget(line)

            text = QtGui.QLabel("Log")
            l.addWidget(text)

            self.log = QtGui.QPlainTextEdit()
            l.addWidget(self.log)
            self.log.setLineWrapMode(0)

        self.closeCallBack=closeCallback
        self.skipCallBack=skipCallBack

    def closeEvent(self, event):
        if not self.skipCallBack and self.closeCallBack!=None:
            if not self.closeCallBack():
                return event.ignore()



    def hide(self):
        super(ProgressWrapper, self).hide()
        super(ProgressWrapper, self).close()


    def setMax(self, v, i=0):
        if self.hasProgressBars:
            self.ps[i].setMaximum(v)

    def setValue(self, v, i=0):
        if self.hasProgressBars:
            self.ps[i].setValue(v)

    def setUntils(self, u, i=0):
        if self.hasProgressBars:
            self.ps[i].setUntils(u)

    def setText(self, t, i=0):
        if self.hasProgressBars:
            self.texts[i].setText(t)

    def setHeader(self, t, i=0):
        if self.hasProgressBars:
            self.setWindowTitle(t)

    def appendToLog(self, t, i=0):
        if self.hasLog:
            self.log.setVisible(True)
            self.log.appendPlainText(t)

    def setStatusColor(self, status, colorName=None, r=255, g=0, b=0):
        if self.hasIndFileGroups:
            for lab in self.statuss[status]:
                setLabelBackground(lab, colorName=colorName, r=r, g=g, b=b)

    def setStatusText(self, status, text):
        if self.hasIndFileGroups:
            for lab in self.statuss[status]:
                lab.setToolTip(text)

    def setCloseCallback(self, closeCallBack=None):
        self.closeCallBack=closeCallBack

    def setSkipCallBack(self, skipCallBack):
        self.skipCallBack=skipCallBack

    def setMaxu(self, v, i=0):
        self.setMax(v, i)
        QtGui.QApplication.processEvents();

    def setValueu(self, v, i=0):
        self.setValue(v, i)
        QtGui.QApplication.processEvents()

    def setUntilsu(self, u, i=0):
        self.setUntils(u, i)
        QtGui.QApplication.processEvents()

    def setTextu(self, t, i=0):
        self.setText(t, i)
        QtGui.QApplication.processEvents();

    def setHeaderu(self, t, i=0):
        self.setHeader(t, i)
        QtGui.QApplication.processEvents();

    def appendToLogu(self, t, i=0):
        self.appendToLog(t, i)
        QtGui.QApplication.processEvents();

    def setStatusColoru(self, status, colorName=None, r=255, g=0, b=0):
        self.setStatusColor(status, colorName=colorName, r=r, g=g, b=b)
        QtGui.QApplication.processEvents();

    def setStatusTextu(self, status, text):
        self.setStatusText(status, text)
        QtGui.QApplication.processEvents();


    def getMaxSetter(self, i=0):
        return lambda x: self.setMaxu(x, i)

    def getValueSetter(self, i=0):
        return lambda x: self.setValueu(x, i)

    def getUntilsSetter(self, i=0):
        return lambda x: self.setUntilsu(x, i)

    def getTextSetter(self, i=0):
        return lambda x: self.setTextu(x, i)

    def getHeaderSetter(self, i=0):
        return lambda x: self.setHeaderu(x, i)

    def getAppendToLog(self, i=0):
        return lambda x: self.appendToLogu(x, i)

    def getSetStatusColorSetter(self):
        return self.setStatusColoru

    def getSetStatusTextSetter(self):
        return self.setStatusText

    def getCallingFunctionu(self, x, i=0):
        x = x.lower()
        if x == "settext" or x == "text":
            return self.getTextSetter(i)
        elif x == "setmax" or x == "max":
            return self.getMaxSetter(i)
        elif x == "setvalue" or x == "value":
            return self.getValueSetter(i)
        elif x == "setuntils" or x == "untils":
            return self.getUntilsSetter(i)
        elif x == "setheader" or x == "header":
            return self.getHeaderSetter(i)
        elif x == "appendtolog" or x == "log":
            return self.getAppendToLog(i)
        elif x == "setstatuscolor" or x == "statuscolor":
            return self.getSetStatusColorSetter()
        elif x == "setstatustext" or x == "statustext":
            return self.getSetStatusTextSetter()

    def getCallingFunction(self, i=0):
        return lambda x, i=i: self.getCallingFunctionu(x, i)


def callBack(a="hello", b="world", qt=None):
    if QtGui.QMessageBox.question(qt, a+b, "Close", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)==QtGui.QMessageBox.Yes:
        print "closing"
        return True
    else:
        print "not closing"
        return False

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)

    # create a ProgressWrapper with 3 progress bars, a log and 3 groups with individual files in each group
    pw = ProgressWrapper(3, showProgressBars=True, showLog=False, showIndProgress=True,
                         indGroups={"blanks": ["blank_1", "blank_2", "blank_3"],
                                    "PH1": ["PH1_%d"%i for i in range(1, 28)],
                                    "TRI5": ["TRI5_1", "TRI5_2", "TRI5_3", "TRI5_4"],
                                    "test1":["Test1", "Test2"],
                                    "test2":["Test2", "Test4"]})
    pw.setCloseCallback(closeCallBack=CallBackMethod(_target=callBack, b="Joe", qt=pw).getRunMethod())
    pw.show()

    pw.getCallingFunction(0)("text")("<p align='right' style='background-color: red;'>first element</p>\n\nfirst elementfirst elementfirst elementfirst elementfirst elementfirst element")
    pw.getCallingFunction(1)("text")("second element")
    pw.getCallingFunction(2)("text")("third element")

    # set a maximum of 10 on the second progress bar and set its current value to 5 (50%)
    pw.getCallingFunction(0)("max")(1.)
    pw.setUntils({.25: "#E3E9E2", .35: "#D3D9D2", .65: "#52626D", .70: "#6C7B8B", .75: "#92A1AB", .9: "#CACDC9", 1.: "darkorange"}, 0)
    pw.getCallingFunction(0)("value")(1.)
    pw.getCallingFunction(1)("max")(10)
    pw.setUntils({.25: "#E3E9E2", .35: "#D3D9D2", .65: "#52626D", .70: "#6C7B8B", .75: "#92A1AB", .9: "#CACDC9", 1.: "darkorange"}, 1)
    pw.getCallingFunction(1)("value")(3)
    pw.getCallingFunction(2)("max")(10)
    pw.setUntils({.25: "#E3E9E2", .35: "#D3D9D2", .65: "#52626D", .70: "#6C7B8B", .75: "#92A1AB", .9: "#CACDC9", 1.: "darkorange"}, 2)
    pw.getCallingFunction(2)("value")(5)

    pw.getCallingFunction()("appendtolog")("some log entry")
    pw.getCallingFunction()("appendtolog")("another log entry")
    pw.getCallingFunction()("appendtolog")("-----------------")
    pw.getCallingFunction()("appendtolog")("")
    pw.getCallingFunction()("appendtolog")("a severe exception occurred")

    # set color of file blank_2 to green and its status to finished
    # note: this will affect two groups in total since they share this respective file
    pw.setStatusColoru("blank_2", colorName="olivedrab")
    pw.setStatusTextu("blank_2", text="File: %s\nStatus: %s" % ("blank_1", "finished"))

    # set color of file blank_1 to red and status to failed
    pw.getSetStatusColorSetter()("blank_1", colorName="firebrick")
    pw.getSetStatusTextSetter()("blank_1", text="File: %s\nStatus: %s" % ("blank_1", "failed"))

    # set color of file PH1_4 to orange and status to running
    pw.getCallingFunction()("statuscolor")("PH1_4", "Orange")
    pw.getCallingFunction()("statustext")("PH1_4", text="File: %s\nStatus: %s" % ("PH1_4", "running"))

    sys.exit(app.exec_())
    
    
    
    
    
    
    
    
    
    
    
    
    
    