import sys
from PyQt4 import QtGui

class QScrollableMessageBox(QtGui.QDialog):
    def __init__(self, parent=None, text="", title="message", width=200, height=100):
        super(QScrollableMessageBox, self).__init__()
        self.setModal(True)
        self.setWindowTitle(title)

        self.setMinimumHeight(height)
        self.setMinimumWidth(width)

        l = QtGui.QGridLayout(self)

        k = QtGui.QWidget()
        l.addWidget(k)
        o = QtGui.QGridLayout(k)
        o.setContentsMargins(10, 10, 10, 10)

        self.scroll = QtGui.QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.res=QtGui.QPlainTextEdit(text)
        self.res.setFont(QtGui.QFont("Consolas"))
        self.res.setWordWrapMode(0)
        self.scroll.setWidget(self.res)
        o.addWidget(self.scroll, 2, 1)

        self.okBut=QtGui.QPushButton("OK")
        self.okBut.clicked.connect(self.hide)
        o.addWidget(self.okBut, 3,1)


if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)

    pw = QScrollableMessageBox(parent=None, text="texttexttexttexttexttexttexttexttexttexttexttexttexttexttext texttexttexttexttexttexttexttexttexttexttext\n111", title="")
    pw.show()

    sys.exit(app.exec_())








