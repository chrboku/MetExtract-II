import sys
from PySide6 import QtGui, QtWidgets


class QScrollableMessageBox(QtWidgets.QDialog):
    def __init__(self, parent=None, text="", title="message", width=200, height=100):
        super(QScrollableMessageBox, self).__init__()
        self.setModal(True)
        self.setWindowTitle(title)

        self.setMinimumHeight(height)
        self.setMinimumWidth(width)

        l = QtWidgets.QGridLayout(self)

        k = QtWidgets.QWidget()
        l.addWidget(k)
        o = QtWidgets.QGridLayout(k)
        o.setContentsMargins(10, 10, 10, 10)

        self.scroll = QtWidgets.QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.res = QtWidgets.QPlainTextEdit(text)
        self.res.setFont(QtGui.QFont("Consolas"))
        self.res.setWordWrapMode(QtGui.QTextOption.WrapMode.NoWrap)
        self.scroll.setWidget(self.res)
        o.addWidget(self.scroll, 2, 1)

        self.okBut = QtWidgets.QPushButton("OK")
        self.okBut.clicked.connect(self.hide)
        o.addWidget(self.okBut, 3, 1)


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)

    pw = QScrollableMessageBox(
        parent=None,
        text="texttexttexttexttexttexttexttexttexttexttexttexttexttexttext texttexttexttexttexttexttexttexttexttexttext\n111",
        title="",
    )
    pw.show()

    sys.exit(app.exec())
