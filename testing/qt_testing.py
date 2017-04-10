import sys, os
from PyQt4 import QtGui, uic

qtCreatorFile = r"..\analysis\biscotti_gui.ui"  # Enter file here.

Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class test_Gui(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self):
        # Setup calls
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        paramtable = QtGui.QGridLayout()
        for _ in range(4):
            paramtable.addWidget(QtGui.QTextBrowser())
        qmain = QtGui.QWidget()
        qmain.setLayout(paramtable)
        QtGui.QMainWindow.setCentralWidget(self, qmain)
        # ok now i have a silly list. Lets get the damn index
        paramtable.itemAt(1).widget().setText('test') #OK Now i can access it via the list thing



if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = test_Gui()
    window.show()
    sys.exit(app.exec_())
