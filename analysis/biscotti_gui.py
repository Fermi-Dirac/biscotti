import sys, os
from PyQt4 import QtGui, uic
from biscotti.classes import qecalc
import logging
# Logging level by default
logger = logging.getLogger(__name__)
loglevel = logging.DEBUG
logger.setLevel(loglevel)
console_handler = logging.StreamHandler()
console_handler.setLevel(loglevel)
formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)
# End logging.

# from matplotlib.figure import Figure
# from matplotlib.backends.backend_qt4agg import (
#     FigureCanvasQTAgg as FigureCanvas,
#     NavigationToolbar2QT as NavigationToolbar)

qtCreatorFile = r"biscotti_gui.ui"  # Enter file here.

Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)
# Widgets:
# selection_table    shows the list of things to pick from
# details_table     highlights details of your selection
# load_dir_button     Loads data into the UI
# overview_tab      The first tab opened by default
## overview_grid
### ov_text
### ov_table
### ov_plot1_layout
#### Xaxis_selector_2
#### Yaxis_selector_2
#### plot_label_2
#### plot_selector_2
#### plot_window_2

class Biscotti_Gui(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self):
        # Setup calls
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)

        # Events to listen for
        self.load_dir_button.clicked.connect(self.load_from_dir)
        self.selection_table.cellClicked.connect(self.cell_clicked)

        # Main data structures
        self.last_rootpath = r'D:\Users\Chris\Documents\SivaLab\Type II Superlattices\2016 MDA Type 2 SL\Ab-Initio\Bulk'
        self.all_dft_files = []  # lets currently assume its unsorted
        self.loaded_dft_files = []
        self.selected_indices = []

    def load_from_dir(self):
        # Popup and ask user where to search
        rootpath = QtGui.QFileDialog.getExistingDirectory(caption='Select calculation input and output root folder',
                                                          directory=self.last_rootpath)
        self.last_rootpath = rootpath

        # Get file extensions from checkbox list
        ext_list = []
        for index in range(self.file_ext_boxes.count()):
            item = self.file_ext_boxes.item(index)
            if item.checkState():
                ext_list.append(item.text()[1:])

        # Populate based on the new rootpath
        alloutpaths = []
        for folder in os.walk(rootpath):
            alloutpaths.extend([folder[0] + os.sep + file for file in os.listdir(folder[0]) if
                                any(file.endswith(ext) for ext in ext_list)])
        self.selection_table.setRowCount(len(alloutpaths))
        for i, path in enumerate(alloutpaths):
            folder, file = os.path.split(path)
            folder_name = os.path.split(folder)[1]
            ext = os.path.splitext(path)[1]
            self.selection_table.setItem(i, 0, QtGui.QTableWidgetItem(folder_name))
            self.selection_table.setItem(i, 1, QtGui.QTableWidgetItem(file))
            self.selection_table.setItem(i, 2, QtGui.QTableWidgetItem(ext))
        self.all_dft_files = list(alloutpaths)

        # Now we update the widgets
        if len(self.all_dft_files) > 0:
            all_params = list(qecalc.QECalcOut(outpath=self.all_dft_files[0]).calc_overview_dict().keys())
            default_params = ['Final Free Energy (eV/atom)','Number of Atoms', 'Job Complete?', 'Calc time (hr)'
                              ,'Final Pressure(kbar)']
            default_indices = [all_params.index(default) for default in default_params]
            for i, param in enumerate(default_indices):
                cmbox = self.param_table.itemAt(i).itemAt(0).widget()
                for option in all_params:
                    cmbox.addItem(option)
                # cmbox = QtGui.QComboBox()
                cmbox.setCurrentIndex(i)

    def cell_clicked(self):
        # A user has selected 1 or more calculation files.
        self.selected_indices = [row.row() for row in self.selection_table.selectedIndexes()]
        # Update the details table

        self.details_table.setRowCount(len(self.selected_indices))
        if __name__ == '__main__':
            for i, row in enumerate(self.selected_indices):
                self.details_table.setItem(i, 0, QtGui.QTableWidgetItem(self.all_dft_files[row]))

            # Populate based on the current tab
            # Assume overview tab
            # Populate summary table
            # Get list of parameters for the table

            # Recalculate each time, i guess


    def addmpl(self, fig):
        self.canvas = FigureCanvas(fig)
        self.mplvl.addWidget(self.canvas)
        self.canvas.draw()

    def prep_plot(self, plot_name, x_axis, y_axis):
        pass


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = Biscotti_Gui()
    window.show()
    sys.exit(app.exec_())