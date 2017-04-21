import sys, os
from PyQt4 import QtGui, uic
from biscotti.classes import qecalc
from biscotti.functions.base import scan_folder_for_calcs, list_folder
from collections import OrderedDict as odict
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
        self.last_rootpath = r'D:\Users\Chris\Documents\SivaLab\Type II Superlattices\2016 MDA Type 2 SL\Ab-Initio'
        self.all_dft_folders = []  # lets currently assume its unsorted
        self.loaded_dft_folders = [] # List of folders already calculated and present in 'calc objects dict'. Probably not needed
        self.selected_indices = [] # These indices were selected by the user
        self.calc_objects_dict = {} # Key: Full path to calc folder, Values: List of calc objects
        headers = ['ID', 'Filename', 'Folder', 'Title', 'Calc Type',
                    'Final Energy (Ry)', 'Last electron step dE (Ry)','Last ion step dE (Ry)','Final Free Energy (Ry)',
                       'Final Free Energy (eV)',
                       'Final Free Energy (eV/atom)',
                       'Number of Atoms',
                       'Cutoff (Ry)',
                       'Total # of K-points',
                       'Job Complete?',
                       'Calc time (hr)',
                       'Start Date-Time',
                       'End Date-Time',
                       'Final Pressure(kbar)',
                       'Initial Volume (A^3)',
                       'Final Volume (A^3)',
                    'Final Pressure (kbar)'
                       ]
        self.summary_data = (headers, []) # Header and summary data (2d);
        self.summary_dict = odict()

    def load_from_dir(self):
        # Clear old variables
        self.calc_objects_dict = {}

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
        ext_dict = {}
        self.all_dft_folders = []
        for folder in os.walk(rootpath):
            exts = [os.path.splitext(file)[1] for file in os.listdir(folder[0]) if any(file.endswith(ext) for ext in ext_list)]
            if len(exts) > 0:
                self.all_dft_folders.append(folder[0])
                ext_dict[folder[0]] = exts
                logger.debug("Adding " + str(exts) + " to " + str(folder[0]))
            else:
                logger.debug("No calc files found in " + folder[0])
        logger.info("Have found " + str(len(self.all_dft_folders)) + " calc folders.")

        self.selection_table.setRowCount(len(self.all_dft_folders))
        for i, calc_folder in enumerate(self.all_dft_folders):
            # We iterate over the keys of the dictionary
            pathlist = calc_folder.split(os.sep)
            self.selection_table.setItem(i,0, QtGui.QTableWidgetItem(pathlist[-2]))
            self.selection_table.setItem(i, 1, QtGui.QTableWidgetItem(pathlist[-1]))
            exts_in_calcfolder = str(ext_dict[calc_folder]).replace("'","")
            self.selection_table.setItem(i, 2, QtGui.QTableWidgetItem(exts_in_calcfolder))
        # self.all_dft_folders =list(allcalcs_dict.keys())

        # Now we update the widgets
        if len(self.all_dft_folders) > 0:
            all_params = self.summary_data[0] # Headers
            default_params = ['Final Free Energy (eV/atom)', 'Last ion step dE (Ry)','Number of Atoms',
                              'Job Complete?', 'Calc time (hr)','Final Pressure(kbar)', 'Final Volume (A^3)', 'Cutoff (Ry)', ]
            default_indices = [all_params.index(default) for default in default_params]
            logger.debug("Default indices are " + str(default_indices))
            for i, param_index in enumerate(default_indices):
                cmbox = self.param_table.itemAt(i).itemAt(0).widget()
                value_label = self.param_table.itemAt(i).itemAt(1).widget().setText('No file loaded')
                for option in all_params:
                    cmbox.addItem(option)
                # cmbox = QtGui.QComboBox()
                cmbox.setCurrentIndex(param_index)

    def cell_clicked(self):
        # A user has selected 1 or more calculation files.
        self.selected_indices = [row.row() for row in self.selection_table.selectedIndexes()]

        # Get the list of qecalc objects
        calc_folders = [self.all_dft_folders[index] for index in self.selected_indices]
        these_calcs = []
        for calcpath in calc_folders:
            if calcpath not in self.calc_objects_dict:
                # this not yet loaded
                # 1. scan for calcs available
                logger.debug("This path not previously scanned! Now scanning " + calcpath + " for QE calculations")
                calcs_found = scan_folder_for_calcs(calcpath)
                listf = list_folder(calcpath, ['.out'])
                logger.debug(str(calcs_found) + " and " + str(listf) + " and, " + str(list_folder(calcpath)))
                logger.debug("Found " + str(len(calcs_found)) + " suitable file extensions with calculations")
                # 2. create related objects and place in the list
                calclist = []
                for calc, file in calcs_found:
                    if calc == 'pw.x input':
                        calclist.append(qecalc.QECalcIn.import_from_file(file))
                    elif calc == 'pw.x output':
                        calclist.append(qecalc.QECalcOut.import_from_file(file))
                # 3 add to calc objects dictionary
                self.calc_objects_dict[calcpath] = calclist
                logger.debug("Adding " + str(len(calclist)) + " calcs to dictionary for path " + calcpath)

            these_calcs.extend(self.calc_objects_dict[calcpath])
            logger.debug(str(len(these_calcs)) + " calcs selected")

        # Update the details table
        # Update headers
        self.update_qtable_headers(self.details_table, self.summary_data[0])
        # Update rows based on selection
        for calc in these_calcs:
            logger.debug("Type is " + str(type(calc)))
        detail_data = [list(calc.calc_overview_dict().values()) for calc in these_calcs if type(calc) is qecalc.QECalcOut]
        logger.debug("Now adding " + str(len(detail_data)) + " rows to the details table")
        self.update_qtable_entires(self.details_table, detail_data)

            # Populate based on the current tab
            # Assume overview tab
            # Populate summary table
            # Get list of parameters for the table

            # Recalculate each time, i guess

    def update_qtable_headers(self, q_table_widget: QtGui.QTableWidget, newheaders = None):
        if newheaders is None:
            newheaders = self.summary_data[0]
        q_table_widget.setColumnCount(len(newheaders))
        for i, column in enumerate(self.summary_data[0]):
            q_table_widget.setHorizontalHeaderItem(i, QtGui.QTableWidgetItem(column))

    def update_qtable_entires(self, q_table_widget: QtGui.QTableWidget, data : list):
        q_table_widget.setRowCount(len(data))
        for i, row in enumerate(data):
            for j, entry in enumerate(row):
                q_table_widget.setItem(i, j, QtGui.QTableWidgetItem(str(entry)))

    # def addmpl(self, fig):
    #     self.canvas = FigureCanvas(fig)
    #     self.mplvl.addWidget(self.canvas)
    #     self.canvas.draw()

    def prep_plot(self, plot_name, x_axis, y_axis):
        pass


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = Biscotti_Gui()
    window.show()
    sys.exit(app.exec_())