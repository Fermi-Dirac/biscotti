import sys, os
from PyQt4 import QtGui, uic
from PyQt4.QtCore import Qt
from biscotti.classes import qecalc
import pyqtgraph as pg
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

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

qtCreatorFile = r"pw.x_analysis.ui"  # Enter file here.

Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)
# Widgets:
# filename_label    Qlabel
# load_button       QPushButton
# param_table       QTabelWidget
# plot_1            PlotWidget
# unit_combo_1      QComboBox
# free_energy_1     QCheckBox
# delta_1           QCheckbox
# plot_2            PlotWidget
# etc. for plot 2

class PWX_Analysis(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self):
        # Setup calls
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)

        # Events to listen for
        self.load_button.clicked.connect(self.load_calc)
        self.unit_combo_1.currentIndexChanged.connect(self.make_ion_plot)
        self.free_energy_1.clicked.connect(self.make_ion_plot)
        self.delta_1.clicked.connect(self.make_ion_plot)
        self.unit_combo_2.currentIndexChanged.connect(self.make_e_plot)
        self.free_energy_2.clicked.connect(self.make_e_plot)
        self.delta_2.clicked.connect(self.make_e_plot)

        # Main data structures
        self.last_rootpath = None
        self.labelstyle = {'font-size' : '16pt'}

        # init plots
        self.setup_plots()

    def load_calc(self):
        # Load the path
        calcpath = QtGui.QFileDialog.getOpenFileName(caption='Select calculation output file',
                                               directory=self.last_rootpath)
        self.last_rootpath = os.path.dirname(calcpath)
        # run import from file
        logger.debug("Now loading " + calcpath)
        self.setCursor(Qt.WaitCursor)
        self.calc = qecalc.QECalcOut.import_from_file(calcpath)
        logger.debug("Loaded " + self.calc.name)
        self.unsetCursor()

        # Set the label
        maxlabel = 50
        if len(calcpath) > maxlabel :
            label = calcpath[-maxlabel:]
        else:
            label = calcpath
        self.filename_label.setText(label)

        # Populate the table
        entires = []
        headers = []
        for key, value in self.calc.calc_overview_dict().items():
            headers.append(key)
            entires.append([value])
        self.param_table.setRowCount(len(headers))
        for i, header in enumerate(headers):
            self.param_table.setVerticalHeaderItem(i, QtGui.QTableWidgetItem(header))
        self.update_qtable_entires(self.param_table, entires)

        # Add plots
        self.make_ion_plot()
        self.make_e_plot()
        # Get plot details

        # free_energy= True
        # unit = 'eV/atom'
        #
        # # Prepare arrays
        # energies_plotted = self.calc.scale_energies(free_energy, unit)
        # energies = [ionstep[-1] for ionstep in energies_plotted]
        #
        # # Plot them
        # self.plot_1.clear()
        # self.plot_2.clear()
        # self.plot_1.plot(energies, symbol = 'o', size=100, pen='k')
        # i = 0
        # pens = ['r', 'g', 'b', 'c', 'm', 'k']
        # for ionstep in energies_plotted:
        #     steps = [step+i for step in range(len(ionstep))]
        #     thispen = pens[i % len(pens)]
        #     self.plot_2.plot(steps, ionstep, linestyle = '-', symbol='o', pen=thispen, symbolPen=thispen)
        #     i += len(ionstep)

    def setup_plots(self, units = 'eV', fontsize = 16):
        pass

    def update_qtable_entires(self, q_table_widget: QtGui.QTableWidget, data : list):
        q_table_widget.setRowCount(len(data))
        for i, row in enumerate(data):
            for j, entry in enumerate(row):
                q_table_widget.setItem(i, j, QtGui.QTableWidgetItem(str(entry)))

    def make_ion_plot(self):
        # Read from drop down and checkboxes
        unit = str(self.unit_combo_1.currentText())
        free_energy = self.free_energy_1.isChecked()
        delta = self.delta_1.isChecked()
        # Axis labels
        ylabel = ('Free' if free_energy else 'Total') + ' energy' + (' delta ' if delta else '')
        labelstyle = self.labelstyle
        self.plot_1.setLabel('left', ylabel, units=unit, **labelstyle)
        self.plot_1.setLabel('bottom', 'Step #', **labelstyle)
        # Prepare arrays
        energies_plotted = self.calc.scale_energies(free_energy, unit)
        energies = [ionstep[-1] for ionstep in energies_plotted]
        if delta:
            delta_e = [0]
            for i in range(len(energies)-1):
                delta_e.append(energies[i] - energies[i+1])
            energies = delta_e
        # Plot them
        self.plot_1.clear()
        self.plot_1.plot(energies, symbol='o', size=100, pen='k')

    def make_e_plot(self):
        unit = str(self.unit_combo_2.currentText())
        free_energy = self.free_energy_2.isChecked()
        delta = self.delta_2.isChecked()
        # Axis labels
        ylabel = ('Free' if free_energy else 'Total') + ' energy' + (' delta ' if delta else '')
        labelstyle = self.labelstyle
        self.plot_2.setLabel('left', ylabel, units=unit, **labelstyle)
        self.plot_2.setLabel('bottom', 'Step #', **labelstyle)
        # Prepare arrays
        energies_plotted = self.calc.scale_energies(free_energy, unit)

        # Plot them
        self.plot_2.clear()
        i = 0
        pens = ['r', 'g', 'b', 'c', 'm', 'k']
        for ionstep in energies_plotted:
            steps = [step + i for step in range(len(ionstep))]
            thispen = pens[i % len(pens)]
            self.plot_2.plot(steps, ionstep, linestyle='-', symbol='o', pen=thispen, symbolPen=thispen)
            i += len(ionstep)

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = PWX_Analysis()
    window.show()
    sys.exit(app.exec_())
