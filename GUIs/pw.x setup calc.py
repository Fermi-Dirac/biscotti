import sys, os, pickle
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

qtCreatorFile = r"pw.x setup calc.ui"  # Enter file here.

Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)
# Widgets:
# commands_layout       QVBoxLayout
# # defaults_label       QLabel
# # defaults_list        QListWidget
# # description          QTextBrowser
# # command_val_layout
# ## command_cmbox      QComboBox
# ## command_label      QLabel
# ## val_cmbox          QComboBox
# ## val_label          QLabel
# text_file_layout
# # filename_layout
# ## filename_lineedit  QLineEdit
# ## folder_button      Qbutton
# file_display          QPlainTextEdit

class PWX_Setup(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self):
        # Setup calls
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)

        # Events to listen for
        self.defaults_list.itemClicked.connect(self.load_selected_default)
        self.command_cmbox.currentIndexChanged.connect(self.change_command)
        self.val_cmbox.currentIndexChanged.connect(self.change_val)
        self.load_file_button.clicked.connect(self.load_file_dialog)
        self.file_display.cursorPositionChanged.connect(self.cursor_changed)

        # Main data structures
        self.last_rootpath = None
        self.defaults_path = os.path.abspath('../input_templates/pw.x_templates')
        self.input_text = ''
        self.input_list = []
        self.namelist_dict = pickle.load(open('..\data\pw.x namelist_dict.p', 'rb'))
        # Usage: namelist_dict[namelist] = [list of possible commands]
        self.command_dict = pickle.load(open('..\data\pw.x command_dict.p', 'rb'))
        # Usage: command_dict[command] = (default, options or a list, description string)
        self.val_comments = pickle.load(open('..\data\pw.x val_comments.p', 'rb'))
        # Usage: val_comments[value] = comment

        # Pre-load data
        self.populate_default_list()
        self.load_selected_default()
        self.command_list = []
        for key, list in self.namelist_dict.items():
            self.command_list.append("  " + key)
            self.command_list.extend(list)
        self.command_cmbox.addItems(self.command_list)

    def populate_default_list(self):
        for filename in os.listdir(self.defaults_path):
            self.defaults_list.addItem(QtGui.QListWidgetItem(filename))
        self.defaults_list.setCurrentRow(0)

    def load_selected_default(self):
        filename = self.defaults_list.currentItem().text()
        path = os.path.join(self.defaults_path, filename)
        self.filename_lineedit.setText(path)
        with open(path, 'r') as fileobj:
            self.input_text = fileobj.read()
            self.input_list = self.input_text.split('\n')
            self.file_display.setPlainText(self.input_text)

    def load_file_dialog(self):
        calcpath = QtGui.QFileDialog.getOpenFileName(caption='Select calculation in file',
                                               directory=self.last_rootpath)
        self.last_rootpath = calcpath
        self.load_file(calcpath)

    def load_file(self, path):
        with open(path, 'r') as fileobj:
            self.input_text = fileobj.read()
            self.input_list = self.input_text.split('\n')
            self.file_display.setPlainText(self.input_text)

    def change_command(self):
        command = self.command_cmbox.currentText()
        self.val_cmbox.clear()
        if command in self.command_dict:
            default, options, description = self.command_dict[command]
            self.description.setText(description)
            if type(options) is list:
                self.val_cmbox.addItems(options)
            elif type(options) is bool:
                if default.lower() == '.true.':
                    other = '.false.'
                else:
                    other = '.true.'
                self.val_cmbox.addItems([default, other])
            else:
                self.val_cmbox.addItem(default)
            self.val_cmbox.setCurrentIndex(0)
        else:
            self.description.setText('Cannot find description for command ' + command)

    def change_val(self):
        command = self.command_cmbox.currentText()
        value = self.val_cmbox.currentText()
        comment = ''
        if command in self.val_comments:
            if value in self.val_comments[command]:
                comment = self.val_comments[command][value]
        self.val_comment_label.setText(comment)

    def cursor_changed(self):
        cursor_pos = self.file_display.textCursor().position()
        linestart = 0
        lineend = len(self.input_text)
        start = 0 if cursor_pos < 200 else cursor_pos-200
        end = cursor_pos+200 if cursor_pos + 200 < lineend else lineend
        for index, char in enumerate(self.input_text[cursor_pos:end]):
            if char == '\n':
                lineend = cursor_pos + index
                break
        for index, char in enumerate(self.input_text[start:cursor_pos][::-1]):
            if char == '\n':
                linestart = cursor_pos - index
                break
        selected_line = self.input_text[linestart:lineend]
        linesplit = selected_line.split('=')
        if len(linesplit) > 1:
            command = linesplit[0].strip()
            value = linesplit[1].split(',')[0].strip()
            print(command,value)
            if command in self.command_list:
                self.command_cmbox.setCurrentIndex(self.command_list.index(command))

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = PWX_Setup()
    window.show()
    sys.exit(app.exec_())
