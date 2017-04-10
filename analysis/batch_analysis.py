from classes import qecalc
from tkinter import filedialog
rootpath = filedialog.askdirectory()
qecalc.makeSummaryFile(rootpath)