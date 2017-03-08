from classes import qecalc
from tkinter import filedialog
calcpath = filedialog.askopenfilename()

calc = qecalc.QECalcOut.import_from_file(calcpath)
print(calc.calc_overview_string(transpose=True))
