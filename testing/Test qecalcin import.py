import os
from biscotti.classes import qecalc
testfile = r'D:\Users\Chris\Documents\SivaLab\Type II Superlattices\2016 MDA Type 2 SL\Ab-Initio\Superlattice\Batch020 - Clusters\T2SL_[2,2,3,5],strain=83.3A-20.2A-SCF\T2SL_[2,2,3,5],strain=83.3A-20.2A-SCF.in'
print(os.path.basename(testfile))
calc = qecalc.QECalcIn.import_from_file(testfile)
print(calc.control['title'], calc.name)