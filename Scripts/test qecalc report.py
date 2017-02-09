from classes import qecalc
from classes import atoms

demofile = r'D:\Users\Chris\Documents\SivaLab\2016 MDA Type 2 SL\Ab-Initio\Bulk\InAsSb_Ideal\Cutoff-kpts\kpts-6\cutoff=40.9090909091\InAsSb_Ideal_Bulk_PAW.in.out'
atoms.AtomicStructure.from_QEOutput(demofile)

# democalc = qecalc.QECalcOut.import_from_file(demofile)