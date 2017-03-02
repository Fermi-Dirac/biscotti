from biscotti.classes import qecalc

bugfile = r'D:\Users\Chris\Documents\SivaLab\Type II Superlattices\2016 MDA Type 2 SL\Ab-Initio\Superlattice\Batch012 - BigCalcs\InAsSb_1,1,3,3_relax2\InAsSb_1,1,3,3_fast_relax.in'
bugin = qecalc.QECalcIn.import_from_file(bugfile)
print(bugin.pseudopots)