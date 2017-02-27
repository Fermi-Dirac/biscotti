from biscotti.classes import qecalc

bugfile = r'D:\Users\Chris\Documents\SivaLab\Type II Superlattices\2016 MDA Type 2 SL\Ab-Initio\Learning\post_proc\Si-2\si.scf.in'

bugin = qecalc.QECalcIn.import_from_file(bugfile)

bugfile2 = r'D:\Users\Chris\Documents\SivaLab\Type II Superlattices\2016 MDA Type 2 SL\Ab-Initio\Learning\post_proc\Si-2\si.scf.in.out'

bugout = qecalc.QECalcOut.import_from_file(bugfile2, inpath=bugfile)