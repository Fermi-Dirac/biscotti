from biscotti.classes import qecalc

bugfile = r'D:\Users\Chris\Documents\SivaLab\Type II Superlattices\2016 MDA Type 2 SL\Ab-Initio\Superlattice\Batch010 - Test Cluster\InAsSb_1,1,3,3_relax\InAsSb_10-3Strain_base.in.out'
bugout = qecalc.QECalcOut.import_from_file(bugfile)
print(bugout.qecalcin.structure)
print(bugout.finalstructure)