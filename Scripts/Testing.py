from classes import qecalc
from classes import atoms
outfile = r'D:\Users\Chris\Documents\SivaLab\2016 MDA Type 2 SL\Ab-Initio\Superlattice\Batch006 - Small, Ideal SL\InAs_1,1,3\InAs_Best_Cubic.in expanded + InAsSb_Best_Cubic_Strained.in expanded.in.out'
outobj = qecalc.QECalcOut.import_from_file(outfile)
for ionstep in outobj.relax_list:
    print ("This ion step had " + str(len(ionstep)) + " electron steps")