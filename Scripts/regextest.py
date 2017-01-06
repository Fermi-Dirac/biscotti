import re
import logging

regex = r'CELL_PARAMETERS'
regex2 = r'ATOMIC_SPECIES'
QEinfile = r'D:\Users\Chris\Documents\SivaLab\2016 MDA Type 2 SL\Ab-Initio\Superlattice\Batch006 - Small, Ideal SL\InAs_1,1,9\InAs_Best_Cubic.in expanded + InAsSb_Best_Cubic_Strained.in expanded.in'
with open(QEinfile, 'r') as fileobj:
    filestring = fileobj.read()
    startoflattice = re.search(regex, filestring).end()
    print (re.search(regex, filestring).group(0))
    endoflat = re.search(regex2, filestring).start()
    print (filestring[startoflattice:endoflat])