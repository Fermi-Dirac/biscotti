from Main import Biscotti
from Main import Parsers

import numpy as np
InAsPurefile = r'D:\Documents\SivaLab\2016 MDA Type 2 SL\Ab-Initio\Superlattice\Bulk Pieces\InAs_Best_Cubic.in'
InAsSbPurefile = r'D:\Documents\SivaLab\2016 MDA Type 2 SL\Ab-Initio\Superlattice\Bulk Pieces\InAsSb_Best_Cubic_Strained.in'

InAsDict = Parsers.parseConfig(InAsPurefile)
for key in InAsDict:
    print (key)
    for subkey in InAsDict[key]:
        print (" " + str(subkey))
