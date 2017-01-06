from classes import qecalc, atoms
import os
import numpy as np
bulkpieces = r'D:\Users\Chris\Documents\SivaLab\2016 MDA Type 2 SL\Ab-Initio\Superlattice\Bulk Pieces'
InAsPure = r'InAs_Best_Cubic.in'
InAsSbPure = r'InAsSb_Best_Cubic.in'

AsStruct = atoms.AtomicStructure.from_QEinput(os.path.join(bulkpieces,InAsPure))
SbStruct = atoms.AtomicStructure.from_QEinput(os.path.join(bulkpieces,InAsSbPure))

# With experimental strain in this SL has n,m layers of each piece..:
layers_As = 10
layers_Sb = 3

newlattice = []
for lat1v, lat2v in zip(AsStruct.lattice, SbStruct.lattice):
    newlattice.append((layers_As*lat1v + layers_Sb*lat2v )/(layers_As + layers_Sb))
print (newlattice)

AsStruct.latticeA = newlattice[0]
AsStruct.latticeB = newlattice[1]
AsStruct.latticeC = newlattice[2]
SbStruct.latticeA = newlattice[0]
SbStruct.latticeB = newlattice[1]
SbStruct.latticeC = newlattice[2]
# So ugly...
AsStruct.name = 'InAs_10-3_strained_SL'
SbStruct.name = 'InAsSb_X=0.5_10-3_strained_SL'

infile = qecalc.QECalcIn.import_from_file(os.path.join(bulkpieces, InAsPure))
infile.structure = AsStruct
infile.name = AsStruct.name
infile.write_to_file(bulkpieces)

infile.structure = SbStruct
infile.name = SbStruct.name
infile.write_to_file(bulkpieces)