from classes import atoms as Biscotti
from classes import qecalc as Parsers
import os
from os.path import join as joinpath
import numpy as np

bulkpeices = r'D:\Users\Chris\Documents\SivaLab\2016 MDA Type 2 SL\Ab-Initio\Superlattice\Bulk Pieces'
InAsPurefile = joinpath(bulkpeices, r'InAs_10-3_strained_SL.in')
InAsSbPurefile = joinpath(bulkpeices, r'InAsSb_X=0.5_10-3_strained_SL.in')

newbatch = r'Batch007 - Faster'
newpath = r'D:\Users\Chris\Documents\SivaLab\2016 MDA Type 2 SL\Ab-Initio\Superlattice' + os.path.sep + newbatch
layersWithSb = 3
layerswithInAs_List = range(3,10)
Xlength = 1
Ywidth = 1

InAs = Biscotti.AtomicStructure.from_QEinput(InAsPurefile)
InAsSb = Biscotti.AtomicStructure.from_QEinput(InAsSbPurefile)
SbSlab = InAsSb.supercell([Xlength, Ywidth, layersWithSb])

# Remove Central Sb Atom
SbVacancy = False
Swap = False
# Remove a single atom from the middle of the InSb Section

if SbVacancy or Swap:
    SbOnly = []
    AsOnly = []
    bestindex = 0
    lastdistance = 999
    for i, atom in enumerate(SbSlab.atomsdir):
        if atom.species == 'Sb':
            print(atom)
            newdistance = np.linalg.norm(atom.position - np.array([0.5, 0.5, 0.5]))
            print("Distance to center is " + str(newdistance))
            if newdistance < lastdistance:
                lastdistance = newdistance
                bestindex = i
    AsIndex = 0
    lastdistance = 999
    for i, atom in enumerate(SbSlab.atomscart):
        if atom.species == 'As':
            print(atom)
            newdistance = np.lingalg.norm(atom.position - SbSlab.atomscart[bestindex].position)
            if newdistance < lastdistance:
                lastdistance = newdistance
                AsIndex = i
    if Swap:
        print("Now swapping In and As positions")
        OrgSb = np.array(SbSlab.atoms[bestindex].position)
        OrgAs = np.array(SbSlab.atoms[AsIndex].position)
        SbSlab.atoms[bestindex].position = OrgAs
        SbSlab.atoms[AsIndex].position = OrgAs
        SbSlab.cart_to_direct(SbSlab.atomscart)
    if SbVacancy:
        print ("Now removing Sb at " + str(SbSlab.atomsdir[bestindex].position), bestindex, lastdistance, SbSlab.totalatoms())
        SbSlab.atomsarray.pop(bestindex)
        SbSlab.atomsdir.pop(bestindex)
    print ("Now we have " + str(SbSlab.totalatoms()))
    print (len(SbSlab.atomsdir))

for layersOfInAs in layerswithInAs_List:
    InAsSlab = InAs.supercell([Xlength, Ywidth, layersOfInAs])
    calcStructure = InAsSlab.merge_with(SbSlab)
    folder = newpath + os.path.sep + "InAs_" + str(Xlength) + "," + str(Ywidth) + "," + "%d" % layersOfInAs + os.path.sep
    if not os.path.exists(folder):
        os.makedirs(folder)
    newfile =  folder + os.path.sep + calcStructure.name + ".in"
    calcStructure.write_vasp(newpath, "POSCAR.InAs_" + str(Xlength) + "," + str(Ywidth) + "," + "%d" % layersOfInAs + ".vasp")
    print("new file created! " + newfile)
    # This next line of code is kinda shit and needs replaced with the new class ocne it works
    # TODO read that line above
    Parsers.changeStructure(InAsSbPurefile, newfile, calcStructure, AutoKpts= 20)