from classes import atoms
from classes import qecalc
import os
from os.path import join as joinpath
import numpy as np

bulkpeices = r'D:\Users\Chris\Documents\SivaLab\2016 MDA Type 2 SL\Ab-Initio\Superlattice\Bulk Pieces'
InAsPurefile = joinpath(bulkpeices, r'InAs_10-3_strained_SL.in')
InAsSbPurefile = joinpath(bulkpeices, r'InAsSb_X=0.5_10-3_strained_SL.in')

newbatch = r'Batch007 - Faster'
newpath = r'D:\Users\Chris\Documents\SivaLab\2016 MDA Type 2 SL\Ab-Initio\Superlattice' + os.path.sep + newbatch
layersWithSb = 3
layerswithInAs_List = [3]
Xlength = 1
Ywidth = 1

InAs = atoms.AtomicStructure.from_QEinput(InAsPurefile)
InAsSb = atoms.AtomicStructure.from_QEinput(InAsSbPurefile)
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

InAsSlab = InAs.supercell([Xlength, Ywidth, layerswithInAs_List[0]])
calcStructure = InAsSlab.merge_with(SbSlab)
calcStructure.name = '[InAs]Sb_1x1x3-3_T2SSL'

# Load calculation settings form the InAs file
calcin = qecalc.QECalcIn.import_from_file(InAsSbPurefile)
calcin.structure = calcStructure
calcin.name = calcStructure.name
calcin.control['title'] = calcStructure.name
calcin.control['prefix'] = 'InAsT2SL'
calcin.control['calculation'] = 'relax'

for ecut_rho in np.linspace(40*2, 40*6, 9):
    calcin.system['ecutrho'] = ecut_rho
    folder = joinpath(newpath, "ecutrho=" + str(ecut_rho))
    if not os.path.exists(folder):
        os.makedirs(folder)
    calcin.write_to_file(folder, calcin.name + "ecutrho=" + str(ecut_rho) + ".in")

del calcin.system['ecutrho']

for smear in np.linspace(0.01, 0.0001, 9):
    calcin.system['degauss'] = smear
    folder = joinpath(newpath, "smear="+str(smear))
    if not os.path.exists(folder):
        os.makedirs(folder)
    calcin.write_to_file(folder,calcin.name + "_smear=" + str(smear) + ".in")

calcin.system['degauss'] = 0.001

for dE in np.linspace(1e-3, 1e-4, 5):
    calcin.control['etot_conv_thr'] = dE
    suffix = "ion_dE=" + str(dE)
    folder = joinpath(newpath, suffix)
    if not os.path.exists(folder):
        os.makedirs(folder)
    calcin.write_to_file(folder, calcin.name + "_" + suffix + ".in")

calcin.control['etot_conv_thr'] = 1e-3

calcin.control['calculation'] = 'vc-relax'
folder = joinpath(newpath, 'vc-relax')
if not os.path.exists(folder):
    os.makedirs(folder)
calcin.write_to_file(joinpath(folder), calcin.name + "_vc-relax.in")
