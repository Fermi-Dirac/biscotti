import pymatgen as pmg
import biscotti.classes.atoms as atoms

def biscotti_to_pymatgen(bis_struct : atoms.AtomicStructure):
    coords = []
    species = []
    for atom in bis_struct.atomsdir:
        coords.append(list(atom.position))
        species.append(atom.species)
    lattice = pmg.Lattice(bis_struct.lattice)
    pmg_struct = pmg.Structure(lattice, species, coords)
    return pmg_struct

def pymatgen_Structure_to_AtomicStructure():
    pass

