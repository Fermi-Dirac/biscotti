try:
    import pymatgen as pmg
    pymatgen_avail = True
except ImportError:
    print("Error, pymatgen not found. associated functions disabled")
    pymatgen_avail = False

import biscotti.classes.atoms as atoms

def struct_biscotti_to_pymatgen(bis_struct : atoms.AtomicStructure):
    coords = []
    species = []
    for atom in bis_struct.atomsdir:
        coords.append(list(atom.position))
        species.append(atom.species)
    lattice = pmg.Lattice(bis_struct.lattice)
    pmg_struct = pmg.Structure(lattice, species, coords)
    return pmg_struct

def struct_pymatgen_to_biscotti():
    pass

