from biscotti.classes import atoms
file_primitive = r'D:\Users\Chris\Documents\SivaLab\Type II Superlattices\2016 MDA Type 2 SL\Ab-Initio\Superlattice\Bulk Pieces\InAs_prim.atoms'
prim1 = atoms.AtomicStructure.from_QEinput(file_primitive)
prim2 = atoms.AtomicStructure.from_QEinput(file_primitive)
print(prim1)
prim2.atoms[1].species = 'Sb'
print(prim2)
stacked = prim1.stack_with(prim2)
print(stacked)
stacked.write_vasp(r'D:\Users\Chris\Documents\SivaLab\Type II Superlattices\2016 MDA Type 2 SL\Ab-Initio\Superlattice\Bulk Pieces')