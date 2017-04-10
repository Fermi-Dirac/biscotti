from os.path import join as joinpath
from biscotti.classes import atoms

bulkpeices = r'D:\Users\Chris\Documents\SivaLab\Type II Superlattices\2016 MDA Type 2 SL\Ab-Initio\Superlattice\Bulk Pieces'
file_InAsprim = joinpath(bulkpeices, 'InAs_prim.atoms')
file_InSbprim = joinpath(bulkpeices, 'InAsSb_prim_sorta.atoms')
InAsprim = atoms.AtomicStructure.from_QEinput(file_InAsprim)
InSbprim = atoms.AtomicStructure.from_QEinput(file_InSbprim)
print(InAsprim)
InAsprim.strain(0.006548)
print(InAsprim)
print(InSbprim)
InSbprim.strain(-0.02695)
print(InSbprim)
sillySL = InAsprim.supercell([1,1,27]).stack_with(InSbprim.supercell([1,1,6]))
print(sillySL)
sillySL.write_vasp(bulkpeices, 'POSCAR sillySL.vasp')