from classes import qecalc
import matplotlib.pyplot as plt
demofile = r'D:\Users\Chris\Documents\SivaLab\Type II Superlattices\2016 MDA Type 2 SL\Ab-Initio\Bulk\InAsSb_Ideal\Cutoff-kpts\kpts-6\cutoff=40.9090909091\InAsSb_Ideal_Bulk_PAW.in.out'

democalc = qecalc.QECalcOut.import_from_file(demofile)
fig = plt.figure(figsize=(22,11))
ion_axes = fig.add_subplot(2,2,1)
e_axes = fig.add_subplot(2,2,3)
table1_axes = fig.add_subplot(3,2,2)
table2_axes = fig.add_subplot(3,2,4)
table3_axes = fig.add_subplot(3,2,6)

democalc.plot_ion_energies(ion_axes, free_energy=True, unit='eV/atom')
democalc.plot_e_energies(e_axes, free_energy=True, unit='eV/atom')

headers1 = ['Title', 'Start Date-Time', 'Calc time (hr)', 'Job Complete?']
headers2 = ['Calc Type', 'Cutoff (Ry)', 'Total # of K-points', 'Pseudspots used']
headers3 = ['Final Free Energy (eV/atom)', 'Last electron step dE (Ry)', 'Final Presure(kbar)']
democalc.plot_summary_table(table1_axes,headers1)
democalc.plot_summary_table(table2_axes, headers2)
democalc.plot_summary_table(table3_axes, headers3)

# for key, value in newtable.properties().items():
#     print(str(key) + '\t' + str(value))
plt.show()