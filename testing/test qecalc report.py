from classes import qecalc
import matplotlib.pyplot as plt
demofile = r'D:\Users\Chris\Documents\SivaLab\2016 MDA Type 2 SL\Ab-Initio\Bulk\InAsSb_Ideal\Cutoff-kpts\kpts-6\cutoff=40.9090909091\InAsSb_Ideal_Bulk_PAW.in.out'

democalc = qecalc.QECalcOut.import_from_file(demofile)
fig = plt.figure()

ion_axes = fig.add_subplot(2,1,1) # 2 rows, 2 columns, 1st plot
relaxfigure = democalc.plot_ion_energies(ion_axes, True, unit='eV/atom')

e_axes = fig.add_subplot(2,1,2)
ion = democalc.plot_e_energies(e_axes, free_energy=True, unit='eV/atom')

pl_size = fig.get_size_inches()
scale=1.5
fig.set_size_inches(pl_size[0]*scale, pl_size[1]*scale)
fig.savefig('relaxation.png', dpi=100)

scale=2
piechart = democalc.calctime.pie_charts()
pl_size = piechart.get_size_inches()
piechart.set_size_inches(pl_size[0]*scale, pl_size[1]*scale)
piechart.savefig('time.png', dpi=100)

plt.show()

