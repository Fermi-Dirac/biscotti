from Main.Biscotti import *

basefile = r'D:\Documents\SivaLab\2016 MDA Type 2 SL\Ab-Initio\Bulk\InAsSb_Ideal\Cutoff-kpts\InAsSb_Ideal_Bulk_PAW.in'
kpointslist = ['4 4 4 0 0 0', '5 5 5 0 0 0', '6 6 6 0 0 0', '8 8 8 0 0 0', '10 10 10 0 0 0']
for kpts in kpointslist:
    for cutoff in np.linspace(25,60,12):
        makenewcalc(basefile, kpts,cutoff)
