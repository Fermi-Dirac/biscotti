from biscotti.classes.qecalc import QECalcOut
import numpy as np

# finds the bands at the gamma point for a given file using the 'bands' dictonary
outpath = r'\\cluster\cbuurma\espresso\Superlattice\Batch019 - Long vc-relax2\T2SL_[1,1,3,13],strain=83.3A-20.2A\T2SL_[1,1,3,13],strain=83.3A-20.2A.in.out'
thiscalc = QECalcOut.import_from_file(outpath=outpath)
kpts = list(thiscalc.bands[-1])
print(kpts)
gamma_bands = thiscalc.bands[-1][(0,0,0)] - thiscalc.fermi_levels[-1]  #offset by fermi level
print(thiscalc.fermi_levels)
idx = np.searchsorted(gamma_bands, 0, side='left')
print(gamma_bands[idx-4:idx+4])