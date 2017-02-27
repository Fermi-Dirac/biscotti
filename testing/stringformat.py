class bored:
    def wtf(self, nameofthing = 'lawl'):
        slurm1 = (
            "#!/bin/env python3\n"
            "#\n"
            "#SBATCH --job-name=%(nameofthing)s\n")
        slurm2 = slurm1 % locals()
        print(slurm2)

lawl = bored()
lawl.wtf()

from biscotti.classes import qecalc

democalc = qecalc.QECalcIn()

democalc.write_slurm_jobscript(folder=r'D:\Users\Chris\Documents\SivaLab\Python\biscotti\testing')