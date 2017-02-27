#!/bin/env python3
#
#SBATCH --job-name=Default
#SBATCH --output=slurmout.txt

##SBATCH --ntasks=16
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=MaxMemPerNode
import subprocess as subpr
import datetime as dt
import os

from biscotti.classes import qecalc
from biscotti.reporting import email

# Email at start
body = 'Your QE calculation Default began on ' + str(dt.datetime.now()) + '\nThe full execution path is: \n' + os.path.abspath('')
email.send_mail('', 'QE calculation Default has started', body)
# Begin pw.x call
subpr.call('mpirun -np 16 pw.x -i None > None.out', shell=True)
# Email at end
body_end = 'Your QE calculation Default ended on ' + str(dt.datetime.now()) + '\nThe full execution path is: \n' + os.path.abspath('')
calcout = qecalc.QECalcOut.import_from_file('None.out', 'None')
calcout.make_report(reportname='Default report.png')
email.send_mail('', 'QE calculation Default has ended', body_end, ['Default report.png'])
subpr.call('python script complete!', shell=True

