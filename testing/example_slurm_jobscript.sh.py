#!/bin/env python3
#
#SBATCH --job-name=%(jobname)s
#SBATCH --output=%(slurmout)s
#
#SBATCH --ntasks=%(num_cores)s
#SBATCH --time=%(time)s
#SBATCH --mem-per-cpu=%(memory)
#
# Setup Imports
# basic imports
import subprocess as subpr
import datetime as dt
import os
# biscotti imports
from biscotti.classes import qecalc
from biscotti.reporting import email
# Email at start'
body_start = 'Your QE calculation %(jobname)s began on ' + str(dt.datetime.now()) + '\n' + "full path is: \n" + os.path.abspath('')
email.send_mail('%(email_addr)s', 'QE calculation %(jobname)s has started', body_start)

# Begin pw.x
subpr.call('mpirun -np %(num_cores)s pw.x -i %(infile)s > %(infile)s.out', shell=True)
# Email at end
body_end = 'Your QE calculation %(jobname)s ended on ' + str(dt.datetime.now()) + '\n' + "full path is: \n" + os.path.abspath('')
calcout = qecalc.QECalcOut.import_from_file('%(infile)s.out', '%(infile)s')
calcout.make_report(reportname='%(jobname)s report.png')
email.send_mail('%(email_addr)s', 'QE calculation %(jobname)s has ended', body_end, ['%(jobname)s report.png'])
