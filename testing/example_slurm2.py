import subprocess as subpr
import datetime as dt
import os

from biscotti.classes import qecalc
from biscotti.reporting import email

# Email at start
body = 'Your QE calculation [InAs]Sb_1x1x3-3_T2SSL began on ' + str(dt.datetime.now()) + '\nThe full execution path is: \n' + os.path.abspath('')
email.send_mail('cbuurma@sivananthanlabs.us', 'QE calculation [InAs]Sb_1x1x3-3_T2SSL has started', body)
# Begin pw.x call
subpr.call('mpirun -np 16 pw.x -i [InAs]Sb_1x1x3-3_T2SSLecutrho=80.0.in > [InAs]Sb_1x1x3-3_T2SSLecutrho=80.0.in.out', shell=True)
# Email at end
body_end = 'Your QE calculation [InAs]Sb_1x1x3-3_T2SSL ended on ' + str(dt.datetime.now()) + '\nThe full execution path is: \n' + os.path.abspath('')
calcout = qecalc.QECalcOut.import_from_file('[InAs]Sb_1x1x3-3_T2SSLecutrho=80.0.in.out', '[InAs]Sb_1x1x3-3_T2SSLecutrho=80.0.in')
calcout.make_report(reportname='[InAs]Sb_1x1x3-3_T2SSL report.png')
email.send_mail('cbuurma@sivananthanlabs.us', 'QE calculation [InAs]Sb_1x1x3-3_T2SSL has ended', body_end, ['[InAs]Sb_1x1x3-3_T2SSL report.png'])
subpr.call('python script complete!', shell=True)
