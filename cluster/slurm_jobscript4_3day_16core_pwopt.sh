#!/bin/env python3
#
#SBATCH --job-name=opt16cre
#SBATCH --output=slurmout.txt
#SBATCH --partition=cluster
#
#SBATCH --ntasks=16
#SBATCH --time=71:30:00
#SBATCH --mem-per-cpu=MaxMemPerNode

import subprocess as subpr
import datetime as dt
import os
import sys
try:
    from biscotti.classes import qecalc
    has_qecalc = True
except ImportError:
    print("Cannot load Biscotti")
    has_qecalc = False
try:
    from biscotti.reporting import email
    has_email = True
    email_addr = 'cbuurma@sivananthanlabs.us'
except ImportError:
    print("Cannot load email code")
    has_email = False

if len(sys.argv) > 1:
    input_file = sys.argv[1]
else:
    input_file = [file for file in os.listdir() if os.path.splitext(file)[1] == '.in'][0] # Get first .in file
# Email at start
if has_email:
    body = 'Your QE calculation' + input_file + ' began on ' + str(dt.datetime.now()) + \
           '\nThe full execution path is: \n' + os.path.abspath('')
    email.send_mail(email_addr, 'Starting QE calculation: ' + input_file, body)
    print("Email sent!")

# Begin pw.x call
subpr.call('mpirun -np 16 --map-by core --bind-to core pw_opt.x -i ' + input_file + ' > ' + input_file + '.out', shell=True)

# Email at end
if has_email:
    body_end = 'Your QE calculation' + input_file + ' ended on ' + str(dt.datetime.now()) + \
               '\nThe full execution path is: \n' + os.path.abspath('')
    if has_qecalc:
        calcout = qecalc.QECalcOut.import_from_file(input_file + '.out', input_file)
        body_end += '\n\n' + calcout.calc_overview_string(transpose=True)
        try: # Just in case no Matplotlib got loaded, but qecalc loaded ok
            calcout.make_report(reportname=input_file + ' report.png')
            email.send_mail(email_addr, 'QE calculation ' + input_file + ' has ended', body_end,
                            [input_file + ' report.png'])
        except Exception:
            email.send_mail(email_addr, 'QE calculation: ' + input_file + ' has ended', body_end)
    else:
        email.send_mail(email_addr, 'QE calculation ' + input_file + ' has ended', body_end)
    print("Email sent!")
subpr.call('Python script complete!', shell=True)
