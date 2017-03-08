#!/bin/env python3
#
#SBATCH --job-name=nd1nk64
#SBATCH --output=slurmout.txt
#SBATCH --partition=bigmem
#
#SBATCH --ntasks=64
#SBATCH --time=5:30:00
#SBATCH --mem-per-cpu=MaxMemPerNode

import subprocess as subpr
import datetime as dt
import os
import sys
pw_x_flags = %(pw_x_flags)s
send_email = %(send_email)s
email_addr = %(email_addr)s
attach_report = %(attach_report)s

try:
    from biscotti.classes import qecalc
except ImportError:
    print("Cannot load qecalc")
    attach_report = False
try:
    from biscotti.reporting import email
except ImportError:
    print("Cannot load email code")
    send_email = False

if len(sys.argv) > 1:
    input_file = sys.argv[1]
else:
    input_file = [file for file in os.listdir() if os.path.splitext(file)[1] == '.in'][0] # Get first .in file
# Email at start
if send_email:
    body = 'Your QE calculation' + input_file + ' began on ' + str(dt.datetime.now()) + \
           '\nThe full execution path is: \n' + os.path.abspath('')
    email.send_mail(email_addr, 'Starting QE calculation ' + input_file, body)
    print("Email sent!")

# Begin pw.x call
subpr.call('mpirun -np %(num_cores)s --map-by core --bind-to core pw.x ' + pw_x_flags + ' -i ' + input_file + ' > ' + input_file + '.out', shell=True)

# Email at end
if send_email:
    body_end = 'Your QE calculation' + input_file + ' ended on ' + str(dt.datetime.now()) + \
               '\nThe full execution path is: \n' + os.path.abspath('')
    if attach_report:
        calcout = qecalc.QECalcOut.import_from_file(input_file + '.out', input_file)
        body_end += '\n\n' + calcout.calc_overview_string(transpose=True)
        try: # Just in case no Matplotlib got loaded, but qecalc loaded ok
            calcout.make_report(reportname=input_file + ' report.png')
            email.send_mail(email_addr, 'QE calculation ' + input_file + ' has ended', body_end,
                            [input_file + ' report.png'])
        except Exception:
            email.send_mail(email_addr, 'QE calculation ' + input_file + ' has ended', body_end)
    else:
        email.send_mail(email_addr, 'QE calculation ' + input_file + ' has ended', body_end)
    print("Email sent!")
print("Slurm job complete!")