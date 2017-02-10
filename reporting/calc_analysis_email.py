# command line usage syntax
# $python calc_analysis_email.py address outputfile inputfile

import smtplib
from os.path import basename
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.utils import COMMASPACE, formatdate
import sys
from classes import qecalc
import os
# --------------
# Prepare Report
# --------------
rootpath = os.path.abspath('') # Get calling directory
try:
    outpath = sys.argv[2]
except Exception:
    outpath =

qecalc.QECalcOut.import_from_file()

# --------------
# Prepare email
# --------------
fromaddr = 'epirlab@gmail.com'
try:
    toaddrs  = sys.argv[1] #Note, this can also be a list
except Exception:
    toaddrs = 'cbuurma@sivananthanlabs.us'

subject = 'Analysis Report of QE calculation: '
body = 'Quantum Espresso calculation "'
relaxplot = r'D:\Users\Chris\Documents\SivaLab\Python\biscotti\testing\relaxation.png'
attachments = []# [relaxplot]

def send_mail(send_from, send_to, subject, text, files=None,
              server="127.0.0.1"):
    #assert isinstance(send_to, list)
    username = 'epirlab@gmail.com' # I know this is terribly unsafe
    password = '31a4!tBzP!xb'

    msg = MIMEMultipart()
    msg['From'] = send_from
    msg['To'] = COMMASPACE.join(send_to)
    msg['Date'] = formatdate(localtime=True)
    msg['Subject'] = subject

    msg.attach(MIMEText(text))

    for f in files or []:
        with open(f, "rb") as fil:
            part = MIMEApplication(
                fil.read(),
                Name=basename(f)
            )
            part['Content-Disposition'] = 'attachment; filename="%s"' % basename(f)
            msg.attach(part)

    smtp = smtplib.SMTP(server)
    smtp.starttls()
    smtp.login(username, password)
    smtp.sendmail(send_from, send_to, msg.as_string())
    # smtp.close()
    smtp.quit()

send_mail(fromaddr, toaddrs, subject, body, attachments, server='smtp.gmail.com:587')

# # Now send that email!
# message = 'Subject: %s\n\n%s' % (subject, body)
#
#
# #now connect and login
# server = smtplib.SMTP('smtp.gmail.com:587')
# server.starttls()
# server.login(username, password)
# server.sendmail(fromaddr, toaddrs, message)
# server.quit()
# print("Email sent!")
