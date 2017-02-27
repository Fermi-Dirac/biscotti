import smtplib
from os.path import basename, exists
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.utils import formatdate

def send_mail(send_to : str, subject : str, text : str, files=None,
              server="smtp.gmail.com:587"):
    if type(send_to) is list:
        send_to = ", ".join(send_to) #check for multiple recipients
    username = 'epirlab@gmail.com' # I know this is terribly unsafe
    send_from = username
    password = '31a4!tBzP!xb'

    msg = MIMEMultipart()
    msg['From'] = send_from
    msg['To'] = send_to
    msg['Date'] = formatdate(localtime=True)
    msg['Subject'] = subject
    msg.attach(MIMEText(text))

    for f in files or []:
        if exists(f):
            with open(f, "rb") as fil:
                part = MIMEApplication(fil.read(), Name=basename(f))
                part['Content-Disposition'] = 'attachment; filename="%s"' % basename(f)
                msg.attach(part)
        else:
            print("Error, file: " + f + " does not exist.")

    smtp = smtplib.SMTP(server)
    smtp.starttls()
    smtp.login(username, password)
    smtp.sendmail(send_from, send_to, msg.as_string())
    # smtp.close()
    smtp.quit()

if __name__ == '__main__':
    import sys
    send_mail(sys.argv[1], sys.argv[2], sys.argv[3])