# Email syntax:
# $python pymail.py address subject body

import smtplib
import sys
fromaddr = 'epirlab@gmail.com'
toaddrs  = sys.argv[1] #Note, this can also be a list
subject =sys.argv[2]

body = ''
for i in range(3,len(sys.argv)):
    body = body + sys.argv[i] + " "

message = 'Subject: %s\n\n%s' % (subject, body)

username = 'epirlab@gmail.com'
password = '31a4!tBzP!xb'

#now connect and login
server = smtplib.SMTP('smtp.gmail.com:587')
server.starttls()
server.login(username,password)
server.sendmail(fromaddr, toaddrs, message)
server.quit()
print("Email sent!")
