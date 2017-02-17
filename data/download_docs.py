import requests

pwx_url = r'http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PW.html'
ppx_url = r'http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PP.html'
ppx_text = requests.get(ppx_url).text