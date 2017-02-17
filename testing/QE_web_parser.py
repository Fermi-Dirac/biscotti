import parsel
import requests

useweb = False
if useweb:
    ppx_url = r'http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PP.html'
    ppx_text = requests.get(ppx_url).text
else:
    with open('../data/pp.x input description.html', mode='r') as fileobj:
        ppx_text = fileobj.read()

sel = parsel.Selector(text=ppx_text)
print(sel.css('title::text').extract()[0])
print(sel.css('h3::text'))