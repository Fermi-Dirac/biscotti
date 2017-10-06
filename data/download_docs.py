# I think this is not needed, i should totally use the def files:
# http://web.mit.edu/espresso_v6.1/amd64_ubuntu1404/qe-6.1/PP/Doc/INPUT_PROJWFC.def

import requests
from bs4 import BeautifulSoup
import re
from collections import OrderedDict as odict
import pickle

# These will be populated
namelist_dict = odict({'None' : []}) # Usage: namelist_dict[namelist] = [list of possible commands]
command_dict = odict() # Usage: command_dict[command] = (default, options or a list, description string)
val_comments = odict() # Usage: val_comments[value] = comment

pwx_url = r'http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PW.html'
pwx_text = requests.get(pwx_url).text
soup = BeautifulSoup(pwx_text, 'html.parser')
tags = soup.find_all('a')


key_list = ['&CONTROL', '&SYSTEM', '&ELECTRONS', '&IONS', '&CELL', 'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POONTS', 'CELL_PARAMETERS', 'CONSTRAINTS', 'OCCUPATIONS', 'ATOMIC_FORCES']
last_key = 'None'
for tag in tags:
    if tag.text.strip() == '':
        continue
    if tag.text in key_list:
        last_key = tag.text
        namelist_dict[last_key] = []
    elif tag.text == 'Back to Top':
        break
    else:
        namelist_dict[last_key].append(tag.text)

# for key in namelist_dict:
#     print(key + '\n' + str(namelist_dict[key]))

tags2 = soup.find_all('table')

def table_tag_to_data(tag_text):
    regex_dict = {'command' : r'\n.+\n[A-Z]',
                 'data_type' : r'\n[A-Z]+\n',
                 'default' : r'Default:\n.*\n'}
    default = ''
    dtype = 'unknown'
    options = []
    description = '\n'.join([line.strip() for line in tag_text.split('\n') if line.strip() is not ''])
    command = ''
    val_comments_dict = {}
    for key, regex in regex_dict.items():
        match = re.findall(regex, tag_text)
        if match:
            match_text = match[0]
            if key == 'command':
                command = match_text[:-1].strip()
            if key == 'data_type':
                dtype = match_text.strip()
            if key == 'default':
                default = match_text.split(":")[1].strip()
        if dtype == 'INTEGER':
            options = int
        elif dtype == 'LOGICAL':
            options = bool
        elif dtype == 'REAL':
            options = float
        elif dtype == 'CHARACTER':
            matches = re.findall(r"'.*'", tag_text)
            if len(matches) < 2:
                options = str
            else:
                options = []
                for match in matches:
                    if match[1:-1] not in options:
                        options.append(match[1:-1])
                print(options)
                for option in options:
                    regex = option + r".*:[\s\w\.\,\(\)\-\"\:\;]+\n\n"
                    print("checking option: " + option + " with regex: " + regex)
                    matches = re.findall(regex, tag_text)
                    print(len(matches))
                    if matches:
                        comment = matches[0].split(':')[1].replace('\n', '')
                        val_comments_dict[option] = comment
    command_entry = (default, options, description)
    return command, command_entry, val_comments_dict


for tag in tags2[3:]:
    command, command_entry, val_comments_dict = table_tag_to_data(tag.text)
    command_dict[command] = command_entry
    val_comments[command] = odict(val_comments_dict)

namelist_dict.__delitem__('None')
pickle.dump(namelist_dict, open('pw.x namelist_dict.p', 'wb'))
pickle.dump(command_dict, open('pw.x command_dict.p', 'wb'))
pickle.dump(val_comments, open('pw.x val_comments.p', 'wb'))