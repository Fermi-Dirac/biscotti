import logging
logging.basicConfig(level=logging.DEBUG, format=' %(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

import re

pwxfile = r'D:\Documents\SivaLab\2016 MDA Type 2 SL\Ab-Initio\pw.x description.txt'

with open(pwxfile, 'r') as pwxfileobj:
    filelist = pwxfileobj.readlines()
    logger.debug("File read complete!")
    allflags = {}
    defaultflags = {}
    comments = {}
    namelist = ''
    variable = ''
    regexlist = [r'Namelist: \&[A-Z]+' , r'[a-zA-Z\s]+CHARACTER', r'[a-zA-Z\s]+INTEGER', r'[a-zA-Z\s]+LOGICAL', '[Back to Top]', 'Default:[\sA-Za-z0-9\']+']
    for i, line in enumerate(filelist):
        if line is '':
            pass
        else:
            for regex in regexlist:
                match = re.search(regex, line)
                if match : # is a match a regex
                    if 'Namelist' in regex:
                        logger.debug("Found new namelist! " + match.string)
                        namelist = match.string
                        allflags[namelist] = {}
                        defaultflags[namelist] = {}
                    if 'CHARACTER' in regex:
                        variable = match.string.split()[0]
                    if 'Default' in regex:
                        allflags[namelist][variable] = [match.string.split()[1]]
                        defaultflags[namelist][variable] = match.string.split()[1]
                        # TODO I was here last. ugg this sucks
