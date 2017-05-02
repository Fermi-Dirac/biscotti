# This contains a list of essential 'base' programs.
import os, sys

import logging
# Setup a logger
logger = logging.getLogger(__name__)
loglevel = logging.INFO  # <---- Logger Level
logger.setLevel(loglevel)
console_handler = logging.StreamHandler()
console_handler.setLevel(loglevel)
console_handler.setFormatter(logging.Formatter('%(name)s - %(levelname)s - %(message)s'))
logger.addHandler(console_handler)
# End logging.

def list_folder(rootpath : str, extensions = None):
    if not os.path.isdir(rootpath):
        rootpath = os.path.basename(rootpath)
    if extensions is None:
        files = os.listdir(rootpath)
    else:
        files = [file for file in os.listdir(rootpath) if any([file.endswith(ext) for ext in extensions])]
    return files

def scan_folder_for_calcs(rootpath):
    calcs_found = [] # List of tuples, first is the calc type, second is the full path to the file
    # Right now, this is lame and hopes the file extensions match the calc types. Lame I know
    exts = ['.in', '.out', 'slurmout.txt', '.sh']
    files = list_folder(rootpath, exts)
    for file in files:
        if file.endswith('.in'):
            calcs_found.append(('pw.x input', rootpath + os.sep + file))
        elif file.endswith('.out'):
            calcs_found.append(('pw.x output', rootpath + os.sep + file))
        elif file.endswith('slurmout.txt'):
            calcs_found.append(('slurm output', rootpath + os.sep + file))
        elif file.endswith('.sh'):
            calcs_found.append(('jobscript',rootpath + os.sep + file))
    return calcs_found