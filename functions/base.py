# This contains a list of essential 'base' programs.
import os, sys
from biscotti import setup_logger

# Setup a logger
logger = setup_logger(__name__)
# End logging.

def list_folder(rootpath : str, extensions = None, fullpaths=False, folders_only=False):
    """
    This function only lists the rootpath folder and its immediate contents
    :param rootpath:
    :param extensions:
    :return:
    """
    if not os.path.isdir(rootpath):
        rootpath = os.path.basename(rootpath)
    if extensions is None and folders_only is False:
        files = os.listdir(rootpath)
    else:
        if type(extensions) is not list:
            extensions = [extensions]
        if folders_only:
            files = [file for file in os.listdir(rootpath) if os.path.isdir(os.path.join(rootpath, file))]
        else:
            files = [file for file in os.listdir(rootpath) if any([file.endswith(ext) for ext in extensions])]
    if fullpaths:
        files = [rootpath + os.sep + file for file in files]
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
            calcs_found.append(('jobscript', rootpath + os.sep + file))
    return calcs_found

def list_folder_deep(rootpath : str, extensions = None):
    """
    This function works as List_folder except it searches the entire tree.
    :param rootpath:
    :param extensions:
    :return:
    """
    files = []
    if extensions is None:
        for folder in os.walk(rootpath):
            files.extend([folder[0] + os.sep + file for file in os.listdir(folder[0])])
    else:
        for folder in os.walk(rootpath):
            files.extend([folder[0] + os.sep + file for file in os.listdir(folder[0]) if
                                any(file.endswith(ext) for ext in extensions)])
    return files