# This class is a wrapper for the input/output from the pp.x (pp.exe) function of Quantum Espresso
# Please see "http://www.quantum-espresso.org/wp-content/uploads/Doc/pp_user_guide.pdf"
import os, sys, re
import numpy as np
from collections import OrderedDict as odict
import datetime as dt
import hashlib
BUF_SIZE = 65536

sys.path.append(r"D:\Users\Chris\Documents\SivaLab\Python")

import logging
# Logging level by default
logger = logging.getLogger(__name__)
loglevel = logging.INFO
logger.setLevel(loglevel)

# Handler
console_handler = logging.StreamHandler()
console_handler.setLevel(loglevel)

#formatter
formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)

# add handler
logger.addHandler(console_handler)

class PostProcIn(object):
    def __init__(self, name='Default', inputpp=None, plot=None):
        self.name = name
        if inputpp is None:
            self.inputpp = odict()
        else:
            self.inputpp = odict(inputpp)
        if plot is None:
            self.plot = odict()
        else:
            self.plot = odict(plot)
        pass

    inputpp_flags_string = r"prefix | outdir | filplot | plot_num | spin_component | spin_component | sample_bias | kpoint | kband | lsign | spin_component | emin | emax | spin_component | spin_component | spin_component"
    inputpp_flags = inputpp_flags_string.split("|")
    plot_flags_string = r'nfile | filepp | weight | iflag | output_format | fileout | interpolation | e1 | x0 | nx | e1 | e2 | x0 | nx | ny | e1 | e2 | e3 | x0 | nx | ny | nz | radius | nx | ny'
    plot_flags = plot_flags_string.split("|")
    def write_to_file(self, folder=None, filename = None):
        pass

    @staticmethod
    def import_from_file(path):

class PostProcOut(object):
    def __init__(self):