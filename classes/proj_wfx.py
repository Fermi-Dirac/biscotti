# This class is a wrapper for the input/output from the projwfc.x function of Quantum Espresso
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

class ProjWfcIn(object):
    pass

class ProjWfcOut(object):
    pass