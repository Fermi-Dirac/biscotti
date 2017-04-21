# This class is a wrapper for the input/output from the pp.x (pp.exe) function of Quantum Espresso
# Please see "http://www.quantum-espresso.org/wp-content/uploads/Doc/pp_user_guide.pdf"
import os, sys, re
# import numpy as np
from collections import OrderedDict as odict
import datetime as dt
# import hashlib
# BUF_SIZE = 65536

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

class PostProcIn(object):
    def __init__(self, name=None, inputpp=None, plot=None):
        """
        Class designed to handle input files for pp.x post processing function from quantum espresso.
        From the docs:
        Purpose of pp.x: data analysis and plotting.

        The code performs two steps:

        (1) reads the output produced by pw.x, extracts and calculates
            the desired quantity/quantities (rho, V, ...)

        (2) writes the desired quantity to file in a suitable format for
            various types of plotting and various plotting programs

        The input data of this program is read from standard input
        or from file and has the following format:

        NAMELIST &INPUTPP
           containing the variables for step (1), followed by

        NAMELIST &PLOT
           containing the variables for step (2)

        The two steps can be performed independently. In order to perform
        only step (2), leave namelist &INPUTPP blank. In order to perform
        only step (1), do not specify namelist &PLOT

        Intermediate results from step 1 can be saved to disk (see
        variable filplot in &INPUTPP) and later read in step 2.
        Since the file with intermediate results is formatted, it
        can be safely transferred to a different machine. This
        also allows plotting of a linear combination (for instance,
        charge differences) by saving two intermediate files and
        combining them (see variables weight and filepp in &PLOT)

        All output quantities are in ATOMIC (RYDBERG) UNITS unless
        otherwise explicitly specified.
        :param name:
        :param inputpp:
        :param plot:
        """

        if inputpp is None:
            self.inputpp = odict()
        else:
            self.inputpp = odict(inputpp)
        if plot is None:
            self.plot = odict()
        else:
            self.plot = odict(plot)
        if name is None:
            if 'prefix' in self.inputpp:
                self.name = self.inputpp['prefix']
            else:
                self.name = 'Default'
        else:
            self.name = name

    # inputpp_flags_string = r"prefix | outdir | filplot | plot_num | spin_component | spin_component | sample_bias | kpoint | kband | lsign | spin_component | emin | emax | spin_component | spin_component | spin_component"
    # inputpp_flags = inputpp_flags_string.split("|")
    # plot_flags_string = r'nfile | filepp | weight | iflag | output_format | fileout | interpolation | e1 | x0 | nx | e1 | e2 | x0 | nx | ny | e1 | e2 | e3 | x0 | nx | ny | nz | radius | nx | ny'
    # plot_flags = plot_flags_string.split("|")
    def write_to_file(self, folder=None, filename = None):
        if folder is None:
            folder = os.getcwd()
        else:
            if not os.path.exists(folder):
                os.makedirs(folder)
        if filename is None:
            filename = self.name.replace(' ', '_') + ".ppin"
        # Begin file write
        with open(folder + os.path.sep + filename, 'w') as newfile:
            newfile.write("! Quantum Espresso post processing (pp.x) input file.\n"
                          + '! Generated by Biscotti on ' + str(dt.datetime.now()))
            # Change flags for consistency
            self.check_consistency(fix = True)
            # Write namelists
            namelistdict = {'INPUTPP' : self.inputpp, 'PLOT' : self.plot}
            for namelistkey in namelistdict:
                newfile.write(' &' + str(namelistkey) + '\n')
                for key in namelistdict[namelistkey]:
                    logger.debug("Now on key" + str(key))
                    param = "  " + '{message: <{width}}'.format(message = str(key), width = 13) # keeps nice spacing
                    value = namelistdict[namelistkey][key]
                    logger.debug("Value is " + str(value) + " of type " + str(type(value)))
                    if type(value) is str:
                        value = "'" + value + "'"
                    elif type(value) is bool:
                        if value : # is true
                            value = ".true."
                        else:
                            value = ".false."
                    # TODO support for comments
                    comment = None
                    if comment is not None:
                        newline = param + " = " + str(value) + ", ! " + str(comment) + '\n'
                    else:
                        newline = param + " = " + str(value) + ",\n"
                    newfile.write(newline)
                newfile.write('/ \n\n')
            newfile.write('\n')
        logger.info('Successfully written pp.x input file :' + filename + ' to folder:' + folder)

    @staticmethod
    def import_from_file(path):
        pass

    def check_consistency(self, fix = False):
        inputpp_req_flags = ['prefix', 'filplot', 'plot_num']
        inputpp_req_flags_defaults = dict(zip(inputpp_req_flags, ['default', 'charge_density.ppout', 0]))
        plot_req_flags = ['iflag', 'output_format']
        plot_req_flags_defaults = dict(zip(plot_req_flags, [0, 0]))  # Default is 1D plot of spherical average for gnuplot

        for flag in inputpp_req_flags:
            if flag not in self.inputpp:
                logger.warning("Warning, pp.x flag '" + flag + "' was missing")
                if fix:
                    self.inputpp[flag] = inputpp_req_flags_defaults[flag]

        for flag in plot_req_flags:
            if flag not in self.plot:
                logger.warning("Warning, pp.x flag '" + flag + "' was missing")
                if fix:
                    self.inputpp[flag] = plot_req_flags_defaults[flag]

class PostProcOut(object):
    def __init__(self):
        pass