# This class is a wrapper for the input/output from the pp.x (pp.exe) function of Quantum Espresso
# Please see "http://www.quantum-espresso.org/wp-content/uploads/Doc/pp_user_guide.pdf"
import os, sys, re
# import numpy as np
from collections import OrderedDict as odict
import datetime as dt
# import hashlib
# BUF_SIZE = 65536

# Setup a logger
import logging
from biscotti.functions import base
from biscotti import setup_logger
logger = setup_logger(__name__, logging.INFO)

PLOT_TYPES = {'1D spherical': 0,
              '1D': 1,
              '2D': 2,
              '3D': 3,
              '2D polar' : 4}

OUTPUT_FORMATS = {'gnuplot1D': 0,
                  'contour.x': 1,
                  'plotrho': 2,
                  'XCRYSDEN2D': 3,
                  'gOpenMol': 4,
                  'XCRYSDEN3D': 5,
                  'gaussian cube': 6,
                  'gnuplot2D': 7}

class PostProcIn(object):
    def __init__(self, name=None, inputpp=None, plot=None):
        """
        :param name: Quick name for this post processing calc
        :param inputpp: this is the ordered dictionary of the input namelist
        :param plot: this is the ordered dictionary of the plot namelist
        """

        if inputpp is not None:
            self.inputpp = odict(inputpp)
        self.inputpp = inputpp
        if plot is not None:
            plot = odict(plot)
        self.plot = plot  # Plot is None if we're not plotting just extracting
        if name is None:
            if 'prefix' in self.inputpp:
                self.name = self.inputpp['prefix']
            else:
                self.name = 'Default'
        else:
            self.name = name
        # Check for prefix, filplot and plot_num. Set defaults otherwise

    def write_to_file(self, folder=None, filename = None):
        if folder is None:
            folder = os.getcwd()
        else:
            if not os.path.exists(folder):
                os.makedirs(folder)
        if filename is None:
            filename = re.sub(r'[\s\(\)]', "_", self.name)
        # Begin file write
        with open(folder + os.path.sep + filename, 'w') as newfile:
            newfile.write("! Quantum Espresso post processing (pp.x) input file.\n"
                          + '! Generated by Biscotti on ' + str(dt.datetime.now()) + '\n\n')
            # Change flags for consistency
            self.check_consistency(fix = True)
            if self.inputpp is None:
                self.inputpp = odict() # QE requires the inputpp card to exist but be empty if not used.
                namelistdict = odict([('INPUTPP',self.inputpp), ('PLOT', self.plot)])
            if self.plot is None:
                namelistdict = {'INPUTPP': self.inputpp} # Don't write the plot card if not needed
            # Write namelists
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

    def check_consistency(self, fix=False):
        inputpp_req_flags = ['prefix', 'filplot', 'plot_num']
        inputpp_req_flags_defaults = dict(zip(inputpp_req_flags, ['default', 'charge_density.ppout', 0]))
        plot_req_flags = ['iflag', 'output_format']
        plot_req_flags_defaults = dict(
            zip(plot_req_flags, [0, 0]))  # Default is 1D plot of spherical average for gnuplot

        if self.inputpp is not None:
            for flag in inputpp_req_flags:
                if flag not in self.inputpp:
                    logger.warning("Warning, pp.x flag '" + flag + "' was missing")
                    if fix:
                        self.inputpp[flag] = inputpp_req_flags_defaults[flag]

        if self.plot is not None:
            for flag in plot_req_flags:
                if flag not in self.plot:
                    logger.warning("Warning, pp.x flag '" + flag + "' was missing")
                    if fix:
                        self.inputpp[flag] = plot_req_flags_defaults[flag]

    def setup_inputpp_defaults(self, prefix=None, outdir=None, filplot=None, plot_num=1):
        """
        Creates the first few entries in the inputpp odict
        :param prefix: The short prefix for your PWSCF calculation. Default is pwscf
        :param outdir: The output directory of the data, defualt is current directory
        :param filplot: name of the output file with the requested data. Default is cool
        :param plot_num: Which type of quantity you want extracted.
        :return:
        """
        if prefix is None:
            prefix = 'pwscf'
        self.inputpp['prefix'] = prefix
        if outdir is not None:
            self.inputpp['outdir'] = outdir
            # Otherwise don't worry. it will be current dir
        if filplot is None:
            filplot = self.name + ".plot{0:d}.pp.out".format(plot_num)
        self.inputpp['filplot'] = filplot
        self.inputpp['plot_num'] = plot_num

    def extract_ILDOS_quantities(self, emin, emax=None, spin_comp=0, prefix=None, outdir=None, filplot=None):
        """
        Plot 10
        :param emin:
        :param emax:
        :param spin_comp:
        :return:
        """
        if len(self.inputpp) > 3:
            logger.warning("Warning, over-writing existing inputpp namelist with new.")
        self.inputpp = odict()
        self.setup_inputpp_defaults(prefix, outdir, filplot, plot_num=10)
        self.inputpp['emin'] = float(emin)
        if emax is not None:
            self.inputpp['emax'] = float(emax)
        self.inputpp['spin_component'] = int(spin_comp)
        logger.info("pp.x object " + self.name + " is now configured for extracting Integrated Local Density of States.")

    def setup_plot_defaults(self, quantity_file=None, plot_type='1D spherical', output_format='gnuplot1D', fileout=None, intepolation='fourier'):
        if quantity_file is None:
            if 'filplot' not in self.inputpp:
                logger.error("Must set quantity_file when filplot is not specified in inputpp")
            logger.debug("Using filplot as the input quantity file")
        else:
            self.plot['filepp(1)'] = quantity_file
        if plot_type in PLOT_TYPES:
            self.plot['iflag'] = PLOT_TYPES[plot_type]
        elif type(plot_type) is int:
            self.plot['iflag'] = plot_type
            # plot_type = # Reverse dict lookup? wtf
        else:
            self.plot['iflag'] = 0
            plot_type = '1D spherical'
        if output_format in OUTPUT_FORMATS:
            self.plot['output_format'] = OUTPUT_FORMATS[output_format]
        elif output_format == 'gnuplot':
            logger.warning("Must select gnuplot1D or gnuplot3D. Defaulting to plot type")
            if self.plot['iflag'] < 2:
                self.plot['output_format'] = OUTPUT_FORMATS['gnuplot1D']
            else:
                self.plot['output_format'] = OUTPUT_FORMATS['gnuplot2D']
        elif output_format == 'XCRYSDEN':
            logger.warning("Must select XCRYSDEN2D or XCRYSDEN3D. Defaulting to plot type")
            if self.plot['iflag'] == 2:
                self.plot['output_format'] = OUTPUT_FORMATS['XCRYSDEN2D']
            else:
                self.plot['output_format'] = OUTPUT_FORMATS['XCRYSDEN3D']
        elif type(output_format) is int:
            self.plot['output_format'] = output_format
        else:
            self.plot['output_format'] = 0

        if fileout is not None:
            self.plot['fileout'] = fileout
        self.plot['interpolation'] = intepolation

    def plot_1D(self, plot_line : list, plot_origin:list, spherical_avg=False, num_points=10,
                quantity_file=None, fileout=None, interpolation='fourier'):
        """
        This is for plot type 0 and 1
        :param plot_line: 3D vector which determines the plotting line (in alat units)
        :param plot_origin: 3D vector, origin of the line (in alat units)
        :param spherical_avg: Boolean to see if we're doing a spherical average or not
        :param num_points: number of points in the line:
                        rho(i) = rho( x0 + e1 * (i-1)/(nx-1) ), i=1, nx
        :param fileout: Output file, standard out is default
        :param interpolation: 'fourier' or 'bspline'
        :return:
        """
        if self.plot is not None and len(self.plot) > 2:
            logger.warning("Warning, over-writing existing plot namelist with new.")
        self.plot = odict()
        if spherical_avg is True:
            plot_type = '1D spherical'
        else:
            plot_type = '1D'
        self.setup_plot_defaults(quantity_file, plot_type, 'gnuplot1D', fileout, interpolation)
        for i, coord in enumerate(plot_line):
            self.plot['e1({0:d})'.format(i+1)] = coord
        for i, coord in enumerate(plot_origin):
            self.plot['x0({0:d})'.format(i+1)] = coord  # Lists start at 0, but not in Quantum Espresso baby!
        self.plot['nx'] = num_points

    def plot_2D(self, plot_line1: list, plot_line2: list, plot_origin:list, polar=False, num_x_pts=10, num_y_pts=10,
                quantity_file=None, output_format='contour.x', fileout=None, interpolation='fourier'):
        """
        This method sets the object up for a 2D plot.
        :param plot_line1: 3D vectors which determine the plotting plane (in alat units)
        :param plot_line2: e1 and e2 must be orthogonal
        :param plot_origin: 3D vector, origin of the plane (in alat units)
        :param polar: if true, plot will be a 2D polar plot on a sphere
        :param num_x_pts: Number of points in the plane
        :param num_y_pts: Number of points in the plane
                        rho(i,j) = rho( x0 + e1 * (i-1)/(nx-1)+ e2 * (j-1)/(ny-1) )
                        i=1,nx ; j=1,ny
        :param output_format: contour.x, plotrho, gnuplot or XCRYSDEN
        :param fileout:
        :param interpolation:
        :return:
        """
        if self.plot is not None and len(self.plot) > 2:
            logger.warning("Warning, over-writing existing plot namelist with new.")
        self.plot = odict()
        if polar is True:
            plot_type = '2D polar'
        else:
            plot_type = '2D'
        self.setup_plot_defaults(quantity_file, plot_type, output_format, fileout, interpolation)
        # Lists don't start at 1 for Quantum Espresso!?
        for i, coord in enumerate(plot_line1):
            self.plot['e1({0:d})'.format(i+1)] = coord
        for i, coord in enumerate(plot_line2):
            self.plot['e2({0:d})'.format(i+1)] = coord
        for i, coord in enumerate(plot_origin):
            self.plot['x0({0:d})'.format(i+1)] = coord
        self.plot['nx'] = int(num_x_pts)
        self.plot['ny'] = int(num_y_pts)

class PostProcOut(object):
    def __init__(self, inpath=None):
        if inpath is None:
            calcs = base.scan_folder_for_calcs(os.getcwd())
            for calc_type, path in calcs:
                if calc_type == 'pp.x input':
                    self.inpath = path
                    break
        else:
            self.inpath = inpath
        self.ppin = PostProcIn.import_from_file(self.inpath)