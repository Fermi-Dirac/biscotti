"""
This class is a wrapper for the input/output from pw.x (pw.exe) function of Quantum Espresso
This is the major workhorse of DFT in QE
"""
# Dependancies
import numpy as np
# Inter-package dependancies
from biscotti.classes import atoms as Biscotti
from biscotti.classes.calctime import CalcTime
# built-in libs
import re
import os
import sys
# sys.path.append(r"D:\Users\Chris\Documents\SivaLab\Python")
import datetime as dt
from glob import glob
from collections import OrderedDict as odict
import hashlib
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

try:
    import matplotlib.pyplot as plt
    has_mpl = True
    plt.style.use('fivethirtyeight')
except ImportError:
    logger.error("Cannot load matplotlib! Plotting disabled")
    has_mpl = False

# Constants
BUF_SIZE = 65536  # Hashing buffer size
default_pseudos = {"In" : [114.818, "In.pbe-dn-kjpaw_psl.0.2.2.UPF"],
                   "As" : [74.9220,  "As.pbe-n-kjpaw_psl.0.2.UPF"],
                   "Sb" : [121.6700,  "Sb.pbe-n-kjpaw_psl.0.3.1.UPF"]} # Default pseudopotentials
cardlists = ['CELL_PARAMETERS', 'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS', 'CONSTRAINTS', 'OCCUPATIONS', 'ATOMIC_FORCES']
# TODO Would be nice to have this kind of stuff in some sort of CONSTANTS library

class QECalcIn(object):
    """
    This class object is for a pw.x Quantum Espresso input file
    it is organized by a few dictionaries which control the various cards for a QE input file
    """

    def __init__(self, name = 'Default', control = None, system = None, electrons = None, ions = None, cell = None, structure = Biscotti.AtomicStructure(), pseudopots = None, kpts = None):
        """

        :param name: Name of this QE calc object
        :param control: ordered dictionary of the control namelist
        :param system: ordered dictionary of the system namelist
        :param electrons: ordered dictionary of the ELECTRONS namelist
        :param ions: ordered dictionary of the IONS namelist
        :param cell: ordered dictionary of the CELL namelist
        :param structure: atomic structure object to be submitted
        :param pseudopots: ordered dictionary of the pseudopotentials to be used
        :param kpts: string, or array representing the  kpts to be used.
        """
        self.name = name  # Simple string to name the calc. Probably the same as control[title]
        self.structure = structure  # Atomic structure data type in this calculation. This covers Atomic Species and Atomic Positions,and Cell Paraemters
        if control is None:
            self.control = odict()
        else:
            self.control = odict(control) # This is a dictionary of all non-default settings in the Control namelist
            if 'title' in self.control:
                self.name = self.control['title'] # Change name to title if we found it.
        if system is None:
            self.system = {'nat' : self.structure.totalatoms(), 'ntyp' : 1}
        else:
            self.system = odict(system) # Another such dictionary
        if electrons is None:
            self.electrons = odict()
        else:
            self.electrons = odict(electrons) # ditto
        if ions is None:
            self.ions = odict()
        else:
            self.ions = odict(ions)
        if cell is None:
            self.cell = odict()
        else:
            self.cell = odict(cell)
        if pseudopots is None:
            self.pseudopots = odict()
        else:
            self.pseudopots = odict(pseudopots)
        if kpts is None:
            self.kpts = ['automatic', np.array([1, 1, 1, 0, 0, 0])] # kpts by Type and then Thing which is what goes underneath. Likely a list, or a list of lists, or the string Gamma
        else:
            self.kpts = kpts

        self._defaultflags = {}
        self._allflags = {}
        self.namelistdict = odict([('CONTROL', self.control),
                                   ('SYSTEM', self.system),
                                   ('ELECTRONS', self.electrons),
                                   ('IONS',self.ions),
                                   ('CELL', self.cell)])

    def remove_default_flags(self):
        for namelistkey in self.namelistdict:
            for key in self.namelistdict[namelistkey]:
                if self._defaultflags[namelist][key] == namelist[key]:
                    # this is a default flag and can be removed safely
                    logger.debug("Found default flag, now popping: " + str(namelist[key]))
                    namelist.pop(key)

    def validate_calc(self, fix = False):
        # Checks the allflags dictionary to determine if the flags set do indeed exist for this calculation
        valid = True # innocent until proven guilty
        for namelistkey in self.namelistdict:
            for key in self.namelistdict[namelistkey]:
                if self.namelistdict[namelistkey][key] in _allflags[namelistkey][key]:
                    pass # validated input
                else:
                    valid = False
                    logger.debug("Found invalid flag: " + str(key) + " cannot be set to " + str(namelistdict[namelistkey]))
                    if fix:
                        self.namelistdict[namelistkey].pop(key)
                        logger.debug("Now popping invalid flag")
                    else:
                        break
        return valid

    def write_to_file(self, folder = None, filename = None):
        # Writes this QECalc object to a pw.x QE input file with the .in extension
        if folder is None:
            folder = os.getcwd()
        else:
            if not os.path.exists(folder):
                os.makedirs(folder)
        if filename is None:
            filename = self.name.replace(' ', '_') + ".in"
        # Begin file write
        with open(folder + os.path.sep + filename, 'w') as newfile:
            newfile.write("! Quantum Espresso input file for SCP pw.x DFT.\n"
                          + "! Generated by Biscotti on " + str(dt.datetime.now()) + '\n')
            # 1st, we change some flags for consistancy
            self.system['nat'] = self.structure.totalatoms()
            self.system['ntyp'] = len(self.pseudopots)
            for namelistkey in self.namelistdict:
                newfile.write(' &' + str(namelistkey) + '\n')
                for key in self.namelistdict[namelistkey]:
                    logger.debug("Now on key" + str(key))
                    param = "  " + '{message: <{width}}'.format(message = str(key), width = 13)
                    value = self.namelistdict[namelistkey][key]
                    logger.debug("Value is " + str(value) + " of type " + str(type(value)))
                    if type(value) is str:
                        value = "'" + value + "'"
                    elif type(value) is bool:
                        if value : # is true
                            value = ".true."
                        else:
                            value = ".false."
                    newline = param + " = " + str(value) + ",\n"

                    # TODO add support for comments
                    comment = None
                    if comment is not None:
                        newline = param + " = " + str(value) + ", ! " + comment + '\n'
                    newfile.write(newline)
                newfile.write('/ \n')
            # Now to do the cards since all namelists are done.
            # Write pseudopotentials
            newfile.write("ATOMIC_SPECIES\n")
            atomspec = set([atom.species for atom in self.structure.atomsdir])
            for spec in atomspec:
                if spec in self.pseudopots:
                    # Then this pseudopot was listed
                    newfile.write(str(spec) + "  " + str(self.pseudopots[spec][0]) + "  " + self.pseudopots[spec][1] + '\n')
                else:
                    # This species did not have a pseudopot listed
                    logger.error("Warning, Pseudopotential not specified for " + spec)
                    newfile.write(str(spec) + "  " + str(default_pseudos[spec][0]) + "  " + default_pseudos[spec][1] + '\n')
            # Cell Parameters
            newfile.write('CELL_PARAMETERS angstrom\n' \
                             + '  '.join([str(val) for val in self.structure.latticeA]) + '\n' \
                             + '  '.join([str(val) for val in self.structure.latticeB]) + '\n' \
                             + '  '.join([str(val) for val in self.structure.latticeC]) + '\n')

            # Write atomic positions
            newfile.write('ATOMIC_POSITIONS crystal\n' + '\n'.join( [atom.species
                                                                    + "  " + str(atom.x)
                                                                    + "  " + str(atom.y)
                                                                    + "  " + str(atom.z) for atom in self.structure.atomsdir]) \
                         + '\n')
            kptsstring = ''
            if type(self.kpts[1]) is str:
                kptsstring = self.kpts[1]
            if type(self.kpts[1]) is np.array:
                for val in self.kpts[1]:
                    kptsstring +=  "%d" % val + " "
            newfile.write("K_POINTS " + str(self.kpts[0]) + '\n' + kptsstring)
            logger.info('File writing complete!')
        return newfile

    @staticmethod
    def import_from_file(path: str):
        # This imports a QEcalcin object from an input file
        logger.info("Now parsing input file" + path)
        # Current no support for comments!!
        namelistdict = odict()

        # First populate the dictionaries
        if not os.path.exists(path):
            logger.error("Error! file does not exist! checked :" + path)
            return QECalcIn()
        with open(path, 'r') as fileobj:
            filestring = fileobj.read()
        with open(path, 'r') as fileobj:
            filelist = fileobj.readlines() # legacy, irritating i know. Could maybe replace with a .split('\n')?
        namelist = None

        cardlists_short = [card[:8] for card in cardlists] # truncate to first 8 char for quick compare
        for line in filelist:
            logger.debug("Now on line >" + line.strip() + "<")
            if line.strip() == "":
                firstchar = '!'
                continue # Empty line
            else:
                firstchar = line.strip()[0]
            if firstchar is '!':
                continue # commented line
            elif firstchar is '&':
                namelist = line.strip()[1:].upper() # skip the '&' symbol, demand uppercase convention for QE files
                logger.debug("New namelist: " + namelist)
                namelistdict[namelist] = odict()
            elif firstchar is '/':
                key = None
            elif line.strip()[0:8] in cardlists_short:
                logger.info("Found all namelist dictionaries, now looking for cards")
                break # Stop looking for dictionaries
            else:
                logger.debug("Found a key value pair at line " + line.strip())
                # found a key-value pair
                # syntax is KEY = VALUE , ! comment
                key = line.split('=')[0].strip()
                value = line.split('=')[1].split(',')[0].strip()
                logger.debug("Untyped value is: " + str(value) + " with type " + str(type(value)))
                if value[0] == "'" :
                    value = value[1:-1] # strip off the ' symbol
                elif value[0] == ".":
                    if value == '.true.':
                        value = True
                    else:
                        value = False
                elif '.' in value or 'e' in value: # must be float
                    for exp in ["D","d","E"]:
                        value = value.replace(exp, "e")
                    value = float(value)
                else: # must be int
                    value = int(value)  # must be an integer
                logger.debug("New key value pair is: (" + key + " , " + str(value) + ") of type " + str(type(value)))
                namelistdict[namelist][key] = value
        logger.debug("All input file Dictionaries found!")
        # Now dictionaries are populated, need kpts and atomic species and pseduots
        # Find Pseudopotentials
        regex_pseudo = r'ATOMIC_SPECIES.*(\n[ ]*[A-Z][a-z][0-9\. ].+)+'
        pseudopots = odict()
        pseudo_match = re.search(regex_pseudo, filestring)
        if not pseudo_match:
            logger.error("No pseudo-potentials found!")
        else:
            pseudo_strings = pseudo_match.group().split('\n')[1:] # strip first entry which is literal ATOMIC SPECIES
            for pstr in pseudo_strings:
                linesplit = pstr.split()
                pseudopots[linesplit[0]] = [linesplit[1], linesplit[2]]
            logger.info("Loaded " + str(len(pseudopots)) + ' pseudopotentials')
        # now for kpts
        regex_kpts = r''
        found = False
        kptstring = ''
        kptstype = ''
        for line in filelist:
            if line[0:8] == 'K_POINTS':
                kptstype = line[9:].strip()
                found = True
            elif found:
                kptstring += line
            else:
                pass
        # Check existing namelists
        if 'calculation' in namelistdict['CONTROL']:
            calctype = namelistdict['CONTROL']['calculation']
            if calctype not in ['scf', 'nscf', 'bands', 'relax', 'md', 'vc-relax', 'vc-md']:
                logger.info("flag 'calculation' not set to accepted list. Resetting to 'scf'")
                calctype = 'scf'
                namelistdict['CONTROL']['calculation'] = 'scf'
        else:
            calctype = 'scf'
        control = namelistdict['CONTROL']
        system = namelistdict['SYSTEM']
        electrons = namelistdict['ELECTRONS']
        if calctype in ['relax', 'md', 'vc-relax', 'vc-md']:
            ions = namelistdict['IONS']
        else:
            ions = None
        if calctype in ['vc-relax', 'vc-md']:
            cell = namelistdict['CELL']
        else:
            cell = None
        return QECalcIn(os.path.split(path)[1], control, system, electrons, ions, cell
                        , Biscotti.AtomicStructure.from_QEinput(path)
                        , pseudopots
                        , [kptstype, kptstring])

    def write_slurm_jobscript(self, folder=None, infile = None, slurm_dict = None,
                              send_email=True, email_addr='', report=True):
        """
        Generate a SLURM jobscript file with associated flags. Also supports post-job reporting and email notification
        :param folder:
        :param jobscript_name:
        :param infile:
        :param num_cores:
        :param partition:
        :param email_start:
        :param email_end:
        :param email_addr:
        :param report:
        :return:
        """
        if folder is None:
            folder = os.getcwd()
        else:
            if not os.path.exists(folder):
                os.makedirs(folder)
        if infile is None:
            infile = [file for file in os.listdir() if os.path.splitext(file)[1] == '.in'][0]
        if slurm_dict is None:
            try:
                max_seconds = self.control['max_seconds'] + 60 * 5  # add 5minutes so QE cancels before SLURM
            except:
                max_seconds = 60 * 60  # Default is 1 hour
            time = "%02d:%02d:%02d" % (max_seconds / (60 * 60), (max_seconds / (60)) % 60, max_seconds % 60)
            flags = ['job-name', 'output', 'ntasks', 'time', 'mem-per-cpu']
            values = [self.name[-8:], 'slurmout.txt', '16', time, 'MaxMemPerNode']
            slurm_dict = odict(zip(flags, values))
        # Set SLURM parameters

        slurm_template = "#!/bin/env python3\n"
        for key, value in slurm_dict.items():
            slurm_template += '#SBATCH --' + key + "=" + value + "\n"
        slurm_template += ("import subprocess as subpr\n" + "import datetime as dt\n" + "import os, sys\n")
        email_start_template =(
            "\n\n# Email at start"
            "\nbody = 'Your QE calculation %(jobname)s began on ' + str(dt.datetime.now()) + '\\nThe full execution path is: \\n' + os.path.abspath('')"
            "\nemail.send_mail('%(email_addr)s', 'QE calculation %(jobname)s has started', body)"
        )
        pw_x_template = (
            "\n\n# Begin pw.x call"
            "\nsubpr.call('mpirun -np %(num_cores)s pw.x -i %(infile)s > %(infile)s.out', shell=True)"
        )
        email_end_template = (
            "\n\n# Email at end"
            "\nbody_end = 'Your QE calculation %(jobname)s ended on ' + str(dt.datetime.now()) + '\\nThe full execution path is: \\n' + os.path.abspath('')"
        )
        email_report_template =(
            "\ncalcout = qecalc.QECalcOut.import_from_file('%(infile)s.out', '%(infile)s')"
            "\ncalcout.make_report(reportname='%(jobname)s report.png')"
            "\nemail.send_mail('%(email_addr)s', 'QE calculation %(jobname)s has ended', body_end, ['%(jobname)s report.png'])"
        )
        email_noreport_template = "\nemail.send_mail('%(email_addr)s', 'QE calculation %(jobname)s has ended', body_end)"

        with open(folder + os.path.sep + 'slurm_jobscript.sh', 'w', newline='\n') as fileobj:

            fileobj.write(slurm_template % locals())

            # Only load biscotti if emailing or reporting is on
            if send_email:
                fileobj.write("\nfrom biscotti.reporting import email\n")
            if report:
                fileobj.write("\nfrom biscotti.classes import qecalc")

            # Email at start
            if send_email:
                fileobj.write(email_start_template % locals())

            # pw.x parameters
            fileobj.write(pw_x_template % locals())

            # End email and reporting
            if send_email:
                fileobj.write(email_end_template % locals())
                if report:
                    fileobj.write(email_report_template % locals())
                else:
                    fileobj.write(email_noreport_template % locals())
            fileobj.write("\nsubpr.call('python script complete!', shell=True)")
            fileobj.write("\n\n")

class QECalcOut(object):
    """ This class holds the results of a pw.x QE calculation. This one is very much under development"""
    def __init__(self, outpath = None, inpath = None, refenergies=None, relax_list = None, pressure_list = None, jobstatus = 'unknown' ):
        """

        :param outpath: filepath for the quantum espresso output file
        :param inpath:  filepath to the QE input file which generated this output file. Scans local dir if empty
        :param relax_list: List of ion steps which is itself a list of electron steps as energies.
        :param jobstatus: string indicating how this calculation terminated
        """
        # Things out from a QE calculation:
        # Path Derived Variables
        if outpath is None:
            outfiles = [file for file in os.listdir(os.getcwd()) if file.split('.')[-1] == 'out']
            if outfiles: # is not empty
                outpath = outfiles[0]
            else:
                logger.error("No output files found!")
        self.path = outpath
        md5hash = hashlib.md5()
        with open(outpath, 'rb') as fileobj:
            while True:
                data = fileobj.read(BUF_SIZE)
                if not data:
                    break
                md5hash.update(data)
        self.ID = str("{0}".format(md5hash.hexdigest())) # HASH of Path's file
        self.filename = os.path.basename(outpath)
        self.folder = os.path.dirname(outpath)

        # Calctime derived attributes
        self.calctime = CalcTime.get_time(outpath)
        self.startdt = self.calctime.startdt
        self.enddt = self.calctime.enddt
        self.totalsec = (self.enddt - self.startdt).total_seconds()

        #QECalcIn derived variables
        if inpath is None:
            outdir = os.path.dirname(outpath)
            logger.info("No input file specified, checking outfile folder: " + outdir)
            for file in os.listdir(outdir):
                logger.debug("Maybe this file? " + file )
                if file.split(".")[-1] == 'in':
                    inpath = os.path.join(outdir,file)
                    break
        if inpath is None: # still..
            logger.error("No input files found!")
        else:
            self.qecalcin = QECalcIn.import_from_file(inpath)
            self.initialstructure = self.qecalcin.structure
            self.name = self.qecalcin.name

        # Structurelist derived attributes
        if self.qecalcin.control['calculation'] == 'scf':
            self.structurelist = []
            self.finalstructure = self.initialstructure
        else:
            self.structurelist = Biscotti.AtomicStructure.from_QEOutput(outpath)
            if len(self.structurelist) > 0 : # is not empty or None
                self.finalstructure = self.structurelist[-1]
            else:
                self.finalstructure = self.initialstructure
        self.volumeList = None

        # Reference Energies, These are calculations of 'free' atoms to create total free energy, similar to VASP
        if refenergies is None:
            # These are the energies of these species in 'free' cells, e.g. huge boxes
            refenergies = {'As': -175.65034951, 'In': -410.56656045, 'Sb': -347.3619941658}
        self.refenergy = 0
        for atom in self.initialstructure.atoms:
            if atom.species in refenergies:
                self.refenergy += refenergies[atom.species]
            else:
                logger.error("Missing ref energy for atomic species: " + atom.species)
                pass
        # Energy convergances
        if relax_list is None or not relax_list: # is empty
            logger.error("No relaxation data or energies")
            relax_list = [[0]]
        self.relax_list = relax_list
        self.relax_list_ev = []
        self.relax_list_free = []
        self.relax_list_free_ev= []
        for ionstep in relax_list:
            self.relax_list_ev.append([estep*13.605698066 for estep in ionstep])
            self.relax_list_free.append([estep - self.refenergy for estep in ionstep])
            self.relax_list_free_ev.append([(estep - self.refenergy)*13.605698066 for estep in ionstep])

        # Final energies, stress and pressure
        self.final_energy = relax_list[-1][-1]
        self.formation_energy = self.relax_list_free_ev[-1][-1] # in eV
        self.ion_steps_count = len(self.relax_list)
        self.pressure_list = pressure_list

        # Job Status
        self.jobstatus = jobstatus
        if jobstatus == 'complete':
            self.jobcomplete = True
        else:
            self.jobcomplete = False

    @staticmethod
    def import_from_file(outpath=None, inpath = None):
        # This function imports a QECalcOut object from an output file from a pw.x calculation
        # very much work in progress
        if outpath is None:
            outfiles = [file for file in os.listdir(os.getcwd()) if file.split('.')[-1] == 'out']
            if outfiles: # is not empty
                outpath = outfiles[0]
            else:
                logger.error("No output files found!")
        path = outpath # legacy
        RelaxationSteps = []
        PressureList = []
        volumeList = []
        ionstep = 0
        electronstep = 99
        completestring = 'Unknown'
        logger.info("--Now parsing file " + path)
        verbose = False

        with open(path) as fileobj:
            for line in fileobj:
                regexlist = {'iteration': r'iteration #[\s]+[0-9]+',
                             'total energy': r'[\s]*total energy[\s]*=[\s]*[\-\.0-9]*',
                             'starttime': r'[0-9]+[A-Za-z]{3}[0-9]{4} at[\s]+[0-9\s]+:[0-9\s]+:[\s]*[0-9]+',
                             'stoptime': r'[0-9]+:[0-9\s]+:[0-9\s]+[A-Za-z]{3}[0-9]{4}',
                             'done': r'JOB DONE.',
                             'k-points': r'k points=[\s]*[0-9]*',
                             'exitcode': r'Exit code:[\s0-9]*',
                             'pressure': r'\(kbar\)[\s]*P=[\s]*[0-9\.\-]*',
                             'volume': r'unit-cell volume[\s]*=[\s]*[0-9\.]*'
                             }
                for regex in regexlist:
                    searchresult = re.findall(regexlist[regex], line)
                    if searchresult:  # isn't empty..
                        logger.debug("We have a match! :" + (searchresult[0]))
                        if 'iteration' in regex:
                            newstep = int(searchresult[0].split('#')[1])
                            logger.debug('electron step #' + str(newstep) + " last step was " + str(electronstep))
                            if electronstep > newstep:  # reset to new ionic step
                                ionstep += 1
                                logger.debug("looks like a new ion step" + str(ionstep))
                            electronstep = newstep
                        elif 'total energy' in regex:  # total energy
                            newenergy = float(searchresult[0].split('=')[1])
                            logger.debug("\tthis step energy : " + str(newenergy))
                            logger.debug("Relax steps size: " + str(len(RelaxationSteps)) + " ionstep: " + str(ionstep))
                            if len(RelaxationSteps) < ionstep:
                                RelaxationSteps.append([])
                            RelaxationSteps[ionstep - 1].append(newenergy)  # ionstep starts at 1, index starts at 0
                        elif 'pressure' in regex:
                            logger.debug("Pressure of " + searchresult[0])
                            PressureList.append(float(searchresult[0].split('=')[1]))
                            pass
                        elif 'volume' in regex:
                            volumeList.append(float(searchresult[0].split('=')[1]))
                        elif 'exitcode' in regex:
                            logger.debug("Timeout detected! This calculation did not complete! " + searchresult[0])
                            completestring = searchresult[0].strip()
                            jobComplete = False
                        elif 'done' in regex:
                            jobComplete = True
                            completestring = 'completed'
                            # break  #don't check any more regex if empty
        logger.info("Relaxation Steps length : " + str(len(RelaxationSteps)))
        if verbose:
            print("\n\n--This calculation " + completestring)
            print("number of ionic steps is " + str(ionstep))
            print("Here's the giant Relax steps" + str(RelaxationSteps))
            print("Final Energy: " + str(RelaxationSteps[-1][-1]))
        return QECalcOut(outpath=path, inpath=inpath, relax_list=list(RelaxationSteps), pressure_list = PressureList, jobstatus= completestring)

    def calc_overview_dict(self):
        refenergy = self.refenergy
        numatoms = self.initialstructure.totalatoms()
        if numatoms == 0:
            numatoms = -1
        dE = -1
        if len(self.relax_list[-1]) > 1:
            dE = self.relax_list[-1][-1] - self.relax_list[-1][-2]
        iondE = -1
        if len(self.relax_list) > 1:
            iondE = self.relax_list[-1][-1] - self.relax_list[-2][-1]
        keys = ['ID', 'Filename', 'Folder','Title','Calc Type',
                'Final Energy (Ry)','Last electron step dE (Ry)','Last ion step dE (Ry)','Final Free Energy (Ry)',
                   'Final Free Energy (eV)',
                   'Final Free Energy (eV/atom)',
                   'Number of Atoms',
                   'Cutoff (Ry)',
                   'Total # of K-points',
                   'Job Complete?',
                   'Calc time (hr)',
                   'Start Date-Time',
                   'End Date-Time',
                   'Final Pressure(kbar)',
                   'Initial Volume (A^3)',
                   'Final Volume (A^3)',
                'Final Pressure (kbar)'
                   ]
        values = [self.ID,
                  self.filename,
                  self.folder,
                  self.qecalcin.name,
                  self.qecalcin.control['calculation'],
                  self.final_energy,
                  dE,
                  iondE,
                  self.final_energy - refenergy,
                  (self.final_energy - refenergy)*13.605698066,
                  (self.final_energy-refenergy)*13.605698066/numatoms,
                  numatoms,
                  self.qecalcin.system['ecutwfc'],
                  self.qecalcin.kpts,
                  self.jobstatus,
                  self.totalsec/(60*60),
                  self.startdt,
                  self.enddt,
                  0,
                  self.initialstructure.totalvol(),
                  self.finalstructure.totalvol(),
                  self.pressure_list[-1]
                  ]
        return odict(zip(keys, values))

    def calc_overview_string(self, add_headers = True, delim = '\t', transpose = False):
        returnstring = ''
        summary_dict = self.calc_overview_dict()
        if transpose:
            for key, value in summary_dict.items():
                if add_headers:
                    returnstring += key + ":" + delim
                returnstring += str(value) + '\n\r'
        else:
            if add_headers:
                for key in summary_dict:
                    returnstring += key + delim
                returnstring += '\n'
            for value in summary_dict.values():
                returnstring += str(value) + delim
        return returnstring

    def __str__(self):
        return self.calc_overview_string()

    def scale_energies(self, free_energy=False, unit='Ry'):
        scaledenergies = []

        # scale by unit selected, QE outputs in Rydbergs for some God-forsaken reason
        if unit == 'eV':
            scale = 13.605698066
        elif unit == 'eV/atom':
            scale = 13.605698066/self.initialstructure.totalatoms()
        elif unit == 'Ry/atom':
            scale = 1/self.initialstructure.totalatoms()
        else:
            scale = 1

        # scale Free energy vs Total Energy
        if free_energy:
            refenergy = self.refenergy
        else:
            refenergy = 0

        # generate scaled list
        for ionstep in self.relax_list:
            scaledenergies.append([(estep-refenergy)*scale for estep in ionstep])
        return scaledenergies

    def plot_ion_energies(self, axes=None, free_energy=False, unit='Ry', delta=False):
        if has_mpl is False:
            logger.error("Matplotlib dependency not loaded. Plotting disabled!")
            return None
        if axes is None:
            axes = plt.gca() # this is the 'matplotlib' axes class. NOT actual axes of a plot. Its a 'plot' thingy within a figure
        energies_plotted = self.scale_energies(free_energy, unit)
        ylabel = ('Free' if free_energy else 'Total') + ' Energy (' + ('delta ' if delta else '') + unit + ')'
        # Create Ion Plot
        axes.set_xlabel('Ion Step #')
        axes.set_ylabel(ylabel)
        energies = [ionstep[-1] for ionstep in energies_plotted]
        if delta:
            delta_e = [0]
            for i in range(len(energies)-1):
                delta_e.append(energies[i] - energies[i+1])
            energies = delta_e
        plotout = axes.plot(energies, linestyle = '-', marker='o', color='r')
        return plotout

    def plot_e_energies(self, axes=None, free_energy=False, unit='Ry'):
        if has_mpl is False:
            logger.error("Matplotlib dependency not loaded. Plotting disabled!")
            return None
        if axes is None:
            axes = plt.gca() # this is the 'matplotlib' axes class. NOT actual axes of a plot. Its a 'plot' thingy within a figure
        energies_plotted = self.scale_energies(free_energy, unit)
        ylabel = ('Free' if free_energy else 'Total') + ' Energy (' + unit + ')'
        # Create electron plot
        axes.set_xlabel('Electron Step #')
        axes.set_ylabel(ylabel)
        i=0
        for ionstep in energies_plotted:
            steps = [step+i for step in range(len(ionstep))]
            axes.plot(steps, ionstep, linestyle = '-', marker='o')
            i += len(ionstep)
        return axes

    def plot_summary_table(self, axes=None, headers_list = None):
        if has_mpl is False:
            logger.error("Matplotlib dependency not loaded. Plotting disabled!")
            return None
        # Get dictionary of the summary
        summary_dict = self.calc_overview_dict()
        if headers_list is None:
            headers_list = summary_dict.keys()

        # initialize and fill the table
        values = []
        final_headers = []
        for header in headers_list:
            if header in summary_dict:
                final_headers.append(header)
                newval = summary_dict[header]
                if type(newval) is float:
                    newval = '%0.4g' % newval
                values.append(newval)

        # make 'axes' for matplotlib
        if axes is None:
            axes = plt.gca()
        axes.axis('off')
        table = axes.table(cellText=[values], colLabels=final_headers, loc='center')
        #table.auto_set_font_size(False)
        table.set_fontsize(12)
        table.scale(1,5)
        table_props = table.properties()
        table_cells = table_props['child_artists']
        for cell in table_cells:
            #cell._text.set_fontsize(66)
            cell.set_height(0.4)
            #cell._text.set_color('blue')
            pass
        return table

    def plot_report(self, reportfig = None):
        if has_mpl is False:
            logger.error("Matplotlib dependency not loaded. Plotting disabled!")
            return None
        if reportfig is None:
            reportfig = plt.figure()
        ion_axes = reportfig.add_subplot(2,1,1)
        self.plot_ion_energies(ion_axes, free_energy=True, unit='eV/atom')
        # electron axes
        e_axes = reportfig.add_subplot(2,1,2)
        self.plot_e_energies(e_axes, free_energy=True, unit='eV/atom')
        # Rescale
        pl_size = reportfig.get_size_inches()
        scale = 1.5
        reportfig.set_size_inches(pl_size[0]*scale, pl_size[1]*scale)
        return reportfig

    def make_report(self, outputdir=None, reportname = None):
        if has_mpl is False:
            logger.error("Matplotlib dependency not loaded. Plotting disabled!")
            return None
        if outputdir is None:
            outputdir = os.path.abspath('')
        if reportname is None:
            reportname = self.initialstructure.name + " report.png"
        reportfig = self.plot_report()
        reportfig.savefig(outputdir + os.path.sep + reportname)
        return outputdir + os.path.sep + reportname

def makeSummaryFile(rootpath=None, outputfile=None, delim = '\t') -> object:

    """
    # This function finds all of the .out files from a path and subdirectories and then calls calcoverview-to-string
    # this gives a spreadsheet-like answer to a batch calculation
    :param rootpath: The absolute path to the folder which contains all subfoldeers and calculation files
    :param outputfile: name of the text file where delimited data will be stored
    :param delim: character which is used as a delimiter between fields.
    :return: none
    """
    if rootpath is None:
        rootpath = os.getcwd()
    logger.info("Creating summary file for root path: " + rootpath)
    # alloutputfiles = [y for x in os.walk(rootpath) for y in glob(os.path.join(x[0], '*.out'))] obfuscated code
    alloutputfiles = []
    for folder in os.walk(rootpath):
        for file in glob(os.path.join(folder[0], '*.out')):
            alloutputfiles.append(file)
    if outputfile is None:
        outputfile = os.path.join(rootpath, 'Summary for ' + os.path.split(rootpath)[1] + '.tsv')

    with open(outputfile, 'w') as fileobj:
        for index, calcfile in enumerate(alloutputfiles):
            logger.info("Now on calc file at: " + calcfile)
            calcresult = QECalcOut.import_from_file(calcfile)
            add_headers = index == 0
            writestring = calcresult.calc_overview_string(add_headers= add_headers, delim = delim) + '\n'
            fileobj.write(writestring)

"""
---------------------------------------------------
Begin Legacy Code. :(
---------------------------------------------------
"""

def writeConfig(outputFileString, calculationDict):
    # Legacy before QECalcIn class
    with open (outputFileString, 'w') as fileobj :  # write a new output file, here we go!
        fileobj.write('! Quantum esspresso input file created using Biscotti\n\n')
        namelists = ['CONTROL','SYSTEM','ELECTRONS','IONS','CELL'] # print them in this order
        for namelist in namelists :
            fileobj.write('&'+namelist + '\n')
            for command in calculationDict[namelist]:
                commandstring = "\t" + command + ' = ' + calculationDict[namelist][command][0] + ', !' + calculationDict[namelist][command][1] + '\n' #prints the command, equates to the value, and appends any comments
                print (commandstring)
                fileobj.write(commandstring)
                # end all commands for this namelist
            fileobj.write('/ \n')  # end this namelist
        print ("File created!")

def calcOverview(outputfile, verbose = False, GetDOS = False):
    # Legacy before QECalcOut class
    # Still somewhat useful as it gets more data

    startson_dtformat = '%d%b%Y at %H:%M:%S'
    ends_dtformat = '%H:%M:%S  %d%b%Y'
    # This quickly looks through a QE output file and generates a calculational summary.
    # Outputs:
    # Ionic parts, a list of lists of ionic energies. for a static calc, just a list of length 1 of all electornic steps
    RelaxationSteps = []
    #   Final unit cell pressure
    PressureList = []
    #   Final Free Energy eV, and Energy/atom and total formation energy
    #   atomic structures. Lattice vectors and atomic positions. List of these for each ionic step.
    atomicStructure = []
    #   Final Structure Volume, lattice vector lengths
    volumeList = []
    #   Number of atoms
    #   Calculation Time (seconds?)
    #   Execution Date
    startdt = dt.datetime.today()
    enddt = dt.datetime.today()
    #   Density of States Optional
    #   K-points
    #   Job completion
    # This summary dictionary will have each ionic step as an array element.
    ionstep = 0
    electronstep = 99
    jobComplete = False
    kpts = 0
    cutoff = 0
    if verbose : print ("\n--Now parsing file " + outputfile)
    with open(outputfile) as fileobj:
        for line in fileobj :
            regexlist = {'iteration' : r'iteration #[\s]+[0-9]+',
                         'total energy' : r'[\s]*total energy[\s]*=[\s]*[\-\.0-9]*',
                         'cutoff' : r'kinetic-energy cutoff[\s]*=[\s]*[\-\.0-9]*',
                         'starttime': r'[0-9]+[A-Za-z]{3}[0-9]{4} at[\s]+[0-9\s]+:[0-9\s]+:[\s]*[0-9]+',
                         'stoptime': r'[0-9]+:[0-9\s]+:[0-9\s]+[A-Za-z]{3}[0-9]{4}',
                         'done' : r'JOB DONE.',
                         'Cell Parameters' : r'CELL_PARAMETERS',
                         'k-points' : r'k points=[\s]*[0-9]*',
                         'exitcode' : r'Exit code:[\s0-9]*',
                         'pressure' : r'\(kbar\)[\s]*P=[\s]*[0-9\.\-]*',
                         'volume' : r'unit-cell volume[\s]*=[\s]*[0-9\.]*'
                         }
            for regex in regexlist :
                searchresult = re.findall(regexlist[regex],line)
                if searchresult :# isn't empty..
                    if verbose: print ("We have a match! :" +(searchresult[0]))
                    if 'iteration' in regex:
                        newstep = int(searchresult[0].split('#')[1])
                        if verbose: print ('electron step #' + str(newstep) + " last step was " + str(electronstep))
                        if electronstep > newstep : #reset to new ionic step
                            ionstep += 1
                            if verbose: print ("looks like a new ion step" + str(ionstep))
                            RelaxationSteps.append([])
                        electronstep = newstep
                    elif 'total energy' in regex : #total energy
                        newenergy = float(searchresult[0].split('=')[1])
                        if verbose: print ("\tthis step energy : " + str(newenergy))
                        RelaxationSteps[ionstep-1].append(newenergy)  #ionstep starts at 1, index starts at 0
                    elif 'Cell Parameters' in regex:
                        lattice = []
                        atoms = []
                        #fileobj.readline()  #read the rest of this line
                        for i in range(3):
                            lattice.append(fileobj.readline().split())
                            if verbose : print (lattice)
                        fileobj.readline()
                        fileobj.readline()  #read the empty line and Atomic_Positions
                        while True: # is not null
                            newatom = fileobj.readline().split()
                            if verbose : print ("new atom: " + str(newatom))
                            if not newatom or newatom[0] == 'End':  #if empty
                                break
                            atoms.append(newatom)
                        atomicStructure = [lattice,atoms]
                    elif 'pressure' in regex:
                        if verbose: print("Pressure of " + searchresult[0])
                        PressureList.append(float(searchresult[0].split('=')[1]))
                        pass
                    elif 'volume' in regex:
                        volumeList.append(float(searchresult[0].split('=')[1]))
                    elif 'cutoff' in regex:
                        cutoff = float(searchresult[0].split('=')[1]) # in Ry
                    elif 'k-points' in regex:
                        kpts = float(searchresult[0].split('=')[1]) # total # of kpts
                    elif 'starttime' in regex:
                        dtsplit = searchresult[0].split(' at ')
                        dtstring = dtsplit[0] + ' at ' + dtsplit[1].replace(' ', '0')
                        startdt = dt.datetime.strptime(dtstring, startson_dtformat)
                    elif 'stoptime' in regex:
                        dtsplit = searchresult[0].split('  ')
                        dtstring = dtsplit[0].replace(' ', '0') + '  ' + dtsplit[1]
                        enddt = dt.datetime.strptime(dtstring,ends_dtformat)
                    elif 'exitcode' in regex:
                        if verbose : print ("Timeout detected! This calculation did not complete! "+ searchresult[0])
                        jobComplete = False
                    elif 'done' in regex:
                        jobComplete = True
                # break  #don't check any more regex if empty
    if verbose :
        if jobComplete:
            completestring = "completed"
        else:
            completestring = "did not complete"
        print ("\n\n--This calculation " + completestring)
        print ("The cutoff is " + str(cutoff))
        print ("number of ionic steps is " + str(ionstep))
        print ("Here's the giant Relax steps" + str(RelaxationSteps))
        print("Final Energy: " + str(RelaxationSteps[-1][-1]))
    if not PressureList : # is empty
        PressureList = [0]
    if not volumeList : # is empty
        volumeList = [0]
    return {'RelaxationSteps': RelaxationSteps,
            'AtomicStructure' : atomicStructure,
            'Total # of K-points' : kpts,
            'Cutoff (Ry)' : cutoff,
            'Job Complete?': jobComplete,
            'starttime': startdt,
            'endtime': enddt,
            'pressure': PressureList,
            'volume': volumeList}

def AtomCount(atomicStructure):
    # Legacy since creation of AtomicStructure class
    species = []
    counts = []
    output = []
    if not atomicStructure:  #if empty
        return output
    for atom in atomicStructure[1]:
        if atom[0] in species:
            counts[species.index(atom[0])]+=1
        else:
            species.append(atom[0])
            counts.append(1)
    for index, element in enumerate(species):
        output.append([element,counts[index]])
    return output

def calcOverview_to_string(calcOverviewDict, refenergies = None):
    # Legacy before QECalcOut class
    if refenergies is None:
        refenergies = {'As':-175.65034951, 'In':-410.56656045, 'Sb':-347.3619941658}
    numatoms = 0
    refenergy = 0
    strings = []
    for species in AtomCount(calcOverviewDict['AtomicStructure']):
        refenergy += refenergies[species[0]]*species[1]
        numatoms += species[1]
    if numatoms < 1:
        print('Error, no atoms?!')
        numatoms=1
    try:
        totalenergy = calcOverviewDict['RelaxationSteps'][-1][-1]
    except (Exception) as exc:
        totalenergy = 0
    headers = ['Final Energy (Ry)',
               'Final Free Energy (Ry)',
               'Final Free Energy (eV/atom)',
               'Number of Atoms',
               'Cutoff (Ry)',
               'Total # of K-points',
               'Job Complete?',
               'Calc time (min)',
               'Start Date-Time',
               'End Date-Time',
               'Final Pressure(kbar)',
               'Initial Volume (A^3)',
               'Final Volume (A^3)',
               ]
    values = [totalenergy,
              totalenergy-refenergy,
              (totalenergy - refenergy)*13.605698066/numatoms,
              numatoms,calcOverviewDict['Cutoff (Ry)'],
              calcOverviewDict['Total # of K-points'],
              calcOverviewDict['Job Complete?'],
              (calcOverviewDict['endtime'] - calcOverviewDict['starttime']).total_seconds()/60,
              calcOverviewDict['starttime'].strftime('%m/%d/%Y %H:%M:%S'),
              calcOverviewDict['endtime'].strftime('%m/%d/%Y %H:%M:%S'),
              calcOverviewDict['pressure'][-1],
              calcOverviewDict['volume'][0]/(1.8897**3),
              calcOverviewDict['volume'][-1]/(1.8897**3)
              ]
    for value in values:
        strings.append(str(value))
    return headers,strings

def parseConfig(fname):
    # Legacy before QECalcOut class
    with open(fname,'r') as fileobj:
        filestring = fileobj.readlines()  #change scope to ensure fileobj is closed
    config = [line.split() for line in filestring]  #split over whitespace, returns a list of list.
    namelistsarray = [wordlist[0][1:] for wordlist in config if '&' in ''.join(wordlist) and ''.join(wordlist)[0]=='&']
    # this iterative searches the list for &'s and puts them in the namelist, getting a list of the special words
    namelistsarray.append("None")  #add on a dummy index for useless crap
    Dictionary = {keyword:{} for keyword in namelistsarray}   #create a dictionary using {}, which contains an array of empty dictionries
    namelist = 'None'  #set first dictionary index to junk, the first namelist is junk
    for lines in config:
        lineAsString_NoWhite = ''.join(lines)  #undo the split, no spaces.
        # Setting the dictionary index
        if len(lineAsString_NoWhite)>1:    # catching zero-length lines
            if lineAsString_NoWhite[0] == '/':
                namelist = "None"  # close a block by sending all data to junk index
            if lineAsString_NoWhite[0]=='&':
                namelist = lines[0][1:]  # set proper index if we found the right one
        # Dict index set, lets populate
        if "=" in lines:   # then we found a command parameter
            if '!' not in lines[0]:  # check for commented line
              Dictionary[namelist][lines[0]]=[lines[2]," ".join(lines[4:])]  # Doesn't check for !; commends
        else:    # handle comment'ed commands
            pass
    return Dictionary

def makenewcalc(basefile, kpts, cutoff):
    # Legacy again
    with open(basefile, 'r') as fileobj:
        mainfolder = os.path.split(basefile)[0]
        newpath = os.path.join(mainfolder,'kpts-'+kpts[:2].strip(),'cutoff='+ str(cutoff))
        if not os.path.exists(newpath):
            os.makedirs(newpath)
        newfile = os.path.join(newpath, os.path.split(basefile)[1])
        with open(newfile, 'w') as newfileobj:
            for line in fileobj:
                searchresult = re.findall(r'ecutwfc[\s]+=',line)
                if searchresult: #isn't null
                    newfileobj.write('    ecutwfc         = ' + str(cutoff)+' ,! autogenerated \n')
                else:
                    searchresult = re.findall(r'K_POINTS',line)
                    if searchresult: # isn't null again
                        newfileobj.write(line)
                        newfileobj.write('  ' + kpts)
                        fileobj.readline()
                    else:
                        newfileobj.write(line)

def changeStructure(QEinputfile, newQEfile = None, structure = Biscotti.AtomicStructure(), AutoKpts = False):
    # Legacy, replaced by QECalcIn classes
    # First find the portion of the .in file which has the atomic structure. Namely these portions

    # Requires that the SYSTEM variable A = 1, use angstrums!!!

    # Really hope that ntyp matches

    # This code is really kind of shit and should be replaced by an actual class with methods.
    # Yay! It sort of was.

    with open(QEinputfile, 'r') as fileobj:

        specieslookup = {'In' : [114.8180, 'In.pbe-dn-kjpaw_psl.0.2.2.UPF'],
                         'As' : [74.9220, 'As.pbe-n-kjpaw_psl.0.2.UPF'],
                         'Sb' : [121.6700,'Sb.pbe-n-kjpaw_psl.0.3.1.UPF' ]}

        filestring = fileobj.read()

        regex0 = r'nat[\s]*=[\s0-9]+'
        regex1 = r'CELL_PARAMETERS[\sA-Za-z]+'
        regex2 = r'ATOMIC_SPECIES'
        regex3 = r'ATOMIC_POSITION[\sA-Za-z]+'
        regex4 = r'K_POINTS[\sA-Za-z]+'

        ntypfind = re.search(r'ntype[\s]*=[\s0-9]+', filestring)
        natfind = re.search(regex0, filestring)
        startlat = re.search(regex1, filestring).start()
        logging.debug('regex1 found: ' + str(startlat))
        startatomspec = re.search(regex2, filestring).start()
        startatompos = re.search(regex3, filestring).start()
        startkpts = re.search(regex4, filestring).start()
        logging.info("All regex searching complete!")

        # Lattice Vectors First
        newcellparamstring = 'CELL_PARAMETERS alat\n' \
                             + '  '.join([str(val) for val in structure.latticeA]) + '\n' \
                             + '  '.join([str(val) for val in structure.latticeB]) + '\n' \
                             + '  '.join([str(val) for val in structure.latticeC]) + '\n'
        logging.debug("New cell parameters string is :\n" + newcellparamstring)
        # Now for species
        specieslist = set([atom.species for atom in structure.atomsarray])
        logging.debug("Found all different atomic species :" + str(specieslist))
        newspeciesstring = 'ATOMIC_SPECIES\n' + '\n'.join([species + "  "
                                                           + str(specieslookup[species][0]) + "  "
                                                           + specieslookup[species][1] for species in specieslist])+'\n'
        logging.info("New atomic species string is:\n" + newspeciesstring)

        # And now atom positions

        newatomsstring = 'ATOMIC_POSITIONS crystal\n' + '\n'.join( [atom.species
                                                                    + "  " + str(atom.x)
                                                                    + "  " + str(atom.y)
                                                                    + "  " + str(atom.z) for atom in structure.atomsdir]) \
                         + '\n'
        logging.info('New atom positions string is:\n' + newatomsstring)
        kptsstring = ''
        if AutoKpts:
            # Optionally set a new kpts grid based on the old one
            Gfactor = AutoKpts # use the passed value
            newA = np.ceil(Gfactor / np.linalg.norm(structure.latticeA))
            newB = np.ceil(Gfactor / np.linalg.norm(structure.latticeB))
            newC = np.ceil(Gfactor / np.linalg.norm(structure.latticeC))
            kptsstring = 'K_POINTS automatic\n' + ' '.join(["%d" % newA, "%d" % newB, "%d" % newC, '0', '0', '0'])
            logging.info("New kpts string is:\n" + kptsstring)
        else:
            kptsstring = filestring[re.search("K_POINTS", filestring).start():]
        # Now create new input file after making replacements..
        newfilestring = filestring[:natfind.start()] + "nat = " + "%d" % structure.totalatoms() + filestring[natfind.end():startlat] + newcellparamstring + newspeciesstring + newatomsstring + kptsstring

        logging.info("Update complete! New file is now:\n" + newfilestring)

    if newQEfile is None:
        newQEfile = os.path.split(QEinputfile)[0] + os.path.sep + os.path.split(QEinputfile)[1] + " with " + structure.name
    with open(newQEfile, 'w') as newfile:
        logging.debug("Now writing file " + newQEfile)
        newfile.write(newfilestring)
    return ''

