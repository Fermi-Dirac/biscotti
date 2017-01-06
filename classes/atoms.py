import numpy as np
import re
import os
import logging

# Logging level by default
logger = logging.getLogger(__name__)
# formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
# streamhandler = logging.StreamHandler()
# streamhandler.setFormatter(formatter)
#
# logger.setLevel(logging.INFO)
# logger.addHandler(streamhandler)

logging.basicConfig(level=logging.INFO, format=' %(asctime)s - %(levelname)s - %(message)s')

# End logging
class Atom(object):
    def __init__(self, symbol='X', position = np.zeros(3), velocity = np.zeros(3), mass = 0):
        self.symbol = symbol
        self.species = symbol # The letters of the thing
        self.position = np.array(position) # In cartesean coordinates, angstrums
        self.x = position[0]
        self.y = position[1]
        self.z = position[2]
        self.velocity = velocity # in Angstrum / sec
        self.mass = mass # allows for isotopes

    def __str__(self):
        return (self.symbol + " at " + str(self.position) + " with velocity " + str(self.velocity))

class AtomicStructure(object):
    """
        This class represents an atomic structure for a crystal
    Attributes:
        Name:       Name of the structure
        Lattice:    Vectors a, b, and c for the unit cell
        AtomicPos:  numpy array of atomic positions, organized by species
        atomsin:    This is a list of atomic species, which is a list of atomic positions and a string for the species
    """
    def __init__(self, name = 'Default', A = np.array((1,0,0)), B = np.array((0,1,0)), C = np.array((0,0,1)), atomsarray = None):
        if atomsarray is None:
            atomsarray = []
        self.latticeA = np.array(A)
        self.latticeB = B
        self.latticeC = C
        self.name = name
        self.atomsarray = atomsarray
        self.atomscart = atomsarray
        self.atoms = atomsarray  # Default functionality is to use cartesean!
        self.atomsdir = self.cart_to_direct(atomsarray)

    def __str__(self):
        return ("Structure '" + self.name
                + "' with " + str(self.totalatoms())
                + ' atoms and a total volume of ' + str(self.totalvol()) + ' Å³'
                + " and lattice vectors of: \n" + str(self.latticeA)
                + "\n" + str(self.latticeB)
                + "\n" + str(self.latticeC) + '\n'
                + "  Atomic Positions of Atoms: <- Cartesean and Direct -> \n" +
                '\n'.join([str(self.atomscart[i]) + " direct: " + str(self.atomsdir[i].position) for i in range(len(self.atomscart))])
                + '\n'
                )

    @staticmethod
    def from_QEinput(inputpath):
        regexlist = ['CELL_PARAMETERS', 'ATOMIC_SPECIES', 'ATOMIC_POSITIONS']
        with open(inputpath) as fileobj:
            lineslist = fileobj.readlines()
            speciesdict = {}
            atomsdir = []
            atomscart = []
            lattice = np.zeros([3,3])
            for i, line in enumerate(lineslist):
                for regcount, regex in enumerate(regexlist):
                    searchresult = re.findall(regex, line)
                    if searchresult: # not empty...
                        logging.debug("Found regex match of line: " + line)
                        if regcount == 0:
                            logging.debug("Matching lattice:")
                            for j in range(3):
                                logging.debug(lineslist[i+j+1])
                                lattice[j] = np.array(lineslist[i+j+1].split()).astype(float)
                            logging.info("Lattice vectors found as\n" + '\n'.join([str(lat) for lat in lattice]))
                        elif regcount == 1:
                            logging.info("Getting species masses")
                            offset = 1
                            nextline = lineslist[i+offset]
                            while 'ATOMIC_POSITIONS' not in nextline:
                                species = nextline.split()
                                logging.debug("Next species is " + str(species))
                                speciesdict[species[0]] = species[1]
                                offset += 1
                                nextline = lineslist[i+offset]
                            logging.info("species are " + str(speciesdict))
                        elif regcount == 2:
                            logging.info('Getting atomic positions')
                            offset = 1
                            nextline = lineslist[i+offset]
                            while 'K_POINTS' not in nextline:
                                newatom = nextline.split()
                                logging.debug('Found new atom string ' + str(newatom))
                                logging.debug('Next atom is: ' + str(newatom))
                                atomsdir.append(Atom(newatom[0], np.array(newatom[1:4]).astype(float), mass = speciesdict[newatom[0]]))
                                atomscart.append(Atom(newatom[0], np.dot(np.array(newatom[1:4]).astype(float), lattice), mass = speciesdict[newatom[0]]))
                                offset += 1
                                nextline = lineslist[i+offset]
        return AtomicStructure(os.path.split(inputpath)[1], lattice[0], lattice[1], lattice[2], atomscart)

    @staticmethod
    def from_QEOutput(QEoutputfile):
        # parses a quantum espresso output file and returns the atomic structure
        structs = []  # Array of all structures

        # Parse through the output file to find the atomic structures
        with open(QEoutputfile) as fileobj:
            logger.info("Now importing atomic structures from " + QEoutputfile)
            filestring = fileobj.read()
            # First check if this file has any structures at all:
            regex_relax = r'ATOMIC_POSITIONS'
            regex_vcrelax = r'CELL_PARAMETERS'
            atomsmatch = re.search(regex_relax, filestring)
            if not atomsmatch:
                logger.error("No atoms in this file")
                return AtomicStructure()
            vcrelaxmatch = re.search(regex_vcrelax, filestring)
            lattice = []
            if not vcrelaxmatch:
                logger.info("Fixed unit cell detected")
                starti = re.search(r'crystal axes', filestring).start()
                latticelist = filestring[starti:starti+400].split('\n')
                logger.debug(str(latticelist))
                for i in range(3):
                    lattice.append(np.array(latticelist[i+1].split()[3:6]).astype(float))
                    # string looks like: a(1) = (   6.189404   0.000000   0.000000 )
                logger.debug("Fixed lattice is " + str(lattice))
            else:
                logger.info("Variable unit cell detected")

            # Now find the list of SCF calcs
            regex_start_step = r'\s+Self-consistent Calculation\s+'
            regex_end_step = r'\s+End of self-consistent calculation\s+'
            # for vc-relax and relax, we have a new SCF at each step.
            startmatches = list(re.finditer(regex_start_step, filestring))
            endmatches = list(re.finditer(regex_end_step, filestring))
            logger.debug("Start matches: " + str(len(startmatches)) + " end matches: " + str(len(endmatches)))
            lastendmatch = endmatches.pop() #handling edge case. The last end match may seek to the EOF
            steplist = []
            for index, endmatch in enumerate(endmatches):
                steplist.append(filestring[endmatch.start():startmatches[index+1].start()])
                # Looking between End and the next Start for the analysis proceedures
            if len(startmatches) > len(endmatches) + 1:
                steplist.append(filestring[lastendmatch.start():startmatches[-1].start()])
                logger.info("This calculation has more starts than ends, did not complete!")
            else:
                steplist.append(filestring[lastendmatch.start():])
            logger.info("Found " + str(len(steplist)) + " matches")

            for index, step in enumerate(steplist):
                logger.debug("Starting step "+ str(index))
                if vcrelaxmatch:
                    lattice = []
                    starti = re.search(regex_vcrelax, step).start()
                    latticelist = filestring[starti:starti + 200].split('\n')
                    for i in range(3):
                        lattice.append(np.array(latticelist[i+1].split()).astype(float))
                    logger.info("Lattice is currently: " + str(lattice))
                atomstart = re.search(regex_relax, step).start()
                logger.debug(" now checking for end in : " + step[atomstart:atomstart+100])
                atomsend = re.search('\n\n', step[atomstart:]).start()
                logger.debug("atoms found at indicies: " + str(atomstart) + " to " + str(atomsend))
                atomlist = step[atomstart:atomstart + atomsend].split('\n') # make a list of atom strings
                atomlist.pop(0) #Skip the first which has the header
                atoms = []
                for atomstr in atomlist:
                    logger.debug("Now parsing new atom string: " + atomstr)
                    if atomstr == 'End final coordinates':
                        break
                    newpos = np.array(atomstr.split()[1:]).astype(float)
                    newspec = atomstr.split()[0]
                    logger.debug("-New atom " + newspec + " at " + str(newpos))
                    atoms.append(Atom(newspec, np.dot(newpos, lattice)))
                logger.info("Added " + str(len(atoms)) + " atoms for this structure")
                structs.append(AtomicStructure(os.path.split(QEoutputfile)[1], lattice[0], lattice[1], lattice[2], np.array(atoms)))
                logger.info("Now completing structure " + str(index) + " of " + str(len(steplist)))

            # lineslist = fileobj.readlines()  # sadly need to read all lines anyway. very inefficient!
            # startline = 0
            #  endline = 0
            #
            # for linenumber, line in enumerate(lineslist):
            #     if startline == 0:  # haven't found start yet
            #         searchresult = re.findall(regex, line)
            #         if searchresult:  # is not empty....
            #             logging.debug("Match found at " + str(linenumber))
            #             startline = linenumber
            #     elif endline == 0:  # found the start, but not the end
            #         if line == 'End final coordinates\n':
            #             endline = linenumber
            #             logging.info("found last structure at " + str(linenumber))
            #             structs.append(''.join(lineslist[startline:endline]))
            #             break
            #         if line == '\n' and lineslist[linenumber + 1] == '\n':
            #             # if two consecutive lines are null
            #             logging.info("found end of this structure at " + str(linenumber))
            #             endline = linenumber  # not the end of the file though!
            #     else:
            #         # then we found an atomic structure, but haven't found the last one yet
            #         structs.append(''.join(lineslist[startline:endline]))
            #         startline = 0
            #         endline = 0
        #
        # cellparamstring = ''
        # if len(structs) > 1:
        #     cellparamstring = structs[-1]  # Just use the final structure
        # logging.debug(" found the atoms strings: \n\n" + cellparamstring)
        # cellparamarray = cellparamstring.split('\n')
        # logging.info(cellparamarray)
        # lattice = np.zeros((3, 3))  # create some empty arrys for the lattice vectors
        # atoms = []  # TODO LEGACY not sure if i need to intialize, except for scope reasons.
        # atoms2 = [] # This is the list of atoms from the Atoms class object, not the old thing
        # for i in range(3):
        #     logging.debug(str(i)+ str(lattice[i]) + str(cellparamarray[i + 1]))
        #     lattice[i] = np.array(cellparamarray[i + 1].split()).astype(float)
        #     logging.debug(lattice)
        # for atomstr in cellparamarray[6:-1]:  #only past the 3 lattice vectors, and strip the last empty element
        #     logging.debug(atomstr)
        #     species = atomstr.split()[0]
        #     newatom = np.dot(np.array(atomstr.split()[1:]).astype(float), lattice)  # string to float, then to nparray, then to cartesean from direct
        #     atoms2.append(Atom(species,newatom))
        #     logging.debug("new atom: " + species + " at " + str(newatom))

        return structs

    def totalatoms(self):
        # total = 0
        # for species in self.atoms:
        #     total += len(species[1])
        return len(self.atomsarray)

    def totalvol(self):
        return (0)

    def supercell(self,latticefactors):
        """
        :param latticefactors: List of 3 integers
        :return:
        """
        verbose = False # Debug flag
        logging.info("   >Starting Supercell<\n")
        newLatticeA = self.latticeA * latticefactors[0]
        newLatticeB = self.latticeB * latticefactors[1]
        newLatticeC = self.latticeC * latticefactors[2]
        lattice = np.array([newLatticeA, newLatticeB, newLatticeC])
        newatomsdir = []
        newatomscart = []
        for atom in self.atomsdir:  # for each atom using the Direct lattice
            logging.info("Expanding atom " + str(atom))
            newa = (atom.position[0] / latticefactors[0])
            for a in range(latticefactors[0]):
                newx = newa + (a / latticefactors[0])
                newb = (atom.position[1] / latticefactors[1])
                for b in range(latticefactors[1]):
                    newy = newb + (b/latticefactors[1])
                    newc = (atom.position[2] / latticefactors[2])
                    for c in range(latticefactors[2]):
                        newz = newc + (c/latticefactors[2])
                        newatom = Atom(atom.species, np.array([newx, newy, newz]))
                        logging.debug(" New atom made! : " + str(newatom))
                        newatomsdir.append(newatom)
                        newatomscart.append(Atom(atom.species, np.array(np.dot(np.array([newx, newy, newz]), lattice))))
            logging.info("Atom expanded, next!")

        #have all the data collected! just create a new atomic structure now...
        return AtomicStructure(self.name + ' expanded', newLatticeA, newLatticeB, newLatticeC, newatomscart)

    def direct_to_cart(self, atomsdir = None):
        if atomsdir is None:
            atomsdir = self.atoms
        else:
            self.atoms = atomsdir
        latticematrix = np.array([self.latticeA, self.latticeB, self.latticeC])
        atomscart = []
        for species in atomsdir:
            atomlist = []
            for atom in species[1]:
                atomlist.append(np.array(np.dot(atom,latticematrix)))
            atomscart.append([species[0], atomlist])
        self.atomscart = atomscart
        return atomscart

    def cart_to_direct(self,atomscart = None):
        if atomscart is None:
            atomscart = self.atomscart
        else:
            self.atomscart = atomscart
        latticematrix = np.array([self.latticeA, self.latticeB, self.latticeC])
        invlatmatrix = np.linalg.inv(latticematrix)
        atomsdir = []
        for atom in atomscart:
            newatom = Atom(atom.species, np.array(np.dot(atom.position, invlatmatrix)))
            atomsdir.append(newatom)
            # So take the inverse matrix of the lattice, and dot it on the right with the atom vector
        self.atomsdir = np.array(atomsdir)
        return atomsdir

    def merge_with(self, structure2, translation = None):
        logging.info("-Now starting structure Merge-")
        if translation is None:
            translation = np.array(self.latticeC)
            # By default, glue them end on end via the C vector

        verbose = True
        newatomscart = []

        # Make the new lattice vectors first
        newlatmatrix = []
        latmatrix1 = np.array([self.latticeA, self.latticeB, self.latticeC])
        latmatrix2 = np.array([structure2.latticeA, structure2.latticeB, structure2.latticeC])
        for latvect1, latvect2 in zip(latmatrix1, latmatrix2):
            projection = latvect2 * (1+ np.dot(latvect2,translation)/np.dot(latvect2,latvect2))
            newlatvec = []
            for latcoord1, latcoord2 in zip(latvect1, projection):
                newlatvec.append(latcoord1 if latcoord1 > latcoord2 else latcoord2)
            newlatmatrix.append(np.array(newlatvec))
        logging.debug("New lattice matrix is " + str(newlatmatrix))

        # Now to get all of the atoms together
        logging.debug("Atoms in self: ")

        newatomscart = []
        for atom in self.atomsarray:
            newatomscart.append(atom)
            logging.debug(" " + str(atom.species) + str(atom.position))
        logging.debug("Atoms in 2nd :")
        for atom in structure2.atomsarray:
            # Create new Atom instance!
            newatomscart.append(Atom(atom.species, atom.position + translation))
            logging.debug(" " + str(atom.species) + str(atom.position) +" shifted by "+ str(translation))

        # newatomscart = np.array(list(self.atomsarray) + list(np.array(structure2.atomsarray) + np.array(translation)))
        # tht doesn't work for some reason.... wtf

        logging.debug("Cartesean atoms now total of " + str(len(newatomscart)))

        logging.debug("New cart atoms list is ")
        for atom in newatomscart:
            logging.debug(str(atom))

        mergedstruct = AtomicStructure(self.name + " + " + structure2.name, newlatmatrix[0], newlatmatrix[1], newlatmatrix[2], newatomscart)
        mergedstruct.cart_to_direct(newatomscart)
        return mergedstruct

    def check_collisions(self, tolerance = 0.01):
        collisions = False
        for (i, atom) in enumerate(self.atomscart):
            print(atom)
            for j, atom2 in enumerate(self.atomscart):
                if i != j and np.linalg.norm(atom.position - atom2.position) < tolerance:
                    collisions = True
                    logging.error("Collision!")
                    logging.error("See atoms at " + str(atom2))
                    break
            if collisions:
                break
        logging.info("Collision check complete!")
        return collisions

    def QE_Input_String(self):
        atomstrings = '\n'.join([atom.species + " " + str(atom.x) + " " + str(atom.y) + " " + str(atom.z) for atom in self.atomsdir])
        latticestring = 'CELL_PARAMETERS alat\n' + '  '.join([str(x) for x in self.latticeA]) + '\n' + '  '.join([str(x) for x in self.latticeB]) + '\n' + '  '.join([str(x) for x in self.latticeC]) + '\n'
        return "Total Atoms: nat = " + str(self.totalatoms()) + '\n' + latticestring + "\n\n" + atomstrings

    def write_xyz(self, xyzfolder = None, xyzfilename = None):
        if xyzfolder is None:
            xyzfolder = os.getcwd()
        else:
            if not os.path.exists(xyzfolder):
                os.makedirs(xyzfolder)  # Need to catch errors here
        if xyzfilename is None:
            xyzfilename = self.name + ".xyz"
        with open(xyzfolder + os.path.sep + xyzfilename, 'w') as newfile:
            newfile.write(str(self.totalatoms())+'\n')
            newfile.write(self.name + "generated by Biscotti\n")
            for atom in self.atomscart:
                newfile.write(atom.species + "  " + "  ".join(["%.6f" % val for val in atom.position])+'\n')
        return True

    def write_vasp(self, folder = None, filename = None):
        if folder is None:
            folder = os.getcwd()
        else:
            if not os.path.exists(folder):
                os.makedirs(folder)  # Need to catch errors here
        if filename is None:
            filename = "POSCAR." + self.name + ".vasp"
        with open(folder + os.path.sep + filename, 'w') as newfile:
            newfile.write(str(self.name + " generated by Biscotti\n"))
            newfile.write("1\n")
            for lattice in [self.latticeA, self.latticeB, self.latticeC]:
                for val in lattice:
                    newfile.write("%.6f" % val + " ")
                newfile.write("\n")
            speclist = [atom.species for atom in self.atomsdir]
            speciesdict = {spec:speclist.count(spec) for spec in speclist}
            for species in speciesdict:
                newfile.write(species + " ")
            newfile.write('\n')
            for species in speciesdict:
                newfile.write("%d" % speciesdict[species] + " ")
            newfile.write("\nDirect\n")
            for species in speciesdict:
                for atom in self.atomsdir:
                    if species == atom.species:
                        newfile.write("  ".join(["%.6f" % val for val in atom.position]) + " " + atom.species + '\n')
            logger.debug('Vasp file written!')
            return True

def makenewcalc(basefile, kpts, cutoff):
    """

    :param basefile: File path to make a new calc from
    :param kpts: the kpts as a string to  be used
    :param cutoff: the cutoff in Ry for this calculation
    :return: none
    """
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