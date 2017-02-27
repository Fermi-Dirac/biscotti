import numpy as np
import re
import os
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

# logging.basicConfig(level=logging.DEBUG, format='%(name)s - %(levelname)s - %(message)s')
# End logging.

class Atom(object):
    def __init__(self, symbol='X', position = np.zeros(3), velocity = np.zeros(3), mass = 0):
        """

        :param symbol: Chemical symbol for this atomic species
        :param position: its position in cartesean coordinates. Default usage is Ansgtrom
        :param velocity: velocity in Angstrum per second
        :param mass: rest mass in Kg, allows for isotopes
        """
        self.symbol = symbol
        self.species = symbol # The letters of the thing
        self.position = np.array(position) # In cartesean coordinates, angstroms
        self.x, self.y, self.z = position # for quick reference
        self.velocity = velocity # in Angstrom / sec
        self.mass = mass # allows for isotopes

    def __str__(self):
        returnstring = self.symbol + " at " + str(self.position)
        if not np.isclose(self.velocity, np.zeros(3)).all():
            returnstring += " with velocity " + str(self.velocity) + " A/s"
        # if self.mass != 0:
        #     returnstring += " and isotope mass " + str(self.mass) + "Kg"
        return returnstring

class AtomicStructure(object):
    def __init__(self, name = 'Default', A = np.array((1,0,0)), B = np.array((0,1,0)), C = np.array((0,0,1)), atomsarray = None):
        """

        :param name: Name for this atomic structure
        :param A: Lattice vector 'a'
        :param B: Lattice vector 'b'
        :param C: Lattice vector 'c'
        :param atomsarray: Array of Atom objects whose positions are in cartesean coordinates
        """
        if atomsarray is None:
            # sentinal value
            atomsarray = []
        self.latticeA = np.array(A)
        self.latticeB = np.array(B)
        self.latticeC = np.array(C)
        self.lattice = [self.latticeA, self.latticeB, self.latticeC] # this is the lattice matrix
        self.name = name
        self.atomsarray = atomsarray
        self.atomscart = atomsarray # allowing for specific usage of cartesean
        self.atoms = atomsarray  # Default functionality is to use cartesean! Legacy
        self.atomsdir = self.cart_to_direct(atomsarray) # Direct lattice atomic positions.

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
        # this extracts the atomic structures from the pw.x input file
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
                        logger.debug("Found regex match of line: " + line)
                        if regcount == 0:
                            logger.debug("Matching lattice:")
                            for j in range(3):
                                logger.debug(lineslist[i+j+1])
                                lattice[j] = np.array(lineslist[i+j+1].split()).astype(float)
                            logger.info("Lattice vectors found as\n" + '\n'.join([str(lat) for lat in lattice]))
                        elif regcount == 1:
                            logger.info("Getting species masses")
                            offset = 1
                            nextline = lineslist[i+offset]
                            while 'ATOMIC_POSITIONS' not in nextline:
                                species = nextline.split()
                                logger.debug("Next species is " + str(species))
                                speciesdict[species[0]] = species[1]
                                offset += 1
                                nextline = lineslist[i+offset]
                            logger.info("species are " + str(speciesdict))
                        elif regcount == 2:
                            logger.info('Getting atomic positions')
                            offset = 1
                            nextline = lineslist[i+offset]
                            while 'K_POINTS' not in nextline:
                                newatom = nextline.split()
                                logger.debug('Found new atom string ' + str(newatom))
                                logger.debug('Next atom is: ' + str(newatom))
                                atomsdir.append(Atom(newatom[0], np.array(newatom[1:4]).astype(float), mass = speciesdict[newatom[0]]))
                                atomscart.append(Atom(newatom[0], np.dot(np.array(newatom[1:4]).astype(float), lattice), mass = speciesdict[newatom[0]]))
                                offset += 1
                                nextline = lineslist[i+offset]
        return AtomicStructure(os.path.split(inputpath)[1], lattice[0], lattice[1], lattice[2], atomscart)

    @staticmethod
    def from_QEOutput(QEoutputfile):
        # parses a pw.x quantum espresso output file and returns the atomic structure
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
                return []
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

            # Now to find the steplist of Atomic Positions and Lattices
            if vcrelaxmatch:
                regex_cellparam = r'CELL_PARAMETERS.+(\n[ ]+[0-9\. ]+){3}'
                # Literal CELL_PARAMETERS, then anything followed by newline, then 3 sets of 3 sets of floats and spaces
                cellparam_matches = list(re.finditer(regex_cellparam, filestring))
                logger.debug("Found " + str(len(cellparam_matches)) + " lattice matches")

            regex_atompos = r'ATOMIC_POSITIONS[ A-Za-z\(\)]+(\n[A-Z][a-z][ ]+[0-9\. ]+)+'
            # Literal ATOMIC_POSITION then anything including parathensis, followed by one cap, one lower case, and then floats
            atompos_matches = list(re.finditer(regex_atompos, filestring))
            logger.debug("Found "+ str(len(atompos_matches)) + " atom set matches")
            logger.debug("The first atom set is " + str(atompos_matches[0].group(0)))

            # Iterate over all atompos_matches
            for index, atompos_match in enumerate(atompos_matches):
                if vcrelaxmatch:
                    # add lattice if variable cell
                    lattice = []
                    latticelist = cellparam_matches[index].group(0).split('\n')
                    for i in range(3):
                        lattice.append(np.array(latticelist[i+1].split()).astype(float))
                logger.debug("Lattice is currently: " + str(lattice))
                atomlist = atompos_match.group(0).split('\n')
                atomlist.pop(0) # the first element is 'ATOMIC_POSITIONS' stuff
                atoms = []
                for atomstr in atomlist:
                    logger.debug("Now parsing new atom string: " + atomstr)
                    if atomstr == 'End final coordinates':
                        break
                    newpos = np.array(atomstr.split()[1:]).astype(float)
                    newspec = atomstr.split()[0]
                    newatom = Atom(newspec, np.dot(newpos, lattice))
                    logger.debug("\tNew atom " + str(newatom))
                    atoms.append(newatom)
                logger.debug("Added " + str(len(atoms)) + " atoms for this structure")
                structs.append(AtomicStructure(os.path.split(QEoutputfile)[1], lattice[0], lattice[1], lattice[2], np.array(atoms)))
                logger.info("Now completing structure " + str(index+1) + " of " + str(len(atompos_matches)))
        logger.info(">>\tImport of structure from QE Output file complete!\t<<")
        return structs

    def totalatoms(self):
        return len(self.atomsarray)

    def totalvol(self):
        return (0)

    def supercell(self,latticefactors):
        """
        :param latticefactors: List of 3 integers
        :return: atomic structure grown in those directions
        """
        verbose = False # Debug flag
        logger.info("   >Starting Supercell<\n")
        newLatticeA = self.latticeA * latticefactors[0]
        newLatticeB = self.latticeB * latticefactors[1]
        newLatticeC = self.latticeC * latticefactors[2]
        lattice = np.array([newLatticeA, newLatticeB, newLatticeC])
        newatomsdir = []
        newatomscart = []
        for atom in self.atomsdir:  # for each atom using the Direct lattice
            logger.info("Expanding atom " + str(atom))
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
                        logger.debug(" New atom made! : " + str(newatom))
                        newatomsdir.append(newatom)
                        newatomscart.append(Atom(atom.species, np.array(np.dot(np.array([newx, newy, newz]), lattice))))
            logger.info("Atom expanded, next!")

        #have all the data collected! just create a new atomic structure now...
        return AtomicStructure(self.name + ' expanded', newLatticeA, newLatticeB, newLatticeC, newatomscart)

    def direct_to_cart(self, atomsdir = None):
        # Takes the direct lattice atoms and produces the cartesean ones
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
        # takes the cartesean atom positions and produces direct lattice ones
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
        # This code merges two structures together via a translation vector.
        # By default, it arranges them 'on top' via the C dimension
        logger.info("-Now starting structure Merge-")
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
        logger.debug("New lattice matrix is " + str(newlatmatrix))

        # Now to get all of the atoms together
        logger.debug("Atoms in self: ")

        newatomscart = []
        for atom in self.atomsarray:
            newatomscart.append(atom)
            logger.debug(" " + str(atom.species) + str(atom.position))
        logger.debug("Atoms in 2nd :")
        for atom in structure2.atomsarray:
            # Create new Atom instance!
            newatomscart.append(Atom(atom.species, atom.position + translation))
            logger.debug(" " + str(atom.species) + str(atom.position) +" shifted by "+ str(translation))

        # newatomscart = np.array(list(self.atomsarray) + list(np.array(structure2.atomsarray) + np.array(translation)))
        # tht doesn't work for some reason.... wtf

        logger.debug("Cartesean atoms now total of " + str(len(newatomscart)))

        logger.debug("New cart atoms list is ")
        for atom in newatomscart:
            logger.debug(str(atom))

        mergedstruct = AtomicStructure(self.name + " + " + structure2.name, newlatmatrix[0], newlatmatrix[1], newlatmatrix[2], newatomscart)
        mergedstruct.cart_to_direct(newatomscart)
        return mergedstruct

    def check_collisions(self, tolerance = 0.01):
        # Checks for atomic colisions in a structure via a tolerence factor
        # good for bug checking
        collisions = False
        for (i, atom) in enumerate(self.atomscart):
            print(atom)
            for j, atom2 in enumerate(self.atomscart):
                if i != j and np.linalg.norm(atom.position - atom2.position) < tolerance:
                    collisions = True
                    logger.error("Collision!")
                    logger.error("See atoms at " + str(atom2))
                    break
            if collisions:
                break
        logger.info("Collision check complete!")
        return collisions

    def QE_Input_String(self):
        atomstrings = '\n'.join([atom.species + " " + str(atom.x) + " " + str(atom.y) + " " + str(atom.z) for atom in self.atomsdir])
        latticestring = 'CELL_PARAMETERS alat\n' + '  '.join([str(x) for x in self.latticeA]) + '\n' + '  '.join([str(x) for x in self.latticeB]) + '\n' + '  '.join([str(x) for x in self.latticeC]) + '\n'
        return "Total Atoms: nat = " + str(self.totalatoms()) + '\n' + latticestring + "\n\n" + atomstrings

    def write_xyz(self, xyzfolder = None, xyzfilename = None):
        # Writes a .xyz file for plotting in Vesta and the like
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
        # Writes a POSCAR-type file for visualiztion in VESTA and the like
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
    # TODO Legacy function now replaced by the QECALC class
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