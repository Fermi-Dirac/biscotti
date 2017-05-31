import numpy as np
import re
import os
import logging
# Logging level by default
from biscotti import setup_logger
logger = setup_logger(__name__, logging.WARNING)

class Atom(object):
    def __init__(self, species='X', position = np.zeros(3), velocity = np.zeros(3), mass = 0):
        """

        :param symbol: Chemical symbol for this atomic species
        :param position: its position in cartesean coordinates. Default usage is Ansgtrom
        :param velocity: velocity in Angstrum per second
        :param mass: rest mass in Kg, allows for isotopes
        """
        self.species = species # The letters of the thing
        self.position = np.array(position) # In cartesean coordinates, angstroms
        self.x, self.y, self.z = position # for quick reference
        self.velocity = velocity # in Angstrom / sec
        self.mass = mass # allows for isotopes

    def __str__(self):
        returnstring = self.species + " at " + str(self.position)
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
        self.lattice = np.array([self.latticeA, self.latticeB, self.latticeC]) # this is the lattice matrix
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
        logger.info("Now importing atomic structures from " + inputpath)
        with open(inputpath, 'r') as fileobj:
            filestring = fileobj.read()
            # Now to find the steplist of Atomic Positions and Lattices
            regex_cellparam = r'CELL_PARAMETERS.+(\n[ ]*[0-9\-\. ]+){3}'
            # Literal CELL_PARAMETERS, then anything followed by newline, then 3 sets of 3 sets of floats and spaces
            cellparam_match = re.search(regex_cellparam, filestring)
            if cellparam_match:
                logger.debug("Found CELL_PARAMETERS card")
                lattice = []
                latticelist = cellparam_match.group().split('\n')
                for i in range(3):
                    lattice.append(np.array(latticelist[i + 1].split()).astype(float))
                logger.debug("Lattice is: " + str(lattice))
            else:
                logger.error("No CELL_PARAMETERS card")
                lattice = np.identity(3) #TODO account for situations ibrav != 0
            atoms = []
            regex_atompos = r'ATOMIC_POSITIONS[ A-Za-z\(\)]+(\n[A-Z][a-z][ ]+[0-9\-\. ]+)+'
            # Literal ATOMIC_POSITION then anything including parathensis, followed by one cap, one lower case, and then floats
            atompos_match = re.search(regex_atompos, filestring)
            if not atompos_match:
                logger.error("No atoms found?")
            else:
                atomlist = atompos_match.group(0).split('\n')
                atomlist.pop(0)  # the first element is 'ATOMIC_POSITIONS' stuff
                for atomstr in atomlist:
                    logger.debug("Now parsing new atom string: " + atomstr)
                    if atomstr == 'End final coordinates':
                        break
                    newpos = np.array(atomstr.split()[1:]).astype(float)
                    newspec = atomstr.split()[0]
                    newatom = Atom(newspec, np.dot(newpos, lattice)) # Project from Crystal to Cartesean
                    logger.debug("\tNew atom " + str(newatom))
                    atoms.append(newatom)
                logger.debug("Added " + str(len(atoms)) + " atoms for this structure")
        logger.info(">>\tImport of structure from QE Input file complete!\t<<")
        return AtomicStructure(os.path.split(inputpath)[1], lattice[0], lattice[1], lattice[2], np.array(atoms))

    @staticmethod
    def from_QEOutput(QEoutputfile):
        # parses a pw.x quantum espresso output file and returns the atomic structures
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
                alat = float(re.search(r'lattice parameter \(alat\).+=[ ]+[0-9\.]+', filestring).group(0).split('=')[1].strip())
                alat = alat*0.529177249  # Convert from stupid Bohr's to ansgtroms
                # Ok thats nuts, but i'm finding the lattice parameter alat line, then stripping out the number.
                starti = re.search(r'crystal axes', filestring).start()
                latticelist = filestring[starti:starti+400].split('\n')
                logger.debug(str(latticelist))
                for i in range(3):
                    lattice.append(np.array(latticelist[i+1].split()[3:6]).astype(float) * alat)
                    # string looks like: a(1) = (   6.189404   0.000000   0.000000 )
                    # scaling by 'alat'
                logger.debug("Fixed lattice is " + str(lattice))
            else:
                logger.info("Variable unit cell detected")

            # Now to find the steplist of Atomic Positions and Lattices
            if vcrelaxmatch:
                regex_cellparam = r'CELL_PARAMETERS.+(\n[ ]+[0-9\-\. ]+){3}'
                # Literal CELL_PARAMETERS, then anything followed by newline, then 3 sets of 3 sets of floats and spaces
                cellparam_matches = list(re.finditer(regex_cellparam, filestring))
                logger.debug("Found " + str(len(cellparam_matches)) + " lattice matches")

            regex_atompos = r'ATOMIC_POSITIONS[ A-Za-z\(\)]+(\n[A-Z][a-z][ ]+[0-9\-\. ]+)+'
            # Literal ATOMIC_POSITION then anything including parathensis, followed by one cap, one lower case, and then floats
            atompos_matches = list(re.finditer(regex_atompos, filestring))
            if len(atompos_matches) < 1:
                logger.error("No atoms found?")
            else:
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
                logger.debug("Now completing structure " + str(index+1) + " of " + str(len(atompos_matches)))
        logger.info("Import of structure from QE Output file complete! " + str(len(structs)) + " structures imported")
        return structs

    def totalatoms(self):
        return len(self.atomsarray)

    def totalvol(self):
        return abs(np.dot(np.cross(self.latticeA, self.latticeB), self.latticeC))

    def supercell(self, latticefactors):
        """
        :param latticefactors: List of 3 integers
        :return: atomic structure grown in those directions
        """
        logger.info("Expanding Supercell")
        newLatticeA = self.latticeA * latticefactors[0]
        newLatticeB = self.latticeB * latticefactors[1]
        newLatticeC = self.latticeC * latticefactors[2]
        lattice = np.array([newLatticeA, newLatticeB, newLatticeC])
        newatomsdir = []
        newatomscart = []
        for atom in self.atomsdir:  # for each atom using the Direct lattice
            logger.debug("Expanding atom " + str(atom))
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
            logger.debug("Atom expanded, next!")

        #have all the data collected! just create a new atomic structure now...
        logger.info("Supercell created")
        return AtomicStructure(self.name + ' expanded', newLatticeA, newLatticeB, newLatticeC, newatomscart)

    def direct_to_cart(self, atomsdir = None):
        # Takes the direct lattice atoms and produces the cartesean ones
        if atomsdir is None:
            atomsdir = self.atomsdir
        latticematrix = np.array([self.latticeA, self.latticeB, self.latticeC])
        atomscart = []
        for atom in atomsdir:
            newatom = Atom(atom.species, np.array(np.dot(atom.position,latticematrix)))
            atomscart.append(newatom)
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
        """
        This code merges two structures together via a translation vector.
        By default, it arranges them 'on top' via the C dimension
        :param structure2: structure to merge with this one
        :param translation: translation vector between the two. Default is 'end to end' via C
        :return:
        """
        logger.info("-Now starting structure Merge-")
        if translation is None:
            translation = np.array(self.latticeC)
            # By default, glue them end on end via the C vector
        newatomscart = []

        # Make the new lattice vectors first
        newlatmatrix = []
        latmatrix1 = np.array([self.latticeA, self.latticeB, self.latticeC])
        latmatrix2 = np.array([structure2.latticeA, structure2.latticeB, structure2.latticeC])
        logger.debug("Old lattice 1: " + str(latmatrix1))
        logger.debug("Old lattice 2: " + str(latmatrix2))
        for latvect1, latvect2 in zip(latmatrix1, latmatrix2):
            unit_vec = latvect1/np.linalg.norm(latvect1)
            newlatmatrix.append(unit_vec * np.dot(translation, unit_vec) + latvect2)
        #     projection = latvect2 * (1+ np.dot(latvect2,translation) / np.dot(latvect2,latvect2))
        #     newlatvec = []
        #     for latcoord1, latcoord2 in zip(latvect1, projection):
        #         newlatvec.append(latcoord1 if latcoord1 > latcoord2 else latcoord2)
        #     newlatmatrix.append(np.array(newlatvec))

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
        logger.debug("Cartesean atoms now total of " + str(len(newatomscart)))
        logger.debug("New cart atoms list is ")
        for atom in newatomscart:
            logger.debug(str(atom))

        mergedstruct = AtomicStructure(self.name + " + " + structure2.name, newlatmatrix[0], newlatmatrix[1], newlatmatrix[2], newatomscart)
        mergedstruct.cart_to_direct(newatomscart)
        logger.info("Merger complete!")
        return mergedstruct

    def stack_with(self, struct2, lattice_vec_index = 2):
        """
        Stacks an atomic structure end to end with another.
        Both structures are stacked via the lattice vector index
        Good for combining supercells of different material structures ocne they've been strained.
        :param structure2: The second structure to glue to the first
        :param lattice_vec_index: 0, 1 or 2. The A B or C lattice vector to stack them by
        :return: the stacked structure
        """

        logger.info("Now stacking " + self.name + " with " + struct2.name + " by vector index " + str(lattice_vec_index))
        newlattice = np.identity(3)
        newlattice[lattice_vec_index] = self.lattice[lattice_vec_index] + struct2.lattice[lattice_vec_index] # this better be the numpy add!
        logger.debug("New stacked vector is " + str(newlattice[lattice_vec_index]))
        for i in [1, 2]: # For the other two vectors, pick the longer one. Really hoping that are aligned and correlated.
            lat_i = lattice_vec_index-i
            if all(np.isclose(self.lattice[lat_i], struct2.lattice[lat_i])):
                newlattice[lat_i] = self.lattice[lat_i]
            else:
                if np.linalg.norm(self.lattice[lat_i]) > np.linalg.norm(struct2.lattice[lat_i]):
                    newlattice[lat_i] = self.lattice[lat_i]
                else:
                    newlattice[lat_i] = struct2.lattice[lat_i]
        newatoms = list(self.atomscart)
        newatoms.extend([Atom(atom.species, atom.position + self.lattice[lattice_vec_index]) for atom in struct2.atomscart])
        logger.debug("Now stacking " + str(len(self.atomscart)) + " atoms with " + str(len(struct2.atomscart)) + " new atoms")
        stackedstruct = AtomicStructure(self.name + " + " + struct2.name, newlattice[0], newlattice[1], newlattice[2], newatoms)
        stackedstruct.cart_to_direct() #likely not needed
        return stackedstruct

    def strain(self, trans_strain_pct = 0, poisson_ratio=0.35):
        """
        Should be able to use a matrix
        Direction right now is always Z
        Assumes a cube?
        :param poisson_ratio:
        :return:
        """
        # Compute axial strain
        axial_strain_pct = 1 - (1+trans_strain_pct)**(poisson_ratio)
        logger.debug("Axial strain is " + str(axial_strain_pct))

        # # Find the bounding rectangular cuboid
        # logger.debug(self.lattice)
        # xmax, ymax, zmax = np.amax(self.lattice, axis=0)
        # logger.debug("box is " + str([xmax, ymax, zmax]))
        #
        # # Double check for no reason
        # deltaz = -zmax * (1- (1 + trans_strain_pct)**-poisson_ratio)
        # logger.debug("delta z: " + str(deltaz))
        # new_Z = zmax + deltaz
        # new_Z_2 = zmax * (1 + axial_strain_pct)
        # logger.info("New Z: " + str(new_Z) + " New Z again: " + str(new_Z_2))

        # Make new lattice vectors
        for latvec in self.lattice:
            latvec[0] = latvec[0] * (1+ trans_strain_pct)
            latvec[1] *= (1 + trans_strain_pct)
            latvec[2] *= (1 + axial_strain_pct)

        self.latticeA, self.latticeB, self.latticeC = self.lattice
        self.direct_to_cart()

        logger.debug("New lattice is now " + str(self.lattice))
        return self.lattice

        #First, strain the X and Y dimensions by strain percent
        # oldlattice = np.array(self.lattice)
        # for latvec in self.lattice:
        #     latvec[2] = latvec[2]*(1 + trstrain_percent)

        # Compute the new length

        # new_xy = (InAs_len + Sb_len) / (InAs_len / InAs_A + Sb_len / Sb_A)
        # logger.debug("Weighted average is: " + str(new_xy))
        # change_z_InAs = -InAs_A * (1 - (1 + (new_xy - InAs_A) / InAs_A) ** -poisson_ratio)
        # logger.debug("New Z stretch/compression for InAs will be: " + str(change_z_InAs))
        # change_z_Sb = -Sb_A * (1 - ((1 + (new_xy - Sb_A) / Sb_A) ** -poisson_ratio))
        # logger.debug("New Z stretch/compression for InAsSb will be: " + str(change_z_Sb))
        # logger.debug("Now adjusting x,y from " + str([InAs_A, Sb_A]) + " to " + str(new_xy))
        # logger.debug("Now adjusting z from " + str([InAs_A, Sb_A]) + " to " + str(
        #     [InAs_A + change_z_InAs, Sb_A + change_z_Sb]))
        #
        # InAs_struct.latticeA = np.array([new_xy, 0, 0])
        # InAsSb_struct.latticeA = np.array([new_xy, 0, 0])
        # InAs_struct.latticeB = np.array([0, new_xy, 0])
        # InAsSb_struct.latticeB = np.array([0, new_xy, 0])
        # InAs_struct.latticeC = np.array([0, 0, InAs_A + change_z_InAs])
        # InAsSb_struct.latticeC = np.array([0, 0, Sb_A + change_z_Sb])
        # InAs_struct.lattice = [InAs_struct.latticeA, InAs_struct.latticeB, InAs_struct.latticeC]
        # InAsSb_struct.lattice = [InAsSb_struct.latticeA, InAsSb_struct.latticeB, InAsSb_struct.latticeC]
        # logger.info("New lattice for InAs is now " + str(InAs_struct.lattice))
        # logger.info("New lattice for InAsSb is now " + str(InAsSb_struct.lattice))

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
            newfile.write(str(self.name + ". Generated by Biscotti\n"))
            newfile.write("1\n")
            for lattice in [self.latticeA, self.latticeB, self.latticeC]:
                for val in lattice:
                    newfile.write("%.9f" % val + " ")
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
                        newfile.write("  ".join(["%.9f" % val for val in atom.position]) + " " + atom.species + '\n')
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