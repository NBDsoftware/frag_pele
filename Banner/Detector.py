import prody
import logging

# Getting the name of the module for the log system
logger = logging.getLogger(__name__)


class Detector:
    def __init__(self, pdb, threshold, atom1, atom2, atom3, atom4, lig_chain="L"):
        self.pdb = pdb
        self.threshold = threshold
        self.lig_chain = lig_chain
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4

    def read_pdb(self):
        pdb = prody.parsePDB(self.pdb)
        return pdb

    def get_ligand(self):
        pdb = self.read_pdb()
        ligand = pdb.select("chain {}".format(self.lig_chain))
        if ligand is None:
            logger.critical("Wrong chain selected!")
        elif ligand.ishetero:
            return ligand
        else:
            logger.critical("The selected chain does not contain heteroatoms!")

    def select_atoms(self):
        ligand = self.get_ligand()
        atoms = ligand.select("name {} {} {} {}".format(self.atom1, self.atom2, self.atom3, self.atom4))
        return atoms

    def read_dihedral(self):
        atoms = self.select_atoms()
        dihedral_atoms = []
        for atom in atoms:
            dihedral_atoms.append(atom)
        return prody.calcDihedral(dihedral_atoms[0], dihedral_atoms[1], dihedral_atoms[2], dihedral_atoms[3])

    def check_threshold_dihedral(self):
        value = self.read_dihedral()
        if abs(self.threshold) > value:
            return True
        else:
            return False

