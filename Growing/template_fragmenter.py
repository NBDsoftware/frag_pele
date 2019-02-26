import sys
import re
import os
import logging

# Getting the name of the module for the log system
logger = logging.getLogger(__name__)

# Definition of reggex patterns
HEADER_OPLS2005 = "* LIGAND DATABASE FILE (OPLS2005)\n*\n"
RESX_PATTERN_OPLS2005 = "(\w{4,})\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)"
# Example: LYS      22    21     38     51     110
# (Template name) (Num of atoms) (NumElements in BOND) (NumElements in THET) (NumElements in PHI)
ATOM_PATTERN_OPLS2005 = "\s+(\d+)\s+(\d+)\s+(\w)\s+(\w*)\s+(\w{4,})\s+(\d*)\s+(-?[0-9]*\.[-]?[0-9]*)\s+(-?[0-9]*\.[0-9]*)\s+(-?[0-9]*\.[0-9]*)"
NBON_PATTERN_OPLS2005 = "\s+(\d+)\s+(\d+\.\d{4})\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
BOND_PATTERN_OPLS2005 = "\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)"
WRITE_NBON_PATTERN_OPLS2005 = " {:5d}   {:3.4f}   {:3.4f}  {: 3.6f}   {:3.4f}   {:3.4f}   {:3.9f}  {: 3.9f}\n"
WRITE_BOND_PATTERN_OPLS2005 = " {:5d} {:5d}   {:5.3f} {: 2.3f}\n"


# Classes definitions

class Atom:
    def __init__(self, atom_id, parent_id, location, atom_type, pdb_atom_name, unknown, x_zmatrix=0, y_zmatrix=0,
                 z_zmatrix=0, sigma=0, epsilon=0, charge=0, radnpSGB=0, radnpType=0, sgbnpGamma=0, sgbnpType=0,
                 is_linker=False, is_fragment=False):
        self.atom_id = atom_id
        self.parent_id = parent_id
        self.location = location
        self.atom_type = atom_type
        self.pdb_atom_name = pdb_atom_name
        self.unknown = unknown
        self.x_zmatrix = x_zmatrix
        self.y_zmatrix = y_zmatrix
        self.z_zmatrix = z_zmatrix
        self.sigma = sigma
        self.epsilon = epsilon
        self.charge = charge
        self.radnpSGB = radnpSGB
        self.radnpType = radnpType
        self.sgbnpGamma = sgbnpGamma
        self.sgbnpType = sgbnpType
        self.is_fragment = is_fragment
        self.is_linker = is_linker


