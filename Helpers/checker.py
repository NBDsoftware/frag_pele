import sys
import logging
import glob
import prody
import Growing.add_fragment_from_pdbs as addfr


# Getting the name of the module for the log system
logger = logging.getLogger(__name__)


def check_duplicated_pdbatomnames(pdb_content):
    pdb_atom_names_list = []
    for line in pdb_content:
        if line.startswith("HETATM") and line.split()[3] != "HOH":
            pdb_atom_name = line[12:15]
            pdb_atom_names_list.append(pdb_atom_name)
    set_to_check = set(pdb_atom_names_list)
    list_to_check = sorted(list(set_to_check))
    sorted_list_names = sorted(pdb_atom_names_list)
    if list_to_check != sorted_list_names:
        print(list_to_check)
        print(sorted_list_names)
        sys.exit("SOME REPEATED PDB ATOM NAMES of the ligand IN PDB FILES!!")


def check_and_fix_pdbatomnames(pdb_file):
    with open(pdb_file) as pdb:
        content = pdb.readlines()
        check_duplicated_pdbatomnames(content)
        for i, line in enumerate(content):
            if line.startswith("HETATM"):
                atom_name = line[12:15]
                if atom_name.strip().startswith("G"):
                    new_atom_name = line[77:78]+atom_name.strip()
                    line_to_list = list(line)
                    line_to_list[12:15] = new_atom_name
                    line_to_list = "".join(line_to_list)
                    content[i] = line_to_list
        check_duplicated_pdbatomnames(content)
        new_pdb = "".join(content)
    with open(pdb_file, "w") as writepdb:
        writepdb.write("{}".format(new_pdb))


def check_if_atom_exists_in_ligand(pdb_file, atom_name):
    ligand = addfr.extract_heteroatoms_pdbs(pdb_file, create_file=False, ligand_chain="L", get_ligand=True)
    atom = ligand.select("name {}".format(atom_name))
    try:
        print("PDB ATOM NAME selected: {} in {}.".format(atom.getNames(), pdb_file))
    except Exception as e:
        logger.critical(e)
        sys.exit("Check if the selected atom '{}' exists in '{}'".format(atom_name, pdb_file))
