import prody
import os
import pdb
import logging
logger = logging.getLogger(__name__)

class ProdyMerging():
    def __init__(self):
        pass

    def _extract_heavy_atoms(self,
                             pdb_atom_names,
                             lists_of_bioatoms):
        heavy_atoms = []
        # Select the heavy atoms for each list that we will want to bond together in further steps
        for atom_name, list_of_bioatoms in zip(pdb_atom_names, lists_of_bioatoms):
            atom_heavy = pdb_joiner.select_atoms_from_list(atom_name, list_of_bioatoms)
            heavy_atoms.append(atom_heavy)
        return heavy_atoms

    def _extract_hydrogens(self,
                           pdb_atom_names,
                           lists_of_bioatoms,
                           list_of_pdbs,
                           h_core=None,
                           h_frag=None,
                           c_chain='L',
                           f_chain='L',
                           c_resnum=None,
                           f_resnum=None):
        hydrogens = []
        selected_hydrogens = [h_core, h_frag]
        chains = [c_chain, f_chain]
        resnums = [c_resnum, f_resnum]
        for atom_name, pdb, list_of_bioatoms, sel_h, chain, resnum in zip(pdb_atom_names, list_of_pdbs,
                                                                          lists_of_bioatoms, selected_hydrogens, chains,
                                                                          resnums):
            complex = prody.parsePDB(pdb)
            # Select name of the H atoms bonded to this heavy atom (the place where we will grow)
            atom_name_hydrogens = pdb_joiner.get_H_bonded_to_grow(atom_name, complex, sel_h, chain=chain, resnum=resnum)
            # Select this hydrogen atoms
            atom_hydrogen = pdb_joiner.select_atoms_from_list(atom_name_hydrogens, list_of_bioatoms)
            hydrogens.append(atom_hydrogen)
        return hydrogens

    def _join_structures(self,
                         core_bond,
                         fragment_bond,
                         core_structure,
                         fragment_structure,
                         pdb_complex,
                         pdb_fragment,
                         chain_complex,
                         chain_fragment,
                         output_path,
                         only_grow=False,
                         core_resnum=None):
        name_to_replace_core = core_bond[1].name
        name_to_replace_fragment = fragment_bond[0].name
        if only_grow:
            return 0, name_to_replace_core, name_to_replace_fragment
        if RDKIT:

            if core_resnum:
                atoms_to_delete_core = tree_detector.main(pdb_complex, (core_bond[0].name, name_to_replace_core),
                                                          chain=chain_complex, resnum=core_resnum)
            else:
                atoms_to_delete_core = tree_detector.main(pdb_complex, (core_bond[0].name, name_to_replace_core),
                                                          chain=chain_complex)

            atoms_to_delete_fragment = tree_detector.main(pdb_fragment,
                                                          (fragment_bond[1].name, name_to_replace_fragment),
                                                          chain=chain_fragment)
        else:
            print("WARNING: YOU CAN NOT REPLACE HEAVY ATOMS FOR HYDROGENS WITHOUT RDKIT!")
        fragment_rdkit = rdkit.Chem.MolFromPDBFile(pdb_fragment, removeHs=False)
        bond_type = detect_bond_type(fragment_rdkit, fragment_bond[0], fragment_bond[1])
        print(bond_type)
        if RDKIT:
            if core_bond[1].element != "H":
                atom_replaced_idx = replace_heavy_by_hydrogen(core_bond[1], core_structure)
                new_coords, new_dist = correct_hydrogen_position(hydrogen_atom=core_structure[atom_replaced_idx],
                                                                 atom_to_bond_with=core_bond[0],
                                                                 structure=core_structure[atom_replaced_idx])
                core_structure[atom_replaced_idx].setCoords(new_coords)
                names_to_keep = list(
                    set(core_structure.getNames()) ^ set(
                        atoms_to_delete_core))  # Compare two sets and get the common items
                names_to_keep.remove(name_to_replace_core)
                core_structure = core_structure.select("name {}".format(" ".join(names_to_keep)))
                prody.writePDB(os.path.join(output_path, "{}.pdb".format(core_structure.getResnames()[0])),
                               core_structure)  # Overwrite the initial structure
                name_to_replace_core = core_structure[atom_replaced_idx].getName()
                core_bond[1].coord = new_coords
            if fragment_bond[0].element != "H":
                atom_replaced_idx = replace_heavy_by_hydrogen(fragment_bond[0], fragment_structure)
                new_coords, new_dist = correct_hydrogen_position(hydrogen_atom=fragment_structure[atom_replaced_idx],
                                                                 atom_to_bond_with=fragment_bond[1],
                                                                 structure=fragment_structure[atom_replaced_idx])
                fragment_structure[atom_replaced_idx].setCoords(new_coords)
                names_to_keep = list(
                    set(fragment_structure.getNames()) ^ set(
                        atoms_to_delete_fragment))  # Compare two sets and get the common items
                names_to_keep.remove(name_to_replace_fragment)
                fragment_structure = fragment_structure.select("name {}".format(" ".join(names_to_keep)))
                prody.writePDB(os.path.join(output_path, "{}.pdb".format(fragment_structure.getResnames()[0])),
                               fragment_structure)  # Overwrite the initial structure
                name_to_replace_fragment = fragment_structure[atom_replaced_idx].getName()
                fragment_bond[0].coord = new_coords
        bio_list = \
        from_pdb_to_bioatomlist([os.path.join(output_path, "{}".format(fragment_structure.getResnames()[0]))])[
            0]  # Its a list, so we keep only the unique element that is inside
        # Superimpose atoms of the fragment to the core bond
        pdb_joiner.superimpose(core_bond, fragment_bond, bio_list)
        # Get the new coords and change them in prody
        transform_coords_from_bio2prody(fragment_structure, bio_list)
        # Now, we have to remove the hydrogens of the binding
        h_atom_names = [name_to_replace_core, name_to_replace_fragment]
        # Correcting linking distance
        new_coords, new_distance = correct_bonding_distance(atom_reference=core_bond[0],
                                                            atom_to_correct=fragment_bond[1],
                                                            reference_structure=core_structure,
                                                            movil_structure=fragment_structure,
                                                            bond_type=bond_type)
        if bond_type == "double":
            hydrogen_to_delete = pdb_joiner.get_H_bonded_to_atom(core_bond[0].name, core_structure,
                                                                 bond_dist=new_distance,
                                                                 banned_hydrogen=name_to_replace_core,
                                                                 chain=chain_complex)
            names_to_keep = list(core_structure.getNames())
            names_to_keep.remove(hydrogen_to_delete)
            core_structure = core_structure.select("name {}".format(" ".join(names_to_keep)))

        if bond_type == "triple":
            for n in range(2):
                hydrogen_to_delete = pdb_joiner.get_H_bonded_to_atom(core_bond[0].name, core_structure,
                                                                     bond_dist=new_distance,
                                                                     banned_hydrogen=name_to_replace_core,
                                                                     chain=chain_complex)
                names_to_keep = list(core_structure.getNames())
                names_to_keep.remove(hydrogen_to_delete)
                core_structure = core_structure.select("name {}".format(" ".join(names_to_keep)))

        fragment_structure.setCoords(new_coords)
        merged_structure = bond(h_atom_names, [core_structure, fragment_structure])
        return merged_structure, name_to_replace_core, name_to_replace_fragment, new_distance  # new_distance to modify the clash treshold
        # in check_collisions()

    def _check_collisions(self,
                          merged_structure,
                          bond,
                          theta,
                          theta_interval,
                          core_bond,
                          list_of_atoms,
                          fragment_bond,
                          core_structure,
                          fragment_structure,
                          pdb_complex,
                          pdb_fragment,
                          chain_complex,
                          chain_fragment,
                          output_path,
                          threshold_clash=None,
                          only_grow=False):
        core_resname = bond[0].get_parent().get_resname()
        frag_resname = bond[1].get_parent().get_resname()
        print(core_resname, frag_resname)
        if core_resname is frag_resname:
            logger.critical("The resname of the core and the fragment is the same. Please, change one of both")
        print("resname {} and within {} of resname {}".format(core_resname,
                                                              threshold_clash,
                                                              frag_resname))
        check_possible_collision = merged_structure.select("resname {} and within {} of resname {}".format(core_resname,
                                                                                                           threshold_clash,
                                                                                                           frag_resname))
        # This list only should have the atom of the fragment that will be bonded to the core, so if not we will have to
        # solve it
        if len(check_possible_collision.getNames()) > 1 or check_possible_collision.getNames()[0] != bond[0].name:
            print("We have a collision between atoms of the fragment {} and the core {} within {} A!"
                  " Rotating the fragment to solve it...".format(frag_resname, core_resname, threshold_clash))
            theta = theta + theta_interval
            if theta >= math.pi * 2:
                print("Not possible solution, decreasing the angle of rotation...")
            else:
                rotated_structure = rotation_thought_axis(bond, theta, core_bond, list_of_atoms, fragment_bond,
                                                          core_structure,
                                                          fragment_structure, pdb_complex, pdb_fragment, chain_complex,
                                                          chain_fragment,
                                                          output_path=output_path, only_grow=only_grow)
                if debug:
                    print(theta)
                    prody.writePDB("testing_{}.pdb".format(theta), rotated_structure[0])
                recall = check_collision(merged_structure=rotated_structure[0], bond=bond, theta=theta,
                                         theta_interval=theta_interval,
                                         core_bond=core_bond, list_of_atoms=list_of_atoms, fragment_bond=fragment_bond,
                                         core_structure=core_structure,
                                         fragment_structure=fragment_structure, pdb_complex=pdb_complex,
                                         pdb_fragment=pdb_fragment,
                                         chain_complex=chain_complex, chain_fragment=chain_fragment,
                                         output_path=output_path, threshold_clash=threshold_clash, only_grow=only_grow,
                                         debug=debug)
                return recall
        else:
            return merged_structure

    def _extract_and_change_atomnames(self,
                                      molecule,
                                      selected_resname,
                                      core_resname,
                                      rename=False):
        assert selected_resname != core_resname, "core and fragment residue name must be different"
        fragment = molecule.select("resname {}".format(selected_resname))
        core = molecule.select("resname {}".format(core_resname))
        core_atom_names = [atom.getName() for atom in core]
        fragment_atom_names = [atom.getName() for atom in fragment]
        names_dictionary = {}
        for n, atom_name in enumerate(fragment_atom_names):
            if rename:
                names_dictionary[atom_name] = "G{}".format(n)
            else:
                # If the atomname is repited
                if atom_name in core_atom_names:
                    initial_atom_name = atom_name
                    while atom_name in core_atom_names:
                        atom_name_digit = re.findall('\d+', atom_name)[0]
                        new_atom_name_digit = int(atom_name_digit) + 1
                        atom_name = atom_name.replace(atom_name_digit, str(new_atom_name_digit))
                    final_atom_name = atom_name
                    core_atom_names.append(final_atom_name)
                    names_dictionary[initial_atom_name] = final_atom_name
                else:
                    names_dictionary[atom_name] = atom_name
        for atom in molecule:
            if atom.getResname() == selected_resname:
                if atom.getName() in names_dictionary:
                    atom.setName(names_dictionary[atom.getName()])
        return molecule, names_dictionary

    def _reduce_molecule_size(self,
                              molecule,
                              residue,
                              lambda_in):
        if lambda_in >= 0 and lambda_in <= 1:
            selection = molecule.select("resname {}".format(residue))
            centroid = compute_centroid(selection)
            for atom in selection:
                atom_coords = atom.getCoords()
                new_coords = move_atom_along_vector(atom_coords, centroid, lambda_in)
                atom.setCoords(new_coords)
        else:
            logger.critical("Sorry, reduce_molecule_size() needs a lambda value between 0 and 1!")

    def _extract_atoms_pdbs(self,
                            pdb,
                            create_file=True,
                            chain='L',
                            resnum=None,
                            get_atoms=False,
                            output_folder="."):
        if not resnum:
            selection = self.pdb_parser_ligand(pdb, chain)
        else:
            selection = self.pdb_parser_residue(pdb, chain, resnum)
        if selection is None:
            raise TypeError(
                "The selection can not be found. Selection for {}: chain {} and resnum {}".format(pdb, chain, resnum))
            # Check if the ligand has H
        self.check_protonation(selection)
        # Save the ligand in a PDB (the name of the file is the name of the residue)
        name = selection.getResnames()[0]
        if create_file:
            prody.writePDB(os.path.join(output_folder, name), selection)
            print("The selection of {} has been extracted and saved in '{}.pdb'".format(pdb, os.path.join(output_folder,
                                                                                                          name)))
        if get_atoms is True:
            return selection
        else:
            return name

    def pdb_parser_ligand(self,
                          pdb_file,
                          ligand_chain="L"):

        pdb = prody.parsePDB(pdb_file)
        ligand = pdb.select("chain {}".format(ligand_chain))
        if ligand is None:
            logger.critical("Wrong chain selected!")
        elif ligand.ishetero:
            return ligand
        else:
            logger.critical("The selected chain does not contain heteroatoms!")

    def pdb_parser_residue(self,
                           pdb_file,
                           res_chain,
                           res_num):

        pdb = prody.parsePDB(pdb_file)
        residue = pdb.select("chain {} and resnum {}".format(res_chain, res_num))
        if residue is None:
            logger.critical("Wrong selection!")
        else:
            return residue

    def check_protonation(self,
                          selection):

        try:
            if not selection.select("hydrogen"):
                logger.critical("We have not detected Hydrogens in your ligand. Please, add them before starting.")
        except AttributeError:
            raise AttributeError("Check ligand and core are in the L chain. Otherwise specify their chain using the flags \
    -cc chain_core -fc chain_frag. Fragment & Core must always be in the same chain.")
