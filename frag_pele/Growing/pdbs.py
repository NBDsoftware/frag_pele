
class PDBHandler():
    """
    From a core and fragment pdb files, given the heavy atom
    names that we want to connect, add the fragment to the core
    structure.
    We get three pdb files:
        (1) The ligand core of the complex isolated. Will be used
            in further steps to generate the template of the initial
            structure.
        (2) The ligand completed with the core and the fragment added.
        (3) The pdb file that will be used to initialise PELE simulations.
            Here we have the core structure with the fragment added, but
            this fragment has been size-reduced in order to get small bond
            lengths between its atoms.
    """
    def __init__(self, complex_pdb, pdb_fragment, pdb_atom_core_name, pdb_atom_fragment_name):
        from rdkit import Chem
        from rdkit.Chem import Draw
        from rdkit.Chem.Draw import rdMolDraw2D
        from rdkit.Chem import AllChem
        from rdkit.Chem import rdchem

        """

        """
        # Check that ligand names are not repeated
        # Get PDB components
        pdb_ligand = self._get_pdb_components(complex_pdb)
        # Save PDB
        self._write_pdb(pdb_ligand)
        # Generate RDKit molecules
        fragment_mol, ligand_mol = self._create_rdkit_molecules(pdb_ligand, pdb_fragment)
        # Join Molecules
        self.grown_ligand = self._merge_ligand_and_fragment(fragment_mol, ligand_mol, pdb_atom_core_name, pdb_atom_fragment_name)


    def _create_rdkit_molecules(self, pdb_ligand, pdb_fragment):
        """

        Parameters
        ----------
        complex_pdb : PDB file with protein-ligand complex.
        pdb_fragment : PDB file with only the ligand.

        Returns
        -------
        core_mol : RDKit molecule of initial seed compound.
        framgent_mol : RDKit molecule of fragment.
        """
        from rdkit import Chem
        fragment_mol = Chem.MolFromPDBFile(pdb_fragment)
        ligand_mol = Chem.MolFromPDBFile("ligand.pdb")
        return fragment_mol, ligand_mol

    def _get_pdb_components(self, complex_pdb):
        """
        Split a protein-ligand pdb into protein and ligand components

        """
        import prody

        ligand = prody.parsePDB(complex_pdb, chain='Z')
        return ligand

    def _merge_ligand_and_fragment(self, fragment_mol, ligand_mol, pdb_atom_core_name, pdb_atom_fragment_name):
        """
        Merge ligand and fragment an create bond between
        user-specified atoms.
        Parameters
        ----------
        fragment_mol
        ligand_mol

        Returns
        -------

        """
        from rdkit import Chem
        unlinked_molecules = Chem.CombineMols(fragment_mol, ligand_mol)
        mw = Chem.RWMol(unlinked_molecules)
        for atom in mw.GetAtoms():
            if atom.GetPDBResidueInfo().GetName().split()[0] == pdb_atom_core_name and atom.GetPDBResidueInfo().GetChainId() != 'Z':
                ligand_atom_idx = atom.GetIdx()
            elif atom.GetPDBResidueInfo().GetName().split()[0] == pdb_atom_fragment_name and atom.GetPDBResidueInfo().GetChainId() == 'Z':
                fragment_atom_idx = atom.GetIdx()
        mw.AddBond(ligand_atom_idx, fragment_atom_idx)
        Chem.SanitizeMol(mw)
        grown_ligand = mw.GetMol()
        return grown_ligand

    def _write_pdb(self, pdb_ligand):
        """
        Writes PDB file from extracted ligand with prody.

        Parameters
        ----------
        pdb_ligand: prody extracted ligand

        """
        import prody

        output_pdb_name = f"ligand.pdb"
        prody.writePDB(f"{output_pdb_name}", pdb_ligand)
        print(f"wrote {output_pdb_name}")


