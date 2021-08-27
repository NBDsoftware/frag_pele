
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
        fragment_mol = Chem.MolFromPDBFile(pdb_fragment)
        ligand_mol = Chem.MolFromPDBFile(pdb_ligand)

        return fragment_mol, ligand_mol

    def _get_pdb_components(self, complex_pdb):
        """
        Split a protein-ligand pdb into protein and ligand components

        """
        import prody

        pdb = prody.parsePDB(complex_pdb)
        # protein = pdb.select('protein')
        ligand = pdb.select('not protein and not water')
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
        unlinked_molecules = Chem.CombineMols(fragment_mol, ligand_mol)
        mw = Chem.RWMol(unlinked_molecules)

        for atom in mw.GetAtoms():
            if atom.GetPDBResidueInfo() == pdb_atom_core_name:
                ligand_atom_idx = atom.GetIdx()
            elif atom.GetPDBResidueInfo() == pdb_atom_fragment_name:
                fragment_atom_idx = atom.GetIdx()

        mw.ADDBond(ligand_atom_idx, fragment_atom_idx)
        Chem.SanitizeMol(mw)
        smiles = Chem.MolToSmiles(mw)
        grown_ligand = Chem.MolFromSmiles(smiles)
        return grown_ligand


