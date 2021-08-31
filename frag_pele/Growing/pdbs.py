
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
        self._reduce_molecule_size(self.grown_ligand,0.5)
        breakpoint()

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
        from rdkit.Chem import AllChem

        unlinked_molecules = Chem.CombineMols(fragment_mol, ligand_mol)
        mw = Chem.RWMol(unlinked_molecules)
        for atom in mw.GetAtoms():
            if atom.GetPDBResidueInfo().GetName().split()[0] == pdb_atom_core_name and atom.GetPDBResidueInfo().GetChainId() != 'Z':
                ligand_atom_idx = atom.GetIdx()
                ligand_mol.GetAtomWithIdx(atom.GetIdx()).SetAtomicNum(0)
            elif atom.GetPDBResidueInfo().GetName().split()[0] == pdb_atom_fragment_name and atom.GetPDBResidueInfo().GetChainId() == 'Z':
                fragment_atom_idx = atom.GetIdx()
        mw.AddBond(ligand_atom_idx, fragment_atom_idx)
        Chem.SanitizeMol(mw)
        smiles = Chem.MolToSmiles(mw)
        grown_ligand = Chem.MolFromSmiles(smiles)
        AllChem.Compute2DCoords(grown_ligand)
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

    def _reduce_molecule_size(self, molecule, lambda_in=0.5):
        """
        This function performs a reduction of the size of a given residue of a ProDy molecule object.

        Parameters
        ----------
        molecule: ProDy molecule object.
        residue: Resname of the residue of the molecule that we want to reduce. string
        lambda_in: proportion of reduction of the size that we want to apply to the selected residue (between 0 and

        """
        from rdkit.Geometry import Point3D

        if lambda_in >= 0 and lambda_in <= 1:
            centroid = self._compute_centroid(molecule)
            for atom in molecule.GetAtoms():
                atom_coords = self._get_coords(molecule)[atom.GetIdx()]
                new_coords = self._move_atom_along_vector(atom_coords, centroid, lambda_in)
                conf = molecule.GetConformer()
                x, y, z = new_coords[0], new_coords[1], new_coords[2]
                conf.SetAtomPosition(atom.GetIdx(), Point3D(x, y, z))

    def _move_atom_along_vector(self, initial_coord, final_coord, position_proportion):
        """
        Given two points (atom coordinates), this function moves the initial point a distance of "length of the vector
        formed by the two coordinates" * "position_proportion" on the vector's direction.

        Parameters
        ----------
        initial_coord: initial 3D coordinates (X, Y, Z). numpy.ndarray
        final_coord: final 3D coordinates (X, Y, Z). numpy.ndarray
        position_proportion: proportion of movement that we would like to apply on the initial atom on the vector's
        direction. float(generally between 0 and 1)

        """
        vector = final_coord - initial_coord
        new_coords = initial_coord + (position_proportion * vector)
        return new_coords

    def _get_coords(self, grown_ligand):
        """
        Gets coordinates of grown ligand and returns it into an numpy array.
        Parameters
        ----------
        grown_ligand

        Returns
        -------
        coords

        """
        conf = grown_ligand.GetConformer()
        coords = conf.GetPositions()
        return coords

    def _compute_centroid(self, molecule):
        """
        Given a ProDy molecule, the function extract the coordinates of their atoms and compute the centroid of the
        molecule.
        :param molecule: ProDy molecule object.
        :return: centroid of the molecule, tuple(X,Y,Z).
        """
        import numpy as np

        coords = self._get_coords(molecule)
        x = []
        y = []
        z = []
        for coord in coords:
            x.append(float(coord[0]))
            y.append(float(coord[1]))
            z.append(float(coord[2]))
        centroid = (np.mean(x), np.mean(y), np.mean(z))
        return centroid

