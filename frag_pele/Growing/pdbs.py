
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
    def __init__(self):
        """

        """
        pass
