import os.path
import os
import frag_pele.constants as c
from frag_pele.Errors.custom_errors import WrongComplexName


class ParametersBuilder():
    def __init__(
            self,
            serie_file,
            protocol=False,
            test=False):
        """
        It builds the parameters for FrAG PELE, according to the argument that are supplied and
        return the corresponding parameters.

        Parameters
        ----------
        serie_file : str
           Name of the tabular file which must contain the instructions required to perform several
           successive growings, using different fragments or different growing positions.
        protocol
           Select growing protocol. Choose between: 'SofcoreLike', 'SpreadingHcharge', 'AllLinear'.
           SofcoreLike: Charges initially set to 0. They are added in the mid GS. Then, they grow exponentially
           or linearly (depending on your settings).
           SpreadingHcharge: reimplementation of FragPELE1.0.0 methodology.
           AllLinear: All FF parameters are linearly and equally incremented in each GS.
        test : bool
           Run test config

        """
        self._original_dir = os.path.abspath(os.getcwd())
        self._list_of_instructions = self._read_instructions_from_file(serie_file)
        # self.dict_traceback
        self._ID = self._extract_id(self._list_of_instructions)
        self._core_atom, self._h_core, self._fragment_atom, self._h_frag = self._extract_linker_atoms(
            self._list_of_instructions)
        self._steps, self._pele_eq_steps, self._temp = self._define_protocol_params(protocol, test)


    def _define_protocol_params(self, protocol, test):
        """
        Checks protocol type and assigns steps, pele equilibration steps and temperature parameters.

        Parameters
        ----------
        protocol

        Returns
        -------

        """
        if protocol == "HT":
            # iteration = 1 In code we have 'iteration' but we pass to grow_fragment function 'iterations', is this useless or wrong? Ask Carles.
            steps = 3
            pele_eq_steps = 10
            temp = c.TEMPERATURE

        elif test:
            # iteration = 1
            steps = 1
            pele_eq_steps = 1
            temp = 1000000

        else:
            steps = c.STEPS
            pele_eq_steps = c.PELE_EQ_STEPS
            temp = c.TEMPERATURE

        return steps, pele_eq_steps, temp

    def _read_instructions_from_file(self, file):
        """
        From the serie file reads the instructions.
        Parameters
        ----------
        file

        Returns
        -------
        list_of_instructions

        """
        list_of_instructions = []
        with open(file) as sf:
            instructions = sf.readlines()
            for line in instructions:
                if 0 < len(line.split()) <= 3:
                    # Read from file the required information
                    try:
                        fragment_pdb = line.split()[0]
                        core_atom = line.split()[1]
                        fragment_atom = line.split()[2]
                        ID = "{}{}{}".format(os.path.splitext(fragment_pdb)[0], core_atom, fragment_atom)
                        task = (fragment_pdb, core_atom, fragment_atom, ID)
                        list_of_instructions.append(task)
                    except IndexError:
                        logger.critical(
                            "Check that the serie file {} contains: 'PDB_fragment_file'\t'PDB_core_atom_name'\t'PDB_fragment_atom_name' ".format(
                                file))
                elif len(line.split()) > 3:
                    growing_counter = len(line.split()) / 3
                    successive_tasks = []
                    for i in range(int(growing_counter)):
                        try:
                            fragment_pdb = line.split()[i * 3]
                            core_atom = line.split()[(i * 3) + 1]
                            if "*" in core_atom:
                                # Get fragment number and remove part of string
                                fragment_number = re.findall(r'[*]\d[*]', core_atom)[0].strip("*")
                                core_atom = core_atom.replace("*{}*".format(fragment_number), "")
                            else:
                                fragment_number = None
                            fragment_atom = line.split()[(i * 3) + 2]
                            ID = "{}{}{}".format(os.path.splitext(fragment_pdb)[0], core_atom, fragment_atom)
                            task = (fragment_pdb, core_atom, fragment_atom, ID, fragment_number)
                            successive_tasks.append(task)
                        except IndexError:
                            logger.critical(
                                "Check that the serie file {} contains: 'PDB_fragment_file'\t'PDB_core_atom_name'\t'PDB_fragment_atom_name' ".format(
                                    file))
                    list_of_instructions.append(successive_tasks)
        return list_of_instructions


    def _extract_id(self, path):
        """
        Extracts ID of the complex. Raises and error in the complex path is not correctly specified.

        Parameters
        ----------
        path

        Returns
        -------
        id

        """
        try:
            id = os.path.split("/")[-1]
            return id
        except AssertionError:
            raise WrongComplexName(f"Incorrect format of complex path: " +
                                   f"{path}." +
                                   f"It should be /example/of/path")


    def _are_hydrogens_user_defined(self, instruction):
        """
        Checks if the user has specified the hydrogens in the instructions of the serie file.

        Parameters
        ----------
        self
        instruction

        Returns
        -------
        bool

        """
        if "-" in instruction[1] or "-" in instruction[2]:
            try:
                heavy_core = instruction[1].split("-")[0]
                hydrogen_core = instruction[1].split("-")[1]
                heavy_fragment = instruction[2].split("-")[0]
                hydrogen_fragment = instruction[2].split("-")[1]
                return heavy_core, hydrogen_core, heavy_fragment, hydrogen_fragment
            except IndexError:
                raise IndexError(
                    f"Wrong growing direction in {instruction}. Both, fragment and core atoms must include one or two atoms."
                    "Ex: frag.pdb  C1  C2 |or| frag.pdb  C1-H1  C2-H2")
        else:
            return False


    def _extract_linker_atoms(self, instruction_list):
        """
        Extracts linker atoms from the serie file.

        Parameters
        ----------
        self

        Returns
        -------

        """
        for instruction in instruction_list:
            if type(instruction) == list:
                #Sequential Growing. Still checking how this works.
                pass
            elif type(instruction) != list:
                if self._are_hydrogens_user_defined(instruction):
                    core_atom, h_core, fragment_atom, h_frag = instruction[0], instruction[1], instruction[2], instruction[3]
                    return core_atom, h_core, fragment_atom, h_frag

                else:
                    core_atom, fragment_atom = instruction[1], instruction[2]
                    return core_atom, None, fragment_atom, None