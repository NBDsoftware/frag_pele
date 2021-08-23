import os.path
import frag_pele.constants as c
from frag_pele.Errors.custom_errors import WrongComplexName


class Parametrizer:
    """

    """

    def __init__(
            self,
            complex_pdb,
            serie_file,
            iterations=c.GROWING_STEPS,
            criteria=c.SELECTION_CRITERIA,
            plop_path=c.PLOP_PATH,
            sch_path=c.SCHRODINGER,
            pele_dir=c.PATH_TO_PELE,
            contrl=c.CONTROL_TEMPLATE,
            license=c.PATH_TO_LICENSE,
            resfold=c.RESULTS_FOLDER,
            report=c.REPORT_NAME,
            traject=c.TRAJECTORY_NAME,
            pdbout=c.PDBS_OUTPUT_FOLDER,
            cpus=c.CPUS,
            distcont=c.DISTANCE_COUNTER,
            threshold=c.CONTACT_THRESHOLD,
            epsilon=c.EPSILON,
            condition=c.CONDITION,
            metricweights=c.METRIC_WEIGHTS,
            nclusters=c.NUM_CLUSTERS,
            pele_eq_steps=c.PELE_EQ_STEPS,
            restart=False,
            min_overlap=c.MIN_OVERLAP,
            max_overlap=c.MAX_OVERLAP,
            c_chain="L",
            f_chain="L",
            steps=c.STEPS,
            temperature=c.TEMPERATURE,
            seed=c.SEED,
            rotamers=c.ROTRES,
            banned=c.BANNED_DIHEDRALS_ATOMS,
            limit=c.BANNED_ANGLE_THRESHOLD,
            mae=False,
            rename=None,
            threshold_clash=None,
            steering=c.STEERING,
            translation_high=c.TRANSLATION_HIGH,
            rotation_high=c.ROTATION_HIGH,
            translation_low=c.TRANSLATOIN_LOW,
            rotation_low=c.ROTATION_LOW,
            explorative=False,
            radius_box=c.RADIUS_BOX,
            sampling_control=None,
            data=c.PATH_TO_PELE_DATA,
            documents=c.PATH_TO_PELE_DOCUMENTS,
            only_prepare=False,
            only_grow=False,
            no_check=False,
            debug=False,
            protocol=False,
            test=False,
            cov_res=None,
            dist_constraint=None,
            contraint_core=False,
            dih_constr=None,
            growing_protocol='SoftcoreLike',
            strat_growing_from=0.0,
            min_grow=0.01,
            min_sampling=0.1,
            force_field='OPLS2005',
            dih_to_constraint=None,
            srun=True,
            keep_templates=False):
        """
        Initializes parametrization to obtain FragPELE instructions and proceed with fragment growing.

        Parameters
        ----------
        complex_pdb
        serie_file
        iterations
        criteria
        plop_path
        sch_path
        pele_dir
        contrl
        license
        resfold
        report
        traject
        pdbout
        cpus
        distcont
        threshold
        epsilon
        condition
        metricweights
        nclusters
        pele_eq_steps
        restart
        min_overlap
        max_overlap
        c_chain
        f_chain
        steps
        temperature
        seed
        rotamers
        banned
        limit
        mae
        rename
        threshold_clash
        steering
        translation_high
        rotation_high
        translation_low
        rotation_low
        explorative
        radius_box
        sampling_control
        data
        documents
        only_prepare
        only_grow
        no_check
        debug
        protocol
        test
        cov_res
        dist_constraint
        contraint_core
        dih_constr
        growing_protocol
        strat_growing_from
        min_grow
        min_sampling
        force_field
        dih_to_constraint
        srun
        keep_templates
        """
        self.original_dir = os.path.abspath(os.getcwd())
        self._list_of_instructions = self._read_instructions_form_file(serie_file)
        # self.dict_traceback
        self._ID = self._extract_id(path)
        self._core_atom, self._h_core, self._fragment_atom, self._h_frag = self._extract_linker_atoms(
            self._list_of_instructions)
        self._steps, self._pele_eq_steps, self._temp = self._define_protocol_params(protocol, test)

    def _define_protocol_params(self, protocol, test):
        """
        Checks protocol type and assigns steps, pele quilibration steps and temperature parameters.

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

    def _read_instructions_form_file(self, file):
        """
        From the serie file reads the instructions.
        Parameters
        ----------
        file

        Returns
        -------
        listo_of_instructions

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
            id = path.split("/")[-1]
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