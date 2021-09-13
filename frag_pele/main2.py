import frag_pele.constants as c
from frag_pele.args_parser import ArgsParser
import os


def parse_args(args):
    """
    Command line patser.
    Parameters
    ----------
    args : list[str]
       The list of strings to parse. Default is an empty list.
       If empty, input arguments are retrieve from the command-line.

    Returns
    -------
    parsed_args : 

    """
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description='FrAG is a Fragment-based ligand growing software which performs automatically the '
                    'addition of several fragments to a core structure of the ligand in a protein-ligand complex'
    )

    required_named = parser.add_argument_group('required named arguments')

    # FrAG related arguments
    required_named.add_argument("-cp", "--complex_pdb", required=True,
                                help="""Path to the PDB file which must contain a protein-ligand complex. Its ligand will be
                            used as the core structure. Remember to rename the ligand chain with a different character in
                            order to detect it.""")
    required_named.add_argument("-sef", "--serie_file", required=True,
                                help=""" Name of the tabular file which must contain the instructions required to perform several
                            successive growings, using different fragments or different growing positions.

                            To do simple growings:
                            col1   col2    col3
                            To do successive growings:
                            (col1   col2    col3) x n_growings

                            Where col1 is the path to the PDB file of the fragment that will be added to the core structure
                            (name the chain of the fragment "L" by default).
                            Col2 is a string with the PDB atom name of the heavy atom of the core (the ligand contained in
                            the complex) where you would like to start the growing and create a new bond with the fragment.
                            And col3 is a string with the PDB atom name of the heavy atom of the fragment that will be used
                            to perform the bonding with the core.
                            """)
    parser.add_argument("--debug", action="store_true", help="Run Frag without launching PELE simulation")
    parser.add_argument("-nc", "--no_check", action="store_true", help="Don't perform the environment variables check")
    parser.add_argument("-x", "--growing_steps", type=int, default=c.GROWING_STEPS,
                        help="""Number of Growing Steps (GS). By default = {}.""".format(c.GROWING_STEPS))
    parser.add_argument("-cr", "--criteria", default=c.SELECTION_CRITERIA,
                        help="""Name of the column of the report file to select the structures that will spawn in the 
                            next GS. Additionally, this parameter will be the selection criteria to extract the best 
                            structure after completing the growing. By default = {}.""".format(c.SELECTION_CRITERIA))
    parser.add_argument("-rst", "--restart", action="store_true",
                        help="If set FrAG will continue from the last GS detected. If all GS are finished it will"
                             "restart in the equilibration phase.")
    parser.add_argument("-cc", "--c_chain", default="L", help="Chain name of the core. By default = 'L'")
    parser.add_argument("-fc", "--f_chain", default="L", help="Chain name of the fragment. By default = 'L'")
    parser.add_argument("-tc", "--clash_thr", default=None,
                        help="Threshold distance that would to classify intramolecular"
                             "clashes. If None value set, it will use the distance of the "
                             "bond between the fragment and the core.")
    parser.add_argument("-sc", "--sampling_control", default=None, help="If set, templatized control file to use in the"
                                                                        " sampling simulation.")
    parser.add_argument("-op", "--only_prepare", action="store_true", help="If set, all files to run growing are"
                                                                           " prepared, it stops before running PELE.")
    parser.add_argument("-og", "--only_grow", action="store_true", help="If set, it runs all growings of folders "
                                                                        "already prepared.")
    parser.add_argument("-cov", "--cov_res", default=None,
                        help="Set to do growing onto protein residues. Example of selection: "
                             "'A:145' (chain A and resnum 145).")

    parser.add_argument("-dist", "--dist_const", default=None, nargs="+",
                        help="Atom pairs links id to constraint and equilibrum distance."
                             " p.Ex: 'L:1:_C3_' 'A:145:_C1_' 2.5")

    parser.add_argument("-pro", "--protocol", default="SoftcoreLike",
                        choices=['SpreadHcharge', 'SoftcoreLike', 'AllLinear'],
                        help="Select growing protocol. Choose between: 'SofcoreLike', 'SpreadingHcharge', 'AllLinear'. "
                             "SofcoreLike: Charges initially set to 0. They are added in the mid GS. Then, "
                             "they grow exponentially or linearly (depending on your settings). "
                             "SpreadingHcharge: reimplementation of FragPELE1.0.0 methodology. "
                             "AllLinear: All FF parameters are linearly and equally incremented in each GS. ")
    parser.add_argument("-stf", "--st_from", type=float, default=0.0,
                        help="Lamnda value to start the growing of a fragment from. F.ex: if you"
                             " set a 9 GS simulation, setting this value to 0.3 your fragment "
                             "growing will start from the third step, but second GS. (30%% of the size and FFp).")

    parser.add_argument("--keep_templates", action="store_true", help="If set, templates will not be overwritten "
                                                                      " in templates_generated folder.")

    # Plop related arguments
    parser.add_argument("-pl", "--plop_path", default=c.PLOP_PATH,
                        help="Absolute path to PlopRotTemp.py. By default = {}".format(c.PLOP_PATH))
    parser.add_argument("-sp", "--sch_path", default=c.SCHRODINGER,
                        help="""Absolute path to Schrodinger's directory. 
                            By default = {}""".format(c.SCHRODINGER))
    parser.add_argument("-rot", "--rotamers", default=c.ROTRES, type=int,
                        help="""Rotamers threshold used in the rotamers' library. 
                                By default = {}""".format(c.ROTRES))

    # PELE configuration arguments
    parser.add_argument("-d", "--pele_dir", default=c.PATH_TO_PELE,
                        help="Complete path to Pele_serial. "
                             "By default = {}".format(c.PATH_TO_PELE))
    parser.add_argument("-c", "--contrl", default=c.CONTROL_TEMPLATE,
                        help="Path to PELE's control file templatized. By default = {}".format(c.CONTROL_TEMPLATE))
    parser.add_argument("-l", "--license", default=c.PATH_TO_LICENSE,
                        help="Absolute path to PELE's licenses folder. "
                             " By default = {}".format(c.PATH_TO_LICENSE))
    parser.add_argument("-r", "--resfold", default=c.RESULTS_FOLDER,
                        help="Name for PELE's results folder. By default = {}".format(c.RESULTS_FOLDER))
    parser.add_argument("-rp", "--report", default=c.REPORT_NAME,
                        help="Suffix name of the report file from PELE. By default = {}".format(c.REPORT_NAME))
    parser.add_argument("-tj", "--traject", default=c.TRAJECTORY_NAME,
                        help="Suffix name of the trajectory file from PELE. By default = {}".format(c.TRAJECTORY_NAME))
    parser.add_argument("-cs", "--cpus", type=int, default=c.CPUS,
                        help="Number of cores (computational power) to paralellize PELE's simulations."
                             "By default = {}".format(c.CPUS))
    parser.add_argument("-stp", "--steps", type=int, default=c.STEPS,
                        help="""Number of simulation steps inside each GS. By default = {}""".format(c.STEPS))
    parser.add_argument("-es", "--pele_eq_steps", default=c.PELE_EQ_STEPS,
                        help="Number of PELE steps in equilibration. By default = {}".format(c.PELE_EQ_STEPS))
    parser.add_argument("-miov", "--min_overlap", default=c.MIN_OVERLAP,
                        help="Minimum value of overlapping factor used in the control_file of PELE. "
                             "By default = {}".format(c.MIN_OVERLAP))
    parser.add_argument("-maov", "--max_overlap", default=c.MAX_OVERLAP,
                        help="Maximum value of overlapping factor used in the control_file of PELE."
                             " By default = {}".format(c.MAX_OVERLAP))
    parser.add_argument("-tmp", "--temperature", default=c.TEMPERATURE,
                        help="Temperature value to add in the control file. If the temperature is high more steps of "
                             "PELE will be accepted when applying the Metropolis Criteria. "
                             "By default = {}".format(c.TEMPERATURE))
    parser.add_argument("-sd", "--seed", default=c.SEED,
                        help="Seed to get the random numbers that will be used in the PELE simulation"
                             "By default = {}".format(c.SEED))
    parser.add_argument("-st", "--steering", default=c.STEERING,
                        help="Steering that will be used in the PELE simulation"
                             "By default = {}".format(c.STEERING))
    parser.add_argument("-trh", "--translation_high", default=c.TRANSLATION_HIGH,
                        help="Translation that will be used in the PELE simulation to displace the ligand. Low value."
                             "By default = {}".format(c.TRANSLATION_HIGH))
    parser.add_argument("-roth", "--rotation_high", default=c.ROTATION_HIGH,
                        help="Rad of rotation that will be used in the PELE simulation to perturb the ligand. "
                             "Low value. By default = {}".format(c.ROTATION_HIGH))
    parser.add_argument("-trl", "--translation_low", default=c.TRANSLATION_LOW,
                        help="Translation that will be used in the PELE simulation to displace the ligand. Low value."
                             "By default = {}".format(c.TRANSLATION_LOW))
    parser.add_argument("-rotl", "--rotation_low", default=c.ROTATION_LOW,
                        help="Rad of rotation that will be used in the PELE simulation to perturb the ligand. "
                             "Low value. By default = {}".format(c.ROTATION_LOW))
    parser.add_argument("-rad", "--radius_box", default=c.RADIUS_BOX,
                        help="Size of the radius to define the box in the PELE simulation where the ligand will be"
                             "perturbed. By default = {}".format(c.RADIUS_BOX))
    parser.add_argument("-dat", "--data", default=c.PATH_TO_PELE_DATA,
                        help="Path to PELE Data folder.")
    parser.add_argument("-doc", "--documents", default=c.PATH_TO_PELE_DOCUMENTS,
                        help="Path to PELE Documents folder.")
    parser.add_argument("-sr", "--srun", default=True,
                        help="If true it runs PELE with srun command, else it will run it with mpirun.")
    parser.add_argument("-core", "--constraint_core", action="store_true",
                        help="Set true to apply core constraints at template level."
                             " These atoms will be skipped from rotamers library"
                             " (only PELE minimization will be applyed on them).")
    parser.add_argument("-dhc", "--dih_constr", default=None,
                        help="Spring constant to apply dihedrals constraints in PELE conf."
                             " This constraint will not be applyied only until the 1/2 GS"
                             " to the dihedrals formed by 4 atoms of the fragment.")
    parser.add_argument("-dhl", "--dihedrals_list", default=[], type=str, nargs='+', action='append',
                        help="List of lists of atom names (in PELE format) to constraint dihedrals between them ."
                             "F.ex: '_C16 _C17 _C18 _C19'")
    parser.add_argument("-ming", "--min_grow", default=0.01,
                        help="Minimum RMS for Minimization in PELE conf that will we applied "
                             "  only until the 1/2 GS. Lowest the value stronger minimization.")
    parser.add_argument("-mins", "--min_sampling", default=0.1,
                        help="Minimum RMS for Minimization in PELE conf that will we applied "
                             "  only after the 1/2 GS. Lowest the value stronger minimization.")
    parser.add_argument("-ff", "--force_field", default="OPLS2005", choices=['OPLS2005', 'OFF'],
                        help="Atomic ForceField that will be used in PELE simulation. Choose between: "
                             "'OPLS2005', "
                             "'OFF': OpenForceField "
                             " By default: OPLS2005")

    # Clustering related arguments
    parser.add_argument("-dis", "--distcont", default=c.DISTANCE_COUNTER,
                        help="""Distance used to determine which amino acids are in contact with the ligand to generate 
                            different clusters of structures to initialize the next GS.
                            By default = {}""".format(c.DISTANCE_COUNTER))
    parser.add_argument("-ct", "--threshold", default=c.CONTACT_THRESHOLD,
                        help="Threshold distance used in the clustering. By default = {}".format(c.CONTACT_THRESHOLD))
    parser.add_argument("-e", "--epsilon", default=c.EPSILON,
                        help="An epsilon fraction of processors are distributed proportionally to the value of a metric,"
                             "and the rest are inverselyProportional distributed. A param n can be specified to only "
                             "consider the n clusters with best metric. By default = {}".format(c.EPSILON))
    parser.add_argument("-cn", "--condition", default=c.CONDITION,
                        help="Selects wether to take into account maximum or minimum values in epsilon related spawning,"
                             "values are min or max. By default = {}".format(c.CONDITION))
    parser.add_argument("-mw", "--metricweights", default=c.METRICS_WEIGHTS,
                        help="Selects how to distribute the weights of the cluster according to its metric, "
                             "two options: linear (proportional to metric) or Boltzmann weigths (proportional "
                             "to exp(-metric/T). Needs to define the temperature T. "
                             "By default = {}".format(c.METRICS_WEIGHTS))
    parser.add_argument("-ncl", "--nclusters", default=c.NUM_CLUSTERS,
                        help="Number of initial structures that we want to use in each new GS. "
                             "By default = {}".format(c.NUM_CLUSTERS))

    parser.add_argument("-pdbf", "--pdbout", default=c.PDBS_OUTPUT_FOLDER,
                        help="Folder where PDBs selected to spawn in the next GS will be stored."
                             "By default = {}".format(c.PDBS_OUTPUT_FOLDER))
    parser.add_argument("-ban", "--banned", default=c.BANNED_DIHEDRALS_ATOMS, type=str, nargs='+', action='append',
                        help="List of lists of quartets of PDB atom names that form the banned dihedrals."
                             "By default = {}".format(c.BANNED_DIHEDRALS_ATOMS))
    parser.add_argument("-lim", "--limit", default=c.BANNED_ANGLE_THRESHOLD, type=int,
                        help="Maximum degrees that can accept the banned dihedral."
                             "By default = {}".format(c.BANNED_ANGLE_THRESHOLD))

    # Protocol argument
    parser.add_argument("-HT", "--highthroughput", action="store_true",
                        help="Run frag pele high-throughput mode")

    parser.add_argument("-EX", "--explorative", action="store_true",
                        help="Run frag pele explorative mode: sampling simulation with high movement of the ligand.")

    parser.add_argument("--test", action="store_true", help="run test config")

    # Output format option
    parser.add_argument("--mae", action="store_true",
                        help="Retrieve .mae files intead of pdbs")

    # Others
    parser.add_argument("--rename", action="store_true",
                        help="Avoid core renaming")

    parsed_args = parser.parse_args(args)
    return parsed_args


def run_fragpele(parameters):
    """
    Main function to prepare and launch the simulation.
    1) Create variables
    2) Create folders and define paths
    3) Prepare ligand and growing results.
    4) Launch fragment growing simulations.
    5) Clustering.
    6) Equilibration.

    Parameters
    ----------
    args

    Returns
    -------

    """
    from frag_pele.parametrizer import ParametersBuilder
    from frag_pele.Growing.reduction import FragmentReduction
    from frag_pele.Helpers.path import PathHandler
    from frag_pele.Growing.pdbs import PDBHandler
    from frag_pele.Utils.rdkit_wrapper import RDKitWrapper
    from logging.config import fileConfig
    import logging
    import time
    import sys
    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "AdaptivePELE_repo"))

    # Calling configuration file for log system
    FilePath = os.path.abspath(__file__)
    PackagePath = os.path.dirname(FilePath)
    LogPath = os.path.join(PackagePath, c.CONFIG_PATH)
    fileConfig(LogPath)

    # Getting the name of the module for the log system
    logger = logging.getLogger(__name__)

    # Start timer
    start_time = time.time()

    # Initialize parameters
    sim_parameters = ParametersBuilder()
    PDBreader = PDBHandler()

    ligand_pdb = PDBreader._get_pdb_components(parameters.complex_pdb)
    instructions = sim_parameters._read_instructions_from_file(parameters.serie_file)
    fragment_pdb, core_atom, h_core, fragment_atom, h_frag =  sim_parameters._extract_linker_atoms(instructions)
    fragment_mol, ligand_mol = PDBreader._create_rdkit_molecules(ligand_pdb, fragment_pdb)

    # Preparation
    #    a) Fragment reduction parameters
    Reduction = FragmentReduction()
    lam_initial = Reduction._define_lam(parameters.iterations, parameters.start_growing_from)
    print(f"Reducing fragment size to the {lam_initial * 100} %")

    #    b) Fragment Merging and reduction
    grown_ligand = PDBreader._merge_ligand_and_fragment(fragment_mol,ligand_mol,core_atom,fragment_atom)

    #    c) Create Templates
    #    d) Get Templates
    #    e) Set box center from ligand COM
    #    f) Correct Templates

    # Growing

    # Clustering

    # Equilibration

    # Compute Time



def main(args):
    """
    Reads the command-line arguments and runs FragPELE.
    """
    parameters = ArgsParser(args)
    run_fragpele(parameters)

if __name__ == "__main__":
    import sys
    args = parse_args(sys.argv[1:])
    parameters = ArgsParser(args)
    main(args)