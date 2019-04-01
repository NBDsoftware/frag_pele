import argparse


def parse_arguments():
    """
            Parse user arguments

            Output: list with all the user arguments
        """
    # All the docstrings are very provisional and some of them are old, they would be changed in further steps!!
    parser = argparse.ArgumentParser(description="""Given a trajectory writes all MODELS in different files.""")
    required_named = parser.add_argument_group('required named arguments')
    # Growing related arguments
    required_named.add_argument("trajectory",
                                help="""Trajectory path""")

    args = parser.parse_args()

    return args.trajectory


def main(path_to_pdb):
    with open(path_to_pdb) as pdb:
        content = pdb.read()
    pdbs = content.split("ENDMDL")
    for n, pdb in enumerate(pdbs):
        if n != len(pdbs)-1:
            with open("{}_{}".format(path_to_pdb, n), "w") as out:
                out.write(pdb)


if __name__ == '__main__':
    traj = parse_arguments()
    main(traj)

