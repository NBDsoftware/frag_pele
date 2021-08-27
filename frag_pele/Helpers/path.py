
class PathHandler():
    def __init__(self,
                 PackagePath,
                 plop_path,
                 sch_path,
                 complex_pdb,
                 pdbout,
                 force_field,
                 ID,
                 dist_constraint=None,
                 data=None,
                 documents=None,):
        """

        Parameters
        ----------
        PackagePath
        plop_path
        sch_path
        complex_pdb
        working_dir
        pdbout
        force_field
        """
        import os
        self._current_path = current_path = os.path.abspath(".")
        self._plop_relative_path = os.path.join(PackagePath, plop_path)
        self._current_path = os.path.abspath(".")
        self._sch_python = os.path.join(sch_path, "utilities/python")
        self._pdb_basename = self._get_pdb_basename(complex_pdb)
        self._working_dir = os.path.join(current_path, "{}_{}".format(self._pdb_basename, ID))
        if not os.path.exists(self._working_dir):
            os.mkdir(self._working_dir)  # Creating a working directory for each PDB-fragment combination

        self._pdbout_folder = os.path.join(self._working_dir, pdbout)
        self._path_to_templates_generated, self._path_to_templates = self._get_path_to_templates_generated(self._working_dir, force_field)
        self._path_to_lib = os.path.join(self._working_dir, "DataLocal/LigandRotamerLibs")
        self._create_output_folder(self._working_dir)
        self._create_symbolik_links(data,documents,self._working_dir)
        self._const = self._create_constraints(dist_constraint)

    def _create_output_folder(self, working_dir):
        from frag_pele.Helpers import folder_handler
        folder_handler.check_and_create_DataLocal(working_dir=self._working_dir)

    def _get_pdb_basename(self, complex_pdb):
        pdb_basename = complex_pdb.split(".pdb")[0]  # Get the name of the pdb without extension
        if "/" in pdb_basename:
            pdb_basename = pdb_basename.split("/")[-1]  # And if it is a path, get only the name
        return pdb_basename

    def _get_path_to_templates_generated(self, working_dir, force_field):
        from frag_pele.Errors.custom_errors import WrongForceField
        import os

        if force_field == 'OFF':
            path_to_templates_generated = os.path.join(working_dir,
                                                       "DataLocal/Templates/OpenFF/Parsley/templates_generated".format(
                                                           force_field))
            path_to_templates = os.path.join(working_dir, "DataLocal/Templates/OpenFF/Parsley".format(force_field))
        elif force_field == 'OPLS2005':
            path_to_templates_generated = os.path.join(working_dir,
                                                       "DataLocal/Templates/OPLS2005/HeteroAtoms/templates_generated".format(
                                                           force_field))
            path_to_templates = os.path.join(working_dir,
                                             "DataLocal/Templates/OPLS2005/HeteroAtoms".format(force_field))
        else:
            raise WrongComplexName(f"Incorrect force field type. "
                                   f"Force field should be OFF or OPL2005.")

        return path_to_templates_generated, path_to_templates

    def _create_symbolik_links(self, data, documents, working_dir):
        from frag_pele.Helpers import helpers
        if data:
            helpers.create_symlinks(data, os.path.join(working_dir, 'Data'))
        if documents:
            helpers.create_symlinks(documents, os.path.join(working_dir, 'Documents'))

    def _create_constraints(self, dist_constraint):
        if dist_constraint:
            atom1_info, atom2_info, equil_dist = dist_constraint
            const = "\n".join(constraints.retrieve_constraints(complex_pdb, {}, {}, 5, 5, 10,
                                                               atom1_info, atom2_info, equil_dist))
            return const

