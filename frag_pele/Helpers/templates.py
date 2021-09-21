
import os
import frag_pele
from peleffy.topology import Molecule, Topology, RotamerLibrary
from peleffy.forcefield import OpenForceField, OPLS2005ForceField
from peleffy.forcefield.selectors import ForceFieldSelector
from peleffy.template import Impact
from peleffy.utils import get_data_file_path
from frag_pele.Helpers import folder_handler
from peleffy.utils import Logger


class TemplateGenerator:
    def __init__(self):
        pass

    def _get_datalocal(self,sim_parameters):
        folder_handler.check_and_create_DataLocal(working_dir=sim_parameters.output_dir)
        self._get_template_and_rot(sim_parameters)

    def _create_template_path(self, sim_parameters):
        if sim_parameters.templates_generated:
            templ_string = "templates_generated"
        else:
            templ_string = ""
        if sim_parameters.forcefield == 'OPLS2005' and sim_parameters.protein:
            sim_parameters.output_dir = os.path.join(sim_parameters.output_dir,
                                'DataLocal/Templates/OPLS2005/Protein',
                                templ_string, sim_parameters.template_name.lower())
        if 'openff_' in sim_parameters.forcefield and sim_parameters.protein:
            sim_parameters.output_dir = os.path.join(sim_parameters.output_dir,
                                'DataLocal/Templates/OpenFF/Protein',
                                templ_string, aminoacid.lower())
        if sim_parameters.forcefield == 'OPLS2005' and not sim_parameters.protein:
            sim_parameters.output_dir = os.path.join(sim_parameters.output_dir,
                                'DataLocal/Templates/OPLS2005/HeteroAtoms',
                                templ_string, sim_parameters.template_name.lower() + "z")
        if 'openff_' in sim_parameters.forcefield and not sim_parameters.protein:
            sim_parameters.output_dir = os.path.join(sim_parameters.output_dir,
                                'DataLocal/Templates/OpenFF/Parsley',
                                templ_string, sim_parameters.template_name.lower() + "z")
        return sim_parameters.output_dir

    def _get_template_and_rot(self,sim_parameters):

        for pdb in [sim_parameters.initial_ligand_pdb, sim_parameters.grown_ligand]:
            breakpoint()
            pregrow_folder, sim_parameters.initial_ligand_pdb = os.path.split(sim_parameters.initial_ligand_pdb)
            breakpoint()
            if pregrow_folder == '':
                pregrow_folder = './'
            prepared_ligand_name = pdb.split(".pdb")[0]  + ".pdb"
            currdir = os.getcwd()

            # Check if the residue is an amino-acid from the library
            path = os.path.dirname(frag_pele.__file__)
            aa_pdb = os.path.join(path, f"Templates/Aminoacids/{sim_parameters.aminoacid_type}.pdb")
            if not os.path.exists(aa_pdb) or sim_parameters.aminoacid_type == None:
                prepared_ligand_path = os.path.join(pregrow_folder, prepared_ligand_name)
                print(f"{aa_pdb} does not exist, using {prepared_ligand_path} instead")
            else:
                prepared_ligand_path = aa_pdb
            # Check if the output path exist to dont repeat calculations
            if not os.path.exists(prepared_ligand_path):
                breakpoint()
                os.chdir(pregrow_folder)
                #prepare_pdb(pdb_in=sim_parameters.initial_ligand_pdb,
                #            pdb_out=out,
                #            sch_path=sch_path)
                os.chdir(currdir)
            #os.environ['SCHRODINGER'] = sch_path
            template_path = self._create_template_path(sim_parameters)
            if sim_parameters.aminoacid:
                print("Aminoacid template")
                if not sim_parameters.contrained_atoms:
                    contraints = [' CA ', ' C  ', ' N  ']
                else:
                    contraints = sim_parameters.contrained_atoms
                m = Molecule(prepared_ligand_path,
                             core_constraints=contraints,
                             rotamer_resolution=sim_parameters.rot_res)
                print(f"Using PDB from {prepared_ligand_path} to generate the template!")
            else:
                print("Heteroatom template")
                if not sim_parameters.contrained_atoms:
                    m = Molecule(prepared_ligand_path,
                                 rotamer_resolution=sim_parameters.rot_res)
                else:
                    m = Molecule(prepared_ligand_path,
                                 core_constraints=sim_parameters.contrained_atoms,
                                 rotamer_resolution=sim_parameters.rot_res)
            ff_select = ForceFieldSelector()
            ff = ff_select.get_by_name(sim_parameters.forcefield)
            parameters = ff.parameterize(m)
            topology = Topology(m, parameters)
            if 'openff_' in sim_parameters.forcefield:  # Not tested yet
                from peleffy.solvent import OBC2
                solvent = OBC2(topology)
                solvent.to_file(os.path.join(sim_parameters.outdir,
                                             "DataLocal/OBC/ligandParams.txt"))

            impact = Impact(topology)
            impact.to_file(template_path)
            print("Template in {}.".format(template_path))
            rot_path = os.path.join(sim_parameters.outdir,
                                    "DataLocal/LigandRotamerLibs/{}.rot.assign".format(sim_parameters.template_name.upper()))
            rotamer_library = RotamerLibrary(m)
            rotamer_library.to_file(rot_path)
            print("Rotamer library stored in {}".format(rot_path))

