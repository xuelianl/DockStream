import os
import tempfile
import shutil
import multiprocessing
from copy import deepcopy
from typing import Optional, List, Any

import rdkit.Chem as Chem
from pydantic import BaseModel, Field
from typing_extensions import Literal

from dockstream.core.Schrodinger.Glide_docker import Parallelization
from dockstream.core.docker import Docker
from dockstream.core.Dock37.Dock37_result_parser import Dock37ResultParser
from dockstream.utils.enums.logging_enums import LoggingConfigEnum
from dockstream.utils.execute_external.Dock37 import Dock37Executor
# from dockstream.utils.enums.AutodockVina_enums import AutodockVinaExecutablesEnum, AutodockVinaOutputEnum, AutodockResultKeywordsEnum
from dockstream.utils.execute_external.OpenBabel import OpenBabelExecutor
from dockstream.utils.enums.OpenBabel_enums import OpenBabelExecutablesEnum
from dockstream.utils.enums.RDkit_enums import RDkitLigandPreparationEnum
from dockstream.utils.general_utils import gen_temp_file

# from dockstream.utils.translations.molecule_translator import MoleculeTranslator
from dockstream.utils.dockstream_exceptions import DockingRunFailed

_LP = RDkitLigandPreparationEnum()
# _RKA = AutodockResultKeywordsEnum()
# _BEE = OpenBabelExecutablesEnum()
_LE = LoggingConfigEnum()
# _ROE = AutodockVinaOutputEnum()
# _EE = AutodockVinaExecutablesEnum()


class Dock37Parameters(BaseModel):
    check_dir: Optional[str] = None
    dockfiles_indock_dir: Optional[List[str]] = None
    prefix_execution: Optional[str] = None
    binary_location: Optional[str] = None

    # def get(self, key: str) -> Any:
    #     """Temporary method to support nested_get"""
    #     return self.dict()[key]


class Dock37(Docker, BaseModel):
    """Interface to the "UCSF DOCK3.7" backend."""

    backend: Literal["Dock37"] = "Dock37"
    parameters: Dock37Parameters

    _DOCK37_executor: Dock37Executor = None
    # _OpenBabel_executor: OpenBabelExecutor = None

    class Config:
        underscore_attrs_are_private = True

    def __init__(self, **data):
        super().__init__(**data)

    def _initialize_executors(self):
        """Initialize executors and check if they are available."""

        self._DOCK37_executor = Dock37Executor(
            prefix_execution=self.parameters.prefix_execution, 
            binary_location=self.parameters.binary_location
        )
        if not self._DOCK37_executor.is_available():
            raise DockingRunFailed("Cannot initialize DOCK 3.7 docker, as backend is not available - abort.")
        self._logger.log(f"Checked DOCK 3.7 backend availability (prefix_execution={self.parameters.prefix_execution}).",
                         _LE.DEBUG)

    #     # initialize the executor for all "OpenBabel" related calls and also check if it is available
    #     # note, that while there is an "OpenBabel" API (python binding) which we also use, the match to the binary
    #     # options is not trivial; thus, use command-line here
    #     self._OpenBabel_executor = OpenBabelExecutor()
    #     if not self._OpenBabel_executor.is_available():
    #         raise DockingRunFailed(
    #             "Cannot initialize OpenBabel external library, which should be part of the environment - abort.")

    def _get_score_from_conformer(self, conformer):
        return float(conformer.GetProp('dock37_score'))

    def add_molecules(self, molecules: list):
        """This method overrides the parent class, docker.py add_molecules method. This method appends prepared
        ligands to a list for subsequent docking. Note, DOCK 3.7 needs db2 format files.

        :param molecules: A list that is to contain all prepared ligands for subsequent docking
        :type molecules: list
        :raises NotImplementedError: Each backend must override the parent class, docker.py add_molecules method.
            Inability to do so or a bug causing incorrect implementation will raise a NotImplementedError
        """
        # mol_trans = MoleculeTranslator(self.ligands, force_mol_type=_LP.TYPE_DB2)
        # mol_trans.add_molecules(molecules)
        # self.ligands = mol_trans.get_as_rdkit()
        self.ligands = [mol.get_clone() for mol in molecules]
        self._docking_performed = False

    # def _write_molecule_to_pdbqt(self, path, molecule) -> bool:
    #     # generate temporary copy as PDB
    #     temp_pdb = gen_temp_file(suffix=".pdb")
    #     Chem.MolToPDBFile(mol=molecule, filename=temp_pdb)

    #     # Note: In contrast to the target preparation,
    #     # we will use a tree-based flexibility treatment here -
    #     # thus, the option "-xr" is NOT used.
    #     arguments = [temp_pdb,
    #                  _BEE.OBABEL_OUTPUT_FORMAT_PDBQT,
    #                  "".join([_BEE.OBABEL_O, path]),
    #                  _BEE.OBABEL_PARTIALCHARGE, _BEE.OBABEL_PARTIALCHARGE_GASTEIGER]
    #     obabel_execution_result = self._OpenBabel_executor.execute(command=_BEE.OBABEL,
    #                                      arguments=arguments,
    #                                      check=False)
    #     self._logger.log(f"obabel convert ligand.pdb to ligand.pdbqt output: {obabel_execution_result.stdout}.",
    #                      _LE.INFO)
    #     self._logger.log(f"obabel convert ligand.pdb to ligand.pdbqt error: {obabel_execution_result.stderr}.",
    #                      _LE.INFO)

    #     if os.path.exists(temp_pdb):
    #         os.remove(temp_pdb)
    #     else:
    #         print(f"Can not delete the file {temp_pdb} as it doesn't exists")

    #     if os.path.exists(path):
    #         return True
    #     else:
    #         return False

    def _generate_temporary_input_output_files(self, start_indices, sublists):
        # in case singletons are handed over, wrap them in a list for "zipping" later
        if not isinstance(start_indices, list):
            start_indices = [start_indices]
        if not isinstance(sublists, list):
            sublists = [sublists]

        tmp_output_dirs = []
        # tmp_input_paths = []
        # tmp_output_paths = []
        # ligand_identifiers = []
        db2_paths = []
        for start_index, sublist in zip(start_indices, sublists):
            cur_tmp_output_dir = tempfile.mkdtemp(dir='/tmp/xli02/', prefix='dock37')
            cur_tmp_docking_dir = f'{cur_tmp_output_dir}/docking'
            os.mkdir(cur_tmp_docking_dir)
            cur_tmp_input_db2_index = f'{cur_tmp_docking_dir}/split_database_index'
            # cur_tmp_output_sdf = gen_temp_file(prefix=str(start_index), suffix=".sdf", dir=cur_tmp_output_dir)
            with open(cur_tmp_input_db2_index, 'w') as sp_f:
                for ligand in sublist:
                    # write-out the temporary input file
                    if ligand.get_molecule() is None:
                        # if os.path.isdir(cur_tmp_output_dir):
                        #     shutil.rmtree(cur_tmp_output_dir)
                        continue
                    db2_path = ligand.get_molecule()
                    sp_f.write(db2_path + '\n')
                    # self._write_molecule_to_pdbqt(cur_tmp_input_pdbqt, mol)
                    # ligand_identifiers.append(ligand.get_identifier())
                    db2_paths.append(db2_path)
            tmp_output_dirs.append(cur_tmp_output_dir)
            # tmp_input_paths.append(cur_tmp_input_db2_index)
            # tmp_output_paths.append(cur_tmp_output_sdf)
        # return tmp_output_dirs, tmp_input_paths, tmp_output_paths, ligand_identifiers
        return tmp_output_dirs, db2_paths


    def _dock(self, number_cores):

        self._initialize_executors()

        start_indices, sublists = self.get_sublists_for_docking(number_cores=number_cores)
        number_sublists = len(sublists) # 1
        self._logger.log(f"Split ligands into {number_sublists} sublists for docking.", _LE.DEBUG)

        if not os.path.exists(self.parameters.dockfiles_indock_dir[0]):
            raise DockingRunFailed("Specified dockfiles and INDOCK directory to target (receptor) does not exist - abort.")

        sublists_submitted = 0
        slices_per_iteration = min(number_cores, number_sublists)
        while sublists_submitted < len(sublists):
            upper_bound_slice = min((sublists_submitted + slices_per_iteration), len(sublists))
            cur_slice_start_indices = start_indices[sublists_submitted:upper_bound_slice]
            cur_slice_sublists = sublists[sublists_submitted:upper_bound_slice]

            # generate paths and initialize molecules (so that if they fail, this can be covered)
            # tmp_output_dirs, tmp_input_paths, tmp_output_paths, \
            # ligand_identifiers = self._generate_temporary_input_output_files(cur_slice_start_indices,
            #                                                                  cur_slice_sublists)
            tmp_output_dirs, db2_paths = self._generate_temporary_input_output_files(cur_slice_start_indices, cur_slice_sublists)

            # run in parallel; wait for all subjobs to finish before proceeding
            processes = []
            for chunk_index in range(len(tmp_output_dirs)):
                p = multiprocessing.Process(target=self._dock_subjob, args=(tmp_output_dirs[chunk_index],))
                processes.append(p)
                p.start()
            for p in processes:
                p.join()

            # add the number of input sublists rather than the output temporary folders to account for cases where
            # entire sublists failed to produce an input structure
            sublists_submitted += len(cur_slice_sublists)



            # parse the resulting sdf files
            for tmp_output_dir in tmp_output_dirs:
                # add conformations
                tmp_docked_poses = f'{tmp_output_dir}/poses.mol2'
                tmp_docked_unique_scores = f'{tmp_output_dir}/extract_all.sort.uniq.txt'
                tmp_docked_poses_unicon_sdf = f'{tmp_output_dir}/poses.mol2.sdf'
                if not os.path.isfile(tmp_docked_poses) or os.path.getsize(tmp_docked_poses) == 0:
                    continue

                # read pose.mol2 as rdkit mol obj
                # with open(tmp_docked_poses, 'r') as pose_f:
                #     pose_lines = pose_f.readlines()
                # mol_s=[]
                # header = True
                # for pose_line in pose_lines:
                #     if pose_line.startswith("##########                 Name"):
                #         if header:
                #             header = False
                #         else:
                #             mol_s.append(single_mol)
                #         single_mol = [pose_line]
                #     single_mol.append(pose_line)
                # mol_s.append(single_mol)
                # rdkit_mols = []
                # for mol in mol_s:
                #     block = ",".join(mol).replace(',','')
                #     rdkit_mol = Chem.MolFromMol2Block(block, removeHs=False,sanitize=False)
                #     rdkit_mols.append(rdkit_mol)

                # idx_to_score = {}
                # for idx, molecule in enumerate(rdkit_mols):
                #     if molecule is None:
                #         continue
                #     score = self._extract_score_from_Dock37Result(molecule=molecule)
                #     idx_to_score[idx] = score

                # unicon_sdf_path = gen_temp_file(suffix=".sdf")
                os.system(f'/pubhome/xli02/Download/ZBH/unicon/unicon -i {tmp_docked_poses} -o {tmp_docked_poses_unicon_sdf} -v 0')

                if os.path.getsize(tmp_docked_poses_unicon_sdf) == 0:
                    self._logger.log(f"The size of {tmp_docked_poses_unicon_sdf} (converted by {tmp_docked_poses}) is 0.", _LE.WARNING)
                    continue

                for molecule in Chem.SDMolSupplier(tmp_docked_poses_unicon_sdf, removeHs=False):
                    if molecule is None:
                        continue

                    # extract the score from the Dock 3.7 output and update some tags
                    cur_identifier = molecule.GetProp('_Name').split()[0]
                    score = self._extract_score_from_Dock37Result(cur_identifier=cur_identifier, tmp_docked_unique_scores=tmp_docked_unique_scores)
                    molecule.SetProp("_Name", cur_identifier)
                    molecule.SetProp('dock37_score', score)
                    # molecule.ClearProp(_ROE.REMARK_TAG)

                    # add molecule to the appropriate ligand
                    for ligand in self.ligands:
                        if ligand.get_identifier() == cur_identifier:
                            ligand.add_conformer(molecule)
                            break

                # if os.path.exists(unicon_sdf_path):
                #     os.remove(unicon_sdf_path)
                # else:
                #     self._logger.log(f"Can not delete the file {unicon_sdf_path} as it doesn't exists", _LE.WARNING)

            # clean-up
            for path in tmp_output_dirs:
                shutil.rmtree(path)
            shutil.rmtree(os.path.dirname(os.path.dirname(db2_paths[0])))
            self._log_docking_progress(number_done=sublists_submitted, number_total=number_sublists)

        # the conformers are already sorted, but some tags are missing
        # -> <ligand_number>:<enumeration>:<conformer_number>
        for ligand in self.ligands:
            ligand.add_tags_to_conformers()

        # log any docking fails
        self._docking_fail_check()

        # generate docking results as dataframe
        result_parser = Dock37ResultParser(ligands=[ligand.get_clone() for ligand in self.ligands])
        self._df_results = result_parser.as_dataframe()

        # set docking flag
        self._docking_performed = True

    def _extract_score_from_Dock37Result(self, cur_identifier, tmp_docked_unique_scores) -> str:
        with open(tmp_docked_unique_scores, 'r') as score_f:
            score_lines = score_f.readlines()
        result_line = [line for line in score_lines if cur_identifier in line][0]
        parts = result_line.split()
        return parts[-1]

    def _dock_subjob(self, tmp_output_dir):

        # set up arguments list and execute
        # TODO: support "ensemble docking" - currently, only the first entry is used
        # search_space = self.parameters.search_space
        # arguments = [self.parameters.dockfiles_indock_dir[0], self.parameters.check_dir, f'&> {tmp_output_dir}/dock.log']
        arguments = [self.parameters.dockfiles_indock_dir[0], self.parameters.check_dir]


        execution_result = self._DOCK37_executor.execute(command='bash /pubhome/xli02/project/mpro/DockStream/dockstream/core/Dock37/dock37.sh',
                                                      arguments=arguments,
                                                      check=True, location=tmp_output_dir)
        # bash /pubhome/xli02/project/mpro/DockStream/dockstream/core/Dock37/dock37.sh $dockfiles_indock_dir $check_dir
        self._logger.log(f"DOCK 3.7 output: {execution_result.stdout}.",
                         _LE.INFO)
        self._logger.log(f"DOCK 3.7 error: {execution_result.stderr}.",
                         _LE.INFO)
        if execution_result.returncode != 0:
            msg = f"Could not dock with DOCK 3.7, error message: {execution_result.stdout}."
            self._logger.log(msg, _LE.ERROR)
            raise DockingRunFailed(msg)
        
        tmp_docked_poses = f'{tmp_output_dir}/poses.mol2'
        self._delay4file_system(path=tmp_docked_poses)

        # translate the parsed output PDBQT into an SDF
        # arguments = [tmp_pdbqt_docked,
        #              _BEE.OBABLE_INPUTFORMAT_PDBQT,
        #              _BEE.OBABEL_OUTPUT_FORMAT_SDF,
        #              "".join([_BEE.OBABEL_O, output_path_sdf])]
        # obabel_execution_result = self._OpenBabel_executor.execute(command=_BEE.OBABEL,
        #                                  arguments=arguments,
        #                                  check=False)
        # self._logger.log(f"obabel convert docked.pdbqt to docked.sdf output: {obabel_execution_result.stdout}.",
        #                  _LE.INFO)
        # self._logger.log(f"obabel convert docked.pdbqt to docked.sdf error: {obabel_execution_result.stderr}.",
        #                  _LE.INFO)
        # self._delay4file_system(path=output_path_sdf)

    def write_docked_ligands(self, path, mode="all"):
        """This method overrides the parent class, docker.py write_docked_ligands method. This method writes docked
        ligands binding poses and conformers to a file. There is the option to output the best predicted binding pose
        per ligand, the best predicted binding pose per enumeration, or all the predicted binding poses

        :param path: Contains information on results output path
        :type path: string
        :param mode: Determines whether the output contains the best predicted binding pose per ligand, the best
            predicted binding pose per enumeration, or all the predicted binding poses
        :type mode: string, optional, default value is "all". Other possible values are "best_per_ligand" and
            "best_per_enumeration"
        :raises DockingRunFailed Error: This error is raised if the docking run has not been performed
        :raises OpenEye (OE) Fatal Error: This error is raised if the output file was unable to be created. Issues may
            be due to problems with the ligand structure
        :raises ValueError: This error is raised if the ligands are neither RDkit nor OpenEye readable
        """
        self._write_docked_ligands(path, mode, mol_type=_LP.TYPE_RDKIT)
