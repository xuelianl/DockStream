import os
import tempfile
from copy import deepcopy

from typing import Optional, List
from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from typing_extensions import Literal

from dockstream.utils.dockstream_exceptions import LigandPreparationFailed

from dockstream.core.ligand_preparator import LigandPreparator, _LE
from dockstream.utils.enums.RDkit_enums import RDkitLigandPreparationEnum
from dockstream.utils.smiles import to_mol
from dockstream.core.ligand.ligand import Ligand

from dockstream.utils.general_utils import gen_temp_file

_LP = RDkitLigandPreparationEnum()

# This class is inspired / based on code written by
# Peter Schmidtke (https://github.com/Discngine/rdkit_tethered_minimization), which is accessible under the MIT license
# and the blogspot entry called "more-on-constrained-embedding" from the RDkit guys.


class Db2_converterParameters(BaseModel):
    max_conf: Optional[int] = 1000
    sampletp: Optional[bool] = False
    checkstereo: Optional[bool] = False
    useff: Optional[bool] = False
    db2_method: Optional[str] = 'conformator'


class Db2_converter(LigandPreparator, BaseModel):
    """Class that deals with all the preparatory steps needed before actual docking using "rDock" can commence."""

    type: Literal["db2_converter"] = "db2_converter"
    parameters: Db2_converterParameters = Db2_converterParameters()

    # class Config:
    #     underscore_attrs_are_private = True

    def __init__(self, **data):
        super().__init__(**data)

    def _load_references(self):
        references = []
        ref_format = self.align.reference_format.upper()
        for path in self.align.reference_paths:
            if ref_format == _LP.ALIGN_REFERENCE_FORMAT_PDB:
                ref_mol = Chem.MolFromPDBFile(path, sanitize=True)
                ref_mol.SetProp("_Name", os.path.basename(path))
                references.append(ref_mol)
            elif ref_format == _LP.ALIGN_REFERENCE_FORMAT_SDF:
                mol_supplier = Chem.SDMolSupplier(path)
                for mol in mol_supplier:
                    references.append(mol)
            else:
                raise IOError("Specified format not supported.")
        if len(references) == 0:
            raise LigandPreparationFailed("No reference molecules could be loaded with path(s) specified.")
        self._references = references
        self._logger.log(f"Stored {len(references)} reference molecules.", _LE.DEBUG)

    def _smiles_to_molecules(self, ligands: List[Ligand]) -> List[Ligand]:
        for lig in ligands:
            mol = to_mol(lig.get_smile())
            lig.set_molecule(mol)
            lig.set_mol_type(_LP.TYPE_RDKIT)
        return ligands

    def generate3Dcoordinates(self, converged_only=False):
        """Method to generate 3D coordinates, in case the molecules have been built from SMILES."""

        for lig in self.ligands:
            lig.set_molecule(None)
            lig.set_mol_type(None)
        ligand_list = self._smiles_to_molecules(deepcopy(self.ligands))

        failed = 0
        succeeded = 0

        # tmp_dir = tempfile.mkdtemp(dir='/tmp/xli02/')
        # smi_paths = []
        # for idx, lig in enumerate(self.ligands):
        #     smi_path = f'{tmp_dir}/{idx}.smi'
        #     with open(smi_path, 'w') as f:
        #         f.write(lig.get_smile() + " " + lig.get_identifier() + "\n")
        #     smi_paths.append(smi_path)

        if not os.path.exists('/tmp/xli02/'):
            os.mkdir('/tmp/xli02/')
        tmp_dir = tempfile.mkdtemp(dir='/tmp/xli02/')
        all_smi_path = f'{tmp_dir}/all.smi'
        with open(all_smi_path, 'w') as f:
            for lig in self.ligands:
                f.write(lig.get_smile() + " " + lig.get_identifier() + "\n")

        os.system(f"cwd=`pwd`; cd {tmp_dir}; cat {all_smi_path} |parallel -I line -k 'timeout 3600 bash /pubhome/xli02/project/mpro/DockStream/dockstream/core/db2_converter/run_db2_converter.sh line {self.parameters.max_conf} {tmp_dir}/input {self.parameters.db2_method} {self.parameters.checkstereo} {self.parameters.useff} {self.parameters.sampletp}'; cd $cwd")
        # os.system(f"cwd=`pwd`; cd /tmp/xli02/tmporasrb3g/; cat /tmp/xli02/tmporasrb3g/test.smi |parallel -I line -k 'bash /pubhome/xli02/project/mpro/test/test_dock_reinvent/test_db2_converter/run_db2_converter.sh line 1000 False False False'; cd $cwd")


        for idx, lig_obj in enumerate(ligand_list):
            if not os.path.exists(f'{tmp_dir}/input/{lig_obj.get_identifier()}.db2.gz'):
                self._logger.log(f"Could not embed molecule number {lig_obj.get_ligand_number()} (smile: {lig_obj.get_smile()}) - no 3D coordinates generated.",
                                    _LE.DEBUG)
                failed += 1
                continue

            self.ligands[idx] = Ligand(smile=lig_obj.get_smile(),
                                       original_smile=lig_obj.get_original_smile(),
                                       ligand_number=lig_obj.get_ligand_number(),
                                       enumeration=lig_obj.get_enumeration(),
                                       molecule=f'{tmp_dir}/input/{lig_obj.get_identifier()}.db2.gz',
                                       mol_type=lig_obj.get_mol_type(),
                                       name=lig_obj.get_name())
            succeeded += 1

        if failed > 0:
            self._logger.log(f"Of {len(self.ligands)}, {failed} could not be embedded.",
                             _LE.WARNING)
        self._logger.log(f"In total, {succeeded} ligands were successfully embedded (db2).", _LE.DEBUG)

    def align_ligands(self):
        self.ligands = self._align_ligands_with_RDkit_preparator(self.ligands)
