import os
import traceback
import warnings

from Bio import BiopythonWarning, BiopythonParserWarning, BiopythonDeprecationWarning, BiopythonExperimentalWarning
from django.core.management.base import BaseCommand, CommandError
from tqdm import tqdm

from pdbdb.models import PDB
from SNDG.Structure.PDBs import PDBs

warnings.simplefilter('ignore', RuntimeWarning)
warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', BiopythonParserWarning)
warnings.simplefilter('ignore', BiopythonDeprecationWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)

from SNDG.Structure.StructureAnnotator import StructureAnnotator
from pdbdb.models import ResidueSet, PDBResidueSet, Residue, ResidueSetResidue, Property, ResidueSetProperty


class Command(BaseCommand):
    help = 'Process and load PDB pocket information'

    DEFAULT_PDB_PATH = "/data/databases/pdb/"

    def add_arguments(self, parser):
        parser.add_argument('--pdbs_dir', default=os.environ.get("DEFAULT_PDB_PATH", Command.DEFAULT_PDB_PATH))
        parser.add_argument('--pdb', default=None, help="4 letter pdb code, process all in the db if none selected")
        parser.add_argument('--force', action="store_true", help="runs fpockets and overwrites the current result")

    def handle(self, *args, **options):

        qs = PDB.objects.filter(code__iexact=options["pdb"]) if options["pdb"] else PDB.objects.all()
        total = qs.count()


        with tqdm(qs, total=total) as pbar:

            for pdb in pbar:
                pbar.set_description(pdb.code)
                ligand = [x["resname"] for x in
                          Residue.objects.filter(pdb=pdb).exclude(type__in="R").exclude(resname="HOH").values(
                              "resname")]

                sa = StructureAnnotator()
                binding_rs = ResidueSet.objects.get_or_create(name=ResidueSet.BINDING_SITE)[0]
                bindind_data = sa.load_pdb_binding_data(pdb.code)
                evcode = Property.objects.get_or_create(name="evidence_code")[0]
                source = Property.objects.get_or_create(name="source")[0]
                bound_ligand = Property.objects.get_or_create(name=Property.bound_ligand)[0]


                for bind_site in bindind_data:
                    rs = (binding_rs if "BINDING SITE" in bind_site["details"] else
                          ResidueSet.objects.get_or_create(name=bind_site["details"])[0])

                    PDBResidueSet.objects.filter(residue_set=rs,pdb=pdb, name=bind_site["site_id"]).delete()
                    pdbrs = PDBResidueSet.objects.get_or_create(residue_set=rs, pdb=pdb, name=bind_site["site_id"],
                                                                description="" if bind_site["details"] == rs.name else
                                                                bind_site["details"])[0]

                    ResidueSetProperty.objects.get_or_create(pdbresidue_set=pdbrs, property=evcode,
                                                             value_text=bind_site["evidence_code"])
                    ResidueSetProperty.objects.get_or_create(pdbresidue_set=pdbrs, property=source,
                                                             value_text="PDBe")


                    bind_site_ligands = []
                    for residue in bind_site["site_residues"]:
                        if residue["chem_comp_id"] != "HOH":
                            db_res = Residue.objects.get(pdb=pdb, resid=residue["author_residue_number"],
                                                         chain=residue["chain_id"])
                            if db_res.resname in ligand:
                                bind_site_ligands.append(str(db_res))
                            ResidueSetResidue.objects.get_or_create(residue=db_res, pdbresidue_set=pdbrs)

                    for ligand in bind_site_ligands:
                        ResidueSetProperty.objects.get_or_create(pdbresidue_set=pdbrs, property=bound_ligand,
                                                                 value_text=ligand)
