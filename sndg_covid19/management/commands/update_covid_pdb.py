import os
import subprocess as sp
from tqdm import tqdm
import requests
from datetime import datetime
from io import StringIO
import Bio.SeqIO as bpio
import tempfile
from django_tqdm import BaseCommand

from bioseq.models.Biodatabase import Biodatabase
from bioseq.models.Variant import Variant
from pdbdb.models import Residue, ResidueProperty, Property, ResidueSet, PDBResidueSet
from SNDG.Structure.StructureVariant import StructureVariant
from bioseq.models.PDBVariant import PDBVariant

from SNDG.WebServices.PDBsWS import PDBsWS


class Command(BaseCommand):
    """

    """

    DEFAULT_PDB_PATH = "/data/databases/pdb/"

    help = 'Tinny DB for documentation and testing purposes'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dbx_dict = {}

    def add_arguments(self, parser):
        parser.add_argument('--pdb_dir',
                            default=os.environ.get("DEFAULT_PDB_PATH", Command.DEFAULT_PDB_PATH))
        parser.add_argument('--force', action="store_true")

    def remove_obsolete(self, protein):

        for pdbref in list(protein.dbxrefs.all()):
            if pdbref.dbxref.dbname == "PDB":
                pdb = pdbref.dbxref.accession
                if PDBsWS.is_obsolete(pdb):
                    pdbref.delete()
                    self.stderr.write(f"deleting {pdb}")
                yield pdb

    def handle(self, *args, **options):

        bdb = Biodatabase.objects.get(name="COVID19" + Biodatabase.PROT_POSTFIX)
        with self.tqdm(bdb.entries.prefetch_related("dbxrefs").all(), total=bdb.entries.count()) as pbar:
            for protein in pbar:
                pbar.set_description(protein.accession)
                pdbs = self.remove_obsolete(protein)

                if pdbs:
                    prot_path = tempfile.NamedTemporaryFile().name
                    bpio.write(protein.to_seq_record(), prot_path, "fasta")
                    protein_variants = list(Variant.objects.filter(bioentry=protein))
                    pbar_pdbs = self.tqdm(pdbs)
                    for pdb in pbar_pdbs:

                        pbar_pdbs.set_description(pdb)
                        vs = StructureVariant(options["pdb_dir"])
                        vs.load_msa(prot_path, pdb, pdb_chain=None)

                        pvar_variants = self.tqdm(protein_variants)
                        pvar_variants.set_description("processing variants...")
                        for pv in pvar_variants:
                            pos_data = vs.residues_from_pos(pv.pos)
                            for res_data in pos_data:
                                res = Residue.objects.get(pdb__code=pdb.lower(), resid=res_data["resid"],
                                                          chain=res_data["chain"])
                                pdbvar = PDBVariant.objects.get_or_create(variant=pv, residue=res)[0]
                                # prop = Property.objects.get_or_create(name="")[0]
                                # ann = vs.annotate_resid(pdb, res.resid)
                                # if "pockets" in ann:
                                #     PDBResidueSet.objects.filter(residue_set__name=PDBResidueSet.pocket_name,
                                #                                  name=str(ann["pockets"][0]),pdb__code=pdb)
                                #     ann

                        # for _,record in df.iterrows():
                        #     res = Residue.objects.filter(pdb__code=pdb,resid=df)
