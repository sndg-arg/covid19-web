# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import os
import warnings

from Bio import BiopythonWarning, BiopythonParserWarning, BiopythonDeprecationWarning, BiopythonExperimentalWarning

warnings.simplefilter('ignore', RuntimeWarning)
warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', BiopythonParserWarning)
warnings.simplefilter('ignore', BiopythonDeprecationWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)

from pdbdb.models import PDB
from pdbdb.io.FPocket2SQL import FPocket2SQL
from pdbdb.io.PDB2SQL import PDB2SQL
from bioseq.models.Biodatabase import Biodatabase
from bioseq.models.Biosequence import Biosequence
from bioseq.models.Bioentry import Bioentry

from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from Bio.PDB.PDBParser import PDBParser


class PDBIO():

    def __init__(self, pdbs_dir="/data/databases/pdb/divided/",
                 entries_path="/data/databases/pdb/entries.idx",
                 tmp="/tmp/PDBIO"):
        self.pdbs_dir = pdbs_dir
        self.entries_path = entries_path
        self.tmp = tmp

    def init(self):
        self.pdb2sql = PDB2SQL(self.pdbs_dir, self.entries_path)
        self.pdb2sql.load_entries()
        self.fpocket2sql = FPocket2SQL()
        self.fpocket2sql.create_or_get_pocket_properties()
        self.create_or_get_biodb()

    def create_or_get_biodb(self):
        self.biodb = Biodatabase.objects.get_or_create(name="PDB")[0]

    def pdb_path(self, pdb_code):
        return os.path.sep.join([self.pdbs_dir, pdb_code[1:3], "pdb" + pdb_code + ".ent"])

    def create_sequence(self, pdb_code, pdb_path):
        pdb = PDB.objects.get(code=pdb_code)

        struct = PDBParser(PERMISSIVE=1, QUIET=1).get_structure(pdb_code, pdb_path)
        for chain in struct[0].get_chains():

            residues = []
            for residue in chain.get_residues():

                if is_aa(residue, standard=True):
                    # alts = [a.get_altloc() for a in residue.get_atoms() if a.get_altloc()]
                    # if len(alts) > 1 :
                    #     print(pdb_code)
                    #     disordered_select
                    #     print("alternative residue %s from %s was removed from sequence" % (
                    #         str(residue.id), pdb_code
                    #     ))
                    # else:
                    residues.append(residue)



            if residues:
                seq = "".join([seq1(x.resname) for x in residues])
                start = str(residues[0].id[1])
                end = str(residues[-1].id[1])
                seqid = "_".join([pdb_code, chain.id, start, end])
                if not Bioentry.objects.filter(biodatabase=self.biodb, identifier=seqid).exists():
                    be = Bioentry(biodatabase=self.biodb, accession=seqid, identifier=seqid, name=pdb.code)
                    be.save()
                    Biosequence(bioentry=be, seq=seq, length=len(seq)).save()

    def process_pdb(self, pdb_code):
        assert self.pdb2sql, "PDBIO not initialized"
        pdb_code = pdb_code.lower()

        if PDB.objects.filter(code=pdb_code).exists():
            raise Exception("PDB already exists")

        self.create_sequence(pdb_code)

        pdb_path = self.pdb_path(pdb_code)
        if not os.path.exists(pdb_path):
            pdb_path = self.pdb2sql.download(pdb_code)

        self.pdb2sql.create_pdb_entry(pdb_code, pdb_path)
        self.pdb2sql.update_entry_data(pdb_code, pdb_path)
        self.fpocket2sql.load_pdb(pdb_code)
        self.fpocket2sql.run_fpocket(self.tmp)
        self.fpocket2sql.load_pockets()
