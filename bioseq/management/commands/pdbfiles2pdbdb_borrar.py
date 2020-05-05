import math
import os
import warnings

import pandas as pd
from Bio import BiopythonWarning, BiopythonParserWarning, BiopythonDeprecationWarning, BiopythonExperimentalWarning
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa
from django.core.management.base import BaseCommand, CommandError
from django.db import transaction
from tqdm import tqdm

from pdbdb.models import PDB, Residue, Atom

warnings.simplefilter('ignore', RuntimeWarning)
warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', BiopythonParserWarning)
warnings.simplefilter('ignore', BiopythonDeprecationWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)


def iterpdbs(pdbs_dir, pdb_extention=".ent"):
    for index_dir in os.listdir(pdbs_dir):
        if len(index_dir) == 2:
            for x in os.listdir(pdbs_dir + "/" + index_dir):
                if x.endswith(pdb_extention):
                    yield x[3:7], pdbs_dir + "/" + index_dir + "/" + x


def entries_df(entries_path):
    entries_columns = ["IDCODE", "HEADER", "ACCESSIONDATE", "COMPOUND", "SOURCE", "AUTHORS", "RESOLUTION",
                       "EXPERIMENT"]
    return pd.read_table(entries_path, skiprows=[0, 1, 2], sep='\t', names=entries_columns)


class Command(BaseCommand):
    help = 'Loads the pdb files to the database'

    def add_arguments(self, parser):
        parser.add_argument('--pdbs_dir', required=True)
        parser.add_argument('--entries_path', required=True)

    def handle(self, *args, **options):

        pdbs_dir = options['pdbs_dir']
        entries_path = options['entries_path']

        df = entries_df(entries_path)
        pdbs = list(tqdm(iterpdbs(pdbs_dir)))
        # 4zux 42 mer 2lo7("5my5","/data/databases/pdb/divided/my/pdb5my5.ent")

        with tqdm([("4zu4", "/data/databases/pdb/divided/zu/pdb4zu4.ent")]) as pbar:
            for (code, pdb_path) in pbar:
                pbar.set_description(code)
                try:
                    entry = df[df.IDCODE == code.upper()].iloc[0]
                except IndexError:
                    continue

                pdb_model = PDB(code=code, experiment=str(entry.EXPERIMENT))

                resolution = None
                try:
                    resolution = float(entry.RESOLUTION)
                except:
                    if resolution and not math.isnan(resolution):
                        pdb_model.resolution = resolution
                p = PDBParser(PERMISSIVE=True, QUIET=True)
                chains = list(p.get_structure(code, pdb_path)[0].get_chains())
                pdb_model.save()

                try:
                    for chain in tqdm(chains):
                        idx = 0
                        with transaction.atomic():
                            residues = []
                            for residue in chain.get_residues():
                                residue_model = Residue(pdb=pdb_model, chain=chain.id, resid=residue.id[1],
                                                        icode=residue.id[2],
                                                        type="R" if not residue.id[0].strip() else residue.id[
                                                            0].strip(),
                                                        resname=residue.resname, disordered=residue.is_disordered())
                                if is_aa(residue, standard=True):
                                    residue_model.seq_order = idx
                                    idx += 1
                                residues.append(residue_model)
                            Residue.objects.bulk_create(residues)

                        with transaction.atomic():
                            residues = {"_".join([str(x.resid), x.icode, x.resname]): x
                                        for x in Residue.objects.filter(pdb__code=code, chain=chain.id)}
                            atoms = []
                            for residue in chain.get_residues():
                                residue_model = residues["_".join([str(residue.id[1]), residue.id[2], residue.resname])]
                                for atom in residue.get_atoms():
                                    if atom.is_disordered():
                                        for altLoc, a in atom.child_dict.items():
                                            atm = Atom(residue=residue_model, serial=a.serial_number, name=a.id,
                                                       x=float(a.coord[0]), y=float(a.coord[1]),
                                                       z=float(a.coord[2]), altLoc=altLoc,
                                                       occupancy=float(a.occupancy), bfactor=float(a.bfactor),
                                                       element=a.element)
                                            atoms.append(atm)
                                    else:
                                        atm = Atom(residue=residue_model, serial=atom.serial_number, name=atom.id,
                                                   x=float(atom.coord[0]), y=float(atom.coord[1]),
                                                   z=float(atom.coord[2]), altLoc=" ",
                                                   occupancy=float(atom.occupancy), bfactor=float(atom.bfactor),
                                                   element=atom.element)
                                        atoms.append(atm)
                            Atom.objects.bulk_create(sorted(atoms,key=lambda x:x.serial))
                except Exception as ex:
                    raise CommandError(ex)
