import os

from django.core.management.base import BaseCommand, CommandError
from tqdm import tqdm

from pdbdb.io.PDB2SQL import PDB2SQL
from pdbdb.models import PDB
from SNDG.Structure.PDBs import  PDBs


def iterpdbs(pdbs_dir, pdb_extention=".ent"):
    for index_dir in os.listdir(pdbs_dir):
        if len(index_dir) == 2:
            for x in os.listdir(pdbs_dir + "/" + index_dir):
                if x.endswith(pdb_extention):
                    yield x[3:7], pdbs_dir + "/" + index_dir + "/" + x


class Command(BaseCommand):
    help = 'Loads the pdb files to the database'

    def add_arguments(self, parser):
        pdbs = PDBs()
        parser.add_argument('--pdbs_dir', default="data/pdb/")
        parser.add_argument('--entries_path', default="data/pdb/entries.idx")
        parser.add_argument('--only_annotated', action='store_false')
        parser.add_argument('--entries_url', default=pdbs.url_pdb_entries)

    def handle(self, *args, **options):
        pdbs = PDBs(pdb_dir=options['pdbs_dir'])
        pdbs.url_pdb_entries = options["entries_url"]
        if not os.path.exists(options["entries_path"]):
            pdbs.download_pdb_entries()


        pdb2sql = PDB2SQL(options['pdbs_dir'], options['entries_path'])
        pdb2sql.load_entries()
        if options["only_annotated"]:
            self.stderr.write("only_annotated option activated by default")
            from bioseq.models.Dbxref import Dbxref
            pdbs = [(x.accession.lower(),pdbs.pdb_path( x.accession.lower()))
                    for x in Dbxref.objects.filter(dbname="PDB")]
        else:
            pdbs = list(tqdm(iterpdbs(options['pdbs_dir'])))
        # 4zux 42 mer 2lo7("5my5","/data/databases/pdb/divided/my/pdb5my5.ent")
        # ("4zu4", "/data/databases/pdb/divided/zu/pdb4zu4.ent")

        with tqdm(pdbs) as pbar:
            for code,pdb_path in pbar:
                code = code.lower()
                try:
                    pdb_path = pdb2sql.download(code)
                except:
                    self.stderr.write("PDB %s could not be downloaded" % code)
                    continue

                if PDB.objects.filter(code=code).exists():
                    self.stderr.write("PDB %s already exists" % code)
                    continue

                pbar.set_description(code)
                try:
                    pdb2sql.create_pdb_entry(code, pdb_path)
                    pdb2sql.update_entry_data(code, pdb_path)
                except Exception as ex:
                    raise CommandError(ex)
