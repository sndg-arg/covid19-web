import os
import traceback
import warnings

from django.core.management.base import BaseCommand, CommandError

from pdbdb.models import PDB


from pdbdb.models import PDB
from pdbdb.io.PDBIO import PDBIO
from SNDG.Structure.PDBs import  PDBs

class Command(BaseCommand):
    help = 'Imports a PDB'

    def add_arguments(self, parser):
        pdbs = PDBs()
        parser.add_argument('--code', required=True,help="4 letter PDB code")
        parser.add_argument('--tmp', default="data/tmp/load_pdb")
        parser.add_argument('--pdbs_dir', default="/data/databases/pdb/divided/")
        parser.add_argument('--entries_path', default="/data/databases/pdb/entries.idx")
        parser.add_argument('--entries_url', default=pdbs.url_pdb_entries)


    def handle(self, *args, **options):

        pdbs = PDBs()
        pdbs.url_pdb_entries = options["entries_url"]
        if not os.path.exists(options["entries_path"]):
            pdbs.download_pdb_entries()

        pdbio = PDBIO(options['pdbs_dir'] + "/",options['entries_path'],options['tmp'])
        pdbio.init()

        try:
            pdbio.process_pdb(options['code'])

        except IOError as ex:
            traceback.print_exc()
            self.stderr.write("error processing pockets from %s: %s" % (options['code'], str(ex)))
        except Exception as ex:
            traceback.print_exc()
            raise CommandError(ex)
