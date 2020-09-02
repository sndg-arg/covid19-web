import os
import traceback
import warnings

from django.core.management.base import BaseCommand, CommandError

from pdbdb.models import PDB

from pdbdb.models import PDB
from pdbdb.io.PDBIO import PDBIO
from SNDG.Structure.PDBs import PDBs


class Command(BaseCommand):
    help = 'Imports a PDB'

    def add_arguments(self, parser):
        pdbs = PDBs()
        parser.add_argument('--code', required=True, help="4 letter PDB code")
        parser.add_argument('--tmp', default="data/tmp/load_pdb")
        parser.add_argument('--pdbs_dir', default="/data/databases/pdb/")
        parser.add_argument('--entries_path', default="/data/databases/pdb/entries.idx")
        parser.add_argument('--entries_url', default=pdbs.url_pdb_entries)
        parser.add_argument('--force', action="store_true")
        parser.add_argument('--verbose', default=0, choices=[0, 1], type=int)

    def handle(self, *args, **options):

        if options["verbose"] == 1:
            import logging
            logging.basicConfig(level=logging.DEBUG)

        pdbs = PDBs(options["pdbs_dir"])
        pdbs.url_pdb_entries = options["entries_url"]
        if not os.path.exists(options["entries_path"]):
            pdbs.download_pdb_entries()

        pdbio = PDBIO(options['pdbs_dir'] + "/", options['entries_path'], options['tmp'])
        pdbio.init()

        try:
            pdbs.update_pdb(options['code'])
            pdbio.process_pdb(options['code'], force=options['force'],
                              pocket_path=pdbs.pdb_pockets_path(options['code']),
                              pdb_path=pdbs.pdb_path(options['code']))

        except IOError as ex:
            traceback.print_exc()
            self.stderr.write("error processing pockets from %s: %s" % (options['code'], str(ex)))
        except Exception as ex:
            traceback.print_exc()
            raise CommandError(ex)
