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

from pdbdb.io.FPocket2SQL import FPocket2SQL


def iterpdbs(pdbs_dir, pdb_extention=".ent"):
    for index_dir in os.listdir(pdbs_dir):
        if len(index_dir) == 2:
            for x in os.listdir(pdbs_dir + "/" + index_dir):
                if x.endswith(pdb_extention):
                    yield x[3:7], pdbs_dir + "/" + index_dir + "/" + x


class Command(BaseCommand):
    help = 'Process and load PDB pocket information'

    DEFAULT_PDB_PATH = "/data/databases/pdb/"

    def add_arguments(self, parser):
        parser.add_argument('--tmp', default="/tmp/fpocket")
        parser.add_argument('--pdbs_dir', default=os.environ.get("DEFAULT_PDB_PATH", Command.DEFAULT_PDB_PATH))
        parser.add_argument('--pdb', default=None, help="4 letter pdb code, process all in the db if none selected")
        parser.add_argument('--force', action="store_true",help="runs fpockets and overwrites the current result")

    def handle(self, *args, **options):

        tmp = os.path.abspath(options['tmp'])
        if not os.path.exists(tmp):
            os.makedirs(tmp)
        qs = PDB.objects.filter(code=options["pdb"]) if options["pdb"] else PDB.objects.all()
        total = qs.count()
        utils = PDBs(options["pdbs_dir"])

        with tqdm(qs, total=total) as pbar:
            for pdb in pbar:
                pbar.set_description(pdb.code)

                try:
                    fpocket2sql = FPocket2SQL()
                    fpocket2sql.create_or_get_pocket_properties()
                    fpocket2sql.load_pdb(pdb.code)
                    fpocket2sql.run_fpocket(options['tmp'],pdb_path=utils.pdb_path(pdb.code),
                                            pockets_path=utils.pdb_pockets_path(pdb.code),
                                            force=options["force"])
                    fpocket2sql.load_pockets()
                    # res.delete_dir()


                except IOError as ex:
                    traceback.print_exc()
                    self.stderr.write("error processing pockets from %s: %s" % (pdb.code, str(ex)))

                except Exception as ex:
                    traceback.print_exc()
                    raise CommandError(ex)
