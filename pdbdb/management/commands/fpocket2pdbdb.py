import os
import traceback
import warnings

from Bio import BiopythonWarning, BiopythonParserWarning, BiopythonDeprecationWarning, BiopythonExperimentalWarning
from django.core.management.base import BaseCommand, CommandError
from tqdm import tqdm

from pdbdb.models import PDB

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
    help = 'Loads the pdb files to the database'

    def add_arguments(self, parser):
        parser.add_argument('--tmp', default="/tmp/fpocket")

    def handle(self, *args, **options):

        tmp = os.path.abspath(options['tmp'])
        if not os.path.exists(tmp):
            os.makedirs(tmp)

        total = PDB.objects.count()


        with tqdm(PDB.objects.all(), total=total) as pbar:
            for pdb in pbar:
                pbar.set_description(pdb.code)

                try:
                    fpocket2sql = FPocket2SQL()
                    fpocket2sql.create_or_get_pocket_properties()
                    fpocket2sql.load_pdb(pdb.code)
                    fpocket2sql.run_fpocket(options['tmp'])
                    fpocket2sql.load_pockets()
                    # res.delete_dir()

                except IOError as ex:
                    traceback.print_exc()
                    self.stderr.write("error processing pockets from %s: %s" % (pdb.code, str(ex)))

                except Exception as ex:
                    traceback.print_exc()
                    raise CommandError(ex)
