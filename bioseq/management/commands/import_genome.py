import os
import warnings
from tqdm import tqdm
import subprocess as sp
import gzip
from io import StringIO

import Bio.SeqIO as bpio

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from bioseq.io.BioIO import BioIO
from Bio import BiopythonWarning, BiopythonParserWarning, BiopythonDeprecationWarning, BiopythonExperimentalWarning
warnings.simplefilter('ignore', RuntimeWarning)
warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', BiopythonParserWarning)
warnings.simplefilter('ignore', BiopythonDeprecationWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)


class Command(BaseCommand):
    help = 'Loads a genome in the database'

    def add_arguments(self, parser):
        parser.add_argument('--input', '-i', required=True)
        parser.add_argument('--accession', '-a', required=True)
        parser.add_argument('--name', '-n', required=True)
        parser.add_argument('--taxon', '-t', type=int, required=True)

    def handle(self, *args, **options):
        input_file = options['input']
        accession = options['accession']
        name = options['name']
        taxon = options['taxon']
        assert os.path.exists(input_file),"'%s' does not exists" % input_file
        io = BioIO(accession, taxon)

        grep_cmd = 'grep -c "FEATURES *Location/Qualifiers" "%s"' % input_file
        if input_file.endswith(".gz"):
            grep_cmd = 'z' + grep_cmd
            input_file = gzip.open(input_file, "rt")
        total = int(sp.check_output(grep_cmd, shell=True))
        io.create_db()
        io.process_record_list(bpio.parse(input_file, "gb"), total)

        self.stderr.write("genome imported!")
