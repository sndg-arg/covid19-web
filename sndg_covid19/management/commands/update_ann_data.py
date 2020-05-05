import os
import subprocess as sp
from tqdm import tqdm
import requests
from datetime import datetime
from io import StringIO
import Bio.SeqIO as bpio

from django.core.management.base import BaseCommand

from bioseq.models.Biodatabase import Biodatabase
from sndg_covid19.io.CovidIO import CovidIO
from bioseq.io.BioIO import BioIO

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


class Command(BaseCommand):
    """

    """



    help = 'Tinny DB for documentation and testing purposes'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dbx_dict = {}

    def add_arguments(self, parser):
        parser.add_argument('--unip_data_url',
                            default=os.environ.get("DEFAULT_COVID_DATA_URL", Command.DEFAULT_COVID_DATA_URL))


    def handle(self, *args, **options):

        bdb = Biodatabase.objects.get(name="COVID19" + Biodatabase.PROT_POSTFIX)
        for protein in bdb.entries.prefetch_related("dbxrefs").all():
            pdbs = [x for x in protein.dbxrefs.all() if x.dbname == "PDB"]





