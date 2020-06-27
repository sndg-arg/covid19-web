import json
from json import JSONDecodeError

import Bio.SeqIO as bpio
from django.core.management.base import BaseCommand
from tqdm import tqdm

from bioseq.models.Bioentry import Bioentry


class Command(BaseCommand):
    help = 'Loads a genome in the database'

    def add_arguments(self, parser):
        parser.add_argument('--dbname', '-n', required=False)
        parser.add_argument('--query', '-q', required=True)
        parser.add_argument('--format', '-f', default="fasta", choices=['fasta'], )  # , 'gb', 'gff'

    def handle(self, *args, **options):
        dbname = options['dbname']
        query = options['query']
        out_format = options['format']

        try:
            query = json.loads(query)
        except JSONDecodeError:
            query = {"accession": query}
        if dbname:
            query["biodatabase__name"] = dbname

        self.stderr.write(f"quering... {json.dumps(query)}")

        qs = Bioentry.objects.prefetch_related("seq").filter(**query)

        self.stderr.write(f"retreived sequences: {qs.count()}")

        for be in tqdm(qs, file=self.stderr, total=qs.count()):
            r = be.to_seq_record()
            bpio.write(r, self.stdout, out_format)

        self.stderr.write("genome imported!")
