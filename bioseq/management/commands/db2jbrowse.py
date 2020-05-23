import os
from django.core.management.base import BaseCommand, CommandError
from bioseq.io.DB2JBrowse import DB2JBrowse


class Command(BaseCommand):
    help = 'Loads a genome in the database'

    def add_arguments(self, parser):
        parser.add_argument('--dbname', '-n', required=True)
        parser.add_argument('--gene',   action='store_true')
        parser.add_argument('--jbrowse_path', '-j', default="data/jbrowse/")
        parser.add_argument('--tmp', default="data/tmp/")

    def handle(self, *args, **options):
        dbname = options['dbname']
        jbrowse_path = options['jbrowse_path']
        assert os.path.exists(jbrowse_path)
        io = DB2JBrowse(jbrowse_path=os.path.abspath(jbrowse_path),tmp=os.path.abspath(options["tmp"]))
        io.ovewrite = True
        if options["gene"]:
            io.excluded.append("gene")
        io.db2fs(dbname)
        self.stderr.write("genome imported!")
