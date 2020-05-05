import os
from glob import glob

import environ
import pandas as pd
from biosql.io.NCBI2SQL import NCBI2SQL
# from bioresources.models import Assembly
from biosql.models import Biodatabase, Term, Bioentry, Ontology, BioentryQualifierValue, Dbxref, BioentryDbxref
from django.core.management.base import BaseCommand
from django.db import transaction
from tqdm import tqdm


# import logging
# log = logging.getLogger('django.db.backends')
# log.setLevel(logging.DEBUG)
# log.addHandler(logging.StreamHandler())

class Command(BaseCommand):
    help = 'Loads and '

    def __init__(self, stdout=None, stderr=None, no_color=False):
        super().__init__(stdout=stdout, stderr=stderr, no_color=no_color)
        self.ontology = None
        self.dbmap = {}
        self.is_a = None
        self.relmap = {}
        self.ontology_name = ""
        self.obo_path = None

    def add_arguments(self, parser):
        parser.add_argument('--accession', required=True)
        parser.add_argument('--driverDB', default="MySQLdb")
        parser.add_argument('--skipAnn', action='store_false')
        parser.add_argument('--workDir', default='/tmp')
        parser.add_argument('--jbrowse', default=None)
        parser.add_argument('--krona', default=None)

    def handle(self, *args, **options):
        #GCF_002103795.1 test
        ncbi2sql = NCBI2SQL()
        dbhost = env.db()["HOST"]
        dbname = env.db()["NAME"]
        dbpass = env.db()["PASSWORD"]
        dbuser = env.db()["USER"]
        # dbport = env.db()["PORT"]
        # dbtype = "mysql" if "mysql" in env.db()["ENGINE"] else "pg"


        options["workDir"] = os.path.abspath(options["workDir"])
        assert os.path.exists(options["workDir"]), "%s does not exists" % options["workDir"]


        # assemblyqs = Assembly.objects.filter(external_ids__identifier=options["accession"])
        # if not assemblyqs.exists():
        #     self.stderr.write(options["accession"] + " not found")
        #     return 1
        # assembly = assemblyqs.first()
        #ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gbff.gz
        acc = "_".join( options["accession"].split("_")[:2])
        if not Biodatabase.objects.filter(name=acc).exists():
            ncbi2sql.download(options["accession"], options["workDir"])
            genebank = glob(options["workDir"] + "/" + options["accession"] + "*.gbff.gz")[0]

            if "mysql" in env.db()["ENGINE"]:
                dbdriver = "MySQLdb"
            elif "post" in env.db()["ENGINE"]:
                dbdriver = "psycopg2"
            else:
                raise Exception("database not supported")

            ncbi2sql.connect_to_server(dbuser, dbpass, dbname, dbdriver, dbhost)
            ncbi2sql.create_contigs(acc,genebank)
            ncbi2sql.create_proteins(acc)

        faa_path = options["workDir"] + "/" + acc + ".faa"
        ncbi2sql.protein_fasta(acc,faa_path)


        # self.annotate(options["accession"], faa_path)
        # if options["jbrowse"]:
        #     self.jbrowse(options)
        # if options["krona"]:
        #     self.krona(options)

    def krona(self):
        pass

    def jbrowse(self):
        pass

    def annotate(self, accession, faa_path):

        ann_eggnog = faa_path + ".emapper.annotations"
        if not os.path.exists(ann_eggnog):
            "run emmaper"

        cols = "query_name,seed_eggNOG_ortholog,seed_ortholog_evalue,seed_ortholog_score,predicted_gene_name,GO_terms,KEGG_KOs,BiGG_reactions,Annotation_tax_scope,OGs,bestOG|evalue|score,COG cat,eggNOG annot"
        df_egg = pd.read_table(ann_eggnog, comment="#", names=cols.split(","))
        GO_ROOT_TERMS = ["GO:0008150", "GO:0005575", "GO:0003674"]
        go_ontology = Ontology.objects.filter(name=Ontology.GO).get()
        df_egg = df_egg.fillna(value=False)
        for _, row in tqdm(df_egg.iterrows(), total=len(df_egg)):
            with transaction.atomic():
                entry = Bioentry.objects.filter(biodatabase__name=accession + "_prots",
                                                accession=row["query_name"]).get()
                entry.description = row['eggNOG annot']

                dbxref = Dbxref.objects.get_or_create(dbname="eggnog", accession=row["seed_eggNOG_ortholog"])[0]
                BioentryDbxref.objects.get_or_create(bioentry=entry, dbxref=dbxref)

                if row["predicted_gene_name"]:
                    entry.name = row["predicted_gene_name"]
                if row["GO_terms"]:
                    for go in row["GO_terms"].split(","):
                        if go not in GO_ROOT_TERMS:
                            gotermqs = Term.objects.filter(ontology=go_ontology, identifier=go)
                            if gotermqs.exists():
                                go_term = gotermqs.get()
                                BioentryQualifierValue(bioentry=entry, term=go_term).save()
                            else:
                                self.stderr.write("%s term not found" % go)

                if row["KEGG_KOs"]:
                    for ko in row["KEGG_KOs"].split(","):
                        dbxref = Dbxref.objects.get_or_create(dbname="kegg_ko", accession=ko)[0]
                        BioentryDbxref.objects.get_or_create(bioentry=entry, dbxref=dbxref)
