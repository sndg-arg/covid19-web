import os
from bioseq.models.Ontology import Ontology
from bioseq.models.Term import Term, TermRelationship, TermDbxref, TermSynonym
from bioseq.models.Dbxref import DBx

from django.core.management.base import BaseCommand
from django.db import transaction
from goatools.obo_parser import GODag
from tqdm import tqdm

from SNDG.WebServices import download_file
import xmltodict


class Command(BaseCommand):
    DEFAULT_GO_URL = "http://current.geneontology.org/ontology/go.obo"
    DEFAULT_GO_REL_URL = "http://current.geneontology.org/ontology/go-basic.obo"
    DEFAULT_TAX_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
    DEFAULT_CROSS_REF_DBS = "https://www.uniprot.org/database/?query=*&format=rdf&force=true"

    help = 'Loads initial ontology and terms: GO and taxonomy'

    def __init__(self, stdout=None, stderr=None, no_color=False):
        super().__init__(stdout=stdout, stderr=stderr, no_color=no_color)
        self.ontology = None
        self.ontology_name = ""
        self.obo_path = None
        self.is_a = None

    def add_arguments(self, parser):
        parser.add_argument('--obo_path', default="data/tmp/go.obo")
        parser.add_argument('--relationships_obo_path', default="data/tmp/go-basic.obo")
        parser.add_argument('--go_url', default=os.environ.get("DEFAULT_GO_URL", Command.DEFAULT_GO_URL))
        parser.add_argument('--go_basic_url', default=os.environ.get("DEFAULT_GO_REL_URL", Command.DEFAULT_GO_REL_URL))
        parser.add_argument('--tax_url', default=os.environ.get("DEFAULT_TAX_URL", Command.DEFAULT_TAX_URL))

        parser.add_argument('-go', action='store_false', )
        parser.add_argument('-tax', action='store_false', )

    def handle(self, *args, **options):
        Ontology.load_ann_terms()
        Ontology.load_go_base()

        if options["go"]:
            if not os.path.exists(options["obo_path"]):
                download_file(options["go_url"], options["obo_path"])
            if not os.path.exists(options["relationships_obo_path"]):
                download_file(options["go_basic_url"], options["relationships_obo_path"])

            self.ontology = Ontology.objects.get(name=Ontology.GO)
            self.is_a = Ontology.relmap["is_a"]

            self.create_terms(options["obo_path"], "go")
            self.create_relationships(options["relationships_obo_path"], "go")

        if options["tax"]:
            pass

        self.stderr.write("Finished!")

    def cross_reference_dbs(self):
        download_file(Command.DEFAULT_CROSS_REF_DBS, "data/tmp/database-all.rdf")
        with open("data/tmp/database-all.rdf") as h:
            data = xmltodict.parse(h.read(), "utf-8")
            for db in data["rdf:RDF"]["rdf:Description"]:
                DBx.objects.get_or_create(
                    url=db['@rdf:about'],
                    name=db['abbreviation'] if "abbreviation" in db else db['dcterms:identifier'],
                    category=db.get('category', ""),
                    description=db.get('rdfs:label', ""),
                    url_template=db.get('urlTemplate', db['dcterms:identifier'])
                )
        DBx.objects.get_or_create(
            url="www.uniprot.org",
            name="UnipAcc",
            category='Protein annotation databases',
            description="UNIPROT",
            url_template="https://www.uniprot.org/uniprot/%s",
        )

    def create_terms(self, obo_path, shortname, attrs=["def", "synonym", "subset", "alt_id", "dbxref"]):
        self.cache = {}
        ids = []
        with open(obo_path) as h:  # filter the ids because GODag iterates over alt_ids too
            for l in h.readlines():
                if l.startswith("id: " + shortname.upper() + ":"):
                    ids.append(l.split(" ")[1].strip())
        ids = set(ids)

        # PAra SO attrs = ["def", "subset", "dbxref", "alt_id"]

        go_dag = GODag(obo_path, load_obsolete=True,
                       optional_attrs=attrs)

        finished = False
        pbar = iter(tqdm(ids))
        while not finished:
            with transaction.atomic():
                for _ in range(2000):
                    try:
                        go = next(pbar)
                        if go not in go_dag:
                            continue
                        term = go_dag[go]
                        if not Term.objects.filter(ontology=self.ontology, identifier=go).exists():
                            dbTerm = Term(name=term.name,
                                          definition=term.defn if hasattr(term, "defn") else "",
                                          identifier=go,
                                          is_obsolete="T" if term.is_obsolete else "F",
                                          ontology=self.ontology)
                            dbTerm.save()
                            if term.namespace:
                                termdbref = TermDbxref(term=dbTerm, dbxref=Ontology.dbmap[term.namespace], rank=1)
                                termdbref.save()

                            for subset in term.subset:
                                if subset in Ontology.dbmap:
                                    termdbref = TermDbxref(term=dbTerm, dbxref=Ontology.dbmap[subset], rank=1)
                                    termdbref.save()
                            if hasattr(term, "synonym"):
                                for synonym in term.synonym:
                                    TermSynonym.objects.get_or_create(term=dbTerm, synonym=synonym[0][:255])

                            for synonym in term.alt_ids:
                                TermSynonym.objects.get_or_create(term=dbTerm, synonym=synonym[0][:255])

                            self.cache[go] = dbTerm
                        else:
                            self.cache[go] = Term.objects.filter(ontology=self.ontology, identifier=go).get()
                            print("repeated: " + go)
                    except StopIteration:
                        finished = True

    def create_relationships(self, relationships_obo_path, load_obsolete=False):
        """

        :param relationships_obo_path:
        :param shortname:
        :param load_obsolete:   --> (shortname == "so")
        :return:
        """
        go_dag = GODag(relationships_obo_path, optional_attrs=['relationship'], load_obsolete=load_obsolete)

        terms = Term.objects.filter(ontology=self.ontology)

        finished = False
        pbar = iter(tqdm(terms.all(), total=terms.count()))
        while not finished:
            with transaction.atomic():
                for _ in range(2000):
                    try:
                        dbTerm = next(pbar)
                        go = dbTerm.identifier
                        if dbTerm.identifier not in go_dag:
                            continue
                        term = go_dag[dbTerm.identifier]

                        for child in term.children:
                            if go in self.cache and child.id in self.cache:
                                r = TermRelationship(
                                    subject_term=self.cache[go],  # parent
                                    predicate_term=self.is_a,
                                    object_term=self.cache[child.id],  # child
                                    ontology=self.ontology
                                )
                                r.save()

                        for rel, terms in term.relationship.items():
                            for child in terms:
                                if go in self.cache and child.id in self.cache:
                                    r = TermRelationship(
                                        subject_term=self.cache[go],  # parent
                                        predicate_term=Ontology.relmap[rel],
                                        object_term=self.cache[child.id],  # child
                                        ontology=self.ontology
                                    )
                                    r.save()
                    except StopIteration:
                        finished = True
