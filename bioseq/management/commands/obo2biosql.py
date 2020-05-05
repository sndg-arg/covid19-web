from bioseq.models.Ontology import Ontology
from bioseq.models.Term import Term, TermRelationship, TermDbxref, TermSynonym

from django.core.management.base import BaseCommand
from django.db import transaction
from goatools.obo_parser import GODag
from tqdm import tqdm


class Command(BaseCommand):
    help = 'Loads the obo files to the database. Prepared for GO and SO'

    def __init__(self, stdout=None, stderr=None, no_color=False):
        super().__init__(stdout=stdout, stderr=stderr, no_color=no_color)
        self.ontology = None
        self.ontology_name = ""
        self.obo_path = None
        self.is_a = None

    def add_arguments(self, parser):
        parser.add_argument('--obo_path', required=True)
        parser.add_argument('--relationships_obo_path', required=False)
        parser.add_argument('--ontology', required=True, choices=["go", "so"])



    def create_terms(self, shortname):
        self.cache = {}
        ids = []
        with open(self.obo_path) as h:  # filter the ids because GODag iterates over alt_ids too
            for l in h.readlines():
                if l.startswith("id: " + shortname.upper() + ":"):
                    ids.append(l.split(" ")[1].strip())
        ids = set(ids)
        if shortname == "go":
            attrs = ["def", "synonym", "subset", "alt_id", "dbxref"]
        else:
            attrs = ["def", "subset", "dbxref", "alt_id"]
        go_dag = GODag(self.obo_path, load_obsolete=True,
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
                        if not Term.objects.filter(ontology=self.ontology,identifier=go).exists():
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
                            self.cache[go] = Term.objects.filter(ontology=self.ontology,identifier=go).get()
                            print("repeated: " + go)
                    except StopIteration:
                        finished = True



    def create_relationships(self,shortname):
        go_dag = GODag(self.relationships_obo_path, optional_attrs=['relationship'], load_obsolete=(shortname == "so"))

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

    def handle(self, *args, **options):
        if "relationships_obo_path" not in options or not options["relationships_obo_path"]:
            options["relationships_obo_path"] = options["obo_path"]
        self.relationships_obo_path = options["relationships_obo_path"]
        if options["ontology"] == "go":
            self.ontology_name = Ontology.GO
        elif options["ontology"] == "so":
            self.ontology_name = Ontology.SO

        self.obo_path = options["obo_path"]
        self.ontology = Ontology.objects.get(name=self.ontology_name)
        Ontology.load_go_base()
        self.is_a = Ontology.relmap["is_a"]

        self.create_terms(options["ontology"])
        self.create_relationships(options["ontology"])

