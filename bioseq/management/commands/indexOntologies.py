import sys

from django.db import transaction
from django.core.management.base import BaseCommand, CommandError

from bioseq.models.Ontology import Ontology
from bioseq.models.Term import TermRelationship, Term, TermIdx
from bioseq.models.Taxon import Taxon, TaxonName, TaxIdx
from tqdm import tqdm


class Command(BaseCommand):
    help = 'Loads the obo files to the database. Prepared for GO and SO'

    def __init__(self, stdout=None, stderr=None, no_color=False):
        super().__init__(stdout=stdout, stderr=stderr, no_color=no_color)
        self.is_a = Term.objects.get(identifier="is_a")
        self.processed_tax = {1: 1}

    def travel_tax_child(self, parent, iter_adv, txt_acc, buffer, genus, family):

        txt = [x.name for x in parent.names.all()] + [str(parent.ncbi_taxon_id)]
        txt += txt_acc

        genus2 = parent.scientific_name() if parent.node_rank == "genus" else ""
        family2 = parent.scientific_name() if parent.node_rank == "family" else ""

        if parent.ncbi_taxon_id != 1:
            if parent.ncbi_taxon_id not in self.processed_tax:
                self.processed_tax[parent.ncbi_taxon_id] = 1
                buffer["data"][parent.ncbi_taxon_id] = TaxIdx(tax=parent, text=" ".join(txt), genus=genus,
                                                              family=family)

        for t in Taxon.objects.prefetch_related("names", "children").filter(parent_taxon=parent):
            # if t.ncbi_taxon_id in [131567, 2, 1783257, 203682, 203683, 112, 1763524,127]:
            if t.ncbi_taxon_id not in self.processed_tax:
                self.travel_tax_child(t, iter_adv, txt, buffer, genus2, family2)

        if len(buffer["data"]) > 5000:
            TaxIdx.objects.bulk_create(buffer["data"].values())
            buffer["data"] = {}

        iter_adv.update(1)

    # def travel_go_child(self, cache, parent, iter_adv, text_acc, buffer):
    #     if parent.term_id in cache:
    #         return
    #     txt = " ".join([parent.name, parent.identifier, parent.definition] +
    #                    [x.term.name for x in parent.dbxrefs.all()])  # [x.synonym for x in parent.synonyms.all()]
    #     txt += text_acc
    #     for c in TermRelationship.objects.prefetch_related("object_term__dbxrefs__term",
    #                                                        ).filter(subject_term=parent,
    #                                                                 predicate_term=self.is_a):
    #         if c.object_term.identifier not in cache:
    #             self.travel_go_child(cache, c.object_term, iter_adv, txt, buffer)
    #
    #     if parent.term_id in cache:
    #         return
    #     cache[parent.term_id] = 1
    #     buffer["data"][parent.ncbi_taxon_id] = TermIdx(term=parent, text=txt)
    #     iter_adv.update(1)
    #     if len(buffer["data"]) == 5000:
    #         with transaction.atomic():
    #             TermIdx.objects.bulk_create(buffer["data"].values())
    #         buffer["data"] = {}

    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):
        TaxIdx.objects.all().delete()
        root = Taxon.objects.prefetch_related("names").get(ncbi_taxon_id=1)
        count = Taxon.objects.count()
        iter_adv = tqdm(total=count, file=sys.stderr)

        buffer = {"data": {}}

        for c in tqdm(root.children.all(), total=root.children.count()):
            # if c.ncbi_taxon_id in [131567, 2, 1783257, 203682, 203683, 112, 1763524,127]:
            self.travel_tax_child(c, iter_adv, "", buffer,"","")
        TaxIdx.objects.bulk_create(buffer["data"].values())

        # buffer = {"data": []}
        # cache = {x.term_id: 1 for x in self.tqdm(TermIdx.objects.select_related("term").all())}
        # go = Ontology.objects.get(name=Ontology.GO)
        # count = Term.objects.filter(ontology=go).count()
        #
        # # Term.objects.prefetch_related("synonyms").filter(ontology=go).first() "synonyms",
        # iter_adv = self.tqdm(total=count)
        # iter_adv.update(len(cache))
        # gos = list(Term.objects.filter(ontology=go))
        # for c in gos:
        #     iter_adv.set_description(c.identifier)
        #     self.travel_go_child(cache, c, iter_adv, "", buffer)
