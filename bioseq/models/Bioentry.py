# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from collections import defaultdict
from django.db.models import Count

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from django.utils.translation import gettext_lazy as __
from django.db import models
from django.shortcuts import reverse

from ..managers.BioentryManager import BioentryManager
from ..models.Biodatabase import Biodatabase



class Bioentry(models.Model):
    bioentry_id = models.AutoField(primary_key=True)
    biodatabase = models.ForeignKey(Biodatabase, models.CASCADE, "entries")
    taxon = models.ForeignKey('Taxon', models.DO_NOTHING, blank=True, null=True)
    name = models.CharField(max_length=40,default="")
    accession = models.CharField(max_length=128)
    identifier = models.CharField(max_length=40, blank=True, null=True,default="") #TODO should not be nullable
    division = models.CharField(max_length=6, blank=True, null=True)
    description = models.TextField(blank=True, null=True)
    version = models.PositiveSmallIntegerField(default=1, null=True)

    index_updated = models.BooleanField(default=False)

    objects = BioentryManager()

    def __str__(self):
        return self.identifier + " " + self.biodatabase.name

    class Meta:
        managed = True
        db_table = 'bioentry'
        unique_together = (('accession', 'biodatabase', 'version'), ('identifier', 'biodatabase'),)


    def get_absolute_url(self):
        return reverse(
            'bioresources:' + ('protein_view' if self.biodatabase.name.endswith("_prots") else "nucleotide_view"),
            args=[str(self.bioentry_id)])  # TODO: parametrizar la app del link

    def groupedFeatures(self):
        group = defaultdict(lambda: [])
        for f in self.features.all():
            group[f.type_term.identifier].append(f.first_location())

        return dict(group)

    def feature_counts(self):

        return {x["type_term__identifier"]: x["total"] for x in
                self.features.values('type_term__identifier').annotate(total=Count("type_term")) if
                x["type_term__identifier"] != "source"}

        # data = defaultdict(lambda: 0)
        # for f in self.features.all():
        #     data[f.type_term.name] += 1
        # return dict(data)

    def genes(self):
        # from ..models.Seqfeature import Seqfeature
        # beg = Biodatabase.objects.get(name=self.biodatabase.name.replace("_prots", ""))
        # feature = Seqfeature.objects.seqfeature_from_locus_tag(beg.biodatabase_id, self.accession)
        # feature = list(feature)[0]
        if not hasattr(self,"_genes"):
            self._genes= [x.value for x in
                          self.qualifiers.all() if x.term.name in
                          ["gene_symbol", "old_locus_tag", "protein_id", "Alias", "gene","gene_synonym"]]
        return self._genes

    def product_description(self):
        qs = self.qualifiers.filter(term__name="product")
        return qs.value if qs.exists() else None

    def molecular_weight(self):
        qs = self.qualifiers.filter(term__name="molecular_weight")
        return qs.value if qs.exists() else None

    def go_terms(self, database):
        return [x for x in self.qualifiers.all()
                if (x.term.dbxrefs.dbxref.dbname == "go")
                and (x.term.dbxrefs.dbxref.accession == database)]
        # term__dbxrefs__dbxref__accession="goslim_generic",

    def biological_process(self):
        return self.go_terms("biological_process")

    def idx_biological_process(self):
        return [x.term.identifier for x in self.biological_process()]

    def txt_biological_process(self):
        return [x.term.keywords.text for x in self.biological_process()]

    def molecular_function(self):
        return self.go_terms("molecular_function")

    def cellular_component(self):
        return self.go_terms("cellular_component")

    def ftype(self):
        return 40  # "protein"

    def qualifiers_dict(self):
        if not hasattr(self,"_qualifiers_dict"):
            self._qualifiers_dict = {x.term.identifier: x.value for x in self.qualifiers.all()}
        return self._qualifiers_dict

    def __str__(self):
        return "BioEntry('%s')" % self.accession

    def structures(self):
        #TODO fix import order
        from pdbdb.models import PDB
        from ..models.Seqfeature import Seqfeature
        codes = [x.qualifiers["structure_code"] for x in self.features.all() if
                 x.type_term.name == Seqfeature.EXPERIMENTAL_STRUCTURE]
        return PDB.objects.filter(code__in=codes)

    def __repr__(self):
        return str(self)

    def to_seq_record(self, addDesc=True):
        desc = self.description if (addDesc and self.description) else ""

        r = SeqRecord(id=self.accession.replace(" ","_"), name="", description=desc, seq=Seq(self.seq.seq))
        return r


class BioentryDbxref(models.Model):
    bioentry_dbxref_id = models.AutoField(primary_key=True)
    bioentry = models.ForeignKey(Bioentry, models.CASCADE, related_name="dbxrefs")
    dbxref = models.ForeignKey('Dbxref', models.DO_NOTHING, related_name="dbxrefs")
    rank = models.SmallIntegerField(default=1, null=True)

    class Meta:
        managed = True
        db_table = 'bioentry_dbxref'
        unique_together = (('bioentry', 'dbxref', 'rank'),)
        indexes = [
            models.Index(fields=['bioentry',]),
        ]


class BioentryPath(models.Model):
    bioentry_path_id = models.AutoField(primary_key=True)
    object_bioentry = models.ForeignKey(Bioentry, models.DO_NOTHING, related_name="object_bioentry_path")
    subject_bioentry = models.ForeignKey(Bioentry, models.DO_NOTHING, related_name="subject_bioentry_path")  # parent
    term = models.ForeignKey('Term', models.DO_NOTHING)
    distance = models.PositiveIntegerField(blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'bioentry_path'
        unique_together = (('object_bioentry', 'subject_bioentry', 'term', 'distance'),)


class BioentryQualifierValue(models.Model):
    bioentry_qualifiervalue_id = models.AutoField(primary_key=True)
    bioentry = models.ForeignKey(Bioentry, models.CASCADE, related_name="qualifiers")
    term = models.ForeignKey('Term', models.DO_NOTHING)
    value = models.TextField(blank=True, null=True)
    rank = models.IntegerField(default=1, null=True)

    class Meta:
        managed = True
        db_table = 'bioentry_qualifier_value'
        # unique_together = (('bioentry', 'term', 'rank'),)

    def __str__(self):
        return "%s - %s" % (self.value,self.term.identifier)


# class BioentryReference(models.Model):
#     bioentry_relationship_id = models.AutoField(primary_key=True)
#     bioentry = models.ForeignKey(Bioentry, models.CASCADE)
#     reference = models.ForeignKey('Reference', models.DO_NOTHING)
#     start_pos = models.IntegerField(blank=True, null=True)
#     end_pos = models.IntegerField(blank=True, null=True)
#     rank = models.SmallIntegerField(default=1, null=True)
#
#     class Meta:
#         managed = True
#         db_table = 'bioentry_reference'
#         unique_together = (('bioentry', 'reference', 'rank'),)


class BioentryRelationship(models.Model):
    bioentry_relationship_id = models.AutoField(primary_key=True)
    object_bioentry = models.ForeignKey(Bioentry, models.DO_NOTHING, related_name="object_bioentry_relationship")
    subject_bioentry = models.ForeignKey(Bioentry, models.DO_NOTHING,
                                         related_name="subject_bioentry_relationship")  # parent
    term = models.ForeignKey('Term', models.DO_NOTHING)
    rank = models.IntegerField(default=1, null=True)

    class Meta:
        managed = True
        db_table = 'bioentry_relationship'
        unique_together = (('object_bioentry', 'subject_bioentry', 'term'),)
