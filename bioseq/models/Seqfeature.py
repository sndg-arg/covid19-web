# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.db import models
from django.shortcuts import reverse

from ..models.Bioentry import Bioentry
from ..models.Dbxref import Dbxref

from ..managers.SeqfeatureManager import SeqfeatureManager


class Seqfeature(models.Model):
    ENTRY_TYPES = ["CDS", "rRNA", "tRNA", "regulatory", "ncRNA", "mRNA", "repeat"]

    EXPERIMENTAL_STRUCTURE = "pdb"

    seqfeature_id = models.AutoField(primary_key=True)
    bioentry = models.ForeignKey(Bioentry, models.CASCADE, related_name="features")
    type_term = models.ForeignKey('Term', models.DO_NOTHING, related_name="features_of_type")
    source_term = models.ForeignKey('Term', models.DO_NOTHING, related_name="source_of")
    display_name = models.CharField(max_length=64, blank=True, null=True)
    rank = models.PositiveSmallIntegerField(default=1, null=True)

    index_updated = models.BooleanField(default=False)



    objects = SeqfeatureManager()

    def __str__(self):
        ls = list(self.locations.all())
        return "%s:%i-%i %s" % (self.bioentry.name, ls[0].start_pos, ls[-1].end_pos,
                                "|".join([k + ":" + v for k, v in self.qualifiers_dict().items()]))

    class Meta:
        managed = True
        db_table = 'seqfeature'
        indexes = [
            models.Index(fields=['bioentry',]),
            models.Index(fields=['type_term',]),
        ]
        # unique_together = (('bioentry', 'type_term', 'source_term', 'rank'),)

    def strand(self):
        return "+" if self.locations.all()[0].strand > 0 else "-"

    def locus_tag(self):
        qs = [x for x in self.qualifiers.all() if x.term.name == 'locus_tag']

        return qs[0].value if qs else ""

    def genes(self):
        return [x.value for x in
                self.qualifiers.all() if x.term_name in
                ["gene_symbol", "old_locus_tag", "protein_id", "Alias", "gene"]]

    def description(self):
        qs = [x for x in self.qualifiers.all() if x.term.name == 'product']

        return qs[0].value if qs else ""

    def length(self):
        return sum([abs(x.end_pos - x.start_pos) for x in self.locations])

    def subfeatures(self):
        return [x.object_seqfeature for x in self.object_relationships.all()]


    def qualifiers_dict(self):
        if not hasattr(self,"_qualifiers_dict"):
            self._qualifiers_dict = {x.term.identifier: x.value for x in self.qualifiers.all()}
        return self._qualifiers_dict

    def is_pseudo(self):
        # count = Seqfeature.objects.filter(qualifiers__term__identifier="pseudo",
        #                                                 bioentry=self.bioentry,
        #                                                 type_term__identifier="gene").count()
        # [x for x in self.bioentry.features.filter(qualifiers__term__identifier="pseudo").all()][0]

        f = [f for f in self.bioentry.features.filter(type_term__identifier="gene", qualifiers__value=self.locus_tag())]
        if f:
            f = f[0]
            return "pseudo" in f.qualifiers_dict()
        return False

    def first_location(self):
        return self.locations.all()[0]


class SeqfeatureDbxref(models.Model):
    seqfeature_dbxref_id = models.AutoField(primary_key=True)
    seqfeature = models.ForeignKey(Seqfeature, models.DO_NOTHING, related_name="dbxrefs")
    dbxref = models.ForeignKey(Dbxref, models.DO_NOTHING)
    rank = models.SmallIntegerField(default=1, null=True)

    class Meta:
        managed = True
        db_table = 'seqfeature_dbxref'
        unique_together = (('seqfeature', 'dbxref'),)


class SeqfeaturePath(models.Model):
    seqfeature_path_id = models.AutoField(primary_key=True)
    object_seqfeature = models.ForeignKey(Seqfeature, models.DO_NOTHING, related_name="object_paths")
    subject_seqfeature = models.ForeignKey(Seqfeature, models.DO_NOTHING, related_name="subject_paths")  # parent
    term = models.ForeignKey('Term', models.DO_NOTHING)
    distance = models.PositiveIntegerField(blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'seqfeature_path'
        unique_together = (('object_seqfeature', 'subject_seqfeature', 'term', 'distance'),)


class SeqfeatureQualifierValue(models.Model):
    seqfeature_qualifiervalue_id = models.AutoField(primary_key=True)
    seqfeature = models.ForeignKey(Seqfeature, models.CASCADE, related_name="qualifiers")
    term = models.ForeignKey('Term', models.DO_NOTHING)
    rank = models.SmallIntegerField(default=1, null=True)
    value = models.TextField()

    class Meta:
        managed = True
        db_table = 'seqfeature_qualifier_value'
        unique_together = (('seqfeature', 'term', 'rank'),)
        indexes = [
            models.Index(fields=['term']),
        ]

    def __str__(self):
        return str(self.term.name) + ":" + self.value


class SeqfeatureRelationship(models.Model):
    seqfeature_relationship_id = models.AutoField(primary_key=True)
    object_seqfeature = models.ForeignKey(Seqfeature, models.DO_NOTHING, related_name="object_relationships")
    subject_seqfeature = models.ForeignKey(Seqfeature, models.CASCADE, related_name="subject_relationships")  # parent
    term = models.ForeignKey('Term', models.DO_NOTHING)
    rank = models.IntegerField(default=1, null=True)

    class Meta:
        managed = True
        db_table = 'seqfeature_relationship'
        unique_together = (('object_seqfeature', 'subject_seqfeature', 'term'),)
