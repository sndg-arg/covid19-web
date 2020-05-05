# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.db import models
from django.shortcuts import reverse

from .Ontology import Ontology
from .Dbxref import Dbxref

class Term(models.Model):
    term_id = models.AutoField(primary_key=True)
    name = models.TextField(blank=True, null=True)
    definition = models.TextField(blank=True, null=True)
    identifier = models.CharField(max_length=255, blank=True, null=True)
    is_obsolete = models.CharField(max_length=1, blank=True, null=True)
    ontology = models.ForeignKey(Ontology, models.DO_NOTHING,related_name="terms")
    version = models.PositiveSmallIntegerField(default=1, null=True)

    class Meta:
        managed = True
        db_table = 'term'
        unique_together = (('identifier', 'ontology', 'is_obsolete'),)

    def display_name(self):
        return self.name if self.name else self.identifier

    def __str__(self):
        return "%s - %s (%s)" % (self.identifier, self.name, self.ontology.name) if self.name and (
                self.identifier != self.name) else (
                "%s  (%s)" % (self.identifier,  self.ontology.name)
        )

class TermDbxref(models.Model):
    term_dbxref_id = models.AutoField(primary_key=True)
    term = models.ForeignKey(Term, models.CASCADE, related_name="dbxrefs")
    dbxref = models.ForeignKey(Dbxref, models.DO_NOTHING)
    rank = models.SmallIntegerField(default=1, null=True)

    class Meta:
        managed = True
        db_table = 'term_dbxref'
        unique_together = (('term', 'dbxref'),)


class TermPath(models.Model):
    term_path_id = models.AutoField(primary_key=True)
    subject_term = models.ForeignKey(Term, models.DO_NOTHING, related_name="subject_termpaths") # parent
    predicate_term = models.ForeignKey(Term, models.DO_NOTHING, related_name="predicate_termpaths")
    object_term = models.ForeignKey(Term, models.DO_NOTHING, related_name="object_termpaths")
    ontology = models.ForeignKey(Ontology, models.DO_NOTHING)
    distance = models.PositiveIntegerField(blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'term_path'
        unique_together = (('subject_term', 'predicate_term', 'object_term', 'ontology', 'distance'),)


class TermRelationship(models.Model):
    term_relationship_id = models.AutoField(primary_key=True)
    subject_term = models.ForeignKey(Term, models.DO_NOTHING, related_name="subject_termrelationships")  # parent
    predicate_term = models.ForeignKey(Term, models.DO_NOTHING, related_name="predicate_termrelationships")
    object_term = models.ForeignKey(Term, models.DO_NOTHING, related_name="object_termrelationships")  # child
    ontology = models.ForeignKey(Ontology, models.DO_NOTHING)

    class Meta:
        managed = True
        db_table = 'term_relationship'
        unique_together = (('subject_term', 'predicate_term', 'object_term', 'ontology'),)


class TermRelationshipTerm(models.Model):
    term_relationship = models.OneToOneField(TermRelationship, models.DO_NOTHING, primary_key=True)
    term = models.OneToOneField(Term, models.DO_NOTHING, unique=True)

    class Meta:
        managed = True
        db_table = 'term_relationship_term'


class TermSynonym(models.Model):
    term_synonym_id = models.AutoField(primary_key=True)
    synonym = models.CharField(max_length=255)
    term = models.ForeignKey(Term, models.CASCADE, related_name="synonyms")

    class Meta:
        managed = True
        db_table = 'term_synonym'
        unique_together = (('term', 'synonym'),)

class TermIdx(models.Model):
    """
    Created for indexing purposes
    """
    term = models.OneToOneField(Term, models.CASCADE, primary_key=True, db_column="term_id", related_name="keywords")
    text = models.TextField()

    class Meta:
        managed = True
        db_table = 'term_idx'
