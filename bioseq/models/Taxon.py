# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.db import models
from django.shortcuts import reverse


class Taxon(models.Model):
    SCIENTIFIC_NAME = "scientific name"

    taxon_id = models.AutoField(primary_key=True)
    ncbi_taxon_id = models.IntegerField(unique=True, blank=True, null=True)
    parent_taxon = models.ForeignKey("self", models.DO_NOTHING, related_name="children", null=True)
    node_rank = models.CharField(max_length=32, blank=True, null=True)
    genetic_code = models.PositiveIntegerField(blank=True, null=True)
    mito_genetic_code = models.PositiveIntegerField(blank=True, null=True)
    left_value = models.PositiveIntegerField(unique=True, blank=True, null=True)
    right_value = models.PositiveIntegerField(unique=True, blank=True, null=True)

    def scientific_name(self):
        for x in self.names.all():
            if x.name_class == Taxon.SCIENTIFIC_NAME:
                return x.name
        return str(self.ncbi_taxon_id)

    def other_names(self):
        return {x.name_class: x.name for x in self.names.filter(name_class__in=["synonym", "common name"])}

    def parents(self):
        parents = []
        curr = self.parent_taxon
        while curr.taxon_id != 1:
            parents.append(curr)
            curr = curr.parent_taxon

        return parents

    class Meta:
        managed = True
        db_table = 'taxon'

    def __str__(self):
        return "(" + str(self.ncbi_taxon_id) + ")" + self.scientific_name()


class TaxonName(models.Model):
    id = models.AutoField(primary_key=True)
    taxon = models.ForeignKey(Taxon, models.DO_NOTHING, related_name="names")
    name = models.CharField(max_length=255)
    name_class = models.CharField(max_length=32)

    class Meta:
        managed = True
        db_table = 'taxon_name'
        unique_together = (('taxon', 'name', 'name_class'),)


class TaxIdx(models.Model):
    """
    Created for indexing purposes
    """
    tax = models.OneToOneField(Taxon, models.CASCADE, primary_key=True, db_column="tax_id", related_name="keywords")
    text = models.TextField()
    genus = models.CharField(max_length=255,default="")
    family = models.CharField(max_length=255,default="")

    class Meta:
        managed = True
        db_table = 'tax_idx'
