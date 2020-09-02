# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.db import models
from django.shortcuts import reverse

from ..models.Dbxref import Dbxref

class Location(models.Model):
    location_id = models.AutoField(primary_key=True)
    seqfeature = models.ForeignKey('Seqfeature', models.CASCADE, related_name="locations")
    dbxref = models.ForeignKey(Dbxref, models.DO_NOTHING, blank=True, null=True)
    term = models.ForeignKey('Term', models.DO_NOTHING, blank=True, null=True)
    start_pos = models.IntegerField(blank=True, null=True)
    end_pos = models.IntegerField(blank=True, null=True)
    strand = models.IntegerField(default=1)
    rank = models.SmallIntegerField(default=1, null=True)

    class Meta:
        managed = True
        db_table = 'location'
        unique_together = (('seqfeature', 'rank'),)

    def __str__(self):
        return "%i-%i(%s)" % (self.start_pos, self.end_pos, "+" if self.strand == 1 else "-")


class LocationQualifierValue(models.Model):
    location = models.OneToOneField(Location, models.DO_NOTHING, primary_key=True)
    term = models.ForeignKey('Term', models.DO_NOTHING)
    value = models.CharField(max_length=255)
    int_value = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'location_qualifier_value'
        unique_together = (('location', 'term'),)
