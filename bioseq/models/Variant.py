# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models

from bioseq.models.Bioentry import Bioentry


class Variant(models.Model):
    variant_id = models.AutoField(primary_key=True)
    bioentry = models.ForeignKey(Bioentry, models.CASCADE, "variants")
    pos = models.PositiveSmallIntegerField(default=1, null=True)
    ref = models.CharField(max_length=10)

    def __str__(self):
        return f'{self.ref}{self.pos}'

    class Meta:
        managed = True
        db_table = 'variant'
        unique_together = (('bioentry_id', 'pos', 'ref',),)


class Sample(models.Model):
    name = models.CharField(max_length=255)
    country = models.CharField(max_length=255)
    gisaid = models.CharField(max_length=255)
    date = models.DateField()

    class Meta:
        managed = True
        db_table = 'samples'
        unique_together = ('name', "gisaid")


class SampleVariant(models.Model):
    variant = models.ForeignKey(Variant, models.CASCADE, "sample_variants")
    sample = models.ForeignKey(Sample, models.CASCADE, "variants",null=True)
    alt = models.CharField(max_length=10)

    class Meta:
        managed = True
        db_table = 'sample_variant'
        unique_together = (('variant_id', 'sample_id', 'alt',),)

    def __str__(self):
        return f'{self.sample.name}{str(self.variant)}{self.alt}'
