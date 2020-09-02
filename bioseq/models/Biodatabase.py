# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.db import models
from django.shortcuts import reverse
from ..managers.BiodatabaseManager import BiodatabaseManager

class Biodatabase(models.Model):

    PROT_POSTFIX = "_prots"
    RNAS_POSTFIX = "_rnas"


    biodatabase_id = models.AutoField(primary_key=True)
    name = models.CharField(unique=True, max_length=128)
    authority = models.CharField(max_length=128, blank=True, null=True)
    description = models.TextField(blank=True, null=True)

    objects = BiodatabaseManager()

    class Meta:
        managed = True
        db_table = 'biodatabase'

    def get_absolute_url(self):
        return reverse('bioseq:assembly_view', args=[str(self.biodatabase_id)])

