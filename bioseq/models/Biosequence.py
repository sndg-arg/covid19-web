# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.db import models
from django.shortcuts import reverse

from ..models.Bioentry import Bioentry

class Biosequence(models.Model):
    bioentry = models.OneToOneField(Bioentry, models.CASCADE, primary_key=True, related_name="seq")
    version = models.SmallIntegerField(blank=True, default=1, null=True)
    length = models.IntegerField(blank=True, null=True)
    alphabet = models.CharField(max_length=10, blank=True, null=True)
    seq = models.TextField(blank=True, null=True)



    class Meta:
        managed = True
        db_table = 'biosequence'