# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.db import models
from django.db.models.query import QuerySet

class BiodatabaseQuerySet(QuerySet):

    def dbs(self,db_name):
        from bioseq.models.Biodatabase import Biodatabase
        return Biodatabase.objects.filter(name__startswith=db_name)



class BiodatabaseManager(models.Manager):

    def get_query_set(self):
        return BiodatabaseQuerySet(self.model)

    def __getattr__(self, attr, *args):
        # see https://code.djangoproject.com/ticket/15062 for details
        if attr.startswith("_"):
            raise AttributeError
        return getattr(self.get_query_set(), attr, *args)