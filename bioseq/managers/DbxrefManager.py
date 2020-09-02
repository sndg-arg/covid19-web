# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.db import models
from django.db.models.query import QuerySet

class DbxrefQuerySet(QuerySet):

    def latest_version(self, dbname,accession):
        from ..models.Dbxref import Dbxref
        return Dbxref.objects.filter(dbname=dbname,accession=accession).order_by("-version")[0]




class DbxrefManager(models.Manager):

    def get_query_set(self):
        return DbxrefQuerySet(self.model)

    def __getattr__(self, attr, *args):
        # see https://code.djangoproject.com/ticket/15062 for details
        if attr.startswith("_"):
            raise AttributeError
        return getattr(self.get_query_set(), attr, *args)
