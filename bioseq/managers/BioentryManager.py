# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.db import models
from django.db.models.query import QuerySet

class BioentryQuerySet(QuerySet):

    def proteins(self, index_updated=False):
        from bioseq.models.Bioentry import Bioentry
        return Bioentry.objects.prefetch_related("dbxrefs__dbxref", "qualifiers__term__dbxrefs__dbxref", "seq").filter(
            index_updated=index_updated)
            # ,bioentry_id = 3724603,
            # biodatabase_id=388)


class BioentryManager(models.Manager):

    def get_query_set(self):
        return BioentryQuerySet(self.model)

    def __getattr__(self, attr, *args):
        # see https://code.djangoproject.com/ticket/15062 for details
        if attr.startswith("_"):
            raise AttributeError
        return getattr(self.get_query_set(), attr, *args)