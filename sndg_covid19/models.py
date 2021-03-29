# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from collections import defaultdict
from django.db import models
from django.shortcuts import reverse
import zipfile


from sndg_covid19.users.models import User
from django.core.files.storage import FileSystemStorage
from django.conf import settings


upload_storage = FileSystemStorage(location=str(settings.UPLOAD_ROOT), base_url='/uploads')

class ImportJob(models.Model):
    import_job_id = models.AutoField(primary_key=True)
    user = models.ForeignKey(User, models.SET_NULL, "import_jobs", null=True)
    name = models.CharField(max_length=255,blank=False,null=False)
    aln_type = models.CharField(max_length=20, blank=False, null=False, choices=[("genome","genome"),("spike","spike")])
    version = models.PositiveSmallIntegerField(default=1, null=True)
    fasta = models.FileField(upload_to='uploads/',storage=upload_storage,help_text="archivo de secuencias en formato .fasta.zip")
    csv = models.FileField(upload_to='uploads/',storage=upload_storage,help_text="archivo de propiedades .csv")
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)


    status = models.CharField(max_length=255,default="iniciado",blank=False,null=False)
    status_desc = models.TextField(blank=False,null=True)
    debug_status_desc = models.TextField(blank=False,null=True)
    errors = models.TextField(blank=False,null=True)

    # objects = BioentryManager()

    def __str__(self):
        return self.name

    class Meta:
        managed = True
        db_table = 'importjobs'
        unique_together = (('name', 'version'),)
