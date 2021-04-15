from django.contrib import admin
from django.urls import reverse
from django.utils.html import format_html

from sndg_covid19.models import ImportJob


@admin.register(ImportJob)
class ImportJobAdmin(admin.ModelAdmin):
    list_display = (
    'import_job_id', 'user', 'name', 'aln_type', 'version', "status", "status_desc", "fasta", "csv", "created_at")
