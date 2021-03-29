from django.contrib import admin
from django.urls import reverse
from django.utils.html import format_html

from sndg_covid19.models import ImportJob

@admin.register(ImportJob)
class ImportJobAdmin(admin.ModelAdmin):
    pass
