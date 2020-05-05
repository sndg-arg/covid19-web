from django.contrib import admin
# from easy_select2 import select2_modelform


from .models.Biodatabase import Biodatabase
from .models.Bioentry import Bioentry
from .models.Biosequence import Biosequence
from .models.Ontology import Ontology
from .models.Term import Term
from .models.Bioentry import BioentryQualifierValue
from .models.Seqfeature import Seqfeature
from .models.Dbxref import Dbxref
from .models.Seqfeature import SeqfeatureDbxref
from .models.Taxon import Taxon,TaxonName
from .models.Location import Location

admin.site.register(Biosequence)
admin.site.register(Biodatabase)
admin.site.register(Ontology)

# @admin.register(Tool)
# class ToolAdmin(admin.ModelAdmin):
#     list_display = ["name","version","url"]
#     search_fields = ["name"]
#
# admin.site.register(ToolRun)

# BioentryForm = select2_modelform(Bioentry, attrs={'width': '250px'})

@admin.register(Bioentry)
class BiosequenceAdmin(admin.ModelAdmin):
    search_fields = ["name","accession","identifier"]
    # form = BioentryForm
    # fields = ('name','biodatabase','accession','identifier','division','version')
    raw_id_fields = (
        'taxon',
    )
    # readonly_fields = ('taxon_id',)

@admin.register(Term)
class BiosequenceAdmin(admin.ModelAdmin):
    search_fields = ["name","identifier"]


@admin.register(BioentryQualifierValue)
class BioentryQualifierValueAdmin(admin.ModelAdmin):
    raw_id_fields = (
        'bioentry', "term"
    )


@admin.register(Seqfeature)
class SeqfeatureAdmin(admin.ModelAdmin):
    search_fields = ["type_term__name","qualifiers__value"]
    raw_id_fields = (
        'bioentry', "type_term", "source_term"
    )

    def get_search_results(self,request, qs, term,**kwargs):
        return super().get_search_results(request, qs, term,**kwargs)

@admin.register(Location)
class LocationAdmin(admin.ModelAdmin):
    raw_id_fields = (
        'seqfeature', "dbxref", "term"
    )


@admin.register(Dbxref)
class DbxrefAdmin(admin.ModelAdmin):
    search_fields = ["accession"]

@admin.register(SeqfeatureDbxref)
class SeqfeatureDbxrefAdmin(admin.ModelAdmin):
    autocomplete_fields = ["dbxref"]
    raw_id_fields = (
        'seqfeature',
    )

@admin.register(Taxon)
class TaxonAdmin(admin.ModelAdmin):
    search_fields = ["ncbi_taxon_id","node_rank"]
    list_display = ('ncbi_taxon_id','scientific_name','node_rank','genetic_code','parent_taxon')
    raw_id_fields = (
        'parent_taxon',
    )

    def get_queryset(self, request):
        return super(TaxonAdmin, self).get_queryset(request).prefetch_related("names","parent_taxon__names")

@admin.register(TaxonName)
class TaxonNameAdmin(admin.ModelAdmin):

    list_display = ('name','name_class','parent_id')
    search_fields = ["taxon__ncbi_taxon_id","name"]
    raw_id_fields = (
        'taxon',
    )
    def parent_id(self, obj):
        return obj.taxon.ncbi_taxon_id

    def get_queryset(self, request):
        return super(TaxonNameAdmin, self).get_queryset(request).select_related("taxon")