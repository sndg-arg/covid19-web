from django.contrib import admin
from django.urls import reverse
from django.utils.html import format_html

from .models import Property, ResidueSet, PDB, Residue, PDBResidueSet


@admin.register(Property)
class PropertyAdmin(admin.ModelAdmin):
    list_display = ["name","description"]
    search_fields = ["name","description"]

@admin.register(ResidueSet)
class ResidueSetAdmin(admin.ModelAdmin):
    list_display = ["name","description"]
    search_fields = ["name","description"]


class PDBResidueSetInline(admin.TabularInline):

    model = PDBResidueSet
    can_delete = False
    # def has_add_permission(self, request, obj=None):
    #     return False
    #
    # def has_delete_permission(self, request, obj=None):
    #     return False

@admin.register(PDB)
class PDBAdmin(admin.ModelAdmin):
    list_display = ["code","resolution","taxon","residues","ligands","solvent"]
    search_fields = ["code","experiment"]
    raw_id_fields = (
        'taxon',
    )
    def residues(self, obj):
        #https://en.proft.me/2014/10/12/reversing-admin-urls-django/
        return format_html('<a href="{url}?pdb_id={{pdb_id}}">Residues</a>',
                           pdb_id=obj.id,url=reverse('admin:pdbdb_residue_changelist'))

    residues.short_description = "Residues"

    def ligands(self, obj):
        return format_html("<a href='{url}'>{url}</a>", url=obj.code)

    ligands.short_description = "Ligands"

    def solvent(self, obj):
        return format_html("<a href='{url}'>{url}</a>", url=obj.code)

    solvent.short_description = "Solvents"



    # autocomplete_fields = ["dbxref"]
    # inlines = [
    #     PDBResidueSetInline,
    # ]

@admin.register(Residue)
class ResidueAdmin(admin.ModelAdmin):
    # list_display = ["code","resolution","tax"]
    # search_fields = ["code","experiment"]
    raw_id_fields = ("pdb",)
    # autocomplete_fields = ["dbxref"]
    list_display = ["chain","resid","icode","resname"]



    def get_queryset(self, request):

        #pdb = PDB.objects.get(id=int(request.GET["pdbid"]))
        qs = super(ResidueAdmin, self).get_queryset(request)
        #qs.filter(pdb=pdb)
        return qs
