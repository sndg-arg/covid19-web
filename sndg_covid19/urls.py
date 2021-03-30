from django.urls import path

from .views.AssemblyView import assembly_view
from .views.ProteinView import ProteinView,MSAView
from .views.StructureView import StructureView,pdb_download,StructureStaticView
from .views.ProjectStrainsView import ProjectStrainsView
from .views.VariantView import VariantView,pdb_variants_download,InmunovaView,variant_data
from .views.BackOffice.AlnImportView import AlnImportView
from .views.SamplesFromVariantView import SamplesFromVariantView


from django.views.generic import RedirectView

from django.conf import settings

app_name = 'covid'
# if settings.MINCYT_URL:

urlpatterns = [


    path('', assembly_view, name='genome_view'),
    path('<int:pk>', ProteinView, name='gene_view'),
    path('<int:pk>/msa', MSAView, name='msa_view'),
    path('structure/<int:prot_id>/<str:pdbid>', StructureView.as_view(), name='structure_view'),
    path('structure_static/<int:prot_id>/<str:pdbid>', StructureStaticView.as_view(), name='structure_static_view'),
    path('structure/download/<str:pdbid>', pdb_download, name='download_pdb_data'),

    path('project/', ProjectStrainsView.as_view(), name='project_strains'),
    path('variant/<str:gene>/<int:pos>', VariantView.as_view(), name='pos_stats'),
    path('variant/<str:gene>/<int:pos>/sample', SamplesFromVariantView.as_view(), name='variant_samples'),
    path('variant_data/<str:gene>/<int:pos>', variant_data, name='pos_stats_data'),


    # TODO: link en el menu para S
    # Ver muestras original
    # Pagina de combinaciones

    path('variant/', pdb_variants_download, name='download_pdb_variants'),
    path('variant/spike', InmunovaView.as_view(), name='spike_variants'),
    path('backoffice/aln_import/', AlnImportView.as_view(), name='aln_import'),



    # path('static/jbrowse/data/COVID19/trackList2.json', assembly_view, name='genome_view'),
]
