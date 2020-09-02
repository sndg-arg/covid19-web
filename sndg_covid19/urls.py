from django.urls import path

from .views.AssemblyView import assembly_view
from .views.ProteinView import ProteinView,MSAView
from .views.StructureView import StructureView,pdb_download,StructureStaticView
from .views.ProjectStrainsView import ProjectStrainsView
from .views.VariantView import VariantView,pdb_variants_download




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
    path('variant/', pdb_variants_download, name='download_pdb_variants'),
    # path('static/jbrowse/data/COVID19/trackList2.json', assembly_view, name='genome_view'),
]
