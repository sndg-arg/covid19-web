from django.urls import path

from .views.AssemblyView import assembly_view
from .views.ProteinView import ProteinView
from .views.StructureView import StructureView
from .views.ProjectStrainsView import ProjectStrainsView
from .views.VariantView import VariantView




from django.views.generic import RedirectView

from django.conf import settings

app_name = 'covid'
# if settings.MINCYT_URL:

urlpatterns = [


    path('', assembly_view, name='genome_view'),
    path('<int:pk>', ProteinView, name='gene_view'),
    path('structure/<int:prot_id>/<str:pdbid>', StructureView.as_view(), name='structure_view'),

    path('project/', ProjectStrainsView.as_view(), name='project_strains'),
    path('variant/<str:gene>/<int:pos>', VariantView.as_view(), name='pos_stats'),
    # path('static/jbrowse/data/COVID19/trackList2.json', assembly_view, name='genome_view'),
]
