from django.urls import path

from . import views
from .controllers import tools

app_name = 'bioseq'
urlpatterns = [
    path('about', views.AboutView.as_view(), name='index'),
    path('tax/<int:pk>', views.TaxView.as_view(), name='tax_view'),
    path('seq/<int:pk>', views.sequence_view, name='seq_view'),
    path('assembly/<int:pk>', views.assembly_view, name='assembly_view'),
    path('seq/<str:locus_tag>', views.sequence_lt_view, name='sequence_lt_view'),

    path('geneproducts', views.gene_product_list_view, name='gene_product_list'),
    path('proteins/<int:assembly_id>', views.protein_list_view, name='protein_list'),


    path('variant/<int:pk>', views.TaxView.as_view(), name='variant_view'),

    path('blast', tools.blast, name='blast_view'),
    path('blast_result', tools.blast_result, name='blast_result_view'),
    path('msa/', views.AboutView.as_view(), name='msa_view'),
    path('primer/', views.AboutView.as_view(), name='primer_view'),



    path('analysis/<int:pk>', views.TaxView.as_view(), name='analysis_view'), # tree / msa / blast



    # path('pathway/<int:pk>', views.TaxView.as_view(), name='tax_view'),

    # path('publications', views.FilteredPersonListView.as_view(), name='publications'),
    # path('search/', views.BioSearchView.as_view(), name='search_view'),


]