from django.urls import path

from . import views

app_name = 'pdbdb'
urlpatterns = [

    # path('structure/<str:pdbid>', views.StructureView.as_view(), name='structure_view'),
    path('structure_raw/<str:pdbid>', views.structure_raw, name='structure_raw_view'),



]
