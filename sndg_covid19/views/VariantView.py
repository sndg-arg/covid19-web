# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os

# from django.shortcuts import redirect, reverse
from django.views.generic import TemplateView
from config.settings.base import STATICFILES_DIRS
from ..tasks import variant_graphics
from bioseq.models.Variant import Variant

class VariantView(TemplateView):
    # PermissionRequiredMixin permission_required = 'polls.add_choice'
    # login_url = '/login/'
    # redirect_field_name = 'redirect_to'
    template_name = "variant_view.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        gene = context["gene"]
        pos = context["pos"] - 1

        context["gene"] = self.kwargs["gene"]
        context["pos"] = self.kwargs["pos"]


        variant = Variant.objects.get(bioentry__accession=context["gene"],pos=pos)

        context["ref"] = variant.ref
        context["gene_id"] = variant.bioentry_id
        context["gene_desc"] = variant.bioentry.description



        context["fig_avail"] =  os.path.exists( f'{STATICFILES_DIRS[0]}/auto/posfigs/{gene}{pos}.png')



        return context
