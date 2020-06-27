# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
# from django.shortcuts import redirect, reverse

from django.db.models import Q
from django.views.generic import TemplateView

from django.contrib.auth.mixins import LoginRequiredMixin
from django.contrib.auth.mixins import PermissionRequiredMixin

class ProjectStrainsView(LoginRequiredMixin,TemplateView):
    # PermissionRequiredMixin permission_required = 'polls.add_choice'
    # login_url = '/login/'
    # redirect_field_name = 'redirect_to'
    template_name = "project_view.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)



