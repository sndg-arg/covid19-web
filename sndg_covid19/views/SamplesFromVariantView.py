from datetime import datetime

from django.contrib.auth.mixins import LoginRequiredMixin
from django.views.generic import TemplateView

from bioseq.models.Variant import Sample, SampleVariant


class SamplesFromVariantView(LoginRequiredMixin, TemplateView):
    # PermissionRequiredMixin permission_required = 'polls.add_choice'
    # login_url = '/login/'
    # redirect_field_name = 'redirect_to'
    template_name = "samples_from_variant_view.html"

    def get_context_data(self, gene, pos, **kwargs):
        context = super().get_context_data(**kwargs)
        year_month = datetime.strptime(self.request.GET["year_month"],"%Y%m")
        context["year_month"] = self.request.GET["year_month"]
        context["pos"] = pos
        context["gene"] = gene
        vs = list(SampleVariant.objects.values("sample__subdivision", "alt",
                                               "sample__name").filter(
            sample__date__month=year_month.month,sample__date__year=year_month.year,
            variant__bioentry__accession=gene, sample__country="Argentina", variant__pos=pos - 1).order_by(
            "alt").distinct())
        context["samples"] = [
            {"alt": x["alt" ], "subdivision": x["sample__subdivision"], "name": x["sample__name"]} for x in vs]
        return context
