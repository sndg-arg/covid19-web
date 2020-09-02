# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os
import json
# from django.shortcuts import redirect, reverse
from django.http.response import HttpResponse
from django.views.generic import TemplateView
from config.settings.base import STATICFILES_DIRS
from ..tasks import variant_graphics
from bioseq.models.Variant import Variant
from bioseq.models.PDBVariant import PDBVariant
from itertools import groupby


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

        variant = Variant.objects.prefetch_related("pdb_variants__residue__pdb",
                                                   "pdb_variants__residue__residue_sets__pdbresidue_set"
                                                   ).get(bioentry__accession=context["gene"], pos=pos)
        # "pdb_variants__residue__residue_sets__pdbresidue_set"
        context["ref"] = variant.ref
        context["gene_id"] = variant.bioentry_id
        context["gene_desc"] = variant.bioentry.description
        context["residues"] = []
        for pdb_variant in variant.pdb_variants.all():
            context["residues"].append(pdb_variant.residue)
            ann = pdb_variant.ann()
            pdb_variant.residue.ann = [x["desc"] for x in ann]
            pdb_variant.residue.layers = [rs["name"] for rs in ann]

        context["residues"] = [(x, list(y)) for x, y in
                               groupby(sorted(context["residues"], key=lambda x: x.pdb.code), lambda x: x.pdb.code)]
        context["fig_avail"] = os.path.exists(f'{STATICFILES_DIRS[0]}/auto/posfigs/{gene}{pos}.png')

        return context


def pdb_variants_download(request):
    response = HttpResponse(content_type='text/json')
    if "download" in request.GET:
        response['Content-Disposition'] = 'attachment; filename="var2pdb.json"'

    fields = {"variant__pos": "pos", "variant__ref": "ref",
              "residue__chain": "chain", "residue__resid": "resid",
              "residue__pdb__code": "pdb", "variant__bioentry__accession": "gene",
              "residue__residue_sets__pdbresidue_set__name": "residue_set",
              "residue__residue_sets__pdbresidue_set__residue_set__name": "residue_set_type",
              "residue__residue_sets__pdbresidue_set__description": "residue_set_desc"
              }
    qs = {} if not request.GET.get("country", None) else (
        {"variant__sample_variants__sample__country__in": request.GET["country"].split(",")})
    pdb_vars = [{v: x[k] for k, v in fields.items()} for x in
                PDBVariant.objects.filter(**qs).values(*fields.keys()).order_by("variant__bioentry__accession",
                                                                                "-variant__pos")]

    fields = {"variant__pos": "pos", "variant__sample_variants__alt": "alt",
              "residue__pdb__code": "pdb", "variant__bioentry__accession": "gene",
              "variant__sample_variants__sample__country": "country",
              "variant__sample_variants__sample__name":"sample_name"
              }
    country = [{v: x[k] for k, v in fields.items()} for x in
               PDBVariant.objects.filter(**qs).values(*fields.keys()).order_by("variant__bioentry__accession",
                                                                               "-variant__pos")]
    country2 = {gene: {pos: [{"country": country, "alt": alt, "count": len(set(xx["sample_name"] for xx in alts)) } for (country, alt), alts in
                             groupby(sorted(countries, key=lambda x: (x["country"], x["alt"])),
                                     lambda x: (x["country"], x["alt"]))]
                       for pos, countries in groupby(sorted(g2, key=lambda x: x["pos"]), lambda x: x["pos"])}
                for gene, g2 in groupby(country, lambda x: x["gene"])
                }

    data = {gene: [{"pos": pos + 1, "ref": ref, "countries": country2[gene][pos], "residues":
        [{k: residue[k] for k in ["residue_set_type", "residue_set", "chain", "resid", "residue_set_desc"] if
          residue[k]} for residue in residues]
                    }
                   for (pos, ref), residues in
                   groupby(sorted(g2, key=lambda x: x["pos"]), lambda orow: (orow["pos"], orow["ref"]))]
            for gene, g2 in groupby(pdb_vars, lambda row: row["gene"])
            }
    json.dump(data, response, indent=4, sort_keys=True)
    return response
