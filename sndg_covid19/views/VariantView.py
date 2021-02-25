# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.contrib.auth.mixins import LoginRequiredMixin
import os
import json
# from django.shortcuts import redirect, reverse
from collections import defaultdict

from django.http.response import HttpResponse
from django.views.generic import TemplateView

from bioseq.models.Bioentry import Bioentry
from config.settings.base import STATICFILES_DIRS
from ..tasks import variant_graphics
from bioseq.models.Variant import Variant, Sample,SampleVariant
from bioseq.models.PDBVariant import PDBVariant
from itertools import groupby
import pandas as pd
import numpy as np
from datetime import date
from dateutil.relativedelta import relativedelta
from . import latam_countries


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


# LoginRequiredMixin,
class InmunovaView(TemplateView):
    # PermissionRequiredMixin permission_required = 'polls.add_choice'
    # login_url = '/login/'
    # redirect_field_name = 'redirect_to'
    template_name = "inmunova_view.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        # gene = context["gene"]
        be = (Bioentry.objects
              .prefetch_related("seq", "qualifiers__term", "qualifiers__term__dbxrefs__dbxref",
                                # "dbxrefs__dbxref","qualifiers__term__dbxrefs__dbxref",
                                "features__locations", "features__source_term",  # "dbxrefs__dbxref",
                                "features__type_term__ontology", "features__qualifiers__term"))
        be = be.get(accession="S")
        subdivision = self.request.GET.get("subdivision","")
        vars, month_counts = variants(be, "Argentina",subdivision)
        context["variants"] = vars
        context["month_counts"] = month_counts
        context["object"] = be
        context["year_months"] = sorted([year * 100 + month for year, month in get_last_months(date.today(), 11)])
        context["subdivisions"] = [x["subdivision"] for x in Sample.objects.values("subdivision").distinct()]
        context["object"] = be
        context["selectedSub"] = subdivision
        from bioseq.templatetags.bioresources_extras import getattribute
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
              "variant__sample_variants__sample__name": "sample_name"
              }
    country = [{v: x[k] for k, v in fields.items()} for x in
               PDBVariant.objects.filter(**qs).values(*fields.keys()).order_by("variant__bioentry__accession",
                                                                               "-variant__pos")]
    country2 = {gene: {pos: [{"country": country, "alt": alt, "count": len(set(xx["sample_name"] for xx in alts))} for
                             (country, alt), alts in
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


def variants(be, filter_by_country=None,filter_by_subdivision=None):
    data = []
    if filter_by_country:
        vs = list(SampleVariant.objects.values("variant__ref", "alt", "variant__pos", "variant__variant_id",
                                               "sample__country", "variant__bioentry_id",
                                               "sample__name",
                                               "sample__date").filter(
            variant__bioentry=be, sample__country=filter_by_country).distinct())
        if filter_by_subdivision:
            vs = list(SampleVariant.objects.values("variant__ref", "alt", "variant__pos", "variant__variant_id",
                                             "sample__country", "variant__bioentry_id",
                                             "sample__name",
                                             "sample__date").filter(
                variant__bioentry=be, sample__subdivision=filter_by_subdivision,
                sample__country=filter_by_country).distinct())
    else:
        vs = list(SampleVariant.objects.values("variant__ref", "alt", "variant__pos", "variant__variant_id",
                                               "sample__country", "variant__bioentry_id",
                                               "sample__name",
                                               "sample__date").filter(variant__bioentry=be))

    variant_ids = []
    for sample_variant in vs:
        # if ((sample_variant["ref"] != "X")
        # and            (sample_variant["sample_variants__sample__country"] == "Argentina")):
        if (sample_variant["variant__ref"] != "X"):
            # key = "_".join([ str(pdb_var["pos"]), pdb_var["ref"]])
            variant_ids.append(sample_variant["variant__variant_id"])

            record = {"ref": sample_variant["variant__ref"], "pos": sample_variant["variant__pos"] + 1,
                      "alt": sample_variant["alt"],
                      "count": 1.0,
                      "date": sample_variant['sample__date'].year * 100 +
                              sample_variant['sample__date'].month,
                      "cod": sample_variant["sample__name"]}
            if 317 < record["pos"] < 542:
                data.append(record)
    df_raw = pd.DataFrame(data)
    if len(df_raw):
        df = df_raw.pivot_table(index=["pos", "ref", "alt"], columns=["date"], values="count",
                                aggfunc=np.sum).fillna(0).astype(int).sort_values(
            ["pos", "ref", "alt"])
    else:
        df = pd.DataFrame()

    month_counts = {}
    for year, month in get_last_months(date.today(), 12):
        if filter_by_subdivision:
            month_counts[year * 100 + month] = Sample.objects.filter(date__year=year, date__month=month,
                                                    country="Argentina",subdivision=filter_by_subdivision).count()
        else:
            month_counts[year * 100 + month] = Sample.objects.filter(date__year=year, date__month=month,
                                                                     country="Argentina").count()

    fields = {"variant__pos": "pos", "variant__ref": "ref",
              # "residue__chain": "chain", "residue__resid": "resid",
              "residue__pdb__code": "pdb",
              # "residue__residue_sets__pdbresidue_set__name": "residue_set",
              "residue__residue_sets__pdbresidue_set__residue_set__name": "residue_set_type",
              # "residue__residue_sets__pdbresidue_set__description": "residue_set_desc"
              }
    pdb_vars = [{v: x[k] for k, v in fields.items()} for x in
                PDBVariant.objects.filter(variant__bioentry=be, variant__in=variant_ids).values(*fields.keys())]

    pdb_vars_dict = defaultdict(list)
    for pdb_var in pdb_vars:
        key = "_".join([str(pdb_var["pos"] + 1), pdb_var["ref"]])
        pdb_vars_dict[key].append(pdb_var)
    pdb_vars_dict = dict(pdb_vars_dict)
    pdb_vars_dict2 = defaultdict(list)

    for k, v in pdb_vars_dict.items():
        pdbs = []
        sites = []
        for strd in v:
            # agg[strd["pdb"]].append(" ".join([str(strd[x]) for x in ['chain','resid',
            #             'residue_set_type','residue_set','residue_set_desc'] if strd[x] ] ))
            pdbs.append(strd["pdb"])
            if strd["residue_set_type"]:
                # ,'residue_set','residue_set_desc'
                sites.append(" ".join([str(strd[x]) for x in ['residue_set_type']]))
        # pdb_vars_dict[k] = [{"pdb": k, "res": v,"sites":list(set(sites))} for k, v in agg.items()]
        pdb_vars_dict2[k] = {"pdbs": list(set(pdbs)), "sites": list(set(sites))}
    pdb_vars_dict = dict(pdb_vars_dict2)

    variants = []
    if len(df):

        for i, r in df.reset_index().iterrows():
            key = "_".join([str(r["pos"]), r["ref"]])
            v = r.to_dict()
            v["struct"] = pdb_vars_dict.get(key, "")
            for year_month, sample_count in month_counts.items():
                if year_month in v:
                    v[str(
                        year_month) + "_total"] = sample_count  # f'{sample_count/total_month}({sample_count}/{total_month})'

            variants.append(v)

    return variants, month_counts


def get_last_months(start_date, months):
    for i in range(months):
        yield (start_date.year, start_date.month)
        start_date += relativedelta(months=-1)
