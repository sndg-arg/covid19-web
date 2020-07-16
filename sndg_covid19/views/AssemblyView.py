# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db.models import Count
from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioseq.models.Biodatabase import Biodatabase
from bioseq.models.Bioentry import Bioentry, BioentryDbxref
from bioseq.models.Biosequence import Biosequence
from bioseq.models.Variant import SampleVariant, Variant
from . import latam_countries
import pandas as pd
import numpy as np
import json
from collections import defaultdict


# from bioseq.models.Dbxref import Dbxref


def variants_table():
    data = []
    vs = Variant.objects.values("ref", "sample_variants__alt", "pos", "bioentry__accession",
                                "sample_variants__sample__country", "bioentry_id",
                                "sample_variants__sample__name")
    pdb_vars = Variant.objects.values("ref", "pos", "bioentry__accession",
                                      "pdb_variants__residue__resid",
                                      "pdb_variants__residue__chain",
                                      "pdb_variants__residue__pdb__code",
                                      "pdb_variants__residue__icode",
                                      )
    pdb_vars_dict = defaultdict(list)
    for pdb_var in pdb_vars:
        key = "_".join([pdb_var["bioentry__accession"], str(pdb_var["pos"] + 1), pdb_var["ref"]])
        if pdb_var["pdb_variants__residue__pdb__code"]:
            v = "_".join([pdb_var["pdb_variants__residue__pdb__code"], str(pdb_var["pdb_variants__residue__chain"]),
                          str(pdb_var["pdb_variants__residue__resid"])])
            if pdb_var["pdb_variants__residue__icode"].strip():
                key = key + pdb_var["pdb_variants__residue__icode"]
            pdb_vars_dict[key].append(v)
    pdb_vars_dict = dict(pdb_vars_dict)
    for sample_variant in vs:
        if ((sample_variant["ref"] != "X") and
            (sample_variant["sample_variants__sample__country"] in latam_countries)):
            key = "_".join([pdb_var["bioentry__accession"], str(pdb_var["pos"]), pdb_var["ref"]])
            record = {"ref": sample_variant["ref"], "pos": sample_variant["pos"] + 1,
                      "alt": sample_variant["sample_variants__alt"],
                      "name": sample_variant["bioentry__accession"],
                      "bioentry_id": sample_variant["bioentry__accession"],
                      "country": sample_variant["sample_variants__sample__country"],
                      "count": 1.0,
                      "cod": sample_variant["sample_variants__sample__name"]}
            data.append(record)
    df_raw = pd.DataFrame(data)

    df = df_raw.pivot_table(index=["name", "pos", "ref", "alt"], columns=["country"], values="count",
                            aggfunc=np.sum, margins=True, margins_name='Total').fillna(0).astype(int).sort_values(
        ["name", "pos", "ref", "alt"])
    keys = ["_".join([str(x) for x in r.name[:3]]) for _, r in df.iterrows()]
    # con_struct = [(" ".join(pdb_vars_dict[k]) if k in pdb_vars_dict else "") for k in keys]
    con_struct = [("X" if k in pdb_vars_dict else "") for k in keys]

    # df_raw["pos"] = [f'<a href="{reverse("covid:pos_stats", kwargs={"gene": r["name"], "pos": r.pos})}">{r.pos}</a>'
    #                  for _, r in df_raw.iterrows()]


    df["En PDB"] = con_struct
    df.columns.name = "País"
    df.index.names = ["Gen", "Pos", "Ref", "Alt"]
    html = df.to_html(table_id="variants_table",
                      classes="table table-responsive", escape=False)
    return html


def assembly_view(request):
    params = {}
    bdb = Biodatabase.objects.get(name="COVID19")
    contig = bdb.entries.first()
    bdb_ids = [x.biodatabase_id for x in Biodatabase.objects.filter(name__startswith="COVID19_")]

    lengths = {}
    seqs = Biosequence.objects.prefetch_related("bioentry").raw("""
        SELECT s.bioentry_id, s.length
        FROM biosequence s,bioentry b WHERE b.biodatabase_id IN( %i,%i ) AND  b.bioentry_id = s.bioentry_id ;
        """ % tuple(bdb_ids))
    for seq in seqs:
        lengths[seq.bioentry.accession] = seq.length

    features = list(contig.features.exclude(type_term__identifier__in=["source", "gene", "stem_loop"]
                                            ).prefetch_related("source_term", "type_term", "locations",
                                                               "qualifiers__term__dbxrefs__dbxref"))
    features = sorted(features, key=lambda x: x.first_location().start_pos)
    properties = {}
    feature_ids = [
        x.qualifiers_dict()["BioentryId"] for x in features if "BioentryId" in x.qualifiers_dict()]

    bioentries = Bioentry.objects.prefetch_related("qualifiers__term",  # "dbxrefs__dbxref"
                                                   ).filter(biodatabase_id__in=bdb_ids)  #

    dbxss = {x.bioentry_id: [] for x in bioentries}
    for x in BioentryDbxref.objects.prefetch_related("dbxref").filter(bioentry__biodatabase_id__in=bdb_ids,
                                                                      dbxref__dbname="PDB"):
        dbxss[x.bioentry_id].append(x.dbxref.accession)

    variants = {v["variant__bioentry_id"]: v["count"] for v in
                SampleVariant.objects.exclude(alt="X").values("variant__bioentry_id").annotate(
                    count=Count('variant__pos', distinct=True))}

    for bioentry in bioentries:
        # data = bioentry.qualifiers_dict()

        properties[bioentry.bioentry_id] = {
            "structures": len(dbxss[bioentry.bioentry_id]),
            "description": bioentry.description,
            "variants": variants.get(bioentry.bioentry_id, 0)
        }

    for f in features:
        if "BioentryId" in f.qualifiers_dict() and int(f.qualifiers_dict()["BioentryId"]) in properties:
            f.extra_gene_props = properties[int(f.qualifiers_dict()["BioentryId"])]

    params = {"query": "",
              "lengths": lengths, "variants": variants_table(),
              "genes": features, "sidebarleft": {}}
    return render(request, 'genome_view.html', params)
