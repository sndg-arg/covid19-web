# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os

# from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioseq.models.Bioentry import Bioentry
from bioseq.models.Biosequence import Biosequence
from bioseq.models.Dbxref import DBx, Dbxref
from bioseq.models.Seqfeature import Seqfeature
from bioseq.models.Variant import Variant
from bioseq.models.PDBVariant import PDBVariant
from config.settings.base import STATICFILES_DIRS
from pdbdb.models import PDB, Property

# from ..templatetags.bioresources_extras import split
from . import latam_countries
from collections import defaultdict
import pandas as pd
import numpy as np


# "GU280_gp11": "orf10_prot.fasta", ??


def url_map(db_map, dbname, accession):
    if dbname in db_map:
        dbx = db_map[dbname]
        try:
            return dbx.url_template % accession
        except:
            print(dbx.url_template)
            return (dbx.url_template)
    return "www.google.com?q=" + accession


def ProteinView(request, pk):
    be = (Bioentry.objects
          .prefetch_related("seq", "qualifiers__term", "qualifiers__term__dbxrefs__dbxref",
                            # "dbxrefs__dbxref","qualifiers__term__dbxrefs__dbxref",
                            "features__locations", "features__source_term",  # "dbxrefs__dbxref",
                            "features__type_term__ontology", "features__qualifiers__term"))

    dbxss = list(Dbxref.objects.filter(dbxrefs__bioentry_id=pk))
    # .select_related("biodatabase").select_related("taxon")
    be = be.get(bioentry_id=pk)
    be.genes()

    sfid = be.qualifiers_dict()["GFeatureId"]
    feature = Seqfeature.objects.prefetch_related("locations").filter(bioentry__biodatabase__name="COVID19",
                                                                      seqfeature_id=sfid).get()
    # [x.term  for x in be.qualifiers.all() if x.term.ontology.name == "Gene Ontology"][0].dbxrefs.all()

    locations = list(feature.locations.all())
    start = locations[0].start_pos
    end = locations[-1].end_pos

    seq = Biosequence.objects.raw("""
    SELECT bioentry_id, version , length , alphabet ,SUBSTRING( seq,%i,%i ) seq
    FROM biosequence WHERE bioentry_id = %i ;
    """ % (start, end - start, feature.bioentry_id))[0]
    # UnipName
    functions = go_function(be)
    pfeatures = protein_features(be)
    structures = protein_structures(dbxss)
    pdbxrefs = dbxrefs(dbxss)
    prot_variants = variants(be)
    # prot_variants_table = variants_table(be)

    msa = be.accession + ".faa"
    msa_file = f'{STATICFILES_DIRS[0]}/ORFs/{msa}'
    if not os.path.exists(msa_file):
        msa = None
    grouped_features = {k: [{y: z for y, z in v.__dict__.items() if y != '_state'} for v in vs] for k, vs in
                        be.groupedFeatures().items()}

    return render(request, 'gene_view.html', {"grouped_features": grouped_features,
                                              "functions": functions, "assembly": "assembly", "msa": msa,
                                              "variants": prot_variants, "latam_countries": latam_countries,
                                              # "variants_table": prot_variants_table,
                                              "object": be, "feature": feature, "seq": seq, "start": start, "end": end,
                                              "protein_features": pfeatures, "structures": structures,
                                              "dbxrefs": pdbxrefs,
                                              "sidebarleft": 1})


def dbxrefs(dbxss):
    pdbxrefs = []

    dbs = []
    for dbxref in dbxss:
        if dbxref.dbname.lower() not in ["go", "pdb"]:
            dbs.append(dbxref.dbname)
    dbs = set(dbs)
    dbmap = {x.name: x for x in DBx.objects.filter(name__in=dbs)}

    for dbxref in dbxss:
        if dbxref.dbname.lower() not in ["go", "pdb", "pdbsum"]:
            pdbxrefs.append({
                "dbname": dbxref.dbname,
                "accession": dbxref.accession,
                "url": url_map(dbmap, dbxref.dbname, dbxref.accession)
            })

    return pdbxrefs


def go_function(be):
    functions = {"biological_process": [], "molecular_function": [], "cellular_component": []}
    for qual in [x for x in be.qualifiers.all()]:
        # term__dbxrefs__dbxref__accession="goslim_generic"
        for dbxref in qual.term.dbxrefs.all():
            if dbxref.dbxref.dbname == "go":
                if dbxref.dbxref.accession in functions:
                    functions[dbxref.dbxref.accession].append(qual.term)
    return functions


def protein_structures(dbxss):
    structures = []
    pdbs = []
    for dbxref in dbxss:
        if dbxref.dbname.lower() == "pdb":
            pdbs.append(dbxref.accession.lower())
    pdb_map = {pdb.code: pdb for pdb in
               PDB.objects.prefetch_related("residue_sets__properties__property", "residue_sets__residue_set").filter(
                   code__in=set(pdbs))}

    for pdb_code in pdbs:
        druggability = ""
        header = ""
        if pdb_code in pdb_map:
            pdb = pdb_map[pdb_code]
            header = pdb.header
            p = pdb.max_druggability_pocket()

            druggability = p.properties_dict()[Property.druggability] if p else ""

        structure = {
            "code": pdb_code,
            "name": header,
            "type": "PDB",
            "druggability": druggability,
            "info": ""
        }
        structures.append(structure)
    return structures


def protein_features(protein_entry):
    pfeatures = []
    for feature in protein_entry.features.all():
        if feature.type_term.identifier not in 'SO:0001079':
            pfeatures.append(feature)
            if "InterPro" in feature.qualifiers_dict():
                feature.display_name = feature.display_name + " (" + feature.qualifiers_dict()["InterPro"] + ")"

    return pfeatures


def variants(be):
    data = []
    vs = Variant.objects.values("ref", "sample_variants__alt", "pos", "variant_id",
                                "sample_variants__sample__country", "bioentry_id",
                                "sample_variants__sample__name").filter(bioentry=be)

    variant_ids = []
    for sample_variant in vs:
        if ((sample_variant["ref"] != "X") and
            (sample_variant["sample_variants__sample__country"] in latam_countries)):
            # key = "_".join([ str(pdb_var["pos"]), pdb_var["ref"]])
            variant_ids.append(sample_variant["variant_id"])
            record = {"ref": sample_variant["ref"], "pos": sample_variant["pos"] + 1,
                      "alt": sample_variant["sample_variants__alt"],

                      "bioentry_id": be.bioentry_id,
                      "country": sample_variant["sample_variants__sample__country"],
                      "count": 1.0,
                      "cod": sample_variant["sample_variants__sample__name"]}
            data.append(record)
    df_raw = pd.DataFrame(data)

    df = df_raw.pivot_table(index=["pos", "ref", "alt"], columns=["country"], values="count",
                            aggfunc=np.sum).fillna(0).astype(int).sort_values(
        ["pos", "ref", "alt"])

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

    for i, r in df.reset_index().iterrows():
        key = "_".join([str(r["pos"]), r["ref"]])
        v = r.to_dict()
        v["struct"] = pdb_vars_dict.get(key, "")
        variants.append(v)

    return variants
