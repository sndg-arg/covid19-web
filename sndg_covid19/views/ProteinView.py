# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
# from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioseq.models.Bioentry import Bioentry
from bioseq.models.Biodatabase import Biodatabase
from bioseq.models.Seqfeature import Seqfeature
from bioseq.models.Biosequence import Biosequence
from bioseq.models.Dbxref import DBx, Dbxref
from pdbdb.models import PDB
from pdbdb.models import Property

from bioseq.models.Variant import Variant

# from ..templatetags.bioresources_extras import split
from . import latam_countries


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
    msa = be.accession + ".faa"
    grouped_features = {k: [{y: z for y, z in v.__dict__.items() if y != '_state'} for v in vs] for k, vs in
                        be.groupedFeatures().items()}

    return render(request, 'gene_view.html', {"grouped_features": grouped_features,
                                              "functions": functions, "assembly": "assembly", "msa": msa,
                                              "variants": prot_variants,
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


def variants(protein_entry):
    """
    https://pivottable.js.org/examples/
    :param protein_entry:
    :return:
    """
    data = []
    vs = Variant.objects.filter(bioentry=protein_entry).values("ref", "sample_variants__alt", "pos",
                                                               "sample_variants__sample__country",
                                                               "sample_variants__sample__name")
    for sample_variant in vs:
        if ((sample_variant["ref"] != "X") and
            (sample_variant["sample_variants__sample__country"] in latam_countries)):
            record = {"ref": sample_variant["ref"], "pos": sample_variant["pos"] + 1,
                      "alt": sample_variant["sample_variants__alt"],
                      "country": sample_variant["sample_variants__sample__country"],
                      "cod": sample_variant["sample_variants__sample__name"]}
            data.append(record)
    return data
