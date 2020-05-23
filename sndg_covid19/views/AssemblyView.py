# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
# from django.shortcuts import redirect, reverse
from django.shortcuts import render

from bioseq.models.Biodatabase import Biodatabase
from bioseq.models.Biosequence import Biosequence
from bioseq.models.Bioentry import Bioentry,BioentryDbxref
# from bioseq.models.Dbxref import Dbxref



from django.views.decorators.clickjacking import xframe_options_exempt



def assembly_view(request):
    params = {}
    bdb = Biodatabase.objects.get(name="COVID19")
    contig =  bdb.entries.first()
    bdb_ids = [x.biodatabase_id for x in Biodatabase.objects.filter(name__startswith="COVID19_")]

    lengths = {}
    seqs = Biosequence.objects.prefetch_related("bioentry").raw("""
        SELECT s.bioentry_id, s.length
        FROM biosequence s,bioentry b WHERE b.biodatabase_id IN( %i,%i ) AND  b.bioentry_id = s.bioentry_id ;
        """ % tuple(bdb_ids) )
    for seq in seqs:
        lengths[seq.bioentry.accession] = seq.length

    features = list(contig.features.exclude(type_term__identifier__in=["source","gene","stem_loop"]
        ).prefetch_related("source_term","type_term","locations",
                           "qualifiers__term__dbxrefs__dbxref"))
    features = sorted(features, key=lambda x:x.first_location().start_pos)
    properties = {}
    feature_ids = [
        x.qualifiers_dict()["BioentryId"] for x in features if "BioentryId" in x.qualifiers_dict()]

    bioentries = Bioentry.objects.prefetch_related("qualifiers__term", # "dbxrefs__dbxref"
        ).filter(biodatabase_id__in=bdb_ids) #

    dbxss = {x.bioentry_id:[] for x in bioentries}
    for x in BioentryDbxref.objects.prefetch_related("dbxref").filter(bioentry__biodatabase_id__in=bdb_ids,
                                            dbxref__dbname="PDB"):
        dbxss[x.bioentry_id].append(x.dbxref.accession)

    for bioentry in bioentries:
        # data = bioentry.qualifiers_dict()

        properties[bioentry.bioentry_id] = {
            "structures":len( dbxss[bioentry.bioentry_id] ),
            "description": bioentry.description
        }

    for f in features:
        if "BioentryId" in f.qualifiers_dict() and  int(f.qualifiers_dict()["BioentryId"]) in properties:
            f.extra_gene_props = properties[int(f.qualifiers_dict()["BioentryId"])]


    params = {"query": "",
              "lengths": lengths,
              "genes": features, "sidebarleft": {}}

    return render(request, 'genome_view.html', params)
