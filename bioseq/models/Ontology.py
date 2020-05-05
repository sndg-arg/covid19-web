# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.db import models
from django.shortcuts import reverse

from .Dbxref import Dbxref

class Ontology(models.Model):
    GO = "Gene Ontology"
    SO = "Sequence Ontology"
    EC = "Enzyme Commission number"
    SFS = "SeqFeature Sources"
    SFK = 'SeqFeature Keys'
    GRAPH = "Graph"
    ANNTAGS = "Annotation Tags"

    ontology_id = models.AutoField(primary_key=True)
    name = models.CharField(unique=True, max_length=32)
    definition = models.TextField(blank=True, null=True)

    class Meta:
        managed = True
        db_table = 'ontology'

    def __str__(self):
        return "%s  (%i)" % (self.name, self.ontology_id)

    @staticmethod
    def load_ann_terms():
        from .Term import Term
        sfs_ontology = Ontology.objects.get_or_create(name=Ontology.SFS)[0]
        Term.objects.get_or_create(identifier="manual", name="manual", version=1, ontology=sfs_ontology,
                                   definition="added or corrected by a person")
        Term.objects.get_or_create(identifier="bibliography", name="bibliography", version=1, ontology=sfs_ontology,
                                   definition="found in bibliography")
        Term.objects.get_or_create(identifier="experimental", name="experimental", version=1, ontology=sfs_ontology,
                                   definition="the annotation was obtained experimentally")
        Term.objects.get_or_create(identifier="calculated", name="calculated", version=1, ontology=sfs_ontology,
                                   definition="the annotation was obtained using software")
        Term.objects.get_or_create(identifier="other", name="other", version=1, ontology=sfs_ontology,
                                   definition="")
        sfk_ontology = Ontology.objects.get_or_create(name=Ontology.SFK)[0]
        Term.objects.get_or_create(identifier="gene", name="gene", version=1, ontology=sfk_ontology)
        Term.objects.get_or_create(identifier="CDS", name="CDS", version=1, ontology=sfk_ontology)
        Term.objects.get_or_create(identifier="tRNA", name="tRNA", version=1, ontology=sfk_ontology)
        Term.objects.get_or_create(identifier="rRNA", name="rRNA", version=1, ontology=sfk_ontology)
        Term.objects.get_or_create(identifier="ncRNA", name="ncRNA", version=1, ontology=sfk_ontology)
        Term.objects.get_or_create(identifier="exon", name="exon", version=1, ontology=sfk_ontology)
        Term.objects.get_or_create(identifier="mRNA", name="mRNA", version=1, ontology=sfk_ontology)
        Term.objects.get_or_create(identifier="miRNA", name="miRNA", version=1, ontology=sfk_ontology)

        ann_ontology = Ontology.objects.get_or_create(name=Ontology.ANNTAGS)[0]
        Term.objects.get_or_create(identifier="locus_tag", name="locus_tag", version=1, ontology=ann_ontology)
        Term.objects.get_or_create(identifier="gene", name="gene", version=1, ontology=ann_ontology)
        Term.objects.get_or_create(identifier="product", name="product", version=1, ontology=ann_ontology)



    @staticmethod
    def load_go_base():
        from .Term import Term
        graph_ontology = Ontology.objects.get_or_create(name=Ontology.GRAPH, definition="")[0]

        is_a = Term.objects.get_or_create(
                identifier="is_a", name="is_a", version=1, ontology=graph_ontology, definition="")[0]

        ontology = Ontology.objects.get_or_create(name=Ontology.GO)[0]

        Ontology.dbmap = {
            "biological_process": Dbxref.objects.get_or_create(dbname="go", accession="biological_process", version=1)[
                0],
            "molecular_function": Dbxref.objects.get_or_create(dbname="go", accession="molecular_function", version=1)[
                0],
            "cellular_component": Dbxref.objects.get_or_create(dbname="go", accession="cellular_component", version=1)[
                0],
            "biosapiens": Dbxref.objects.get_or_create(dbname="so", accession="biosapiens", version=1)[0],
            "DBVAR": Dbxref.objects.get_or_create(dbname="so", accession="DBVAR", version=1)[0],
            "SOFA": Dbxref.objects.get_or_create(dbname="so", accession="SOFA", version=1)[0],
        }

        for x in """goslim_agr "AGR slim"
        goslim_aspergillus "Aspergillus GO slim"
        goslim_candida "Candida GO slim"
        goslim_chembl "ChEMBL protein targets summary"
        goslim_generic "Generic GO slim"
        goslim_goa "GOA and proteome slim"
        goslim_metagenomics "Metagenomics GO slim"
        goslim_mouse "Mouse GO slim"
        goslim_pir "PIR GO slim"
        goslim_plant "Plant GO slim"
        goslim_pombe "Fission yeast GO slim"
        goslim_synapse "synapse GO slim"
        goslim_virus "Viral GO slim"
        goslim_yeast "Yeast GO slim"
        gosubset_prok "Prokaryotic GO subset"
        virus_checked "Viral overhaul terms" """.split("\n"):
            k = x.strip().split(' ')[0]
            Ontology.dbmap[k] = Dbxref.objects.get_or_create(dbname="go", accession=k, version=1)[0]

        part_of = Term.objects.get_or_create(identifier="part_of", ontology=graph_ontology)[0]
        regulates = Term.objects.get_or_create(identifier="regulates", ontology=graph_ontology)[0]
        negatively_regulates = Term.objects.get_or_create(identifier="negatively_regulates", ontology=graph_ontology)[0]
        positively_regulates = Term.objects.get_or_create(identifier="positively_regulates", ontology=graph_ontology)[0]

        has_quality = Term.objects.get_or_create(identifier="has_quality", name="has_quality", version=1, ontology=
        graph_ontology)[0]
        derives_from = Term.objects.get_or_create(identifier="derives_from", name="derives_from", version=1, ontology=
        graph_ontology)[0]
        has_origin = Term.objects.get_or_create(identifier="has_origin", name="has_origin", version=1, ontology=
        graph_ontology)[0]
        has_part = Term.objects.get_or_create(identifier="has_part", name="has_part", version=1, ontology=
        graph_ontology)[0]
        transcribed_to = \
            Term.objects.get_or_create(identifier="transcribed_to", name="transcribed_to", version=1, ontology=
            graph_ontology)[0]

        variant_of = Term.objects.get_or_create(identifier="variant_of", name="variant_of", version=1, ontology=
        graph_ontology)[0]

        transcribed_from = \
            Term.objects.get_or_create(identifier="transcribed_from", name="transcribed_from", version=1, ontology=
            graph_ontology)[0]

        adjacent_to = Term.objects.get_or_create(identifier="adjacent_to", name="adjacent_to", version=1, ontology=
        graph_ontology)[0]

        member_of = Term.objects.get_or_create(identifier="member_of", name="member_of", version=1, ontology=
        graph_ontology)[0]

        contains = Term.objects.get_or_create(identifier="contains", name="contains", version=1, ontology=
        graph_ontology)[0]
        non_functional_homolog_of = Term.objects.get_or_create(identifier="non_functional_homolog_of",
                                                               name="non_functional_homolog_of", version=1, ontology=
                                                               graph_ontology)[0]

        overlaps = Term.objects.get_or_create(identifier="overlaps", name="overlaps", version=1, ontology=
        graph_ontology)[0]

        guided_by = Term.objects.get_or_create(identifier="guided_by", name="guided_by", version=1, ontology=
        graph_ontology)[0]

        Ontology.relmap = {
            'negatively_regulates': negatively_regulates, 'regulates': regulates,
            'positively_regulates': positively_regulates, 'part_of': part_of,
            "has_quality": has_quality, 'derives_from': derives_from, 'has_origin': has_origin,
            "has_part": has_part, "transcribed_to": transcribed_to, "variant_of": variant_of,
            "transcribed_from": transcribed_from, "adjacent_to": adjacent_to,
            "member_of": member_of, "contains": contains, "non_functional_homolog_of": non_functional_homolog_of,
            "overlaps": overlaps, "guided_by": guided_by,"is_a":is_a
        }