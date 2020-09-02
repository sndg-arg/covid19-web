import os

from django.db.utils import IntegrityError
from tqdm import tqdm

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "sndg.settings")
import django

django.setup()

from goatools.obo_parser import GODag
from biosql.models import Ontology, Term, TermRelationship, Dbxref

ontology = Ontology.objects.get_or_create(name="SeqFeature Sources")[0]
Term.objects.get_or_create(identifier="manual",name="manual",version=1,ontology=ontology,
                           definition="added or corrected by a person")
Term.objects.get_or_create(identifier="bibliography",name="bibliography",version=1,ontology=ontology,
                           definition="found in bibliography")
Term.objects.get_or_create(identifier="experimental",name="experimental",version=1,ontology=ontology,
                           definition="the annotation was obtained experimentally")
Term.objects.get_or_create(identifier="calculated",name="calculated",version=1,ontology=ontology,
                           definition="the annotation was obtained using software")
Term.objects.get_or_create(identifier="other",name="other",version=1,ontology=ontology,
                           definition="")

ontology = Ontology.objects.get_or_create(name="Gene Ontology")
cache = {}

dbmap = {
    "biological_process": Dbxref.objects.get_or_create(dbname="go", accession="biological_process", version=1)[0],
    "molecular_function": Dbxref.objects.get_or_create(dbname="go", accession="molecular_function", version=1)[0],
    "cellular_component": Dbxref.objects.get_or_create(dbname="go", accession="cellular_component", version=1)[0],
}

# for x in """goslim_agr "AGR slim"
# goslim_aspergillus "Aspergillus GO slim"
# goslim_candida "Candida GO slim"
# goslim_chembl "ChEMBL protein targets summary"
# goslim_generic "Generic GO slim"
# goslim_goa "GOA and proteome slim"
# goslim_metagenomics "Metagenomics GO slim"
# goslim_mouse "Mouse GO slim"
# goslim_pir "PIR GO slim"
# goslim_plant "Plant GO slim"
# goslim_pombe "Fission yeast GO slim"
# goslim_synapse "synapse GO slim"
# goslim_virus "Viral GO slim"
# goslim_yeast "Yeast GO slim"
# gosubset_prok "Prokaryotic GO subset"
# virus_checked "Viral overhaul terms" """.split("\n"):
#     k = x.strip().split(' ')[0]
#     dbmap[k] = Dbxref.objects.get_or_create(dbname="go", accession=k, version=1)[0]
#
# go_dag = GODag("/data/databases/go/go.obo", load_obsolete=True,optional_attrs=["def","synonym","subset","alt_id"])
# for go, term in tqdm(go_dag.items()):
#     if Term.objects.filter(ontology=ontology,identifier=go).exists():
#         t = Term.objects.prefetch_related("dbxrefs__dbxref").get(ontology=ontology,identifier=go)
#         intersetc = (set([ x.dbxref.accession for x in  t.dbxrefs.all()]) &
#                      set(["biological_process","molecular_function","cellular_component"]))
#         if not intersetc:
#             print (go)
#             termdbref = TermDbxref(term=t, dbxref=dbmap[term.namespace], rank=1)
#             termdbref.save()

# for go, term in tqdm(go_dag.items()):
#
#     with transaction.atomic():
#         dbTerm = Term(name=term.name,
#                       definition=term.defn,
#                       identifier=go,
#                       is_obsolete= "T" if term.is_obsolete else "F",
#                       ontology=ontology)
#
#         try:
#             dbTerm.save()
#             termdbref = TermDbxref(term=dbTerm, dbxref=dbmap[term.namespace], rank=1)
#             termdbref.save()
#
#
#             for subset in term.subset:
#                 if subset in dbmap:
#                     termdbref = TermDbxref(term=dbTerm, dbxref=dbmap[subset], rank=1)
#                     termdbref.save()
#             for synonym in term.synonym:
#                 syn = TermSynonym(term=dbTerm, synonym=synonym[0][:255])
#                 syn.save()
#
#
#             cache[go] = dbTerm
#         except IntegrityError as ex:
#             print (ex)
#         # dbTerm = Term.objects.get(ontology=ontology,identifier=go)
#
#


# go_dag = GODag("/data/databases/go/go-basic.obo", optional_attrs=['relationship'], load_obsolete=False)
# cache = {x.identifier:x for x in Term.objects.filter(ontology=ontology)}
is_a = Term.objects.get(identifier="is_a")
part_of = Term.objects.get(identifier="part_of")
regulates = Term.objects.get(identifier="regulates")
negatively_regulates = Term.objects.get(identifier="negatively_regulates")
positively_regulates = Term.objects.get(identifier="positively_regulates")
has_quality = Term.objects.get_or_create(identifier="has_quality",name="has_quality",version=1,ontology=
                                         Ontology.objects.get(name="Graph"))[0]
derives_from = Term.objects.get_or_create(identifier="derives_from",name="derives_from",version=1,ontology=
Ontology.objects.get(name="Graph"))[0]
has_origin = Term.objects.get_or_create(identifier="has_origin",name="has_origin",version=1,ontology=
Ontology.objects.get(name="Graph"))[0]
has_part = Term.objects.get_or_create(identifier="has_part",name="has_part",version=1,ontology=
Ontology.objects.get(name="Graph"))[0]
transcribed_to = Term.objects.get_or_create(identifier="transcribed_to",name="transcribed_to",version=1,ontology=
Ontology.objects.get(name="Graph"))[0]

variant_of = Term.objects.get_or_create(identifier="variant_of",name="variant_of",version=1,ontology=
Ontology.objects.get(name="Graph"))[0]

transcribed_from = Term.objects.get_or_create(identifier="transcribed_from",name="transcribed_from",version=1,ontology=
Ontology.objects.get(name="Graph"))[0]

adjacent_to = Term.objects.get_or_create(identifier="adjacent_to",name="adjacent_to",version=1,ontology=
Ontology.objects.get(name="Graph"))[0]

member_of = Term.objects.get_or_create(identifier="member_of",name="member_of",version=1,ontology=
Ontology.objects.get(name="Graph"))[0]

contains = Term.objects.get_or_create(identifier="contains",name="contains",version=1,ontology=
Ontology.objects.get(name="Graph"))[0]
non_functional_homolog_of = Term.objects.get_or_create(identifier="non_functional_homolog_of",
                                                       name="non_functional_homolog_of",version=1,ontology=
Ontology.objects.get(name="Graph"))[0]

overlaps = Term.objects.get_or_create(identifier="overlaps",name="overlaps",version=1,ontology=
Ontology.objects.get(name="Graph"))[0]

guided_by = Term.objects.get_or_create(identifier="guided_by",name="guided_by",version=1,ontology=
Ontology.objects.get(name="Graph"))[0]


relmap = {
    'negatively_regulates': negatively_regulates, 'regulates': regulates,
    'positively_regulates': positively_regulates, 'part_of': part_of,
    "has_quality":has_quality,'derives_from':derives_from,'has_origin':has_origin,
    "has_part":has_part,"transcribed_to":transcribed_to,"variant_of":variant_of,
    "transcribed_from":transcribed_from,"adjacent_to":adjacent_to,
    "member_of":member_of,"contains":contains,"non_functional_homolog_of":non_functional_homolog_of,
    "overlaps":overlaps,"guided_by":guided_by
}
#
# for go, term in tqdm(go_dag.items()):
#     try:
#         with transaction.atomic():
#             for child in term.children:
#                 if go in cache and child.id in cache:
#                     r = TermRelationship(
#                         subject_term=cache[go],  # parent
#                         predicate_term=is_a,
#                         object_term=cache[child.id],  # child
#                         ontology=ontology
#                     )
#                     r.save()
#
#             for rel, terms in term.relationship.items():
#                 for child in terms:
#                     if go in cache and child.id in cache:
#                         r = TermRelationship(
#                             subject_term=cache[go],  # parent
#                             predicate_term=relmap[rel],
#                             object_term=cache[child.id],  # child
#                             ontology=ontology
#                         )
#                         r.save()
#     except IntegrityError as ex:
#         print(ex)
#         with transaction.atomic():
#             for child in term.children:
#                 if go in cache and child.id in cache:
#                     r = TermRelationship(
#                         subject_term=cache[go],  # parent
#                         predicate_term=is_a,
#                         object_term=cache[child.id],  # child
#                         ontology=ontology
#                     )
#                     r.save()


# xref


# pfam_file_path = "/data/databases/xfam/Pfam-A.hmm"
#
# term = ""
# name = ""
# ident = ""
# ver = 1
# description = ""
# is_obsolete = "F"
# ontology = Ontology.objects.get_or_create(name="Pfam")[0]
#
# with open(pfam_file_path) as pfam_handle:
#     for line in tqdm(pfam_handle):
#         if line.strip().startswith("DESC"):
#             description = line.split("DESC")[1].strip()
#
#         elif line.strip().startswith("ACC"):
#             term = line.split("ACC")[1].strip().lower()
#             ver = int(term.split(".")[1])
#             term = term.split(".")[0]
#         elif line.strip().startswith("NAME"):
#             if term:
#                 dbTerm = Term(name=name, definition=description, identifier=term,
#                               is_obsolete=is_obsolete, ontology=ontology,version=ver).save()
#                 term = ""
#                 ident = ""
#                 description = ""
#
#             name = line.split("NAME")[1].strip()


dbmap = {
    "biosapiens": Dbxref.objects.get_or_create(dbname="so", accession="biosapiens", version=1)[0],
    "DBVAR": Dbxref.objects.get_or_create(dbname="so", accession="DBVAR", version=1)[0],
    "SOFA": Dbxref.objects.get_or_create(dbname="so", accession="SOFA", version=1)[0],
}

ontology = Ontology.objects.get_or_create(name="Sequence Ontology")[0]

# go_dag = GODag("/data/databases/so/so.obo", load_obsolete=True,optional_attrs=["def","synonym","subset"])
# names = {}
# for go, term in tqdm(go_dag.items()):
#     with transaction.atomic():
#
#         if term.name in names:
#             continue
#
#         dbTerm = Term.objects.get(name=term.name,  ontology=ontology)
#
#         names[term.name] = 1
#
#         for subset in term.subset:
#             if subset in dbmap:
#                 termdbref = TermDbxref(term=dbTerm, dbxref=dbmap[subset], rank=1)
#                 termdbref.save()
#
#         for synonym in term.synonym:
#             try:
#                 syn = TermSynonym(term=dbTerm, synonym=synonym[0][:255])
#                 syn.save()
#             except:
#                 pass
#
#         for synonym in term.alt_ids:
#             try:
#                 syn = TermSynonym(term=dbTerm, synonym=synonym[0][:255])
#                 syn.save()
#             except:
#                 pass

names = {}
cache = {x.identifier:x for x in Term.objects.filter(ontology=ontology)}
go_dag = GODag("/data/databases/so/so.obo", optional_attrs=['relationship'], load_obsolete=False)
for go, term in tqdm(go_dag.items()):

    if term.name in names:
        continue


    names[term.name] = 1


    for child in term.children:
        if go in cache and child.id in cache:
            r = TermRelationship(
                subject_term=cache[go],  # parent
                predicate_term=is_a,
                object_term=cache[child.id],  # child
                ontology=ontology
            )
            try:
                r.save()
            except IntegrityError:
                pass

    for rel, terms in term.relationship.items():
        for child in terms:
            if go in cache and child.id in cache:
                r = TermRelationship(
                    subject_term=cache[go],  # parent
                    predicate_term=relmap[rel],
                    object_term=cache[child.id],  # child
                    ontology=ontology
                )
                try:
                    r.save()
                except IntegrityError:
                    pass
