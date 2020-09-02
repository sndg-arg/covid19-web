import os
from tqdm import tqdm
from datetime import datetime
from typing import Iterable, Callable
import sys

import Bio.SeqIO as bpio
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from django.db import transaction

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature as BSeqFeature

from bioseq.models.Taxon import Taxon, TaxonName, TaxIdx
from bioseq.models.Ontology import Ontology
from bioseq.models.Term import Term, TermRelationship, TermDbxref, TermIdx
from bioseq.models.Biodatabase import Biodatabase
from bioseq.models.Bioentry import Bioentry, BioentryQualifierValue, BioentryDbxref
from bioseq.models.Biosequence import Biosequence
from bioseq.models.Seqfeature import Seqfeature, SeqfeatureDbxref, SeqfeatureQualifierValue
from bioseq.models.Location import Location
from bioseq.models.Dbxref import Dbxref, DbxrefQualifierValue


class CovidIO:
    @staticmethod
    def load_unip_data(protein_map, covid_biodb, dbx_dict):
        # for db_protein_id, unip_protein_id in protein_map
        for prot in tqdm(covid_biodb.entries.all()):

            if prot.bioentry_id not in protein_map:


                sys.stderr.write(f"{prot.name} without annotation in uniprot\n")
            else:
                unip_seqrecord = protein_map[prot.bioentry_id]
                prot.description = unip_seqrecord.description
                for feature in unip_seqrecord.features:
                    assert int(feature.location.start) < prot.seq.length,[int(feature.location.start) , prot.seq.length]
                    display_name = feature.qualifiers[
                        "description"] if "description" in feature.qualifiers else feature.type
                    source_term = Term.objects.get(ontology__name=Ontology.SFS, name="calculated")
                    ontology = Ontology.objects.get_or_create(name="UniprotFTypes")[0]
                    type_term = Term.objects.get_or_create(ontology=ontology, identifier=feature.type)[0]
                    sfqs = Seqfeature.objects.filter(bioentry=prot, source_term=source_term, type_term=type_term,
                                                     display_name=display_name[:64])


                    sf = Seqfeature(bioentry=prot, source_term=source_term, type_term=type_term,
                                    display_name=display_name[:64])
                    sf.save()
                    Location(start_pos=int(feature.location.start), end_pos=int(feature.location.end), seqfeature=sf,
                             strand=1).save()


                    if "featureId" in feature.qualifiers:
                        type_term = Term.objects.get_or_create(ontology=ontology, identifier="InterPro")[0]
                        SeqfeatureQualifierValue.objects.get_or_create(
                            seqfeature=sf, term=type_term, value=feature.qualifiers["featureId"])

                for dbx in unip_seqrecord.dbxrefs:

                    db = dbx.split(":")[0]
                    identifier = ":".join( dbx.split(":")[1:])

                    if any([dbx.startswith(x) for x in ["GO:", "EC:"]]):
                        identifier = dbx

                        ontology = Ontology.objects.get(name=Ontology.GO if db == "GO" else Ontology.EC)

                        if not Term.objects.filter(ontology=ontology,
                                                   identifier=identifier).exists():  ## TODO check/update version
                            props = dbx_dict[dbx] if dbx in dbx_dict else {}
                            t = Term(ontology=ontology, identifier=identifier, name="")
                            for field in ["GoTerm", "EntryName", "GeneName"]:
                                t.name = props[field] if field in props else t.name

                            if db == "GO" :
                                Bioentry.objects.get(qualifiers__term__ontology__name = Ontology.GO)
                                dbxref = Dbxref.objects.get(dbname="go", accession=props["database"], version=1)
                                TermDbxref(term=t,dbxref=dbxref).save()
                            t.definition = props["Description"] if "Description" in props else ""
                            t.save()

                        term = Term.objects.get(ontology=ontology, identifier=identifier)

                        BioentryQualifierValue.objects.get_or_create(bioentry=prot, term=term, value=identifier)


                        continue

                    if db == "UnipName":
                        BioentryQualifierValue.objects.get_or_create(bioentry=prot, term=Term.objects.get(
                            ontology__name=Ontology.ANNTAGS,identifier="product"), value=identifier)
                    if db == "UnipGene":
                        BioentryQualifierValue.objects.get_or_create(bioentry=prot, term=Term.objects.get(
                            ontology__name=Ontology.ANNTAGS,identifier="gene"), value=identifier)

                    ver =1
                    if "." in identifier and db not in ["ec"]:
                        ver = identifier.split(".")[-1]
                    try:
                        ver = int(ver)
                    except:
                        ver = 1


                    if Dbxref.objects.filter(dbname=db, accession=identifier).exists():
                        dbxs = Dbxref.objects.latest_version(dbname=db, accession=identifier)
                    else:
                        dbxs = Dbxref.objects.get_or_create(dbname=db, accession=identifier, version=ver)[0]

                    BioentryDbxref.objects.get_or_create(bioentry=prot, dbxref=dbxs)

                products = [x.value for x in prot.qualifiers.filter(term__identifier="product")]
                if products:
                    prot.name = products
                prot.save()
