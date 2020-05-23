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



def bulk_save(iterator: Iterable, action: Callable[[BSeqFeature], Seqfeature], bulk_size: int = 1000,stderr=sys.stderr):
    data = []
    for obj in tqdm(iterator, file=stderr):
        data.append(obj)
        if (len(data) % bulk_size) == 0:
            with transaction.atomic():
                for x in data:
                    action(x)
            data = []
    if data:
        with transaction.atomic():
            for x in data:
                action(x)


class BioIO:
    included_feature_qualifiers = ["gene", "locus_tag", "db_xref","gene_synonym"]
    excluded_feature_qualifiers = ["translation", "gene", "locus_tag", "db_xref", "product", ]

    @staticmethod
    def proteome_fasta(genome_name,fp,dbname_postfix=Biodatabase.PROT_POSTFIX):
        if not hasattr(fp, 'write'):
            h = open(fp,"w")
        else:
            h = fp
        qs = Bioentry.objects.prefetch_related("seq").filter(
            biodatabase__name = genome_name + dbname_postfix)
        try:
            for be in qs:
                r = be.to_seq_record()
                desc = "|".join( [ str(x) if x  else "" for x in [be.accession, "__".join(be.genes()), be.description ]  ] )
                r.id = str(be.bioentry_id)
                r.description = desc
                r.name = ""
                bpio.write(r,h,"fasta")
        finally:
            if not hasattr(fp, 'write'):
                h.close()

    def __init__(self, biodb_name: str, ncbi_tax: int):
        self.biodb_name = biodb_name
        self.ncbi_tax = ncbi_tax
        self.stderr = sys.stderr

    def create_db(self):
        self.genomedb = Biodatabase(name=self.biodb_name)
        self.genomedb.save()
        self.proteindb = Biodatabase(name=self.biodb_name + Biodatabase.PROT_POSTFIX)
        self.proteindb.save()
        self.rnadb = Biodatabase(name=self.biodb_name + Biodatabase.RNAS_POSTFIX)
        self.rnadb.save()

        self.sfk_ontology = Ontology.objects.get(name=Ontology.SFK)
        self.sfs_ontology = Ontology.objects.get(name=Ontology.SFS)
        self.ann_ontology = Ontology.objects.get(name=Ontology.ANNTAGS)
        self.feature_term = Term.objects.get_or_create(ontology=self.ann_ontology, identifier="GFeatureId")[0]
        self.bioentry_term = Term.objects.get_or_create(ontology=self.ann_ontology, identifier="BioentryId")[0]



    def process_feature(self, be: Bioentry, feature: BSeqFeature):
        type_term = Term.objects.get_or_create(ontology=self.sfk_ontology, identifier=feature.type)[0]
        source_term = Term.objects.get_or_create(ontology=self.sfk_ontology, identifier="calculated")[0]

        display_name = feature.qualifiers["gene"][0] if "gene" in feature.qualifiers else (
            feature.qualifiers["locus_tag"][0] if "locus_tag" in feature.qualifiers else feature.type)
        locus_tag = feature.qualifiers["locus_tag"][0] if "locus_tag" in feature.qualifiers else ""
        product = feature.qualifiers["product"][0].replace(" ","_") if "product" in feature.qualifiers else ""
        if feature.type == "mat_peptide" and product:
            display_name = display_name + "_" + product
            locus_tag = locus_tag + "_" + product

        sf = Seqfeature(bioentry=be, type_term=type_term, source_term=source_term, display_name=display_name)
        sf.save()

        for key, value in feature.qualifiers.items():
            value = value[0]
            if not (feature.type in ["CDS" , "RNA"] ) or (key in BioIO.included_feature_qualifiers):
                term = Term.objects.get_or_create(identifier=key, ontology=self.ann_ontology)[0]
                if term.identifier == "locus_tag":
                    value = locus_tag
                sfqv = SeqfeatureQualifierValue.objects.create(seqfeature=sf, term=term, value=value)
                sfqv.save()
        # sub_features

        for rank, location in enumerate(feature.location.parts):
            loc = Location(seqfeature=sf, start_pos=location.start, end_pos=location.end,
                           strand=location.strand, rank=rank)
            loc.save()

        if "pseudo" in feature.qualifiers:
            return

        if feature.type in  ["CDS" , "RNA" ,"mat_peptide" ]:
            self.bioentry_from_feature(be, feature, sf)


    def bioentry_from_feature(self, be, feature, sf):
        description = feature.qualifiers["product"][0] if "product" else (
            feature.qualifiers["note"][0] if "note" in feature.qualifiers else "")

        gene = feature.qualifiers["gene"][0] if "gene" in feature.qualifiers else feature.qualifiers["locus_tag"][0]
        locus_tag = feature.qualifiers["locus_tag"][0] if "locus_tag" in feature.qualifiers else gene
        if feature.type == "mat_peptide":
            gene = gene + "_" + feature.qualifiers["product"][0]
            locus_tag = locus_tag + "_" + feature.qualifiers["product"][0].replace(" ","_")

        if feature.type == "CDS":
            seq = feature.qualifiers["translation"][0]
            db = self.proteindb
        elif feature.type == "mat_peptide":
            seq = str(feature.extract(Seq(be.seq.seq)).translate())

            db = self.proteindb
        else:
            db = self.rnadb
            seq = str(feature.extract(Seq(be.seq.seq)))
        be_count = Bioentry.objects.filter(biodatabase=db, accession=locus_tag).count()
        if be_count:
            be_count = Bioentry.objects.filter(biodatabase=db, accession__startswith=locus_tag + "_").count() + 1
            locus_tag = locus_tag + "_" + str(be_count)
        prot = Bioentry.objects.create(biodatabase=db,
                                       description=description, name=gene,identifier=locus_tag,
                                       accession=locus_tag)
        BioentryQualifierValue.objects.create(bioentry=prot, term=self.feature_term, value=sf.seqfeature_id)
        SeqfeatureQualifierValue.objects.create(seqfeature=sf, term=self.bioentry_term, value=prot.bioentry_id)
        Biosequence.objects.create(bioentry=prot, seq=seq, length=len(seq))
        # unip = Dbxref.objects.create(dbname="Uniprot", accession="P9WHW3")
        # BioentryDbxref.objects.create(bioentry=prot, dbxref=unip)
        for key, value in feature.qualifiers.items():
            value = value[0]
            if key not in BioIO.excluded_feature_qualifiers:
                term = Term.objects.get_or_create(identifier=key, ontology=self.ann_ontology)[0]
                BioentryQualifierValue.objects.create(bioentry=prot, term=term,
                                                      value=value)
        return prot

    def process_seqrecord(self, seqrecord):
        be = Bioentry(biodatabase=self.genomedb,
                      taxon=Taxon.objects.get_or_create(ncbi_taxon_id=self.ncbi_tax)[0],
                      name=seqrecord.name, accession=seqrecord.id,
                      identifier=seqrecord.id)
        be.save()

        seq = Biosequence(bioentry=be, seq=str(seqrecord.seq), length=len(seqrecord.seq))
        seq.save()

        bulk_save(seqrecord.features, action=lambda f: self.process_feature(be, f),stderr=self.stderr)

    def process_record_list(self, seq_record_iterator: Iterable[SeqRecord], contig_count: int):
        with tqdm(seq_record_iterator, total=contig_count, file=self.stderr) as pbar:
            for seqrecord in pbar:
                pbar.set_description(seqrecord.id)
                self.process_seqrecord(seqrecord)
