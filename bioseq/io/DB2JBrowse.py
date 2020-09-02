# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import os
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO as bpio

from bioseq.models.Bioentry import Bioentry
from bioseq.models.Location import Location
from bioseq.models.Seqfeature import Seqfeature

from bioseq.models.Biodatabase import Biodatabase
from collections import defaultdict

import subprocess as sp
from tqdm import tqdm
import json

"""
    Template to add click to trackList.json

    "onClick" : {
          "url" : "function(evt) {  debugger;window.parent.location.href =  window.parent.location.href.split('genome/')[0] + '/genome/BR_Mal/gene/' + this.feature.get('locus_tag') + '/'; }",
          "label" : "go to product"
        }
"""

class DB2JBrowse():

    def __init__(self, jbrowse_path, jbrowse_data_path=None, tmp="/tmp/"):
        self.name = ""

        self.fasta_path = tmp + "/jbrowse.fasta"
        self.gff_path = tmp + "/jbrowse.gff"
        self.ovewrite = True
        self.source = "."
        self.jbrowse_path = jbrowse_path
        self.jbrowse_data_path = jbrowse_data_path if jbrowse_data_path else (jbrowse_path + "/data/")
        self.accession = None
        self.excluded = ["source"]
        self.stderr = sys.stderr
        self.previous_lt = ""
        self.previous_gene = ""

    def db2fs(self, bdb_name: str):
        self.accession = bdb_name
        bdb = Biodatabase.objects.get(name=self.accession)
        self.load(bdb.entries.all())  # "features__qualifiers", "features__type_term"

    def load(self, bioentry_iterator):
        assert os.path.exists(os.path.dirname(self.fasta_path))
        assert os.path.exists(os.path.dirname(self.gff_path))

        if not os.path.exists(os.path.dirname(self.fasta_path)) or self.ovewrite:
            with open(self.fasta_path, "w") as h_fasta, open(self.gff_path, "w") as h_gff:
                for bioentry in tqdm(bioentry_iterator, file=self.stderr, total=bioentry_iterator.count()):
                    self.process_bioentry(bioentry, h_fasta, h_gff)
        self.create_files()

    def process_bioentry(self, bioentry: Bioentry, h_fasta, h_gff):
        r = SeqRecord(id=bioentry.accession, name=bioentry.name,
                      seq=Seq(bioentry.seq.seq))
        if bioentry.description:
            r.description = bioentry.description,
        bpio.write(r, h_fasta, "fasta")

        qs = bioentry.features.prefetch_related("qualifiers", "type_term").exclude(
            type_term__identifier__in=self.excluded)

        for feature in tqdm(qs, file=self.stderr, total=qs.count()):
            if feature.type_term.identifier in ["CDS", "tRNA", "rRNA", "miRNA", "soRNA", "sRNA", "miRNA",
                                                "repeat_region", "mat_peptide"]:
                # ,"3'UTR","5'UTR","stem_loop"
                self.process_seqfeature(bioentry.accession, feature, h_gff)

    def process_seqfeature(self, contig_accession, feature: Seqfeature, h_gff):
        loc = sorted(feature.locations.all(), key=lambda x: x.start_pos)
        start = loc[0].start_pos
        end = loc[-1].end_pos
        strand = loc[0].strand

        fid = str(feature.seqfeature_id)
        attributes = ["ID=" + str(fid), "BioentryId="+str(feature.qualifiers_dict()['BioentryId'])]

        qualifiers = {x.term.identifier: x.value for x in feature.qualifiers.all()}
        desc = ""
        if ("note" in qualifiers) and (self.previous_gene != qualifiers["note"]):
            desc = qualifiers["note"]
        elif ("product" in qualifiers) and (self.previous_gene != qualifiers["product"]):
            desc = qualifiers["product"]

        name = ""
        if ("gene" in qualifiers) and (self.previous_gene != qualifiers["gene"]):
            gene = qualifiers["gene"]
            name = gene
        else:
            gene = str(feature.seqfeature_id)

        attributes.append("gene=" + gene)
        self.previous_gene = gene

        lt = ""
        if  ("locus_tag" in qualifiers) and qualifiers["locus_tag"] and (self.previous_lt != qualifiers["locus_tag"]):
            lt = qualifiers["locus_tag"]
            if not name:
                name = lt
            self.previous_lt = qualifiers["locus_tag"]
        else:
            lt = gene

        attributes.append("locus_tag=" + lt)
        attributes.append("Name=" + name)
        if desc and lt != desc:
            attributes.append("description=" + desc)

        attributes_str = ";".join([x.replace(";", "_") for x in attributes])
        r = "\t".join(
            [str(x) for x in [contig_accession, self.source, feature.type_term.identifier, start, end,
                              ".", strand, 0, attributes_str]])
        h_gff.write(r + "\n")

        # if len(loc) > 1:
        #     for sf in loc:
        #         self.process_location(accession, sf, h_gff, fid)

    def process_location(self, bioentry_accession: str, location: Location, h_gff, parent: str):

        attributes = ["Parent=" + parent]
        attributes_str = ";".join(attributes)
        r = "\t".join(
            [str(x) for x in [bioentry_accession, self.source, "exon", location.start_pos, location.end_pos,
                              ".", location.strand, 0, attributes_str]])
        h_gff.write(r + "\n")
        # if sfs:
        #     for sf in sfs:
        #         self.process_subfeature(bioentry_accession, sf, h_gff, fid)

    # def process_subfeature(self, bioentry_accession: str, feature: Seqfeature, h_gff, parent):
    #     loc = feature.locations.all()[0]
    #     sfs = feature.subfeatures()
    #     attributes = ["Parent=" + parent]
    #     fid = str(feature.seqfeature_id)
    #
    #     if "locus_tag" in feature.qualifiers:
    #         attributes.append("locus_tag=" + feature.qualifiers["locus_tag"][0])
    #
    #     attributes.append("BioentryId=" + feature.qualifiers["BioentryId"][0])
    #
    #     if sfs:
    #         attributes.append("ID=" + fid)
    #     attributes_str = ";".join(attributes)
    #     r = "\t".join(
    #         [str(x) for x in [bioentry_accession, self.source, feature.type_term.identifier, loc.start_pos, loc.end_pos,
    #                           ".", loc.strand, 0, attributes_str]])
    #     h_gff.write(r + "\n")
    #     if sfs:
    #         for sf in sfs:
    #             self.process_subfeature(bioentry_accession, sf, h_gff, fid)

    def create_files(self):
        """
        http://jbrowse.org/docs/html_features.html
        """
        sys.stderr.write("prepare-refseqs\n")
        cmd = """docker run -u $(id -u):$(id -g) -v {jbrowse_data_path}:/jbrowse/data/ -v {fasta_path}:/tmp/jbrowse.fasta \
        jbrowse/jbrowse-1.12.0 bin/prepare-refseqs.pl --fasta /tmp/jbrowse.fasta --out data/{acc} --key "Sequence"
        """.format(jbrowse_path=self.jbrowse_path, jbrowse_data_path=self.jbrowse_data_path,
                   acc=self.accession, fasta_path=self.fasta_path)
        sp.call(cmd, shell=True)

        sys.stderr.write("flatfile-to-json.pl\n")
        cmd = """docker run -u $(id -u):$(id -g) -v {jbrowse_data_path}:/jbrowse/data/ -v {gff_path}:/tmp/jbrowse.gff \
         jbrowse/jbrowse-1.12.0 ./bin/flatfile-to-json.pl --gff "/tmp/jbrowse.gff" --out "/jbrowse/data/{acc}" \
          --key "{label}" --trackLabel "{label}" --trackType CanvasFeatures --nameAttributes "gene,locus_tag" \
           --clientConfig '{clientConfig}'
        """.format(jbrowse_path=self.jbrowse_path, jbrowse_data_path=self.jbrowse_data_path,
                   acc=self.accession, gff_path=self.gff_path, label="Genes",
                   clientConfig=json.dumps({"description": "description"}))

        # --subfeatureClasses '{subfeatureClasses}
        # subfeatureClasses=json.dumps(      {"cds": "cds", "trna": "feature1", "trna": "feature2",
        #      "rrna": "feature", "repeat_region": "feature5"             })
        # --className  feature
        # --trackType FeatureTrack
        sp.call(cmd, shell=True)

        sys.stderr.write("generate-names\n")
        cmd = """docker run -u $(id -u):$(id -g) -v {jbrowse_data_path}:/jbrowse/data/  \
         jbrowse/jbrowse-1.12.0 ./bin/generate-names.pl  --out "/jbrowse/data/{acc}" --tracks "{labels}"
        """.format(jbrowse_path=self.jbrowse_path, jbrowse_data_path=self.jbrowse_data_path,
                   acc=self.accession, labels=",".join(["Sequence", "Genes"]))  # --className  feature

        sp.call(cmd, shell=True)
