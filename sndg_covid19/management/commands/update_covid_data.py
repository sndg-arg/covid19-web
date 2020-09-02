import os
import subprocess as sp
from tqdm import tqdm
import requests
from datetime import datetime
from io import StringIO
import Bio.SeqIO as bpio

from django.core.management.base import BaseCommand,CommandError

from bioseq.models.Biodatabase import Biodatabase
from sndg_covid19.io.CovidIO import CovidIO
from bioseq.io.BioIO import BioIO

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


class Command(BaseCommand):
    """

    """

    DEFAULT_COVID_DATA_URL = "https://www.ebi.ac.uk/uniprot/api/covid-19/uniprotkb/stream?compressed=true&download=true&format=json&query=%2A"
    DEFAULT_COVID_FASTA_URL = "https://www.ebi.ac.uk/uniprot/api/covid-19/uniprotkb/stream?compressed=true&download=true&format=fasta&query=%2A"
    DEFAULT_COVID_FASTA = "data/tmp/covid19.fasta"
    DEFAULT_UNIP_COVID_FASTA = "data/tmp/unip_covid19.fasta"
    DEFAULT_COVID_JSON = "data/tmp/unip_covid19.json"
    DEFAULT_BLAST_RESULT = "data/tmp/blast_unip_genome.tbl"
    DEFAULT_UNIP_BIG_PROT = "P0DTD1"

    help = 'Tinny DB for documentation and testing purposes'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dbx_dict = {}

    def add_arguments(self, parser):
        parser.add_argument('--unip_data_url',
                            default=os.environ.get("DEFAULT_COVID_DATA_URL", Command.DEFAULT_COVID_DATA_URL))
        parser.add_argument('--unip_fasta_url',
                            default=os.environ.get("DEFAULT_COVID_FASTA_URL", Command.DEFAULT_COVID_FASTA_URL))
        parser.add_argument('--tmp_db_fasta', default=os.environ.get("COVID_FASTA", Command.DEFAULT_COVID_FASTA))
        parser.add_argument('--tmp_db_json', default=os.environ.get("COVID_JSON", Command.DEFAULT_COVID_JSON))
        parser.add_argument('--tmp_blast_result', default=os.environ.get("COVID_JSON", Command.DEFAULT_BLAST_RESULT))

    def json2seqrecord(self, json_record):
        uid = json_record["primaryAccession"]
        if "recommendedName" in json_record["proteinDescription"]:
            desc = json_record["proteinDescription"]["recommendedName"]["fullName"]["value"]
        else:
            desc = json_record["proteinDescription"]["submissionNames"][0]["fullName"]["value"]
        ecs = ([x["value"] for x in json_record["proteinDescription"]["recommendedName"]["ecNumbers"]]
               if ("recommendedName" in json_record["proteinDescription"] and
                   "ecNumbers" in json_record["proteinDescription"]["recommendedName"]) else [])

        if "contains" in json_record["proteinDescription"]:
            for pd in json_record["proteinDescription"]["contains"]:

                if "ecNumbers" in pd["recommendedName"]:
                    ecs += [x["value"] for x in pd["recommendedName"]["ecNumbers"]]

        r = SeqRecord(id=uid, name="", description=desc, seq=Seq(json_record["sequence"]["value"]))

        if "genes" in json_record:
            for gene in json_record["genes"]:
                if "geneName" in gene:
                    val = gene["geneName"]["value"]
                    dbx = "UnipGene:" + val
                    r.dbxrefs.append(dbx)
                if "synonyms" in json_record["genes"]:
                    for syn in json_record["genes"]["synonyms"]:
                        val = syn["value"]
                        dbx = "UnipGene:" + val
                        r.dbxrefs.append(dbx)

        if "alternativeNames" in json_record["proteinDescription"]:
            for an in json_record["proteinDescription"]["alternativeNames"]:
                if "shortNames" in an:
                    for x in an["shortNames"]:
                        dbx = "UnipName:" + x["value"]
                        r.dbxrefs.append(dbx)

        for ref in json_record["uniProtKBCrossReferences"]:
            dbx = ref["database"] + ":" + ref["id"] if (ref["database"] + ":") not in ref["id"] else ref["id"]
            r.dbxrefs.append(dbx)
            self.dbx_dict[dbx] = {x["key"]: x["value"] for x in ref["properties"]}
            if ref["database"] == "GO":
                gt = self.dbx_dict[dbx]["GoTerm"]
                self.dbx_dict[dbx]["GoTerm"] = ":".join(gt.split(":")[1:])

                self.dbx_dict[dbx]["database"] = gt.split(":")[0]
                if self.dbx_dict[dbx]["database"] == "b":
                    self.dbx_dict[dbx]["database"] = "biological_process"
                if self.dbx_dict[dbx]["database"] == "c":
                    self.dbx_dict[dbx]["database"] = "molecular_function"
                if self.dbx_dict[dbx]["database"] == "f":
                    self.dbx_dict[dbx]["database"] = "cellular_component"

        for ref in (json_record["secondaryAccessions"] if "secondaryAccessions" in json_record else []
                   ) + [json_record["uniProtkbId"], json_record["primaryAccession"]]:
            dbx = "UnipAcc:" + ref.replace(" ", "_")
            r.dbxrefs.append(dbx)

        for ref in ecs:
            dbx = "EC:" + ref
            r.dbxrefs.append(dbx)
            self.dbx_dict[dbx] = [["description", desc]]

        r.dbxrefs = list(set(r.dbxrefs))

        if "features" in json_record:
            for f in json_record["features"]:
                l = FeatureLocation(start=f["location"]["start"]["value"],
                                    end=f["location"]["end"]["value"])

                qual = {"description": f["description"]} if f["description"].replace("-", "").strip() else {}
                if "featureId" in f:
                    qual["featureId"] = f["featureId"]
                seqf = SeqFeature(type=f["type"], location=l, qualifiers=qual)
                r.features.append(seqf)
        return r

    def process_big_prot(self, proteins, prot_id):
        seqs_extra = []
        bigprot = proteins[prot_id]
        for f in bigprot.features:
            if f.type == "chain" and len(f.location) < 5000:
                seq = str(f.extract(bigprot.seq))
                r = SeqRecord(id=prot_id + "_" + f.qualifiers["featureId"], name="",
                              description=f.qualifiers["description"], seq=Seq(seq))
                seqs_extra.append(r)
                for f2 in bigprot.features:

                    if f2.location.start >= f.location.start and f2.location.end <= f.location.end:
                        if f2.type == "chain":
                            r.dbxrefs.append("InterPro:" + f.qualifiers["featureId"])
                        else:
                            f3 = SeqFeature(location=FeatureLocation(start=f2.location.start - f.location.start,
                                                                     end=f2.location.end - f.location.start),
                                            type=f2.type, qualifiers=f2.qualifiers)
                            r.features.append(f3)
                proteins[r.id] = r
        return seqs_extra

    def handle(self, *args, **options):
        data = requests.get(options["unip_data_url"])
        if not data.ok:
            raise CommandError(data.text)

        protein_data = data.json()["results"]
        with open(options["tmp_db_json"], "w") as h:
            import json
            json.dump(protein_data, h)
        proteins = bpio.to_dict([self.json2seqrecord(x) for x in protein_data])

        seqs = self.process_big_prot(proteins, Command.DEFAULT_UNIP_BIG_PROT)

        fasta = options["tmp_db_fasta"]
        unip = Command.DEFAULT_UNIP_COVID_FASTA

        BioIO.proteome_fasta("COVID19", fasta)

        h = StringIO(requests.get(options["unip_fasta_url"]).text)
        bpio.write(list(bpio.parse(h, "fasta")) + seqs, unip, "fasta")

        sp.call(f'makeblastdb -dbtype prot -in {unip}', shell=True)

        cmd = f'blastp -db {unip} -query {fasta} -evalue 1e-6 -outfmt "6 qseqid sseqid pident slen length qcovs"  2>/dev/null'
        blast_result = sp.getoutput(cmd)

        with open(options["tmp_blast_result"], "w") as h:
            h.write(blast_result)

        prot_dict = {}
        for x in blast_result.split("\n"):
            if x.strip():
                db_id, unip_id, ident, slen, length, qcovs = x.strip().split()
                if ((int(length) * 1.0 / int(slen)) > 0.9) and (int(qcovs) > 90):
                    unip_id = unip_id.split("|")[1] if "PRO_" not in  unip_id else unip_id
                    prot_dict[int(db_id)] = proteins[unip_id]


        bdb = Biodatabase.objects.get(name="COVID19" + Biodatabase.PROT_POSTFIX)

        CovidIO.load_unip_data(prot_dict, bdb, self.dbx_dict)
