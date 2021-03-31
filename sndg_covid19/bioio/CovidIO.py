import gzip
import io
import re
import zipfile
from datetime import datetime

import Bio.SeqIO as bpio
import sys
from Bio.Seq import Seq
from django.core import management
from tqdm import tqdm

from bioseq.models.Bioentry import Bioentry, BioentryQualifierValue, BioentryDbxref
from bioseq.models.Dbxref import Dbxref
from bioseq.models.Location import Location
from bioseq.models.Ontology import Ontology
from bioseq.models.Seqfeature import Seqfeature, SeqfeatureQualifierValue
from bioseq.models.Term import Term, TermDbxref
from sndg_covid19.bioio.JobValidationError import JobValidationError
from sndg_covid19.models import ImportJob


class CovidIO:
    valid_subdivitions = [x.strip() for x in """Buenos Aires
                          CABA
                          Catamarca
                          Chaco
                          Chubut
                          Córdoba
                          Corrientes
                          Entre Ríos
                          Formosa
                          Jujuy
                          La Pampa
                          La Rioja
                          Mendoza
                          Misiones
                          Neuquén
                          Río Negro
                          Salta
                          San Juan
                          San Luis
                          Santa Cruz
                          Santa Fe
                          Santiago del Estero
                          Tierra del Fuego
                          Tucumán
                          turismo""".split("\n")]

    @staticmethod
    def process_import_job(job: ImportJob):
        zipf = zipfile.ZipFile(job.fasta.path, 'r')

        with zipf.open(zipf.namelist()[0]) as hzip:
            h = io.TextIOWrapper(io.BytesIO(hzip.read()))
            try:
                seqs = bpio.to_dict(bpio.parse(h, "fasta"))
            except ValueError:
                raise JobValidationError("error en el formato del archivo fasta")
        if job.aln_type == "spike":
            ref_name = "S"

            result = CovidIO.validate_import_spike(seqs, job.csv.path, ref_name)

            if not result["errors"]:
                for r in seqs.values():
                    seqid = r.id
                    if seqid != ref_name:
                        rdata = result["data"][seqid]
                        r.id = f"Spike|hCoV-19/Argentina/{seqid}/{rdata['date'].split('-')[0]}|{rdata['date']}|-|Original|hCoV-19"
                        r.description = rdata["geo"]
                with gzip.open(job.fasta.path.replace(".zip", "_fixed.fasta.gz"), 'wt') as h:
                    bpio.write(seqs.values(), h, "fasta")

                job.status_desc = f"procesando {len(seqs)} secuencias"
                job.status = "processing"
                job.save()
                CovidIO.process_import_spike(job, seqs, ref_name)
            else:
                raise JobValidationError("\n".join(result["errors"]), result["errores_detallados"])
        else:
            ref_name = "MN908947.3"
            management.call_command("process_covid_msa",
                                    reference=ref_name, input_msa=h,
                                    # outdir_msa=f'{settings.ROOT_DIR}/{alnjobid}.msa',
                                    precompute_graphics=False, override=True, remove=False)

    @staticmethod
    def validate_import_spike(seqs, csv, ref_name):
        result = {"errors": [], "data": {}}
        errores_detallados = []
        result["errores_detallados"] = errores_detallados
        if not len(seqs):
            result["errors"].append("no se detecto niguna secuencia en el fasta")
        if len(set([len(seq) for seq in seqs.values()])) != 1:
            result["errors"].append("no todas las secuencias son del mismo largo")
        error_ids = []
        ids_validos = []
        if ref_name not in seqs:
            result["errors"].append(f"la referencia {ref_name} no esta en el alineamiento")
        for seq in seqs.values():
            if re.match("^PAIS-S-[A-Z]\d{4}$", seq.id):
                ids_validos.append(seq.id)
            elif ref_name == seq.id:
                pass
            else:
                errores_detallados.append(f"fasta error id: {seq.id}")
                error_ids.append(seq.id)
        with open(csv) as h:
            lines = h.readlines()
        error_csv_sin_campos = []
        error_csv_codigo = []
        error_csv_fecha = []
        error_csv_geo = []
        csv_valid_cod = []
        for idx, l in enumerate(lines[1:], start=2):
            if not l.strip():
                continue
            vec = l.split(",")
            if len(vec) < 3:
                error_csv_sin_campos.append(str(idx))
                result["errores_detallados"].append(f"la linea {idx} no tiene 3 campos")
            else:
                cod, sample_date, geo = [x.strip() for x in vec[:3]]
                result["data"][cod] = {"geo": geo, "date": sample_date}
                if not re.match("^PAIS-S-[A-Z]\d{4}$", cod):
                    error_csv_codigo.append(f'{idx} {cod}')
                else:
                    csv_valid_cod.append(cod)
                try:
                    datetime.strptime(sample_date, '%Y-%m-%d')
                except ValueError:
                    error_csv_fecha.append(f'{idx} {sample_date}')
                if geo not in CovidIO.valid_subdivitions:
                    error_csv_geo.append(f'{idx} {geo}')
        if error_csv_sin_campos:
            result["errors"].append(f"hay {len(error_csv_sin_campos)} lineas en el csv que no tienen 3 campos")
            for x in error_csv_sin_campos:
                errores_detallados.append(f"error_csv_sin_campos: {x}")
        if error_csv_codigo:
            result["errors"].append(f"hay {len(error_csv_codigo)} lineas que tienen mal el código de muestra")
            for x in error_csv_codigo:
                errores_detallados.append(f"error_csv_codigo: {x}")
        if error_csv_fecha:
            result["errors"].append(f"hay {len(error_csv_fecha)} lineas que tienen mal la fecha")
            for x in error_csv_fecha:
                errores_detallados.append(f"error_csv_fecha: {x}")
        if error_csv_geo:
            result["errors"].append(f"hay {len(error_csv_geo)} lineas que tienen mal la ubicación")
            for x in error_csv_geo:
                errores_detallados.append(f"error_csv_geo: {x}")




        diff = set(ids_validos) - set(csv_valid_cod)
        if diff:
            msg = "hay códigos fasta que estan en el fasta y no estan en el csv"
            result["errors"].append(msg)
            result["errores_detallados"].append(msg + ": " + ",".join(diff))
        diff = set(csv_valid_cod) - set(ids_validos)
        if diff:
            msg = "hay códigos fasta que estan en el csv y que no estan en el fasta"
            result["errors"].append(msg)
            result["errores_detallados"].append(msg + ": " + ",".join(diff))

        return result

    @staticmethod
    def process_import_spike(job, seqs, ref_name):

        aln_fna_path = job.fasta.path.replace(".zip", ".fna")
        aln_faa_path = job.fasta.path.replace(".zip", ".faa")
        with open(aln_fna_path, "w") as h_fna, open(aln_faa_path, "w") as h_faa:
            ref = seqs[ref_name]
            for sseq in seqs.values():
                if sseq.id != ref_name:
                    seq = []
                    for i, nucleotide in enumerate(str(sseq.seq)):
                        seq.append(ref[i] if nucleotide == "-" else nucleotide)
                    seq = "".join(seq)
                    seq_id = sseq.id
                    sseq.seq = Seq(seq)
                    sseq.id = seq_id
                    sseq.description = sseq.description
                    bpio.write(sseq, h_fna, "fasta")
                    sseq.seq = sseq.seq.translate()
                    bpio.write(sseq, h_faa, "fasta")
                else:
                    bpio.write(sseq, h_fna, "fasta")
                    sseq = sseq.translate()
                    sseq.id = ref_name
                    bpio.write(sseq, h_faa, "fasta")
        management.call_command("process_spike_msa",
                                reference=ref_name, input_msa=aln_faa_path,
                                # outdir_msa=f'{settings.ROOT_DIR}/{alnjobid}.msa',
                                precompute_graphics=False)

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
                    assert int(feature.location.start) < prot.seq.length, [int(feature.location.start), prot.seq.length]
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
                    identifier = ":".join(dbx.split(":")[1:])

                    if any([dbx.startswith(x) for x in ["GO:", "EC:"]]):
                        identifier = dbx

                        ontology = Ontology.objects.get(name=Ontology.GO if db == "GO" else Ontology.EC)

                        if not Term.objects.filter(ontology=ontology,
                                                   identifier=identifier).exists():  ## TODO check/update version
                            props = dbx_dict[dbx] if dbx in dbx_dict else {}
                            t = Term(ontology=ontology, identifier=identifier, name="")
                            for field in ["GoTerm", "EntryName", "GeneName"]:
                                t.name = props[field] if field in props else t.name

                            if db == "GO":
                                Bioentry.objects.get(qualifiers__term__ontology__name=Ontology.GO)
                                dbxref = Dbxref.objects.get(dbname="go", accession=props["database"], version=1)
                                TermDbxref(term=t, dbxref=dbxref).save()
                            t.definition = props["Description"] if "Description" in props else ""
                            t.save()

                        term = Term.objects.get(ontology=ontology, identifier=identifier)

                        BioentryQualifierValue.objects.get_or_create(bioentry=prot, term=term, value=identifier)

                        continue

                    if db == "UnipName":
                        BioentryQualifierValue.objects.get_or_create(bioentry=prot, term=Term.objects.get(
                            ontology__name=Ontology.ANNTAGS, identifier="product"), value=identifier)
                    if db == "UnipGene":
                        BioentryQualifierValue.objects.get_or_create(bioentry=prot, term=Term.objects.get(
                            ontology__name=Ontology.ANNTAGS, identifier="gene"), value=identifier)

                    ver = 1
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
