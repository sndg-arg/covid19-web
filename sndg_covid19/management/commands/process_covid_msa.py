from datetime import datetime
import os

import Bio.SeqIO as bpio
from django.core.management.base import BaseCommand, CommandError
from tqdm import tqdm

from bioseq.io.MSAMap import MSAMap
from bioseq.models.Bioentry import Bioentry
from bioseq.models.Variant import Variant, SampleVariant, Sample
from config.settings.base import STATICFILES_DIRS
from sndg_covid19.io import country_from_gisaid
from sndg_covid19.tasks import variant_graphics
from glob import glob
import traceback
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import math


class Command(BaseCommand):
    """

    """

    help = 'Load variant list from msa'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dbx_dict = {}

    def add_arguments(self, parser):
        parser.add_argument("-r", '--reference', default="MN908947.3", help="ref sequence, default='MN908947.3'")
        # parser.add_argument('-a', '--accession', help="If input_msa is a file is required", required=False)
        parser.add_argument("-i", '--input_msa', required=True,
                            help="Could be a fasta msa file or a directory that contains a list of them")
        parser.add_argument("-o", '--outdir_msa', default="data/processed/", help="directory for protein MSAs")
        parser.add_argument("-pc", '--precompute_graphics', help="skip precompute", action="store_false")
        parser.add_argument("--override", help="override msa and graphics. Default: False", action="store_true")

    def create_aa_msa(self, genomic_msa, ref="MN908947.3", aa_msa_dir="data/processed/", override=False):
        covid = Bioentry.objects.get(biodatabase__name="COVID19")  # genome bioentry
        msas = []
        qs = covid.features.filter(type_term__identifier__in=["CDS", "mat_peptide"])
        for feature in tqdm(qs, file=self.stderr, total=qs.count()):
            gene = Bioentry.objects.get(bioentry_id=feature.qualifiers_dict()["BioentryId"]).accession
            msa_file = f'{aa_msa_dir}{gene}_msa.fasta'
            msas.append(msa_file)
            if override or (not os.path.exists(msa_file)) or (os.path.getsize(msa_file) < 100):
                with open(msa_file, "w") as h:
                    for sample in tqdm(genomic_msa.samples(), file=self.stderr):
                        ok = True
                        if sample != ref:

                            for offset in list(range(1000)) + [-x for x in range(1, 1000)]:
                                seq = str(genomic_msa.subseq(ref, feature.first_location().start_pos + offset,
                                                             feature.first_location().end_pos + offset, sample).seq)
                                if "-" in seq:
                                    gap_pos = seq.index("-")
                                    seq = seq[:gap_pos] + "".join(["N" for _ in seq[gap_pos:]])
                                seq = Seq(seq).translate()
                                if len(refseq) > len(seq):
                                    seq_org = seq
                                    seq = seq + Seq("".join(["X" for _ in range( len(refseq) - len(seq))]))

                                if len(refseq) != len(seq):
                                    raise AssertionError(f"protein reference and sample do not have the same length")

                                missmatches = sum([1 if x != y else 0 for x, y in zip(refseq, seq) if y != "X"])
                                ok = (missmatches / (len([1 for x in seq if x != "X"]) + 1)) < 0.2
                                if ok:
                                    break

                        else:
                            start = feature.first_location().start_pos
                            end = feature.first_location().end_pos
                            seq = Seq(covid.seq.seq[start:end].replace("-", "")).translate()

                        if sample == ref:
                            refseq = seq

                        if ok:
                            record = SeqRecord(id=sample, name="", description="", seq=seq)
                            bpio.write(record, h, "fasta")
                        else:
                            self.stderr.write(f'Invalid sequence "{gene}" for sample "{sample}"')
                            self.stderr.write(str(refseq))
                            self.stderr.write(str(seq))

    def handle(self, *args, **options):

        input_msa = options['input_msa']
        assert os.path.exists(input_msa), f'{input_msa} does not exists'
        ref_seq = options["reference"]
        if not os.path.exists(options['outdir_msa']):
            os.makedirs(options['outdir_msa'])
        assert os.path.exists(options['outdir_msa']), f'could not create {options["outdir_msa"]}'

        msa = MSAMap(bpio.to_dict(bpio.parse(input_msa, "fasta")))
        msa.init()

        self.create_aa_msa(msa, ref=ref_seq, aa_msa_dir=options['outdir_msa'], override=options["override"])

        expected_files = set(
            [x.accession + "_msa.fasta" for x in Bioentry.objects.filter(biodatabase__name="COVID19_prots")])
        files = (set(os.listdir(options['outdir_msa'])) & expected_files) - set(['orf1ab_msa.fasta'])
        if len(files) == 0:
            raise CommandError(f'in {input_msa} cannot find any of the following files {" ".join(expected_files)}')
        pbar = tqdm(files, file=self.stderr)
        for msa_file in pbar:
            gene = msa_file.split("_msa.fasta")[0]
            pbar.set_description(f"processing {gene}")
            msa_path = options['outdir_msa'] + "/" + msa_file
            self.process_msa(gene, msa_path, ref_seq, options["precompute_graphics"], options["override"])

        # if os.path.isdir(input_msa):
        #     expected_files = set(
        #         [x.accession + ".faa" for x in Bioentry.objects.filter(biodatabase__name="COVID19_prots")])
        #     files = set(os.listdir(input_msa)) & expected_files
        #     if len(files) == 0:
        #         raise CommandError(f'in {input_msa} cannot find any of the following files {" ".join(expected_files)}')
        #     pbar = tqdm(files, file=self.stderr)
        #     for msa_file in pbar:
        #         gene = msa_file.split(".faa")[0]
        #         pbar.set_description(f"processing {gene}")
        #         ref_seq = options["reference"] if options["reference"] else gene
        #         self.process_msa(gene, input_msa + "/" + msa_file, ref_seq, options["precompute_graphics"])
        # else:
        #     if not options["accession"]:
        #         raise CommandError(f'if input_msa is a file the gene accession cant be empty ')
        #     if not Bioentry.objects.filter(accession=options["accession"]).exists():
        #         raise CommandError(f'{options["accession"]} is not a valid orf')
        #     ref_seq = options["reference"] if options["reference"] else options["accession"]
        #     self.process_msa(options["accession"], input_msa, ref_seq, options["precompute_graphics"])

        self.stderr.write("Finished!")

    def process_msa(self, gene, msa_file, ref_seq, precompute_graphics=True, override=False):
        self.stderr.write(f"reading msa {msa_file}")
        msa_map = {}
        for r in bpio.parse(msa_file, "fasta"):
            r.name = ""
            r.description = ""
            # hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124|2019-12-30
            if r.id != ref_seq:

                try:
                    rid = r.id.replace("hCoV-19/", "")
                    if "PAIS" in rid:
                        code = rid.split("/")[1]
                        gisaid = ""
                        sdate = datetime.strptime("2020", '%Y').date()
                        country = "Argentina"
                    else:
                        code, gisaid, sdate = rid.split("|")
                        country = country_from_gisaid(r.id)
                        sdate = datetime.strptime(sdate.split("_")[0], '%Y-%m-%d').date()
                except Exception:
                    traceback.print_exc(file=self.stderr)
                    err = f'{r.id} does not have the correct format. Ex: hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124|2019-12-30'
                    raise CommandError(err)
                r.id = code
                Sample.objects.get_or_create(name=r.id, date=sdate, gisaid=gisaid, country=country)
            msa_map[r.id] = r
        if ref_seq not in msa_map:
            raise CommandError(f'{gene} not in {msa_file}')
        self.stderr.write(f"processing msa {msa_file}")
        msa = MSAMap(msa_map)
        msa.init()
        be = Bioentry.objects.get(accession=gene)
        seq = be.seq.seq
        Variant.objects.filter(bioentry=be).delete()
        for ref_pos, variant_samples in tqdm(msa.variants(ref_seq).items(), file=self.stderr):
            ref, pos = ref_pos.split("_")
            pos = int(pos)
            if ref != "*":
                if seq[pos] != ref:
                    err = 'sequence reference and alignment are different: '
                    try:
                        aln_pos = msa.pos_seq_msa_map[pos]
                    except IndexError:
                        aln_pos = "?"
                    except KeyError:
                        aln_pos = "?"
                    err += f'alnpos: {aln_pos}  seqpos: {ref_pos} aln_aa: {seq[pos]}  seq_aa:{ref}'
                    raise CommandError(err)

            for alt, samples in variant_samples.items():
                if (alt != ref) and (alt != "X"):
                    for sample_name in samples:
                        variant = Variant.objects.get_or_create(bioentry=be, pos=pos, ref=ref)[0]
                        sample = Sample.objects.get(name=sample_name)
                        SampleVariant.objects.get_or_create(variant=variant, alt=alt, sample=sample)
        variant_positions = set(Variant.objects.filter(bioentry=be).values_list("pos", flat=True))

        if precompute_graphics:
            self.stderr.write(f"pre computing graphics from variant positions")
            for variant_pos in tqdm(sorted(list(variant_positions))):
                fig_path = f'{STATICFILES_DIRS[0]}/auto/posfigs/{gene}{variant_pos}.png'
                if override or (not os.path.exists(fig_path)) or (os.path.getsize(fig_path) < 100):
                    variant_graphics(ref_seq, variant_pos, fig_path, msa_file, msa)
