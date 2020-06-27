import os

import Bio.SeqIO as bpio

from django.core.management.base import BaseCommand
from bioseq.io.MSAMap import MSAMap
from bioseq.models.Variant import Variant, SampleVariant, Sample
from bioseq.models.Bioentry import Bioentry
from tqdm import tqdm
from collections import defaultdict

from sndg_covid19.tasks import variant_graphics
from django.core.management.base import BaseCommand, CommandError
from datetime import datetime

from sndg_covid19.io import country_from_gisaid

class Command(BaseCommand):
    """

    """

    help = 'Load variant list from msa'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dbx_dict = {}

    def add_arguments(self, parser):
        parser.add_argument("-r", '--reference', default=None)
        parser.add_argument('-a', '--accession', required=True)
        parser.add_argument("-i", '--input_msa', required=True)
        parser.add_argument("-pc", '--precompute_graphics', help="skip precompute", action="store_false")

    def handle(self, *args, **options):
        msa_map = {}
        ref_seq = options["reference"] if options["reference"] else options["accession"]
        if not Bioentry.objects.filter(accession=options["accession"]).exists():
            raise CommandError(f'{options["accession"]} is not a valid orf')

        self.stderr.write(f"reading msa {options['input_msa']}")

        for r in bpio.parse(options["input_msa"], "fasta"):
            r.name = ""
            r.description = ""
            # hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124|2019-12-30
            if r.id != ref_seq:

                try:
                    rid = r.id.replace("hCoV-19/", "")
                    code, gisaid, sdate = rid.split("|")
                    country = country_from_gisaid(r.id)
                    sdate = datetime.strptime(sdate, '%Y-%m-%d').date()
                except Exception as ex:
                    err = f'{r.id} does not have the correct format. Ex: hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124|2019-12-30'
                    raise CommandError(err)
                r.id = code
                Sample.objects.get_or_create(name=r.id, date=sdate, gisaid=gisaid, country=country)
            msa_map[r.id] = r

        if options["accession"] not in msa_map:
            raise CommandError(f'{options["accession"]} not in {options["input_msa"]}')

        self.stderr.write(f"processing msa {options['input_msa']}")
        msa = MSAMap(msa_map)
        msa.init()

        be = Bioentry.objects.get(accession=options["accession"])
        seq = be.seq.seq

        Variant.objects.filter(bioentry=be).delete()

        for ref_pos, variant_samples in tqdm(msa.variants(ref_seq).items(), file=self.stderr):
            ref, pos = ref_pos.split("_")
            pos = int(pos)
            if ref != "*":
                if seq[pos] != ref:
                    err = 'sequence reference and alignment are different: '
                    try:
                        aln_pos = msa.pos_seq_msa_map[ref_pos]
                    except:
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

        if options["precompute_graphics"]:
            self.stderr.write(f"pre computing graphics from variant positions")
            for variant_pos in tqdm(sorted(list(variant_positions))):
                variant_graphics(options["accession"], variant_pos, msa)

        self.stderr.write("Finished!")
