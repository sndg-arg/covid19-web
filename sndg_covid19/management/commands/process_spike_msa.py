from datetime import datetime
import os

import Bio.SeqIO as bpio
from django.core.management.base import BaseCommand, CommandError
from tqdm import tqdm

from bioseq.bioio.MSAMap import MSAMap
from bioseq.models.Bioentry import Bioentry
from bioseq.models.Variant import Variant, SampleVariant, Sample
from config.settings.base import STATICFILES_DIRS
from sndg_covid19.bioio import country_from_gisaid
from sndg_covid19.tasks import variant_graphics
from glob import glob
import traceback
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import math
import json


class Command(BaseCommand):
    """

    """

    help = 'Load variant list from msa'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dbx_dict = {}

    def add_arguments(self, parser):
        parser.add_argument("-r", '--reference', default="S", help="ref sequence, default='S'")
        parser.add_argument("-i", '--input_msa', required=True,
                            help="Could be a fasta msa file or a directory that contains a list of them")
        parser.add_argument("-o", '--outdir_msa', default="data/processed/", help="directory for protein MSAs")
        parser.add_argument("-pc", '--precompute_graphics', help="skip precompute", action="store_false")
        # parser.add_argument("--remove", help="when override is ON, samples with this filter are deleted", default=None)

    def handle(self, *args, **options):

        self.process_msa("S", options['input_msa'], options["reference"], options["precompute_graphics"], True)

        self.stderr.write("Finished!")

    def process_msa(self, gene, msa_file, ref_seq, precompute_graphics=True, override=False):
        self.stderr.write(f"reading msa {msa_file}")
        msa_map = {}
        for r in tqdm(bpio.parse(msa_file, "fasta")):
            r.name = ""
            desc = ("" if len(r.description.split(r.id)) == 1 else r.description.split(r.id)[1]).strip()
            # Spike|hCoV-19/Argentina/PAIS-S-A0001/2020 |2020-04-16|SouthAmerica
            # Spike|hCoV-19/Australia/VIC609/2020       |2020-03-27|EPI_ISL_426863|Original|hCoV-19
            if r.id != ref_seq:

                try:
                    rid = r.id.replace("Spike|hCoV-19/", "")
                    # if "PAIS" in rid:
                    #     code = rid.split("/")[1]
                    #     gisaid = ""
                    #     sdate = datetime.strptime("2020-6", '%Y-%m').date()
                    #     country = "Argentina"
                    # else:
                    code, sdate = rid.split("|")[:2]
                    country = code.split("/")[0]
                    sdate = datetime.strptime(sdate.split("_")[0], '%Y-%m-%d').date()
                    gisaid = ""
                except Exception:
                    traceback.print_exc(file=self.stderr)
                    err = f'{r.id} does not have the correct format. Ex: hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124|2019-12-30'
                    raise CommandError(err)
                r.id = code.split("/")[-2]
                sample = Sample.objects.get_or_create(name=r.id, date=sdate, gisaid=gisaid, country=country)[0]
                assert desc
                sample.subdivision = desc.strip()
                sample.save()
            msa_map[r.id] = r
        if ref_seq not in msa_map:
            raise CommandError(f'{gene} not in {msa_file}')
        self.stderr.write(f"processing msa {msa_file}")
        msa = MSAMap(msa_map)
        msa.init()
        be = Bioentry.objects.get(accession=gene)
        seq = be.seq.seq

        # pbar = tqdm(msa.variants(ref_seq).items(), file=self.stderr)
        for ref_pos, variant_samples in msa.variants(ref_seq).items():
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
                        SampleVariant.objects.filter(sample__name=sample_name, variant__pos=pos).delete()
                        variant = Variant.objects.get_or_create(bioentry=be, pos=pos, ref=ref)[0]
                        sample = Sample.objects.get(name=sample_name)
                        SampleVariant.objects.get_or_create(variant=variant, alt=alt, sample=sample)
        variant_positions = set(Variant.objects.filter(bioentry=be).values_list("pos", flat=True))

        if precompute_graphics:
            self.stderr.write(f"pre computing graphics from variant positions")
            for variant_pos in tqdm(sorted(list(variant_positions))):
                fig_path = f'{STATICFILES_DIRS[0]}/auto/posfigs_spike/{gene}{variant_pos}.png'
                if override or (not os.path.exists(fig_path)) or (os.path.getsize(fig_path) < 100):
                    variant_graphics(ref_seq, variant_pos, fig_path, msa_file, msa,idx_date=-4)
