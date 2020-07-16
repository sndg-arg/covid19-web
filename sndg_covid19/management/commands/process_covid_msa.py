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

class Command(BaseCommand):
    """

    """

    help = 'Load variant list from msa'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dbx_dict = {}

    def add_arguments(self, parser):
        parser.add_argument("-r", '--reference', default=None)
        parser.add_argument('-a', '--accession', help="If input_msa is a file is required", required=False)
        parser.add_argument("-i", '--input_msa', required=True,
                            help="Could be a fasta msa file or a directory that contains a list of them")
        parser.add_argument("-pc", '--precompute_graphics', help="skip precompute", action="store_false")

    def handle(self, *args, **options):



        input_msa = options['input_msa']

        if os.path.isdir(input_msa):
            expected_files = set(
                [x.accession + ".faa" for x in Bioentry.objects.filter(biodatabase__name="COVID19_prots")])
            files = set(os.listdir(input_msa)) & expected_files
            if len(files) == 0:
                raise CommandError(f'in {input_msa} cannot find any of the following files {" ".join(expected_files)}')
            pbar = tqdm(files,file=self.stderr)
            for msa_file in pbar:
                gene = msa_file.split(".faa")[0]
                pbar.set_description(f"processing {gene}")
                ref_seq = options["reference"] if options["reference"] else gene
                self.process_msa(gene, input_msa + "/" + msa_file, ref_seq, options["precompute_graphics"])
        else:
            if not options["accession"]:
                raise CommandError(f'if input_msa is a file the gene accession cant be empty ')
            if not Bioentry.objects.filter(accession=options["accession"]).exists():
                raise CommandError(f'{options["accession"]} is not a valid orf')
            ref_seq = options["reference"] if options["reference"] else options["accession"]
            self.process_msa(options["accession"], input_msa, ref_seq, options["precompute_graphics"])

        self.stderr.write("Finished!")

    def process_msa(self, gene, msa_file, ref_seq, precompute_graphics=True):
        self.stderr.write(f"reading msa {msa_file}")
        msa_map = {}
        for r in bpio.parse(msa_file, "fasta"):
            r.name = ""
            r.description = ""
            # hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124|2019-12-30
            if r.id != ref_seq:

                try:
                    rid = r.id.replace("hCoV-19/", "")
                    code, gisaid, sdate = rid.split("|")
                    country = country_from_gisaid(r.id)
                    sdate = datetime.strptime(sdate, '%Y-%m-%d').date()
                except Exception :
                    traceback.print_exc(file=self.stderr)
                    err = f'{r.id} does not have the correct format. Ex: hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124|2019-12-30'
                    raise CommandError(err)
                r.id = code
                Sample.objects.get_or_create(name=r.id, date=sdate, gisaid=gisaid, country=country)
            msa_map[r.id] = r
        if gene not in msa_map:
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
                        aln_pos = msa.pos_seq_msa_map[ref_pos]
                    except IndexError:
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
                variant_graphics(gene, variant_pos, fig_path, msa_file, msa)
