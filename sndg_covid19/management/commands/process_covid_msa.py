import os

import Bio.SeqIO as bpio

from django.core.management.base import BaseCommand
from bioseq.io.MSAMap import MSAMap
from bioseq.models.Variant import Variant, SampleVariant
from bioseq.models.Bioentry import Bioentry
from tqdm import tqdm


class Command(BaseCommand):
    """

    """

    DEFAULT_REFERENCE = "hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124|2019-12-30"

    help = 'Load variant list from msa'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dbx_dict = {}

    def add_arguments(self, parser):
        parser.add_argument("-r", '--reference',
                            default=os.environ.get("DEFAULT_REFERENCE", Command.DEFAULT_REFERENCE))

        parser.add_argument('-a','--accession', required=True)
        parser.add_argument("-i", '--input_msa', required=True)

    def handle(self, *args, **options):
        msa = {}
        self.stderr.write(f"reading msa {options['input_msa']}")
        for r in bpio.parse(options["input_msa"], "fasta"):
            r.id = r.description.replace(" ", "_")
            assert r.id not in msa
            msa[r.id] = r
        self.stderr.write(f"processing msa {options['input_msa']}")
        msa = MSAMap(msa)
        msa.init()
        be = Bioentry.objects.get(accession=options["accession"])
        seq = be.seq.seq
        for ref_pos, variant_samples in tqdm(msa.variants(options["reference"]).items(), file=self.stderr):
            ref, pos = ref_pos.split("_")
            pos = int(pos)
            if ref != "*":
                assert seq[pos] == ref

            variant = Variant.objects.get_or_create(bioentry=be, pos=pos, ref=ref)[0]
            for alt, samples in variant_samples.items():
                if alt != ref:
                    for sample in samples:
                        SampleVariant.objects.get_or_create(variant=variant, alt=alt, name=sample)

        # counts = {}
        # for pos,variant_samples in msa.variants(options["reference"]).items():
        #     data = []
        #     for variant,samples in variant_samples.items():
        #         data.append(variant + f":{len(samples)}" )
        #     counts[pos] = " ".join(data)
        # for k,v in counts.items():
        #     print (k,v)
