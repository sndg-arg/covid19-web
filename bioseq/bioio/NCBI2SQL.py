import gzip
import os
import subprocess as sp

import Bio.SeqIO as bpio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from BioSQL import BioSeqDatabase
from biosql.models import Biodatabase
from tqdm import tqdm


class NCBI2SQL():

    def __init__(self):
        self.server = None

    def download(self, assembly_name, workdir):
        arr = assembly_name.split("_")
        acc = "_".join(arr[:2])
        name = "_".join(arr[2:])
        asspath = "/".join([acc[0:3], acc[4:7], acc[7:10], acc[10:13],
                            acc + "_" + name.replace(" ", "_").replace("#", "_")])

        cmd = 'rsync --recursive --include="*_genomic.gbff.gz" --exclude="*"   rsync://ftp.ncbi.nlm.nih.gov/genomes/all/' + asspath + '/ "' + \
              workdir + '"'
        with open(os.devnull, 'w') as FNULL:
            sp.call(cmd, shell=True,stdout=FNULL)


