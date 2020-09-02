# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import os
import sys
import subprocess as sp

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.SearchIO as bpsio


class Blast():
    docker_img = "biocontainers/blast:v2.2.31_cv2"

    @staticmethod
    def makeblastdb(fasta_path, dbtype="auto",stdout=sys.stdout,stderr=sys.stderr):
        assert os.path.exists(fasta_path)
        fasta_path = os.path.abspath(fasta_path)
        fasta_dir = os.path.dirname(fasta_path)
        fasta_file = os.path.basename(fasta_path)
        if dbtype == "auto":
            if fasta_path.endswith(".fna"):
                dbtype = "nucl"
            if fasta_path.endswith(".fna"):
                dbtype = "aa"
            raise Exception("unknown file extention, we only recongnized .fna or .faa ")

        cmd = "docker run -v {fasta_dir}:/data {docker_img} makeblastdb -dbtype {dbtype} -in {fasta_file} "
        cmd = cmd.format(docker_img=Blast.docker_img, fasta_dir=fasta_dir, dbtype=dbtype, fasta_file=fasta_file)
        sp.call( cmd  ,shell="True",  stdout=stdout, stderr=stderr)
