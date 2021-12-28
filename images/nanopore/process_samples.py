#!/usr/bin/env python3

import os
import sys
import subprocess as sp
from glob import glob
from tqdm import tqdm

if __name__ == "__main__":
    """
    runs this commands for every barcode
    artic guppyplex --min-length 100 --max-length 1000 --directory barcode01 --prefix run_barcode01/clean
    artic minion --medaka --medaka-model r941_min_fast_g303 --normalise 400 --threads 4 \
        --scheme-directory /artic-ncov2019/primer_schemes --read-file run_barcode01/clean_barcode01.fastq nCoV-2019/V3 run_barcode01
    """

    import argparse
    from functools import reduce
    import re
    import Bio.SeqIO as bpio

    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-c", '--sample-data', type=str, help='csv with directory and sample name', required=True)
    parser.add_argument("-i", '--input-dir', type=str, help='input dir for samples', required=True)
    parser.add_argument("-r", '--run-name', type=str, help='run name to create the fasta result', required=True)
    parser.add_argument("-o", '--output-dir', type=str, help='output dir for processed samples', default="./results/")
    parser.add_argument("-min", '--min-length', default=400, help='filter min size')
    parser.add_argument("-max", '--max-length', default=700, help='filter max size')
    parser.add_argument("-t", '--threads', default=4, help='number of threads')
    parser.add_argument('--normalize', default=None, type=int,
                        help="use --normalize x parameter in 'artic minion' command")
    parser.add_argument('--medakka', default="r941_min_fast_g303", help="--medaka-model parameter in 'artic minion'")
    parser.add_argument('--primers_dir', default="/artic-ncov2019/primer_schemes",
                        help="--scheme-directory parameter in 'artic minion'")
    parser.add_argument('--primers', default="nCoV-2019/V3", help="'scheme' parameter in 'artic minion'")

    args = parser.parse_args()

    assert os.path.exists(args.input_dir), f"{args.input_dir} does not exist"
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    assert os.path.exists(args.output_dir), f"can't create {args.output_dir}"
    assert os.path.exists(args.sample_data), f"{args.sample_data} does not exist"

    normalize = f"--normalise {args.normalize}" if args.normalize else ""

    try:
        barcode_header = {}
        with open(args.sample_data) as h:
            for l in h:
                barcode, header = [x.strip() for x in l.strip().split()[:2]]
                if not re.match("^hCoV\-19\/Argentina\/PAIS\-[A-Za-z]\d\d\d\d\/20\d\d$", header):
                    sys.stderr.write(f"bad header format: {header}\n")
                    sys.exit(1)
                barcode_header[barcode] = header
    except:
        sys.stderr.write("error reading CSV\n")
        sys.exit(1)

    samples = {}
    for sample_dir in os.listdir(args.input_dir):

        sample_fastqs = reduce(list.__add__, [glob(os.sep.join([args.input_dir, sample_dir, x]))
                                              for x in ["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"]])
        if sample_fastqs:
            samples[sample_dir] = sample_fastqs

    if not samples:
        sys.stderr.write(f"No data in {args.input_dir}\n")
        sys.exit(1)
    for sample in samples:
        if sample not in barcode_header:
            sys.stderr.write(f"Sample {sample} is not in the CSV file\n")
            sys.exit(1)

    diff = set(barcode_header) - set(samples)
    if diff:
        sys.stderr.write(f"Samples {','.join(diff)} are in the CSV file but not in the input directory\n")
        sys.exit(1)

    pbar = tqdm(samples.items(), file=sys.stderr)
    for sample, fastqs in pbar:
        pbar.set_description(f"processing: {sample}")
        prefix = os.sep.join([args.output_dir, sample, "clean"])
        outfastq = os.sep.join([prefix + "_" + sample + ".fastq"])
        indir = os.sep.join([args.input_dir, sample])
        outdir = os.sep.join([args.output_dir, sample])
        outdir_prefix = os.sep.join([args.output_dir, sample, sample])
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        cmd1 = f'''artic guppyplex --min-length {args.min_length} --max-length {args.max_length} \
                    --directory {indir} --prefix {prefix}'''

        cmd2 = f'''artic minion --medaka --medaka-model {args.medakka}  {normalize} --threads {args.threads} \
                    --scheme-directory {args.primers_dir} --read-file {outfastq} {args.primers} {outdir_prefix}'''
        with open("process.log", "a") as h:
            h.write("running:\n" + cmd1 + "\n")
            try:
                sp.run(cmd1, shell=True, stderr=h, stdout=h)
            except:
                sys.stderr.write("error running:\n")
                sys.stderr.write(cmd1+ "\n")
                continue
            h.write("running:\n" + cmd2 + "\n")
            try:
                sp.run(cmd2, shell=True, stderr=h, stdout=h)

            except:
                sys.stderr.write("error running:\n")
                sys.stderr.write(cmd2 + "\n")
                continue
            try:
                with open(f"{args.output_dir}/{args.run_name}.fasta","a") as h:
                    r = bpio.read(f'{outdir}.consensus.fasta',"fasta")
                    r.name = ""
                    r.description = ""
                    r.id = barcode_header[sample]
                    bpio.write(r,h,"fasta")
            except:
                sys.stderr.write(f"error writing {outdir}.consensus.fasta\n")
