#!/usr/bin/env python3
# Samples 0001, 0005 and 0033 must be processed with :
# process_batch.py -i samples_dir -o ./results/

import argparse


parser = argparse.ArgumentParser(description='test some fasta consensus results')
parser.add_argument('input_dir', help="where samples processed with process_batch.py are")
parser.add_argument('ref', help="reference fasta",nargs="?", default="/ref/MN996528.fna")

args = parser.parse_args()

import Bio.SeqIO as bpio
from glob import glob

for consensus_file in glob(f"{args.input_dir}/FASTA/*.fasta"):
    consensus = bpio.read(consensus_file, "fasta")
    assert set(["A", "C", "G", "T", "N"]) == set(
        list(str(consensus.seq)))  , f"{consensus.id} : error in alphabet {set(list(str(consensus.seq)))}"
ref = bpio.read(args.ref, "fasta")
seq = bpio.read(glob(f"{args.input_dir}/FASTA/*0005*.fasta")[0], "fasta")
assert (len(ref.seq) - len(seq.seq)) == 44, f"44bp insertion not found in PAIS-A0005"
seq = bpio.read(glob(f"{args.input_dir}/FASTA/*0001*.fasta")[0], "fasta")
assert len(seq.seq) == len(ref.seq), f"PAIS-A0001, should not have any gaps"
assert str(seq.seq[28880:28883]) == "AAC", f"PAIS-A0001 28881-28883 mutation not found"
seq = bpio.read(glob(f"{args.input_dir}/FASTA/*0033*.fasta")[0], "fasta")
assert str(seq.seq[14407]) == "T", f"spanning delete mutation in PAIS-A0033 pos 14408 not found"

print("everything ok!")
