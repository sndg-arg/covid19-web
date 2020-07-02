#!/usr/bin/env python3
from tqdm import tqdm
import argparse
import os
import subprocess as sp
import Bio.SeqIO as bpio
from termcolor import colored
from glob import glob

parser = argparse.ArgumentParser(description='')
parser.add_argument("-i", '--fasta', type=str, help='directorio de fastas no alineados')
parser.add_argument("-o", '--dir_salida', type=str, help='Directorio donde se ponen los reportes de salida')
parser.add_argument("-p", '--primers', type=str, help='fasta con primers')
parser.add_argument("-v", '--verbose', action='store_true')
args = parser.parse_args()

"""
1) % de cada ORF con datos (entiéndase: bases nucleotídicas)
2) % NNNs en cada ORF
3) gaps encontrados en los ORF
4) Codones START atrasados
5) Codones de STOP prematuros
6) Detección de mutaciones en los primers que se usan para hacer las qPCR de diagnóstico
"""
din = args.fasta
out = args.dir_salida
primers = args.primers
verbose = args.verbose

assert os.path.exists(din), f"no existe el directorio de entrada {din}"
if not os.path.exists(out):
    os.makedirs(out)
assert os.path.exists(primers), f"no existe el archivo de primers{primers}"


def exec(cmd, verbose=False):
    if verbose:
        print(cmd)
    sp.call(cmd, shell=True)


if os.path.isdir(din):
    pbar = tqdm(glob(f'{din}/*.fasta'))
    sample_fn = lambda item : os.path.basename(item).split(".fasta")[0]
    sample_file = lambda item : fasta_file
else:
    pbar = tqdm(list(bpio.parse(din,"fasta")))
    sample_fn = lambda item : item.id.replace("/","_").replace(" ","_")
    def sample_file(item) :
        fasta_path = f'/tmp/{item.id.replace("/","_").replace(" ","_")}.fasta'
        bpio.write(item,fasta_path,"fasta")
        return fasta_path



for item in pbar:
    fasta_file = sample_file(item)
    sample = sample_fn(item)
    pbar.set_description(fasta_file)

    if not os.path.exists(f"{out}/{sample}"):
        os.makedirs(f"{out}/{sample}")

    cmd = f'cat  /ref/MN996528.fna {fasta_file} > /tmp/seqs.fasta'
    exec(cmd, verbose=verbose)
    cmd = f'mafft --auto /tmp/seqs.fasta > /tmp/aln.fasta 2>>"{out}/{sample}/log.txt"'
    exec(cmd, verbose=verbose)
    cmd = f'python3 /app/script/MSAMap.py -r MN996528.1  -i /tmp/aln.fasta >> "{out}/{sample}/{sample}_variants.vcf" 2>>"{out}/{sample}/log.txt"'
    exec(cmd, verbose=verbose)

    cmd = f'java -jar /app/snpEff/snpEff.jar covid19  "{out}/{sample}/{sample}_variants.vcf" > "{out}/{sample}/{sample}_ann.vcf" 2>>"{out}/{sample}/log.txt"'
    exec(cmd, verbose=verbose)

    with open(f'{out}/{sample}/{sample}_ann.vcf') as h, open(f'{out}/{sample}/{sample}_vars.txt', "w") as hw:
        hw.write("pos ref alt annotation gene hgvs_c hgvs_p\n")
        for line in h:
            if not line.startswith("#"):
                contig, pos, _, ref, alt, _, _, info = line.split("\t")[:8]
                ann_str = info.split("ANN=")[1].split(",")[0]

                (alt, annotation, impact, gene, geneid, feature_type, feature_id, transcript_biotype,
                 rank_div_total, hgvs_c, hgvs_p, c_dna_pos, cds_pos, (aa_pos_aa_len), dist_to_feature,
                 errors) = ann_str.split( "|")
                hw.write(" ".join([pos, ref, alt, annotation, gene, hgvs_c, hgvs_p]) + "\n")

    if primers:
        cmd = f'blastn -task "blastn-short"  -query {primers} -db /ref/MN996528.fna -qcov_hsp_perc 100 -outfmt "6 sseqid sstart send qseqid" > /tmp/{sample}primers_raw.bed 2>/dev/null'
        exec(cmd, verbose=verbose)
        with open(f"/tmp/{sample}primers.bed","w") as h:
            for line in open(f'/tmp/{sample}primers_raw.bed'):
                vec = line.split("\t")
                if int(vec[1]) > int(vec[2]):
                    vec[1],vec[2] = vec[2],vec[1]
                h.write("\t".join(vec))
        cmd = f'bedtools  intersect -a "{out}/{sample}/{sample}_ann.vcf" -b "/tmp/{sample}primers.bed"  -wb > "/tmp/{sample}primers.intersect" 2>>"{out}/{sample}/log.txt"'
        exec(cmd, verbose=verbose)
        with open(f'/tmp/{sample}primers.intersect') as h, open(f'{out}/{sample}/{sample}_primers.txt', "w") as hw:
            hw.write("pos ref alt primer\n")
            for line in h:
                vec = line.strip().split()
                # MN996528.1      3037    .       C       T       225     .       DP=26;VDB=0.343005;SGB=-0.690438;MQSB=1;MQ0F=0;AC=1;AN=1;DP4=0,0,5,12;MQ=60     GT:PL   1:255,0 MN996528.1      3028    3047    p2
                pos, ref, alt, primer,pstart,pend = vec[1], vec[3], vec[4], vec[-1] , vec[-3] , vec[-2]
                hw.write(" ".join([pos, ref, alt, primer,pstart,pend]) + "\n")
    cmd = f"echo '{sample}' >> {out}/primers.txt ; cat '{out}/{sample}/{sample}_primers.txt' >> {out}/primers.txt"
    exec(cmd, verbose=verbose)
