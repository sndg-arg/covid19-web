#!/usr/bin/env python3
from tqdm import tqdm
import argparse
import os
import subprocess as sp
import Bio.SeqIO as bpio
from  termcolor import colored

parser = argparse.ArgumentParser(description='')
parser.add_argument("-i", '--dir_entrada', type=str, help='Directorio del analisis previo')
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
din = args.dir_entrada
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

pbar = tqdm(os.listdir(din))
for sample in pbar:
    pbar.set_description(sample)
    if not os.path.exists(f"{din}/{sample}/{sample}_call_consensus.vcf.gz"):
        print(colored(f"no existe {din}/{sample}/{sample}_call_consensus.vcf.gz","red"))
        continue

    if not os.path.exists(f"{out}/{sample}"):
        os.makedirs(f"{out}/{sample}")

    cmd = f'makeblastdb -dbtype nucl -in {din}/{sample}/{sample}_consensus.fasta.unmasked >/dev/null'
    exec(cmd, verbose=verbose)
    cmd = f'''blastn -query orfs.fna -db {din}/{sample}/{sample}_consensus.fasta.unmasked \
                -outfmt "6 sseqid sstart send qseqid" -qcov_hsp_perc 100  >  {out}/{sample}/{sample}_orfs.bed'''
    exec(cmd, verbose=verbose)
    cmd = f'seqtk subseq {din}/{sample}/{sample}_consensus.fasta {out}/{sample}/{sample}_orfs.bed > {out}/{sample}/{sample}_orfs.fna'
    exec(cmd, verbose=verbose)

    beds = {}
    for x in open(f'{out}/{sample}/{sample}_orfs.bed'):
        seq,start,end,desc = x.strip().split("\t")
        key = f'{seq}:{int(start)+1}-{end}'
        beds[key] = desc

    with open(f'{out}/{sample}/{sample}_Ns.txt',"w") as h:
        h.write( "orf #Ns\n")
        Ns = {}
        for orf in bpio.parse(f'{out}/{sample}/{sample}_orfs.fna', "fasta"):
            seq = str(orf.seq)
            Ns[orf.id] = 0

            for nucl in seq:

                if nucl == "N":
                    Ns[orf.id] += 1
            h.write(f"{beds[orf.id]} {Ns[orf.id]}" + "\n")

    cmd = f'java -jar /app/snpEff/snpEff.jar covid19  {din}/{sample}/{sample}_call_consensus.vcf.gz > {out}/{sample}/{sample}_ann.vcf'
    exec(cmd, verbose=verbose)


    with open(f'{out}/{sample}/{sample}_ann.vcf') as h , open(f'{out}/{sample}/{sample}_vars.txt',"w") as hw:
        hw.write( "pos ref alt annotation gene hgvs_c hgvs_p\n")
        for line in h:
            if not line.startswith("#"):
                contig, pos, _, ref, alt, _, _, info, _, _ = line.split("\t")
                ann_str = info.split("ANN=")[1].split(",")[0]


                (alt, annotation, impact, gene, geneid, feature_type, feature_id, transcript_biotype,
                 rank_div_total, hgvs_c, hgvs_p, c_dna_pos, cds_pos, (aa_pos_aa_len), dist_to_feature,
                 errors) = ann_str.split(
                    "|")
                hw.write( " ".join([ pos, ref, alt, annotation,gene,hgvs_c,hgvs_p]) + "\n")

    cmd = f'blastn -task "blastn-short"  -query {primers} -db /ref/MN996528.fna -qcov_hsp_perc 100 -outfmt "6 sseqid sstart send qseqid" > /tmp/{sample}primers.bed'
    exec(cmd, verbose=verbose)
    cmd = f'bedtools  intersect -a {din}/{sample}/{sample}_call_consensus.vcf.gz -b /tmp/{sample}primers.bed  -wb > /tmp/{sample}primers.intersect'
    exec(cmd, verbose=verbose)
    with open(f'/tmp/{sample}primers.intersect') as h, open(f'{out}/{sample}/{sample}_primers.txt',"w") as hw:
        hw.write("pos ref alt primer\n")
        for line in h:
            vec = line.strip().split()
            # MN996528.1      3037    .       C       T       225     .       DP=26;VDB=0.343005;SGB=-0.690438;MQSB=1;MQ0F=0;AC=1;AN=1;DP4=0,0,5,12;MQ=60     GT:PL   1:255,0 MN996528.1      3028    3047    p2
            pos,ref,alt,primer = vec[1], vec[3],vec[4],vec[-1]
            hw.write(" ".join([pos,ref,alt,primer]) + "\n")
