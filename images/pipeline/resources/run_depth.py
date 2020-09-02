#!/usr/bin/env python3
import os
import sys
import subprocess as sp
from glob import glob
from collections import defaultdict
from tqdm import tqdm
import argparse
import re
from termcolor import colored
import shutil

parser = argparse.ArgumentParser(description='')
parser.add_argument("-i", '--dir_entrada', type=str, help='Directorio con resuldatos de process_batch')
parser.add_argument("-o", '--dir_salida', type=str, help='Directorio de salida')

args = parser.parse_args()

assert os.path.exists(args.dir_entrada), f"{args.dir_entrada} no existe"

if not os.path.exists(args.dir_salida):
    os.makedirs(args.dir_salida)

assert os.path.exists(args.dir_salida), f"no se pudo crear el directorio de salida {args.dir_salida}"

samples = {}
for bam in glob(args.dir_entrada + "/*/*.bam"):
    sample = bam.split("/")[-2]
    samples[sample] = bam

with open(f'{args.dir_salida}/depth.txt', "w") as hw:
    with tqdm(samples.items()) as pbar:
        for sname, sbam in pbar:
            pbar.set_description(f"procesando {sname}")
            if not os.path.exists(f'{args.dir_salida}/{sname}.depth') or (
                os.path.getsize(f'{args.dir_salida}/{sname}.depth') < 100):
                cmd = f'samtools depth -a {sbam} > {args.dir_salida}/{sname}.depth'
                sp.call(cmd, shell=True)
            with open(f'{args.dir_salida}/{sname}.depth') as h:
                count = 0
                acc = 0
                for line in h:
                    count = count + 1
                    depth = int(line.split()[-1])
                    acc = depth + acc
                prom_depth = acc / count
                hw.write(f'{sname}\t{prom_depth:.2f}\n')
