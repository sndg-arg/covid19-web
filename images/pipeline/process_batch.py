#!/usr/bin/env python3
import os
import sys
import subprocess as sp
from glob import glob
from collections import defaultdict
from tqdm import tqdm
import argparse
import re
from  termcolor import colored
import shutil
parser = argparse.ArgumentParser(description='')
parser.add_argument("-i", '--dir_entrada', type=str, help='Directorio con los fastq, un solo par por muestra')
parser.add_argument("-o", '--dir_salida', type=str, help='Directorio donde se pone las salidas por muestra ')
args = parser.parse_args()

assert os.path.exists(args.dir_entrada), f"{args.dir_entrada} no existe"
if not os.path.exists(args.dir_salida):
    os.makedirs(args.dir_salida)
if not os.path.exists(args.dir_salida + "/FASTA" ):
    os.makedirs(args.dir_salida + "/FASTA" )

muestras = defaultdict(list)
for arch in glob(f"{args.dir_entrada}/*.fastq.gz") + glob(f"{args.dir_entrada}/*.fastq"):
    muestra = arch.split("/")[-1][:10]
    if re.match("^PAIS\-[A-Za-z]\d\d\d\d$",muestra):
        muestras[muestra].append(arch)
    else:
        print(colored(f"{muestra} no tiene el nombre de muestra correcto: PAIS-X0000",color="red"))

sys.stderr.write(f"se encontraron {len(muestras)} en {args.dir_entrada}\n")
for k, v in muestras.items():
    assert len(v) == 2, f"error: {k} tiene {len(v)} archivos"


with tqdm(muestras.items()) as pbar:
    for k, v in pbar:
        pbar.set_description(f"procesando {k}")
        if not os.path.exists(args.dir_salida + "/" + k):
            os.makedirs(args.dir_salida + "/" + k)
        cmd = f'process_sample.sh {v[0]} {v[1]} {k} {args.dir_salida} > {args.dir_salida}/{k}/log.txt   2>&1'
        sp.call(cmd, shell=True)
        shutil.copy(f"{args.dir_salida}/{k}/{k}_consensus.fasta",f"{args.dir_salida}/FASTA")
