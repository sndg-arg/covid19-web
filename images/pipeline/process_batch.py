#!/usr/bin/env python3
import os
import subprocess as sp
from glob import glob
from collections import defaultdict
from tqdm import tqdm

import argparse
parser = argparse.ArgumentParser(description='')
parser.add_argument("-i",'--dir_entrada', type=str, help='Directorio con los fastq, un solo par por muestra')
parser.add_argument("-o",'--dir_salida', type=str, help='Directorio donde se pone las salidas por muestra ')
args = parser.parse_args()

assert os.path.exists(args.dir_entrada), f"{args.dir_entrada} no existe"
if not os.path.exists(args.dir_salida):
    os.makedirs(args.dir_salida)

muestras = defaultdict(list)
for arch in glob(f"{args.dir_entrada}/*.fastq.gz"):
    muestra = arch.split("/")[-1][:9]
    muestras[muestra].append(arch)
print(f"se encontraron {len(muestras)} en {args.dir_entrada}")
for k,v in muestras.items():
    assert len(v) == 2, f"error: {k} tiene {len(v)} archivos"

with tqdm(muestras.items()) as pbar:
    for k,v in pbar:
        pbar.set_description(f"procesando {k}")

        cmd = f'process_sample.sh {v[0]} {v[1]} {k} {args.dir_salida}'
        sp.call(cmd,shell=True)
