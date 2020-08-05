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
parser.add_argument("-o", '--dir_salida', type=str, help='Directorio donde se pone las salidas por muestra')
parser.add_argument("-r", '--reference', default="/ref/MN996528.fna", help='archivo fasta de referencia')
parser.add_argument("-a", '--adapters', default="/ref/adapters.fasta", help='archivo fasta de adaptadores a remover')
parser.add_argument( '--cpus', default=4, type=int, help='cantidad de nucleos a utilizar')


parser.add_argument( '--remove_dupl', action="store_true", help='Cuando se activa NO se ejecuta fastuniq')
args = parser.parse_args()

assert os.path.exists(args.dir_entrada), f"{args.dir_entrada} no existe"
assert os.path.exists(args.reference), f"{args.reference} no existe"
dict_file = ".".join(args.reference.split(".")[:-1]) + ".dict"
assert os.path.exists(dict_file), f"{dict_file} no existe, debe crearlo con samtools dict"
assert os.path.exists(args.reference + ".bwt"), f"{args.reference + '.bwt'} no existe, debe crearlo con bwa index"


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
        cmd = f"export CPUS={args.cpus};export ADAPTERS={args.adapters};export REFERENCE={args.reference}"
        cmd = f'{cmd};process_sample.sh {v[0]} {v[1]} {k} {args.dir_salida} > {args.dir_salida}/{k}/log.txt   2>&1'
        if args.remove_dupl:
            cmd = f"export DEDUP=true;{cmd}"
        sp.call(cmd, shell=True)
        try:
            shutil.copy(f"{args.dir_salida}/{k}/{k}_consensus.fasta",f"{args.dir_salida}/FASTA")
        except:
            sys.stderr.write(f"error al procesar {k}\n")
