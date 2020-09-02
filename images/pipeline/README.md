# build image
```shell script
docker build -t sndg/covid19  .
```

# Interactive use
```shell script
docker run -it --rm -v $PWD:/out sndg/covid19 bash
./process_batch.py -i fastq_dir -o output_dir
```

# Process reads pair
```shell script
process_sample.sh /data/test_set/PAIS-A001_S35_L001_R1_001.fastq.gz /data/test_set/PAIS-A001_S35_L001_R2_001.fastq.gz PAIS-A001 ./results
```

# Run test

```shell script
docker run -it --rm -v $PWD:/out sndg/covid19 bash -c 'process_batch.py -i test_samples/fortest/ -o ./results/'
```
```shell script
docker run -it --rm -v $PWD:/out sndg/covid19 bash -c 'test_consensus.py results/'
```
```shell script
docker run -it --rm -v $PWD:/out sndg/covid19 bash -c 'process_assembly.py -i results/FASTA -o results -p test_primers/PrimersGenomaSARS-CoV2Sanger.txt'
```
Se le cree al variant caller a menos que:
* la cobertura sea baja y haya evidencia de una delecion con AD > 10 => va la delecion y no N
* la profundidad sea menor a 10, caso en el cual queda una N o no se tiene en cuenta el indel
Se le cree al variant caller:
* Heterocigosis: si la variante es Het con prop de alt > 60% va alt, si ref > 60% se borra la variante, sino va base ambigua
* Homocigosis y DP >= 10 => va el alternativo sin duda
