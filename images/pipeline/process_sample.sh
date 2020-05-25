SAMPLE_R1=$1
SAMPLE_R2=$2
SAMPLE_NAME=$3
RESULTS="${4:-./}"
REFERENCE="${REFERENCE:-/ref/genome.fna}"
ADAPTERS="${ADAPTERS:-/ref/adapters.fasta}"
CPUS="${CPUS:-4}"

TRIMR="${TRIMR:-0}"
TRIML="${TRIML:-10}"
MINLEN="${MINLEN:-40}"
QPROM="${QPROM:-30}"
NPROM="${QPROM:-20}"

#*************************************************

#IFS='/' tokens=( $SAMPLE_R1 )
#SAMPLE_NAME=${tokens[1]}
DIR=${RESULTS}/${SAMPLE_NAME}
IFS=' '
[ -d $DIR ] || mkdir $DIR

echo "Current sample: ${SAMPLE_NAME}"
echo "R1: $SAMPLE_R1"
echo "R2: $SAMPLE_R2"
echo "Results: ${RESULTS}"
echo "Ref: $REFERENCE"
echo "Adapters: $ADAPTERS"
echo "TRIMR: $TRIMR"
echo "TRIML: $TRIML"
echo "MINLEN: $MINLEN"
echo "QPROM: $QPROM"
echo "NPROM: $NPROM"
echo "CPUS: $CPUS"

fastqc $SAMPLE_R1 $SAMPLE_R2 -o ${RESULTS}/${SAMPLE_NAME}/ -q

#limpiar fastq con fastp o prinseq
echo "Cleaning FASTQ...."
IN1=$SAMPLE_R1
IN2=$SAMPLE_R2
OUT1=${RESULTS}/${SAMPLE_NAME}/1_fastp.fastq
OUT2=${RESULTS}/${SAMPLE_NAME}/2_fastp.fastq
LOG=${RESULTS}/${SAMPLE_NAME}/fastp.log
echo fastp -f $TRIML -t $TRIMR -l $MINLEN -e $QPROM -n $NPROM --thread $CPUS -i $IN1 -I $IN2 -o $OUT1 -O $OUT2
fastp -f $TRIML -t $TRIMR -l $MINLEN -e $QPROM -n $NPROM --thread $CPUS -i $IN1 -I $IN2 -o $OUT1 -O $OUT2 | tee $LOG

rm fastp.html  fastp.json

#Limpiar adaptadores
echo "Removing adapters..."
IN1=$OUT1
IN2=$OUT2
OUT1=${RESULTS}/${SAMPLE_NAME}/1_cutadapt.fastq
OUT2=${RESULTS}/${SAMPLE_NAME}/2_cutadapt.fastq
LOG=${RESULTS}/${SAMPLE_NAME}/fastp.log
echo cutadapt -a file:$ADAPTERS -o $OUT1 -p $OUT2 $IN1 $IN2
cutadapt --cores=$CPUS -a file:${ADAPTERS} -o $OUT1 -p $OUT2 $IN1 $IN2 | tee $LOG
rm $IN1 $IN2

##remover duplicados con bbmap
echo "Removing duplicates..."
IN1=$OUT1
IN2=$OUT2
OUT1=${RESULTS}/${SAMPLE_NAME}/1_dedupped.fastq
OUT2=${RESULTS}/${SAMPLE_NAME}/2_dedupped.fastq
LOG1=${RESULTS}/${SAMPLE_NAME}/1_dedupped.log
LOG2=${RESULTS}/${SAMPLE_NAME}/2_dedupped.log
dedupe.sh threads=$CPUS overwrite=t in=$IN1 out=$OUT1 |& tee $LOG1
dedupe.sh threads=$CPUS overwrite=t in=$IN2 out=$OUT2 |& tee $LOG2
rm $IN1 $IN2

##reparar los fastq para que coincidan la cantidad y orden de los reads
echo "Re-pairing FASTQ..."
IN1=$OUT1
IN2=$OUT2
OUT1=${RESULTS}/${SAMPLE_NAME}/1_repaired.fastq
OUT2=${RESULTS}/${SAMPLE_NAME}/2_repaired.fastq
SINGLE=${RESULTS}/${SAMPLE_NAME}/repaired_single.fastq
LOG=${RESULTS}/${SAMPLE_NAME}/1_dedupped.log
repair.sh threads=$CPUS overwrite=t in1=$IN1 in2=$IN2 out1=$OUT1 out2=$OUT2 outsingle=$SINGLE  |& tee $LOG
rm $IN1 $IN2 $SINGLE

#Mapeo
echo "Mapping..."
IN1=$OUT1
IN2=$OUT2
OUT1=${RESULTS}/${SAMPLE_NAME}/aln.sam
OUT2=${RESULTS}/${SAMPLE_NAME}/aln.bam
LOG=${RESULTS}/${SAMPLE_NAME}/bwa.log
STAT=${RESULTS}/${SAMPLE_NAME}/flagstat.txt
echo bwa mem -t $CPUS $REFERENCE $IN1 $IN2
bwa mem -t $CPUS $REFERENCE $IN1 $IN2 > $OUT1 2> $LOG
samtools sort $OUT1 -o $OUT2
samtools index $OUT2
samtools flagstat $OUT2 > $STAT
rm $OUT1

#low coverage
IN=$OUT2
TMP1=${RESULTS}/${SAMPLE_NAME}/my.sorted.bedgraph
TMP2=${RESULTS}/${SAMPLE_NAME}/uncovered.intevals.bedgraph
bedtools genomecov -ibam $IN -bg > $TMP1
awk '$4<2' $TMP1 > $TMP2
bedtools merge -i $TMP2 > ${RESULTS}/${SAMPLE_NAME}/uncovered.bed
rm $TMP1 $TMP2

#Llamado de variantes
echo "Calling variants..."
IN1=$OUT2
OUT1=${RESULTS}/${SAMPLE_NAME}/mpileup.vcf.gz
OUT2=${RESULTS}/${SAMPLE_NAME}/call.vcf.gz
bcftools mpileup -f $REFERENCE $IN1  -Oz > $OUT1
bcftools call --ploidy 1 -mv -Oz $OUT1 -o $OUT2
rm $OUT1

#Obtención de secuencia consenso
echo "Getting consensus sequence..."
IN1=$OUT2
OUT1=${RESULTS}/${SAMPLE_NAME}/consensus.fasta
bcftools index -f $IN1
bcftools consensus -f $REFERENCE -H A -o $OUT1 $IN1

#Assembly
echo "Getting assembly..."
IN1=${RESULTS}/${SAMPLE_NAME}/1_repaired.fastq
IN2=${RESULTS}/${SAMPLE_NAME}/2_repaired.fastq
OUT2=${RESULTS}/${SAMPLE_NAME}/assembly
megahit -1 $IN1 -2 $IN2 -o $OUT2
cp ${OUT2}/final.contigs.fa ${OUT2}.fna
rm -r $OUT2 # $IN1 $IN1

#Assembly
echo "denovo"
IN1=${RESULTS}/${SAMPLE_NAME}/1_repaired.fastq
IN2=${RESULTS}/${SAMPLE_NAME}/2_repaired.fastq
OUT2=${RESULTS}/${SAMPLE_NAME}/denovo
spades.py --meta -t 8 --pe-1 1 $IN1 --pe-2 1 $IN2 -o $OUT2
cp ${OUT2}/scaffolds.fasta ${OUT2}.fna
rm -r $OUT2 $IN1 $IN2

#Assembly vs ref
IN=${RESULTS}/${SAMPLE_NAME}/denovo.fna
TMP1=${RESULTS}/${SAMPLE_NAME}/nucdifftmp
nucdiff --vcf yes  /ref/genome.fna $IN $TMP1 $SAMPLE_NAME
mv ${TMP1}/results ${RESULTS}/${SAMPLE_NAME}/nucdiff
rm -r $TMP1

echo "DONE!"
