SAMPLE_R1=$1
SAMPLE_R2=$2
SAMPLE_NAME=$3
RESULTS="${4:-./}"
REFERENCE="${REFERENCE:-/ref/MN996528.fna}"
ADAPTERS="${ADAPTERS:-/ref/adapters.fasta}"
CPUS="${CPUS:-4}"

TRIMR="${TRIMR:-5}"
TRIML="${TRIML:-20}"
MINLEN="${MINLEN:-40}"
QPROM="${QPROM:-30}"
NPROM="${QPROM:-20}"

#*************************************************

#IFS='/' tokens=( $SAMPLE_R1 )
#SAMPLE_NAME=${tokens[1]}
DIR=${RESULTS}/${SAMPLE_NAME}
IFS=' '
[ -d $DIR ] || mkdir -p $DIR

echo "Current sample: ${SAMPLE_NAME}" >> ${RESULTS}/${SAMPLE_NAME}/log.txt
echo "R1: $SAMPLE_R1" >> ${RESULTS}/${SAMPLE_NAME}/log.txt
echo "R2: $SAMPLE_R2" >> ${RESULTS}/${SAMPLE_NAME}/log.txt
echo "Results: ${RESULTS}" >> ${RESULTS}/${SAMPLE_NAME}/log.txt
echo "Ref: $REFERENCE" >> ${RESULTS}/${SAMPLE_NAME}/log.txt
echo "Adapters: $ADAPTERS" >> ${RESULTS}/${SAMPLE_NAME}/log.txt
echo "TRIMR: $TRIMR" >> ${RESULTS}/${SAMPLE_NAME}/log.txt
echo "TRIML: $TRIML" >> ${RESULTS}/${SAMPLE_NAME}/log.txt
echo "MINLEN: $MINLEN" >> ${RESULTS}/${SAMPLE_NAME}/log.txt
echo "QPROM: $QPROM" >> ${RESULTS}/${SAMPLE_NAME}/log.txt
echo "NPROM: $NPROM" >> ${RESULTS}/${SAMPLE_NAME}/log.txt
echo "CPUS: $CPUS" >> ${RESULTS}/${SAMPLE_NAME}/log.txt



#limpiar fastq con fastp o prinseq
echo "Cleaning FASTQ...." >> ${RESULTS}/${SAMPLE_NAME}/log.txt
IN1=$SAMPLE_R1
IN2=$SAMPLE_R2
OUT1=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_1_fastp.fastq
OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_2_fastp.fastq
LOG=${RESULTS}/${SAMPLE_NAME}/fastp.log
echo fastp -f $TRIML -t $TRIMR -l $MINLEN -e $QPROM  --thread $CPUS -i $IN1 -I $IN2 -o $OUT1 -O $OUT2 >> ${RESULTS}/${SAMPLE_NAME}/log.txt
fastp -f $TRIML -t $TRIMR -l $MINLEN -e $QPROM -n $NPROM  --thread $CPUS -i $IN1 -I $IN2 -o $OUT1 -O $OUT2 > $LOG  2>&1

rm fastp.html  fastp.json

#Limpiar adaptadores
echo "Removing adapters..." >> ${RESULTS}/${SAMPLE_NAME}/log.txt
IN1=$OUT1
IN2=$OUT2
OUT1=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_1_cutadapt.fastq
OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_2_cutadapt.fastq
LOG=${RESULTS}/${SAMPLE_NAME}/fastp.log
echo cutadapt -a file:$ADAPTERS -o $OUT1 -p $OUT2 $IN1 $IN2 >> ${RESULTS}/${SAMPLE_NAME}/log.txt
cutadapt --cores=$CPUS -a file:${ADAPTERS} -o $OUT1 -p $OUT2 $IN1 $IN2 > $LOG  2>&1
rm $IN1 $IN2

##remover duplicados con bbmap
echo "Removing duplicates..." >> ${RESULTS}/${SAMPLE_NAME}/log.txt
IN1=$OUT1
IN2=$OUT2
OUT1=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_1_dedupped.fastq
OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_2_dedupped.fastq
LOG1=${RESULTS}/${SAMPLE_NAME}/1_dedupped.log
LOG2=${RESULTS}/${SAMPLE_NAME}/2_dedupped.log
dedupe.sh threads=$CPUS absorbmatch=t  overwrite=t in=$IN1 out=$OUT1 > $LOG1  2>&1
dedupe.sh threads=$CPUS absorbmatch=t  overwrite=t in=$IN2 out=$OUT2 > $LOG2  2>&1
rm $IN1 $IN2

##reparar los fastq para que coincidan la cantidad y orden de los reads
echo "Re-pairing FASTQ..." >> ${RESULTS}/${SAMPLE_NAME}/log.txt
IN1=$OUT1
IN2=$OUT2
OUT1=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_1_repaired.fastq
OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_2_repaired.fastq
SINGLE=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_repaired_single.fastq
LOG=${RESULTS}/${SAMPLE_NAME}/1_dedupped.log
repair.sh threads=$CPUS overwrite=t in1=$IN1 in2=$IN2 out1=$OUT1 out2=$OUT2 outsingle=$SINGLE  > $LOG  2>&1
rm $IN1 $IN2 $SINGLE

fastqc $OUT1 $OUT2 -o ${RESULTS}/${SAMPLE_NAME}/ -q


#Mapeo
echo "Mapping..." >> ${RESULTS}/${SAMPLE_NAME}/log.txt
IN1=$OUT1
IN2=$OUT2
OUT1=${RESULTS}/${SAMPLE_NAME}/aln.sam
OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_aln.bam
LOG=${RESULTS}/${SAMPLE_NAME}/bwa.log
STAT=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_flagstat.txt
echo bwa mem -t $CPUS $REFERENCE $IN1 $IN2 >> ${RESULTS}/${SAMPLE_NAME}/log.txt
bwa mem -t $CPUS $REFERENCE $IN1 $IN2 > $OUT1 2> $LOG
samtools sort $OUT1 -o $OUT2
samtools index $OUT2
samtools flagstat $OUT2 > $STAT  2>&1
rm $OUT1

#low coverage
IN=$OUT2
TMP1=${RESULTS}/${SAMPLE_NAME}/my.sorted.bedgraph
TMP2=${RESULTS}/${SAMPLE_NAME}/uncovered.intevals.bedgraph
bedtools genomecov -ibam $IN -bg > $TMP1  2>>${RESULTS}/${SAMPLE_NAME}/log.txt
awk '$4<10' $TMP1 > $TMP2 2>>${RESULTS}/${SAMPLE_NAME}/log.txt
bedtools merge -i $TMP2 > ${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_uncovered.bed
rm $TMP1 $TMP2

# Create masked ref
echo "Creating masked ref..." >> ${RESULTS}/${SAMPLE_NAME}/log.txt
bedtools maskfasta -fi $REFERENCE -bed ${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_uncovered.bed -fo /tmp/masked_ref.fna
bwa index /tmp/masked_ref.fna  >> ${RESULTS}/${SAMPLE_NAME}/log.txt   2>&1
samtools faidx /tmp/masked_ref.fna  >> ${RESULTS}/${SAMPLE_NAME}/log.txt   2>&1

#Llamado de variantes
echo "Calling variants..." >> ${RESULTS}/${SAMPLE_NAME}/log.txt
IN1=$OUT2
OUT1=${RESULTS}/${SAMPLE_NAME}/mpileup.vcf.gz
OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_call.vcf.gz
OUT3=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_call_consensus.vcf.gz
OUT4=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_call_Ns.vcf.gz
bcftools mpileup -f $REFERENCE $IN1  -Oz > $OUT1 2>>${RESULTS}/${SAMPLE_NAME}/log.txt
bcftools call --ploidy 1 -mv -Oz $OUT1 -o $OUT2 2>>${RESULTS}/${SAMPLE_NAME}/log.txt
bcftools view  --output-type b -i  'DP>9' $OUT2 > $OUT3 2>>${RESULTS}/${SAMPLE_NAME}/log.txt
bcftools view  --output-type b -i  'DP<10' $OUT2 > $OUT4 2>>${RESULTS}/${SAMPLE_NAME}/log.txt # guarda variantes con baja cov
rm $OUT1

#Obtención de secuencia consenso
echo "Getting consensus sequence..." >> ${RESULTS}/${SAMPLE_NAME}/log.txt
IN1=$OUT3
OUT1=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_consensus.fasta
bcftools index -f $IN1 2>>${RESULTS}/${SAMPLE_NAME}/log.txt
bcftools consensus -f /tmp/masked_ref.fna -H A -o ${OUT1}.tmp1 $IN1 2>>${RESULTS}/${SAMPLE_NAME}/log.txt
bedtools maskfasta -fi ${OUT1}.tmp1 -bed $OUT4 -fo ${OUT1}.tmp2 # mask variants with low DP

result=$(python3 <<EOF
import Bio.SeqIO as bpio
r = bpio.read('${OUT1}.tmp2','fasta')
r.description = ''
r.name = ''
r.id = 'hCoV-19/Argentina/${SAMPLE_NAME}/2020'
bpio.write(r,'${OUT1}','fasta')
EOF
)
rm $IN1 ${IN1}.csi ${OUT1}.tmp1 ${OUT1}.tmp2

#Assembly
echo "Getting assembly..." >> ${RESULTS}/${SAMPLE_NAME}/log.txt
IN1=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_1_repaired.fastq
IN2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_2_repaired.fastq
OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_assembly
megahit -1 $IN1 -2 $IN2 -o $OUT2 2>>${RESULTS}/${SAMPLE_NAME}/log.txt
python3 -c "import Bio.SeqIO as bpio;bpio.write([x for x in bpio.parse('${OUT2}/final.contigs.fa','fasta') if len(x) > 300  ],'${OUT2}.fna','fasta') "
rm -r $OUT2 # $IN1 $IN1

#Assembly vs ref
IN=${OUT2}.fna
TMP1=${RESULTS}/${SAMPLE_NAME}/nucdifftmp
nucdiff --vcf yes  $REFERENCE $IN $TMP1 $SAMPLE_NAME >> ${RESULTS}/${SAMPLE_NAME}/log.txt   2>&1
mv ${TMP1}/results ${RESULTS}/${SAMPLE_NAME}/nucdiff_denovo_ref
rm -r $TMP1

#Assembly vs consensus
IN=${OUT2}.fna
IN2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_consensus.fasta
TMP1=${RESULTS}/${SAMPLE_NAME}/nucdifftmp
nucdiff --vcf yes  $IN2 $IN $TMP1 $SAMPLE_NAME >> ${RESULTS}/${SAMPLE_NAME}/log.txt   2>&1
mv ${TMP1}/results ${RESULTS}/${SAMPLE_NAME}/nucdiff_denovo_consensus
rm -r $TMP1

echo "DONE!" >> ${RESULTS}/${SAMPLE_NAME}/log.txt
