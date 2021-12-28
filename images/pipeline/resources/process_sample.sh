#!/bin/bash

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
DEDUP="${DEDUP:-false}"

#*************************************************

#IFS='/' tokens=( $SAMPLE_R1 )
#SAMPLE_NAME=${tokens[1]}
DIR=${RESULTS}/${SAMPLE_NAME}
IFS=' '
[ -d "$DIR" ] || mkdir -p "$DIR"

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

echo "CPUS: $CPUS"

### Clean #####################################33
if [ ! -f "${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_1_filtered.fastq" ]; then

if [[ "${SAMPLE_R1}" == *.gz ]]
then
  zcat $SAMPLE_R1 > /tmp/R1.fastq
  zcat $SAMPLE_R2 > /tmp/R2.fastq
else
  cp $SAMPLE_R1 /tmp/R1.fastq
  cp $SAMPLE_R2 /tmp/R2.fastq
fi

if $DEDUP
then
##remover duplicados
echo "Removing duplicates..."
OUT1=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_1_dedupped.fastq
OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_2_dedupped.fastq
IN1=/tmp/R1.fastq
IN2=/tmp/R2.fastq
echo "$IN1" > /tmp/fastuniq.txt
echo "$IN2" >> /tmp/fastuniq.txt
echo fastuniq -o $OUT1 -p $OUT2 -i /tmp/fastuniq.txt
fastuniq -o $OUT1 -p $OUT2 -i /tmp/fastuniq.txt
tot=$(wc -l $IN1| cut -f1 -d' ')
uniq=$(wc -l $OUT1| cut -f1 -d' ')

python -c "print( 'Duplicates removed: ' + str( ( ${tot} - ${uniq} )*1.0/${tot}   ))"

rm $IN1 $IN2
else
  OUT1=/tmp/R1.fastq
  OUT2=/tmp/R2.fastq
fi


#Limpiar adaptadores
#echo "Removing adapters..."
#IN1=$OUT1
#IN2=$OUT2
#OUT1=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_1_cutadapt.fastq
#OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_2_cutadapt.fastq
#echo cutadapt -a file:$ADAPTERS -o "$OUT1" -p "$OUT2" "$IN1" "$IN2"
#cutadapt --cores="$CPUS" -a file:${ADAPTERS} -o "$OUT1" -p "$OUT2" "$IN1" "$IN2" 1>/dev/null
#rm "$IN1" "$IN2"

#limpiar fastq con fastp o prinseq
echo "Cleaning FASTQ...."
IN1=$OUT1
IN2=$OUT2
OUT1=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_1_filtered.fastq
OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_2_filtered.fastq
echo fastp -f $TRIML -t $TRIMR -l $MINLEN -e $QPROM  --thread $CPUS -i $IN1 -I $IN2 -o $OUT1 -O $OUT2  --adapter_fasta ${ADAPTERS}
fastp -f $TRIML -t $TRIMR -l $MINLEN -e $QPROM  --thread $CPUS -i "$IN1" -I "$IN2" -o "$OUT1" -O "$OUT2" --adapter_fasta "$ADAPTERS"



fastqc "$OUT1" "$OUT2" -o "${RESULTS}/${SAMPLE_NAME}/" -q

rm fastp.html  fastp.json "$IN1" "$IN2" ${RESULTS}/${SAMPLE_NAME}/*_fastqc.zip




fi
################################

###########Mapping#############################
if [ ! -f "${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_call_consensus.vcf.gz" ]; then

echo "Mapping..."
IN1=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_1_filtered.fastq
IN2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_2_filtered.fastq
OUT1=${RESULTS}/${SAMPLE_NAME}/aln.sam
OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_aln.bam
echo bwa mem -t "$CPUS" "$REFERENCE" "$IN1" "$IN2"
bwa mem  -R "@RG\tID:group1\tSM:${SAMPLE_NAME}\tPL:illumina\\tLB:COVID19" -t "$CPUS" "$REFERENCE" "$IN1" "$IN2" > "$OUT1"
samtools sort "$OUT1" -o "$OUT2"
samtools index "$OUT2"
rm "$OUT1"

echo "Calling variants..."
IN1=$OUT2
OUT1=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_call.gvcf.gz
OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_call.vcf.gz

/app/gatk/gatk HaplotypeCaller -ERC GVCF -R "$REFERENCE" -ploidy 2 -I "$IN1" --output-mode EMIT_ALL_CONFIDENT_SITES -O "$OUT1"
/app/gatk/gatk GenotypeGVCFs -R "$REFERENCE" -V "$OUT1" -O "$OUT2"

fi

#low coverage & report
IN=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_aln.bam
TMP1=${RESULTS}/${SAMPLE_NAME}/my.sorted.bedgraph
#pileup.sh in=$IN   out=${RESULTS}/${SAMPLE_NAME}/aln_stats.txt  overwrite=true 2> ${RESULTS}/${SAMPLE_NAME}/cov_med.txt
bedtools genomecov -ibam $IN -bga > "$TMP1"
awk '$4<10' "$TMP1" > "${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_uncovered.bed"
rm "$TMP1"

echo "filter variants"
IN=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_call.vcf.gz
OUT=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_call_consensus.vcf.gz
TMP1="/tmp/unspanning.vcf"
TMP2="/tmp/variants.vcf"
#  exclude variants with more than 2 alleles "INFO/AN>2"
bcftools filter  -e '(alt="*") || (INFO/AN>2)' $IN > $TMP1
bcftools filter  -i '(FORMAT/DP>9) && (GT="hom") && (FORMAT/AD[0:1]>FORMAT/AD[0:0])' $IN | bcftools filter -e "INFO/AN>2"  > $TMP2
#Remove spannig deletions
bcftools filter -i 'alt="*"' $IN  | bcftools norm -m -any | \
         bcftools filter -e 'alt="*"'  | bcftools filter -i 'FORMAT/AD[0:1] > 9' | \
         sed  's|0/1:|1/1:|'  | sed  's|0\|1:|1/1:|' | bcftools view -H   >> $TMP2
# VC exceptions
## HET prop > 60% to homo
bcftools filter -i '(GT="het") && (FORMAT/DP>9) && (((FORMAT/AD[0:1]) / FORMAT/DP) >0.6)' $TMP1 | \
  sed  's|0/1:|1/1:|'  | sed  's|0\|1:|1/1:|' | bcftools view -H   >> $TMP2
## remove HET prop < 60% => otherwise they go to the consensus, even if GT=het
#bcftools filter -i ' ( STRLEN(ALT)==STRLEN(REF)) && (GT="het") && (FORMAT/DP>9) && ( ((FORMAT/AD[0:1]) / (FORMAT/DP) ) <0.6) && ((FORMAT/AD[0:0]/FORMAT/DP)<0.6)' $TMP1 | \
#  bcftools view -H   >> $TMP2
### Indel & COV < 10 -> ya se cubren en el primer filtro
#bcftools filter -i "(STRLEN(REF)>1) && ( STRLEN(ALT)!=STRLEN(REF)) && (FORMAT/DP>9) && (FORMAT/AD[0:1]>FORMAT/AD[0:0])" $TMP1 \
#  | bcftools view -H >> $TMP2

bcftools sort -O z $TMP2 > $OUT
bcftools index "$OUT"
rm $IN ${IN}.tbi

# Remove masking regions that overlaps with deletions with good coverare
IN1=$OUT
IN2="${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_uncovered.bed"
OUT2="${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_uncovered2.bed"
bcftools filter -i 'STRLEN(REF)>STRLEN(ALT)' "$IN1" | bedtools subtract -a "${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_uncovered.bed"  -b - > "$OUT2"
rm $IN2

echo "Getting consensus sequence..."
IN1=$OUT
IN2=$OUT2
OUT1="${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_consensus.fasta"
bcftools index -f "$IN1"
bcftools consensus -m "${IN2}" -f "$REFERENCE"  -o "${OUT1}.tmp1" "$IN1"


result=$(python3 <<EOF
import Bio.SeqIO as bpio
from datetime import date
r = bpio.read('${OUT1}.tmp1','fasta')
r.description = ''
r.name = ''
r.id = 'hCoV-19/Argentina/${SAMPLE_NAME}/' + str(date.today().year)
bpio.write(r,'${OUT1}','fasta')
EOF
)
rm "${OUT1}.tmp1"



##### Denovo Assembly #################3

if [ ! -f "${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_assembly.fna" ]; then


echo "Getting assembly..."
IN1=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_1_filtered.fastq
IN2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_2_filtered.fastq
OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_assembly
megahit -1 "$IN1" -2 "$IN2" -o "$OUT2"

result=$(python3 <<EOF
from datetime import date
import Bio.SeqIO as bpio
with open('${OUT2}.fna','w') as h:
  for i,r  in enumerate(bpio.parse('${OUT2}/final.contigs.fa','fasta')):
      if len(r.seq) > 300:
#          r.description = r.des
          r.name = ''
          r.id = 'hCoV-19/Argentina/${SAMPLE_NAME}/' + str(date.today().year + '_' + str(i)
          bpio.write(r,h,'fasta')
EOF
)

rm -r "$OUT2" # $IN1 $IN1
fi


#Assembly vs ref
IN="${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_assembly.fna"
TMP1="${RESULTS}/${SAMPLE_NAME}/nucdifftmp"
nucdiff --vcf yes  "$REFERENCE" "$IN" "$TMP1" "$SAMPLE_NAME"
mv "${TMP1}/results" "${RESULTS}/${SAMPLE_NAME}/nucdiff_denovo_ref"
rm -r "$TMP1"

#Assembly vs consensus
IN="${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_assembly.fna"
IN2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_consensus.fasta
TMP1=${RESULTS}/${SAMPLE_NAME}/nucdifftmp
nucdiff --vcf yes  "$IN2" "$IN" "$TMP1" "$SAMPLE_NAME"
mv "${TMP1}/results" "${RESULTS}/${SAMPLE_NAME}/nucdiff_denovo_consensus"
rm -r "$TMP1"

echo "DONE!"
