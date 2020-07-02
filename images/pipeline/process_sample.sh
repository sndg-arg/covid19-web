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

#limpiar fastq con fastp o prinseq
echo "Cleaning FASTQ...."
IN1=$SAMPLE_R1
IN2=$SAMPLE_R2
OUT1=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_1_fastp.fastq
OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_2_fastp.fastq
echo fastp -f $TRIML -t $TRIMR -l $MINLEN -e $QPROM  --thread $CPUS -i $IN1 -I $IN2 -o $OUT1 -O $OUT2
fastp -f $TRIML -t $TRIMR -l $MINLEN -e $QPROM  --thread $CPUS -i "$IN1" -I "$IN2" -o "$OUT1" -O "$OUT2"

rm fastp.html  fastp.json

#Limpiar adaptadores
echo "Removing adapters..."
IN1=$OUT1
IN2=$OUT2
OUT1=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_1_cutadapt.fastq
OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_2_cutadapt.fastq
echo cutadapt -a file:$ADAPTERS -o "$OUT1" -p "$OUT2" "$IN1" "$IN2"
cutadapt --cores="$CPUS" -a file:${ADAPTERS} -o "$OUT1" -p "$OUT2" "$IN1" "$IN2"
rm "$IN1" "$IN2"

##remover duplicados con bbmap
echo "Removing duplicates..."
IN1=$OUT1
IN2=$OUT2
OUT1=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_1_dedupped.fastq
OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_2_dedupped.fastq

dedupe.sh threads="$CPUS" ac=f absorbmatch=t  overwrite=t in="$IN1" out="$OUT1"
dedupe.sh threads="$CPUS" ac=f absorbmatch=t  overwrite=t in="$IN2" out="$OUT2"
rm $IN1 $IN2

##reparar los fastq para que coincidan la cantidad y orden de los reads
echo "Re-pairing FASTQ..."
IN1=$OUT1
IN2=$OUT2
OUT1=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_1_repaired.fastq
OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_2_repaired.fastq
SINGLE=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_repaired_single.fastq
repair.sh threads="$CPUS" overwrite=t in1="$IN1" in2="$IN2" out1="$OUT1" out2="$OUT2" outsingle="$SINGLE"
rm "$IN1" "$IN2" "$SINGLE"

fastqc "$OUT1" "$OUT2" -o "${RESULTS}/${SAMPLE_NAME}/" -q


#Mapeo
echo "Mapping..."
IN1=$OUT1
IN2=$OUT2
OUT1=${RESULTS}/${SAMPLE_NAME}/aln.sam
OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_aln.bam
STAT=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_flagstat.txt
echo bwa mem -t "$CPUS" "$REFERENCE" "$IN1" "$IN2"
bwa mem -t "$CPUS" "$REFERENCE" "$IN1" "$IN2" > "$OUT1"
samtools sort "$OUT1" -o "$OUT2"
samtools index "$OUT2"
samtools flagstat "$OUT2" > "$STAT"
rm "$OUT1"

#low coverage
IN=$OUT2
TMP1=${RESULTS}/${SAMPLE_NAME}/my.sorted.bedgraph
TMP2=${RESULTS}/${SAMPLE_NAME}/uncovered.intevals.bedgraph
bedtools genomecov -ibam $IN -bga > "$TMP1"
awk '$4<10' "$TMP1" > "$TMP2"
bedtools merge -i "$TMP2" > "${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_uncovered.bed"
rm "$TMP1" "$TMP2"

#Llamado de variantes
OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_aln.bam
echo "Calling variants..."
IN1=$OUT2
OUT1=${RESULTS}/${SAMPLE_NAME}/mpileup.vcf.gz
OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_call.vcf.gz
OUT3=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_call_consensus_tmp.vcf.gz
OUT4=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_call_Ns.vcf.gz
bcftools mpileup -d 1000 -E -f "$REFERENCE" "$IN1" -Oz > "$OUT1"
bcftools call -mv -Oz "$OUT1" -o "$OUT2"
bcftools view  --output-type z -i  'DP>9' "$OUT2" > "$OUT3"
bcftools view  --output-type z -i  'DP<10' "$OUT2" > "$OUT4"  # guarda variantes con baja cov
rm "$OUT1"


# Si cov es baja y hay una delecion se fuerza que la variante sea homocigota
IN=$OUT3
OUT1=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_call_consensus.vcf
OUT2=${OUT1}.gz
OUT3=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_uncovered2.bed
bcftools view "$IN" > "${OUT1}"
bedtools intersect -a "$IN" -b "${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_uncovered.bed" | grep INDEL  |  sed 's|0/1|1/1|' >> "${OUT1}"
echo "bcftools view  -e 'TYPE=\"indel\" && GT=\"het\"' ${OUT1} |  bcftools norm -f $REFERENCE- | bcftools sort  -Oz -o $OUT2 -"
bcftools view  -e 'TYPE="indel" && GT="het"' "${OUT1}" |  bcftools norm -f "$REFERENCE" - | bcftools sort  -Oz -o "$OUT2" -
bcftools index "$OUT2"

zgrep "^#" "$OUT2" > "${RESULTS}/${SAMPLE_NAME}/low_cov_indels.vcf"
zgrep INDEL "$OUT2"  | grep "1/1:" >> "${RESULTS}/${SAMPLE_NAME}/low_cov_indels.vcf"
bedtools subtract -a "${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_uncovered.bed"  -b "${RESULTS}/${SAMPLE_NAME}/low_cov_indels.vcf" > "$OUT3"

rm "${OUT1}"

#Obtención de secuencia consenso
echo "Getting consensus sequence..."
IN1=$OUT2
IN2=$OUT3 # bed
OUT1=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_consensus.fasta
bcftools index -f "$IN1"
bcftools consensus -m ${IN2} -f "$REFERENCE" -H A -o "${OUT1}.tmp1" "$IN1"
bcftools consensus -f "$REFERENCE" -H A -o "${OUT1}.unmasked" "$IN1"

result=$(python3 <<EOF
import Bio.SeqIO as bpio
r = bpio.read('${OUT1}.tmp1','fasta')
r.description = ''
r.name = ''
r.id = 'hCoV-19/Argentina/${SAMPLE_NAME}/2020'
bpio.write(r,'${OUT1}','fasta')
r = bpio.read('${OUT1}.unmasked','fasta')
r.description = ''
r.name = ''
r.id = 'hCoV-19/Argentina/${SAMPLE_NAME}/2020'
bpio.write(r,'${OUT1}.unmasked','fasta')
EOF
)
rm "${IN1}.csi" "${OUT1}.tmp1" "${OUT1}.tmp2" "${OUT4}" "$IN2" "${IN2}.csi"

#Assembly
echo "Getting assembly..."
IN1=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_1_repaired.fastq
IN2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_2_repaired.fastq
OUT2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_assembly
megahit -1 "$IN1" -2 "$IN2" -o "$OUT2"

result=$(python3 <<EOF
import Bio.SeqIO as bpio
with open('${OUT2}.fna','w') as h:
  for i,r  in enumerate(bpio.parse('${OUT2}/final.contigs.fa','fasta')):
      if len(r.seq) > 300:
#          r.description = r.des
          r.name = ''
          r.id = 'hCoV-19/Argentina/${SAMPLE_NAME}/2020_' + str(i)
          bpio.write(r,h,'fasta')
EOF
)

rm -r "$OUT2" # $IN1 $IN1

#Assembly vs ref
IN=${OUT2}.fna
TMP1=${RESULTS}/${SAMPLE_NAME}/nucdifftmp
nucdiff --vcf yes  "$REFERENCE" "$IN" "$TMP1" "$SAMPLE_NAME"
mv "${TMP1}/results" "${RESULTS}/${SAMPLE_NAME}/nucdiff_denovo_ref"
rm -r "$TMP1"

#Assembly vs consensus
IN=${OUT2}.fna
IN2=${RESULTS}/${SAMPLE_NAME}/${SAMPLE_NAME}_consensus.fasta
TMP1=${RESULTS}/${SAMPLE_NAME}/nucdifftmp
nucdiff --vcf yes  "$IN2" "$IN" "$TMP1" "$SAMPLE_NAME"
mv "${TMP1}/results" "${RESULTS}/${SAMPLE_NAME}/nucdiff_denovo_consensus"
rm -r "$TMP1"

echo "DONE!"
