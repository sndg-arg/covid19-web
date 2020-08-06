FROM ubuntu:18.04

RUN apt update && apt install -y wget curl nano build-essential git  libpng-dev libfreetype6 gcc python3-dev python-dev
WORKDIR /nextstrain
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b

RUN git clone https://github.com/davetang/sars_cov_2.git

ENV PATH=/miniconda/bin:${PATH}



RUN curl http://data.nextstrain.org/nextstrain.yml --compressed -o nextstrain.yml
RUN conda env create -f nextstrain.yml
RUN conda init bash
#RUN source ~/anaconda3/etc/profile.d/conda.sh
#RUN conda activate nextstrain

#RUN curl https://raw.githubusercontent.com/creationix/nvm/v0.35.0/install.sh | bash
#ENV PATH=/root/.nvm:${PATH}
#ENV NVM_DIR="/root/.nvm"
#RUN bash -c "source $NVM_DIR/nvm.sh && nvm install v0.12.18"
#RUN nvm use v10.10.0
RUN apt install -y npm nodejs
RUN npm install --global auspice

#https://github.com/conda/conda/issues/9728
RUN  curl -H 'Cache-Control: no-cache' https://repo.anaconda.com/pkgs/r/noarch/repodata.json && \
     cd sars_cov_2 &&  conda env create --file environment.yml

RUN wget https://iweb.dl.sourceforge.net/project/bbmap/BBMap_38.80.tar.gz && tar xfv BBMap_38.80.tar.gz

ENV PATH=/nextstrain/bbmap/:${PATH}

RUN wget https://vorboss.dl.sourceforge.net/project/deconseq/standalone/deconseq-standalone-0.4.3.tar.gz
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.1.6.0/gatk-4.1.6.0.zip && unzip gatk-4.1.6.0.zip
#conda activate sars_cov_2
#RUN conda activate nextstrain
# SRR10948550 SRR11092056 SRR11092057
# wget -O 1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/056/SRR11092056/SRR11092056_1.fastq.gz

#dedupe.sh in=t1.fastq.gz out=dt1.fastq.gz


#wget -O sarscov2.fasta  "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=1798174254&extrafeat=null&conwithfeat=on&hide-cdd=on"
#bwa index -p sarscov2 -a bwtsw sarscov2.fasta >bwa.log 2>&1
#bcftools mpileup -f sarscov2.fasta aln.bam > pileup.vcf
#bcftools call -mv -Oz pileup.vcf -o variants.vcf.gz
# bcftools index variants.vcf.gz
# echo "$SAMPLE" > seq.fna
#cat sarscov2.fasta | bcftools consensus variants.vcf.gz | grep -v '>' >> seq.fna

#https://samtools.github.io/bcftools/howtos/consensus-sequence.html

