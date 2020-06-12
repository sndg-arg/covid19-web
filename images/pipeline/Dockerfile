FROM ubuntu:18.04
RUN apt-get update && apt-get install -y  openjdk-8-jdk libncurses5-dev libbz2-dev lzma liblzma-dev zlib1g-dev \
    cmake libgtk2.0-dev ncbi-blast+ libeigen3-dev unzip wget nano python3-pip

RUN python3 -m pip install --upgrade cutadapt

WORKDIR /app

RUN wget -O bwa-0.7.17.tar.bz2 https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2/download \
    && tar xfv bwa-0.7.17.tar.bz2 && cd bwa-0.7.17 && make && ln -s /app/bwa-0.7.17/bwa /usr/local/bin

RUN cd /app/ && wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
    tar xfv samtools-1.10.tar.bz2 && rm samtools-1.10.tar.bz2 && cd samtools-1.10 && ./configure && make && make install
RUN cd /app/ && wget https://github.com/samtools/bcftools/releases/download/1.10/bcftools-1.10.tar.bz2 && \
    tar xfv bcftools-1.10.tar.bz2 && rm bcftools-1.10.tar.bz2 && cd bcftools-1.10 && ./configure && make && make install
RUN cd /app/ && wget https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools-2.29.2.tar.gz && \
    tar xfv bedtools-2.29.2.tar.gz && rm bedtools-2.29.2.tar.gz && cd bedtools2 && make && make install
RUN cd /app/ && wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && rm fastqc_v0.11.9.zip && \
    chmod 755 /app/FastQC/fastqc && ln -s /app/FastQC/fastqc /usr/local/bin


RUN cd /app/ && wget http://cab.spbu.ru/files/release3.14.1/SPAdes-3.14.1-Linux.tar.gz && \
    tar xfv SPAdes-3.14.1-Linux.tar.gz && rm SPAdes-3.14.1-Linux.tar.gz && \
    ln -s /app/SPAdes-3.14.1-Linux/bin/spades.py /usr/local/bin

RUN mkdir fastp && cd fastp && wget http://opengene.org/fastp/fastp && chmod a+x ./fastp && \
    ln -s /app/fastp/fastp /usr/local/bin
RUN wget https://github.com/voutcn/megahit/releases/download/v1.2.9/MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz && \
    tar zvxf MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz && rm MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz

RUN wget https://ufpr.dl.sourceforge.net/project/bbmap/BBMap_38.84.tar.gz && tar xfv BBMap_38.84.tar.gz && \
    rm BBMap_38.84.tar.gz

#RUN cd /app/ && wget https://downloads.sourceforge.net/project/quast/quast-5.0.2.tar.gz && \
#    tar xfv quast-5.0.2.tar.gz && cd quast-5.0.2 && ./install.sh && ln -s /app/quast-5.0.2/quast.py /usr/local/bin
#RUN wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip && unzip qualimap_v2.2.1.zip && \
#    rm unzip qualimap_v2.2.1.zip

RUN apt install -y mummer
RUN pip3 install NucDiff

WORKDIR /ref

RUN  wget -O /ref/MN996528.fna "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=1802633808&conwithfeat=on&withparts=on&hide-cdd=on"
RUN makeblastdb -dbtype nucl -in /ref/MN996528.fna && bwa index MN996528.fna && samtools faidx MN996528.fna && \
    samtools dict -o MN996528.dict MN996528.fna

RUN pip3 install tqdm termcolor

ENV PATH $PATH:/app/MEGAHIT-1.2.9-Linux-x86_64-static/bin/:/app/bbmap/:/app/script/
COPY process_sample.sh /app/script/
COPY process_batch.py /app/script/
COPY adapters.fasta /ref/

RUN chmod +x  /app/script/process_sample.sh  /app/script/process_batch.py

RUN groupadd covid && useradd -r -u 1000 -g covid covid
USER covid
WORKDIR /out



