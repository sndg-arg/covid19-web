FROM sndg/covid19:org

USER root

COPY resources/process_sample.sh /app/script/
COPY resources/process_batch.py /app/script/
COPY resources/MSAMap.py /app/script/
COPY resources/list_denovo_diffs.py /app/script/
COPY resources/test_consensus.py /app/script/
COPY resources/run_depth.py /app/script/

COPY process_assembly.py /app/script/




RUN chmod +x  /app/script/process_sample.sh  /app/script/process_batch.py /app/script/process_assembly.py \
    /app/script/MSAMap.py /app/script/list_denovo_diffs.py /app/script/test_consensus.py /app/script/run_depth.py

USER covid
WORKDIR /out



