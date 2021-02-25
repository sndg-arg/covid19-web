import Bio.SeqIO as bpio
from datetime import datetime
from collections import defaultdict


sizes = []
data_g = {}
for i,r in enumerate(open("fechas_genomas.csv")):
    if i == 0:
        continue
    sample,sdate,location = [x.strip() for x in r.strip().split(",")]
    sdate2 = None
    for pattern in ['%Y-%m-%d','%d/%m/%Y','%Y-%m']:
        try:
            sdate2 = datetime.strptime(sdate, pattern).date()
        except ValueError:
            continue
    assert sdate2,sdate
    data_g[sample] = {"sdate":sdate2,"location":location}

with open("fixed_genomes.fna","w") as h:
    for r in bpio.parse("./ALN-MutationTool.fasta","fasta"):
        #hCoV-19/Argentina/PAIS-A0001/2020 -> hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124|2019-12-30'
        r.name = ""
        sample = r.id.split("/")[2]

        code = r.id + "|-|" + data_g[sample]["sdate"].strftime("%Y-%m-%d") + "|SouthAmerica"
        r.id = code
        bpio.write(r,h,"fasta")
        sizes.append(len(r.seq))
print([x for x in set(sizes)])


sizes = defaultdict(list)
data_s = {}
for i,r in enumerate(open("fechas_spike.csv")):
    if i == 0:
        continue
    sample,sdate,location = [x.strip() for x in r.strip().split(",")]
    sdate2 = None
    for pattern in ['%Y-%m-%d','%d/%m/%Y','%Y-%m']:
        try:
            sdate2 = datetime.strptime(sdate, pattern).date()
        except ValueError:
            continue
    assert sdate2,sdate
    data_s[sample] = {"sdate":sdate2,"location":location}
with open("fixed_spike.fna","w") as h:
    for r in bpio.parse("./Spike-ALN-hasta20210201.fasta","fasta"):
        #PAIS-S-A0005-> >Spike|hCoV-19/Australia/VIC546/2020|2020-03-24|EPI_ISL_426809|Original|hCoV-19'
        r.name = ""
        sample = r.id

        code = f"Spike|hCoV-19/Argentina/{r.id}/{data_s[sample]['sdate'].strftime('%Y')}|{data_s[sample]['sdate'].strftime('%Y-%m-%d')}|-|Original|hCoV-19"
        r.id = code
        sizes[len(r.seq)].append(str(r.seq))
        if 2685 == len(r.seq):
            bpio.write(r,h,"fasta")


print([x for x in set(sizes)])
print([(x,len(v)) for x,v in sizes.items()])

print(""" Now add ref sequences:
mafft --thread 8 --addfull ref.fna fixed_genomes.fna > fixed_genomes_with_ref.fna
Create spike ref

covid = Bioentry.objects.get(biodatabase__name="COVID19")
f = covid.features.filter(display_name="S")[0]
covid.seq.seq[f.locations.first().start_pos:f.locations.first().end_pos]
bpio.write(SeqRecord(id="S",description=data_g[seq_id]["location"],seq=Seq(covid.seq.seq[f.locations.first().start_pos:f.locations.first().end_pos])),"s.fna","fasta")

mafft --addfull s.fna fixed_spike.fna  > fixed_spike_with_ref.fna

with open("fixed_genomes_with_ref_and_loc.fna","w") as h:
     for x in bpio.parse("fixed_genomes_with_ref.fna","fasta"):
        x.name = ""
        if "/" in x.id:
           sid = x.id.split("|")[0].split("/")[-2]
           x.description = data_g[sid]["location"]
           bpio.write(x,h,"fasta")

from Bio.Seq import Seq
ref = bpio.read("s.fna","fasta").translate()
with open("fixed_spike.faa","w") as h:
     for sseq in bpio.parse("fixed_spike_with_ref.fna","fasta"):
        if sseq.id != "S":
            seq =  []
            for i,aa in enumerate(str(sseq.seq)):
                seq.append(ref[i] if aa == "-" else aa)
            seq = "".join(seq)
            try:
                seq_id = sseq.id.split("/")[2]
            except:
                print (sseq.id)
                raise
            sseq.seq = Seq(seq).translate()
            sseq.description =  data_s[seq_id]["location"]
            bpio.write(sseq,h,"fasta")
        else:
            sseq = sseq.translate()
            sseq.id = "S"
            bpio.write(,h,"fasta")


""")
