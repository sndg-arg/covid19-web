from django.contrib.auth import get_user_model

from config import celery_app

User = get_user_model()

import datetime
from collections import defaultdict

import Bio.SeqIO as bpio
import matplotlib.pyplot as plt
import pandas as pd

from config.settings.base import STATICFILES_DIRS

from math import ceil

from sndg_covid19.views import latam_countries
from bioseq.io.MSAMap import MSAMap

fechas = ["", "Enero", "Febrero", "Marzo", "Abril", "Mayo", "Junio", "Julio"]

from sndg_covid19.io import country_from_gisaid

@celery_app.task()
def variant_graphics(gene: str, pos: int, msamap=None):
    """
    mafft --keeplength --mapout --addfull orf1ab.faa sndg_covid19/static/rawORFs/orf1ab_prot.fasta  > sndg_covid19/static/ORFs/orf1ab_prot.fasta
    orf1ab -> merge of orf1a and orf1b
    :param gene:
    :param pos:
    :return:
    """
    fig_path = f'{STATICFILES_DIRS[0]}/auto/posfigs/{gene}{pos}.png'
    msa_map = {
        "E": "E_prot.fasta",
        "M": "M_prot.fasta",
        "N": "N_prot.fasta",
        "orf1ab": "orf1ab_prot.fasta",
        "NS3": "orf3a_prot.fasta",
        "NS6": "orf6_prot.fasta",
        "NS7a": "orf7a_prot.fasta",
        "NS7b": "orf7b_prot.fasta",
        "NS8": "orf8_prot.fasta",
        "S": "S_prot.fasta"
    }
    msa_map = {k: f'{STATICFILES_DIRS[0]}/ORFs/{v}' for k, v in msa_map.items()}

    msa_file = msa_map[gene]

    msa = bpio.to_dict(bpio.parse(msa_file, "fasta"))

    if not msamap:
        msamap = MSAMap(msa)
        msamap.init()

    aln_pos = msamap.pos_seq_msa_map[gene][pos]

    prom = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))
    latam = []
    for k, v in msa.items():
        if k == gene:
            continue

        sample_date = datetime.datetime.strptime(k.split("|")[-1], '%Y-%m-%d').date()
        month = sample_date.month
        country = country_from_gisaid(k)

        if country in latam_countries:
            latam.append(country)
            prom[month][country][v.seq[aln_pos]] += 1
            # print(v.seq[aln_pos-1:aln_pos+2])

    latam = sorted(list(set(latam)))
    final = []
    for month in prom:
        for country in latam:
            for aa in prom[month][country]:
                r = {"month": month, "country": country, "aa": aa, "count": prom[month][country][aa]}
                final.append(r)

    df = pd.DataFrame(final)

    if len(df) == 0:
        return

    df = df[df.aa != "X"]

    # df = pd.DataFrame(final, columns=["month"] + latam)
    # d1 = df.set_index(['month']).sort_index()
    with_data = list(df.country.unique())
    fig, axs = plt.subplots(ceil(len(with_data) / 2), 2, sharex=True, figsize=(10, 10))
    for idx, country in enumerate(sorted(with_data)):
        dfp = df[df.country == country]
        dfp["prop"] = [x["count"] / sum(list(dfp[dfp.month == x.month]["count"])) * 100 for _, x in dfp.iterrows()]
        aas = list(dfp.aa.unique())
        for month in range(2, 6):
            if month not in list(dfp.month):
                for aa in aas:
                    dfp = dfp.append({"aa": aa, "count": 0, "month": month}, ignore_index=True)

        dfp["month"] = [fechas[f] for f in dfp.month]
        dfp = pd.pivot(dfp, index="month", columns="aa", values="prop")
        dfp.index = [x for i, x in enumerate(fechas[1:]) if i < len(dfp.index.unique())]
        ax = dfp.plot(kind='bar', title=country, stacked=True, ax=axs.flat[idx])
        ax.set_ylabel("%")
        ax.set_xlabel("Mes")
        ax.legend(title="")
    plt.savefig(fig_path)
