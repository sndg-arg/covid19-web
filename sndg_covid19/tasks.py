import traceback

from celery import shared_task
from django.contrib.auth import get_user_model

from config import celery_app
from sndg_covid19.bioio.JobValidationError import JobValidationError
from .bioio.EncodingUtils import EncodingUtils

User = get_user_model()

import datetime
from collections import defaultdict

import Bio.SeqIO as bpio
import matplotlib.pyplot as plt
import pandas as pd

from math import ceil
from sndg_covid19.views import latam_countries
from bioseq.bioio.MSAMap import MSAMap

fechas = ["", "Enero", "Febrero", "Marzo", "Abril", "Mayo", "Junio", "Julio", "Agosto", "Septiembre", "Octubre",
          "Noviembre", "Diciembre"]

from sndg_covid19.models import ImportJob
from sndg_covid19.bioio import country_from_gisaid
from sndg_covid19.bioio.CovidIO import CovidIO
import numpy as np
from dateutil.relativedelta import relativedelta
import sys

def autolabel(ax, y):
    for i, v in enumerate(y):
        txt = " ".join([str(int(x)) for x in v if not np.isnan(x)])
        ax.text(i - 0.1, 50, txt, color='black', fontweight='bold')


@shared_task
def process_msa(alnjobid: int):
    job = ImportJob.objects.get(import_job_id=alnjobid)
    if job.status in ["processing","finished", "error"]:
        return
    try:
        CovidIO.process_import_job(job)
        job.status_desc = ""
        job.status = "finished"
    except JobValidationError as ex:
        job.status = "error"
        job.status_desc = ex.status_text
        job.errors = "\n".join(ex.errors)
    except Exception:
        job.debug_status_desc = traceback.format_exc()
        job.status = "error"
        job.status_desc = "Error desconocido"
        sys.err.write(job.debug_status_desc)
    job.save()


@celery_app.task()
def variant_graphics(gene: str, pos: int, fig_path, msa_file, msamap=None, idx_date=-2):
    """
    mafft --keeplength --mapout --addfull orf1ab.faa sndg_covid19/static/rawORFs/orf1ab_prot.fasta  > sndg_covid19/static/ORFs/orf1ab_prot.fasta
    orf1ab -> merge of orf1a and orf1b
    :param gene:
    :param pos:
    :return:
    """

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

        sample_date = datetime.datetime.strptime(k.split("|")[idx_date], '%Y-%m-%d').date()
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
    cmap = plt.cm.coolwarm

    handles = []
    color_map = {}
    aas = df.aa.unique()
    for idx2, aa in enumerate(aas):
        color = cmap(1 * idx2 / len(aas))
        color_map[aa] = color
        # handle = Patch(facecolor=color, edgecolor=color, label=aa)
        # handles.append(handle)

    # df = pd.DataFrame(final, columns=["month"] + latam)
    # d1 = df.set_index(['month']).sort_index()
    with_data = list(df.country.unique())
    fig, axs = plt.subplots(ceil(len(with_data) / 2), 2, sharex=True, figsize=(10, 10))
    dates = list(get_last_months(datetime.datetime.now(), 11))
    max_month = dates[-1][1] + dates[-1][0] * 100
    min_month = dates[0][1] + dates[0][0] * 100
    for idx, country in enumerate(sorted(with_data)):
        dfp = df[df.country == country]
        dfp["prop"] = [x["count"] / sum(list(dfp[dfp.month == x.month]["count"])) * 100 for _, x in dfp.iterrows()]
        aas = list(dfp.aa.unique())
        for month in range(min_month, max_month):
            if month not in list(dfp.month):
                for aa in aas:
                    dfp = dfp.append({"aa": aa, "count": 0, "month": month}, ignore_index=True)

        dfp["month"] = [fechas[f] for f in list(dfp.month)]
        dfp.country.fillna(country, inplace=True)
        dfp.prop.fillna(0, inplace=True)
        dfp2 = pd.pivot(dfp, index="month", columns="aa", values="prop")
        dfp2 = dfp2.reindex(sorted(dfp2.index, key=lambda x: fechas.index(x)))

        ax = dfp2.plot(kind='bar', title=country, stacked=True, ax=axs.flat[idx],
                       color=[color_map[x] for x in [aa for aa in aas if aa in list(dfp.aa.unique())]])
        ax.set_ylabel("%")
        ax.set_xlabel("Mes")
        # ax.legend(handles=handles, title="")

        dfp2 = pd.pivot(dfp, index="month", columns="aa", values="count")

        dfp2 = dfp2.reindex(sorted(dfp2.index, key=lambda x: fechas.index(x)))

        autolabel(ax, list(dfp2.values))
    plt.savefig(fig_path)


def get_last_months(start_date, months):
    for i in range(months):
        yield (start_date.year, start_date.month)
        start_date += relativedelta(months=-1)
