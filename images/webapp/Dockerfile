FROM python:3.8-slim-buster

ENV PYTHONUNBUFFERED 1
ENV PYTHONDONTWRITEBYTECODE 1

RUN apt-get update \
  # dependencies for building Python packages
  && apt-get install -y build-essential \
  # psycopg2 dependencies
  && apt-get install -y libpq-dev \
  # Translations dependencies
  && apt-get install -y gettext \
  # cleaning up unused files
  && apt-get purge -y --auto-remove -o APT::AutoRemove::RecommendsImportant=false \
  && rm -rf /var/lib/apt/lists/*

RUN apt install -y ncbi-blast+ unzip
WORKDIR /
RUN apt update && apt install -y git wget nano
RUN git clone https://github.com/sndg-arg/covid19-web.git
WORKDIR /covid19-web
RUN pip install -r  requirements/local.txt




