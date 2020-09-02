# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.db import models
from django.shortcuts import reverse
from collections import defaultdict

"""
ALTER TABLE taxon_name ADD COLUMN id INT AUTO_INCREMENT PRIMARY KEY;
ALTER TABLE bioentry_qualifier_value  ADD COLUMN id INT AUTO_INCREMENT PRIMARY KEY;
ALTER TABLE term_synonym  ADD COLUMN id INT AUTO_INCREMENT PRIMARY KEY;
ALTER TABLE term ADD version INT UNSIGNED DEFAULT 1;

"""

# class Comment(models.Model):
#     comment_id = models.AutoField(primary_key=True)
#     bioentry = models.ForeignKey(Bioentry, models.CASCADE)
#     comment_text = models.TextField()
#     rank = models.SmallIntegerField(default=1, null=True)
#
#     class Meta:
#         managed = True
#         db_table = 'comment'
#         unique_together = (('bioentry', 'rank'),)
#
# class Reference(models.Model):
#     reference_id = models.AutoField(primary_key=True)
#     dbxref = models.OneToOneField(Dbxref, models.DO_NOTHING, unique=True, blank=True, null=True)
#     location = models.TextField()
#     title = models.TextField(blank=True, null=True)
#     authors = models.TextField(blank=True, null=True)
#     crc = models.CharField(unique=True, max_length=32, blank=True, null=True)
#
#     class Meta:
#         managed = True
#         db_table = 'reference'


# class Tool(models.Model):
#     name = models.CharField(max_length=120, blank=False)
#     description = models.TextField(default="")
#     version = models.CharField(max_length=64, blank=True, default=1, null=True)
#     url = models.URLField()
#
#     def rtype(self):
#         return "tool"
#
#
# class ToolRun(models.Model):
#     tool = models.ForeignKey(Tool, models.DO_NOTHING, related_name="runs")
#     parameters = models.TextField(blank=False)
#     result_url = models.CharField(max_length=255, blank=True, null=True)
#     result_path = models.FilePathField(path="/data/runs", null=True, recursive=True, allow_folders=True)
#     created_at = models.DateTimeField(auto_now_add=True)
#     executed_at = models.DateTimeField(auto_now=True)


# class ToolRunResult(models.Model):
#     toolrunresult_id = models.AutoField(primary_key=True)
#     run = models.ForeignKey(ToolRun, models.DO_NOTHING, related_name="results")
#     description = models.CharField(max_length=255, blank=True, null=True)
#
#
# class ToolRunResultAttrs(models.Model):
#     runResult = models.ForeignKey(ToolRunResult, models.DO_NOTHING, related_name="attrs")
#     name = models.CharField(max_length=255, blank=True, null=True)

# class AligmentSimple(ToolRunResult):
#     query = models.TextField(blank=False)
#     query_seq = models.TextField(blank=False)
#     query_start = models.IntegerField()
#     query_end = models.IntegerField()
#     query_strand = models.IntegerField(default=1)


#
#
