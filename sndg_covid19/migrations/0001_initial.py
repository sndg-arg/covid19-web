# Generated by Django 3.1.4 on 2021-03-20 22:54

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='ImportJob',
            fields=[
                ('import_job_id', models.AutoField(primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=255)),
                ('aln_type', models.CharField(choices=[('genome', 'genome'), ('spike', 'spike')], max_length=20)),
                ('version', models.PositiveSmallIntegerField(default=1, null=True)),
                ('fasta', models.FileField(upload_to='uploads/')),
                ('csv', models.FileField(upload_to='uploads/')),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('user', models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='import_jobs', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'db_table': 'importjobs',
                'managed': True,
                'unique_together': {('name', 'version')},
            },
        ),
    ]