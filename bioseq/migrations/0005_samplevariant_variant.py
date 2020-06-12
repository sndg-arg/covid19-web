# Generated by Django 3.0.6 on 2020-06-10 21:21

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('bioseq', '0004_auto_20200523_1818'),
    ]

    operations = [
        migrations.CreateModel(
            name='Variant',
            fields=[
                ('variant_id', models.AutoField(primary_key=True, serialize=False)),
                ('pos', models.PositiveSmallIntegerField(default=1, null=True)),
                ('ref', models.CharField(max_length=10)),
                ('bioentry_id', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='variants', to='bioseq.Bioentry')),
            ],
            options={
                'db_table': 'variant',
                'managed': True,
                'unique_together': {('bioentry_id', 'pos', 'ref')},
            },
        ),
        migrations.CreateModel(
            name='SampleVariant',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=255)),
                ('alt', models.CharField(max_length=10)),
                ('variant_id', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='samples', to='bioseq.Variant')),
            ],
            options={
                'db_table': 'sample_variant',
                'managed': True,
                'unique_together': {('variant_id', 'name', 'alt')},
            },
        ),
    ]
