# Generated by Django 3.0.6 on 2020-07-16 01:35

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('pdbdb', '0004_auto_20200715_2226'),
    ]

    operations = [
        migrations.AddField(
            model_name='residueproperty',
            name='value_text',
            field=models.CharField(max_length=255, null=True),
        ),
        migrations.AddField(
            model_name='residuesetproperty',
            name='value_text',
            field=models.CharField(max_length=255, null=True),
        ),
    ]
