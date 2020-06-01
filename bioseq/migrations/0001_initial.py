# Generated by Django 3.0.6 on 2020-05-15 19:19

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Biodatabase',
            fields=[
                ('biodatabase_id', models.AutoField(primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=128, unique=True)),
                ('authority', models.CharField(blank=True, max_length=128, null=True)),
                ('description', models.TextField(blank=True, null=True)),
            ],
            options={
                'db_table': 'biodatabase',
                'managed': True,
            },
        ),
        migrations.CreateModel(
            name='Bioentry',
            fields=[
                ('bioentry_id', models.AutoField(primary_key=True, serialize=False)),
                ('name', models.CharField(default='', max_length=40)),
                ('accession', models.CharField(max_length=128)),
                ('identifier', models.CharField(blank=True, default='', max_length=40, null=True)),
                ('division', models.CharField(blank=True, max_length=6, null=True)),
                ('description', models.TextField(blank=True, null=True)),
                ('version', models.PositiveSmallIntegerField(default=1, null=True)),
                ('index_updated', models.BooleanField(default=False)),
                ('biodatabase', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='entries', to='bioseq.Biodatabase')),
            ],
            options={
                'db_table': 'bioentry',
                'managed': True,
            },
        ),
        migrations.CreateModel(
            name='Dbxref',
            fields=[
                ('dbxref_id', models.AutoField(primary_key=True, serialize=False)),
                ('dbname', models.CharField(max_length=40)),
                ('accession', models.CharField(max_length=128)),
                ('version', models.PositiveSmallIntegerField(default=1, null=True)),
            ],
            options={
                'db_table': 'dbxref',
                'managed': True,
                'unique_together': {('accession', 'dbname', 'version')},
            },
        ),
        migrations.CreateModel(
            name='Location',
            fields=[
                ('location_id', models.AutoField(primary_key=True, serialize=False)),
                ('start_pos', models.IntegerField(blank=True, null=True)),
                ('end_pos', models.IntegerField(blank=True, null=True)),
                ('strand', models.IntegerField(default=1)),
                ('rank', models.SmallIntegerField(default=1, null=True)),
                ('dbxref', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.DO_NOTHING, to='bioseq.Dbxref')),
            ],
            options={
                'db_table': 'location',
                'managed': True,
            },
        ),
        migrations.CreateModel(
            name='Ontology',
            fields=[
                ('ontology_id', models.AutoField(primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=32, unique=True)),
                ('definition', models.TextField(blank=True, null=True)),
            ],
            options={
                'db_table': 'ontology',
                'managed': True,
            },
        ),
        migrations.CreateModel(
            name='Seqfeature',
            fields=[
                ('seqfeature_id', models.AutoField(primary_key=True, serialize=False)),
                ('display_name', models.CharField(blank=True, max_length=64, null=True)),
                ('rank', models.PositiveSmallIntegerField(default=1, null=True)),
                ('index_updated', models.BooleanField(default=False)),
                ('bioentry', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='features', to='bioseq.Bioentry')),
            ],
            options={
                'db_table': 'seqfeature',
                'managed': True,
            },
        ),
        migrations.CreateModel(
            name='Taxon',
            fields=[
                ('taxon_id', models.AutoField(primary_key=True, serialize=False)),
                ('ncbi_taxon_id', models.IntegerField(blank=True, null=True, unique=True)),
                ('node_rank', models.CharField(blank=True, max_length=32, null=True)),
                ('genetic_code', models.PositiveIntegerField(blank=True, null=True)),
                ('mito_genetic_code', models.PositiveIntegerField(blank=True, null=True)),
                ('left_value', models.PositiveIntegerField(blank=True, null=True, unique=True)),
                ('right_value', models.PositiveIntegerField(blank=True, null=True, unique=True)),
                ('parent_taxon', models.ForeignKey(null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='children', to='bioseq.Taxon')),
            ],
            options={
                'db_table': 'taxon',
                'managed': True,
            },
        ),
        migrations.CreateModel(
            name='Term',
            fields=[
                ('term_id', models.AutoField(primary_key=True, serialize=False)),
                ('name', models.TextField(blank=True, null=True)),
                ('definition', models.TextField(blank=True, null=True)),
                ('identifier', models.CharField(blank=True, max_length=255, null=True)),
                ('is_obsolete', models.CharField(blank=True, max_length=1, null=True)),
                ('version', models.PositiveSmallIntegerField(default=1, null=True)),
                ('ontology', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, related_name='terms', to='bioseq.Ontology')),
            ],
            options={
                'db_table': 'term',
                'managed': True,
                'unique_together': {('identifier', 'ontology', 'is_obsolete')},
            },
        ),
        migrations.CreateModel(
            name='TermRelationship',
            fields=[
                ('term_relationship_id', models.AutoField(primary_key=True, serialize=False)),
                ('object_term', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, related_name='object_termrelationships', to='bioseq.Term')),
                ('ontology', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, to='bioseq.Ontology')),
                ('predicate_term', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, related_name='predicate_termrelationships', to='bioseq.Term')),
                ('subject_term', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, related_name='subject_termrelationships', to='bioseq.Term')),
            ],
            options={
                'db_table': 'term_relationship',
                'managed': True,
                'unique_together': {('subject_term', 'predicate_term', 'object_term', 'ontology')},
            },
        ),
        migrations.CreateModel(
            name='Biosequence',
            fields=[
                ('bioentry', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, primary_key=True, related_name='seq', serialize=False, to='bioseq.Bioentry')),
                ('version', models.SmallIntegerField(blank=True, default=1, null=True)),
                ('length', models.IntegerField(blank=True, null=True)),
                ('alphabet', models.CharField(blank=True, max_length=10, null=True)),
                ('seq', models.TextField(blank=True, null=True)),
            ],
            options={
                'db_table': 'biosequence',
                'managed': True,
            },
        ),
        migrations.CreateModel(
            name='DbxrefQualifierValue',
            fields=[
                ('dbxref', models.OneToOneField(on_delete=django.db.models.deletion.DO_NOTHING, primary_key=True, serialize=False, to='bioseq.Dbxref')),
                ('rank', models.SmallIntegerField(default=1, null=True)),
                ('value', models.TextField(blank=True, null=True)),
            ],
            options={
                'db_table': 'dbxref_qualifier_value',
                'managed': True,
            },
        ),
        migrations.CreateModel(
            name='LocationQualifierValue',
            fields=[
                ('location', models.OneToOneField(on_delete=django.db.models.deletion.DO_NOTHING, primary_key=True, serialize=False, to='bioseq.Location')),
                ('value', models.CharField(max_length=255)),
                ('int_value', models.IntegerField(blank=True, null=True)),
            ],
            options={
                'db_table': 'location_qualifier_value',
                'managed': True,
            },
        ),
        migrations.CreateModel(
            name='TaxIdx',
            fields=[
                ('tax', models.OneToOneField(db_column='tax_id', on_delete=django.db.models.deletion.CASCADE, primary_key=True, related_name='keywords', serialize=False, to='bioseq.Taxon')),
                ('text', models.TextField()),
                ('genus', models.CharField(default='', max_length=255)),
                ('family', models.CharField(default='', max_length=255)),
            ],
            options={
                'db_table': 'tax_idx',
                'managed': True,
            },
        ),
        migrations.CreateModel(
            name='TermIdx',
            fields=[
                ('term', models.OneToOneField(db_column='term_id', on_delete=django.db.models.deletion.CASCADE, primary_key=True, related_name='keywords', serialize=False, to='bioseq.Term')),
                ('text', models.TextField()),
            ],
            options={
                'db_table': 'term_idx',
                'managed': True,
            },
        ),
        migrations.CreateModel(
            name='SeqfeatureQualifierValue',
            fields=[
                ('seqfeature_qualifiervalue_id', models.AutoField(primary_key=True, serialize=False)),
                ('rank', models.SmallIntegerField(default=1, null=True)),
                ('value', models.TextField()),
                ('seqfeature', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='qualifiers', to='bioseq.Seqfeature')),
                ('term', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, to='bioseq.Term')),
            ],
            options={
                'db_table': 'seqfeature_qualifier_value',
                'managed': True,
            },
        ),
        migrations.CreateModel(
            name='SeqfeaturePath',
            fields=[
                ('seqfeature_path_id', models.AutoField(primary_key=True, serialize=False)),
                ('distance', models.PositiveIntegerField(blank=True, null=True)),
                ('object_seqfeature', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, related_name='object_paths', to='bioseq.Seqfeature')),
                ('subject_seqfeature', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, related_name='subject_paths', to='bioseq.Seqfeature')),
                ('term', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, to='bioseq.Term')),
            ],
            options={
                'db_table': 'seqfeature_path',
                'managed': True,
            },
        ),
        migrations.CreateModel(
            name='SeqfeatureDbxref',
            fields=[
                ('seqfeature_dbxref_id', models.AutoField(primary_key=True, serialize=False)),
                ('rank', models.SmallIntegerField(default=1, null=True)),
                ('dbxref', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, to='bioseq.Dbxref')),
                ('seqfeature', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, related_name='dbxrefs', to='bioseq.Seqfeature')),
            ],
            options={
                'db_table': 'seqfeature_dbxref',
                'managed': True,
            },
        ),
        migrations.AddField(
            model_name='seqfeature',
            name='source_term',
            field=models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, related_name='source_of', to='bioseq.Term'),
        ),
        migrations.AddField(
            model_name='seqfeature',
            name='type_term',
            field=models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, related_name='features_of_type', to='bioseq.Term'),
        ),
        migrations.AddField(
            model_name='location',
            name='seqfeature',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='locations', to='bioseq.Seqfeature'),
        ),
        migrations.AddField(
            model_name='location',
            name='term',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.DO_NOTHING, to='bioseq.Term'),
        ),
        migrations.CreateModel(
            name='BioentryRelationship',
            fields=[
                ('bioentry_relationship_id', models.AutoField(primary_key=True, serialize=False)),
                ('rank', models.IntegerField(default=1, null=True)),
                ('object_bioentry', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, related_name='object_bioentry_relationship', to='bioseq.Bioentry')),
                ('subject_bioentry', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, related_name='subject_bioentry_relationship', to='bioseq.Bioentry')),
                ('term', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, to='bioseq.Term')),
            ],
            options={
                'db_table': 'bioentry_relationship',
                'managed': True,
            },
        ),
        migrations.CreateModel(
            name='BioentryQualifierValue',
            fields=[
                ('bioentry_qualifiervalue_id', models.AutoField(primary_key=True, serialize=False)),
                ('value', models.TextField(blank=True, null=True)),
                ('rank', models.IntegerField(default=1, null=True)),
                ('bioentry', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='qualifiers', to='bioseq.Bioentry')),
                ('term', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, to='bioseq.Term')),
            ],
            options={
                'db_table': 'bioentry_qualifier_value',
                'managed': True,
            },
        ),
        migrations.CreateModel(
            name='BioentryPath',
            fields=[
                ('bioentry_path_id', models.AutoField(primary_key=True, serialize=False)),
                ('distance', models.PositiveIntegerField(blank=True, null=True)),
                ('object_bioentry', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, related_name='object_bioentry_path', to='bioseq.Bioentry')),
                ('subject_bioentry', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, related_name='subject_bioentry_path', to='bioseq.Bioentry')),
                ('term', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, to='bioseq.Term')),
            ],
            options={
                'db_table': 'bioentry_path',
                'managed': True,
            },
        ),
        migrations.CreateModel(
            name='BioentryDbxref',
            fields=[
                ('bioentry_dbxref_id', models.AutoField(primary_key=True, serialize=False)),
                ('rank', models.SmallIntegerField(default=1, null=True)),
                ('bioentry', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='dbxrefs', to='bioseq.Bioentry')),
                ('dbxref', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, to='bioseq.Dbxref')),
            ],
            options={
                'db_table': 'bioentry_dbxref',
                'managed': True,
            },
        ),
        migrations.AddField(
            model_name='bioentry',
            name='taxon',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.DO_NOTHING, to='bioseq.Taxon'),
        ),
        migrations.CreateModel(
            name='TermSynonym',
            fields=[
                ('term_synonym_id', models.AutoField(primary_key=True, serialize=False)),
                ('synonym', models.CharField(max_length=255)),
                ('term', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='synonyms', to='bioseq.Term')),
            ],
            options={
                'db_table': 'term_synonym',
                'managed': True,
                'unique_together': {('term', 'synonym')},
            },
        ),
        migrations.CreateModel(
            name='TermRelationshipTerm',
            fields=[
                ('term_relationship', models.OneToOneField(on_delete=django.db.models.deletion.DO_NOTHING, primary_key=True, serialize=False, to='bioseq.TermRelationship')),
                ('term', models.OneToOneField(on_delete=django.db.models.deletion.DO_NOTHING, to='bioseq.Term')),
            ],
            options={
                'db_table': 'term_relationship_term',
                'managed': True,
            },
        ),
        migrations.CreateModel(
            name='TermPath',
            fields=[
                ('term_path_id', models.AutoField(primary_key=True, serialize=False)),
                ('distance', models.PositiveIntegerField(blank=True, null=True)),
                ('object_term', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, related_name='object_termpaths', to='bioseq.Term')),
                ('ontology', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, to='bioseq.Ontology')),
                ('predicate_term', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, related_name='predicate_termpaths', to='bioseq.Term')),
                ('subject_term', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, related_name='subject_termpaths', to='bioseq.Term')),
            ],
            options={
                'db_table': 'term_path',
                'managed': True,
                'unique_together': {('subject_term', 'predicate_term', 'object_term', 'ontology', 'distance')},
            },
        ),
        migrations.CreateModel(
            name='TermDbxref',
            fields=[
                ('term_dbxref_id', models.AutoField(primary_key=True, serialize=False)),
                ('rank', models.SmallIntegerField(default=1, null=True)),
                ('dbxref', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, to='bioseq.Dbxref')),
                ('term', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='dbxrefs', to='bioseq.Term')),
            ],
            options={
                'db_table': 'term_dbxref',
                'managed': True,
                'unique_together': {('term', 'dbxref')},
            },
        ),
        migrations.CreateModel(
            name='TaxonName',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=255)),
                ('name_class', models.CharField(max_length=32)),
                ('taxon', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, related_name='names', to='bioseq.Taxon')),
            ],
            options={
                'db_table': 'taxon_name',
                'managed': True,
                'unique_together': {('taxon', 'name', 'name_class')},
            },
        ),
        migrations.CreateModel(
            name='SeqfeatureRelationship',
            fields=[
                ('seqfeature_relationship_id', models.AutoField(primary_key=True, serialize=False)),
                ('rank', models.IntegerField(default=1, null=True)),
                ('object_seqfeature', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, related_name='object_relationships', to='bioseq.Seqfeature')),
                ('subject_seqfeature', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='subject_relationships', to='bioseq.Seqfeature')),
                ('term', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, to='bioseq.Term')),
            ],
            options={
                'db_table': 'seqfeature_relationship',
                'managed': True,
                'unique_together': {('object_seqfeature', 'subject_seqfeature', 'term')},
            },
        ),
        migrations.AddIndex(
            model_name='seqfeaturequalifiervalue',
            index=models.Index(fields=['term'], name='seqfeature__term_id_9085a4_idx'),
        ),
        migrations.AlterUniqueTogether(
            name='seqfeaturequalifiervalue',
            unique_together={('seqfeature', 'term', 'rank')},
        ),
        migrations.AlterUniqueTogether(
            name='seqfeaturepath',
            unique_together={('object_seqfeature', 'subject_seqfeature', 'term', 'distance')},
        ),
        migrations.AlterUniqueTogether(
            name='seqfeaturedbxref',
            unique_together={('seqfeature', 'dbxref')},
        ),
        migrations.AddField(
            model_name='locationqualifiervalue',
            name='term',
            field=models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, to='bioseq.Term'),
        ),
        migrations.AlterUniqueTogether(
            name='location',
            unique_together={('seqfeature', 'rank')},
        ),
        migrations.AddField(
            model_name='dbxrefqualifiervalue',
            name='term',
            field=models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, to='bioseq.Term'),
        ),
        migrations.AlterUniqueTogether(
            name='bioentryrelationship',
            unique_together={('object_bioentry', 'subject_bioentry', 'term')},
        ),
        migrations.AlterUniqueTogether(
            name='bioentrypath',
            unique_together={('object_bioentry', 'subject_bioentry', 'term', 'distance')},
        ),
        migrations.AlterUniqueTogether(
            name='bioentrydbxref',
            unique_together={('bioentry', 'dbxref', 'rank')},
        ),
        migrations.AlterUniqueTogether(
            name='bioentry',
            unique_together={('identifier', 'biodatabase'), ('accession', 'biodatabase', 'version')},
        ),
        migrations.AlterUniqueTogether(
            name='locationqualifiervalue',
            unique_together={('location', 'term')},
        ),
        migrations.AlterUniqueTogether(
            name='dbxrefqualifiervalue',
            unique_together={('dbxref', 'term', 'rank')},
        ),
    ]
