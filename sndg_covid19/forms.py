import django_tables2 as tables
from django import forms
from django.utils.translation import gettext as _
import zipfile
from .models import ImportJob


# class OverflowTemplateColumn(tables.TemplateColumn):
#     def render(self, record, table, value, bound_column, **kwargs):
#         return super(OverflowTemplateColumn, self).render(record, table, value, bound_column, **kwargs)

class AlnImportTable(tables.Table):
    T2     = '''<div><dialog id="favDialog_{{record.import_job_id}}" >Presione ESC para salir <br /><pre>{{record.errors}}</pre></dialog>
                <button onclick="document.getElementById('favDialog_{{record.import_job_id}}').showModal()">
                    Ver errores </button></div>'''

    name = tables.Column()
    aln_type = tables.Column()
    username = tables.Column(accessor="user.username")
    version = tables.Column(orderable=False)
    errors = tables.TemplateColumn(T2)
    status_desc = tables.Column(orderable=False)
    created_at = tables.Column()

    class Meta:
        model = ImportJob
        template_name = "django_tables2/bootstrap.html"

        fields = ("name", "aln_type","status","status_desc", "errors", "version","username", "created_at")
        sequence = ("name", "aln_type","status","status_desc", "errors","version","username",  "created_at")


class AlnImportForm(forms.ModelForm):

    class Meta:
        model = ImportJob
        fields = ["name", "aln_type", "fasta", "csv"]

    # this function will be used for the validation
    def clean(self):

        # data from the form is fetched using super function
        super(AlnImportForm, self).clean()


        fasta = self.cleaned_data.get('fasta')
        csv = self.cleaned_data.get('csv')

        if not fasta.name.endswith("fasta.zip") :
            self._errors['fasta'] = self.error_class([
                'el archivo de secuencias debe ser un fasta comprimido'])
            try:
                the_zip_file = zipfile.ZipFile(fasta.file)
            except zipfile.BadZipFile:
                self._errors['fasta'] = self.error_class([
                    f'{fasta.name} not esta comprimido'])


        if not csv.name.endswith("csv") :
            self._errors['csv'] = self.error_class([
                'el archivo de propiedades debe ser un csv'])



        return self.cleaned_data
        # if not (self.cleaned_data['fasta'] or self.cleaned_data['text']):
        #     raise forms.ValidationError(
        #         _('Invalid value: %(value)s'),
        #         code='invalid',
        #         params={'value': '42'},
        #     )

