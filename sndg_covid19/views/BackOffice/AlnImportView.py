# from django.views.generic import TemplateView
from django.contrib import messages
from django.contrib.auth.mixins import LoginRequiredMixin
from django.db import transaction
from django.utils.translation import gettext as _
from django_tables2 import RequestConfig
from django_tables2 import SingleTableView

from sndg_covid19.forms import AlnImportForm, AlnImportTable
from sndg_covid19.tasks import process_msa


from sndg_covid19.models import ImportJob


class AlnImportView(LoginRequiredMixin,SingleTableView):
    # PermissionRequiredMixin permission_required = 'polls.add_choice'
    # login_url = '/login/'
    # redirect_field_name = 'redirect_to'
    model = ImportJob
    table_class = AlnImportTable
    template_name = "backoffice/alnimport_view.html"

    @property
    def page_size(self):
        return self.kwargs.get('page_size', 10)

    def post(self, request, **kwargs):
        return self.get(request, **kwargs)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        if self.request.method == "POST":
            form = AlnImportForm(self.request.POST, self.request.FILES)
            if form.is_valid():
                form.instance.user = self.request.user
                form.instance.version = ImportJob.objects.filter(name=form.instance.name).version + 1

                form.save()
                transaction.on_commit(lambda: process_msa.delay(form.instance.import_job_id))
                messages.add_message(self.request, messages.INFO, _(f"Procesando {form.instance.name} ... "))
                context["form"] = AlnImportForm(initial={'aln_type': 'spike'})
            else:
                context["form"] = form
        else:
            context["form"] = AlnImportForm(initial={'aln_type': 'spike'})


        qs = ImportJob.objects.all()

        table = AlnImportTable(qs, order_by="-created_at")

        context["in_progress"] = ImportJob.objects.filter(status__in=["iniciado","processing"]).count()
        if context["in_progress"]:
            context["refresh_delay"] = self.request.GET.get("refresh_delay","60")
            messages.add_message(self.request, messages.INFO,
                                 _("Hay %(in_progress)d en proceso... refrescando en 60 seg") % {
                                     "in_progress": context["in_progress"]})

        context["jobs_table"] = RequestConfig(self.request, paginate={"per_page": self.page_size}).configure(table)

        return context
