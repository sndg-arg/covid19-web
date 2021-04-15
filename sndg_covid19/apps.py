from django.apps import AppConfig
from django.db.models import Q
from django.utils.translation import gettext_lazy as _




class CovidConfig(AppConfig):
    name = "sndg_covid19"
    verbose_name = _("covid19")

    def ready(self):
        from sndg_covid19.models import ImportJob
        ImportJob.objects.filter(
            ~Q(status__in=["finished","error"])
        ).update(status="error", status_desc="operacion interrumpida",
                 debug_status_desc="inconcluso al inicio de la app")
