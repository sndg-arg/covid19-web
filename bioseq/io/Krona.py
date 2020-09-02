

sqs = Sample.objects.filter(ncbi_tax__isnull=False).prefetch_related("ncbi_tax")
sqs = Assembly.objects.filter(ncbi_tax__isnull=False ).prefetch_related("ncbi_tax")
sqs = Barcode.objects.filter(ncbi_tax__isnull=False).prefetch_related("ncbi_tax")

with open("/tmp/assembly_tax.tbl","w") as h:
    for x in sqs:
        h.write("\t".join([x.name,str(x.ncbi_tax.ncbi_taxon_id)]) + "\n")
/opt/KronaTools-2.7/scripts/ImportTaxonomy.pl /tmp/assembly_tax.tbl
firefox taxonomy.krona.html


# ktImportText # /opt/KronaTools-2.7/scripts/ImportText.pl
# 2	Fats	Saturated fat
# 3	Fats	Unsaturated fat	Monounsaturated fat
# 3	Fats	Unsaturated fat	Polyunsaturated fat
# 13	Carbohydrates	Sugars
# 4	Carbohydrates	Dietary fiber
# 21	Carbohydrates
# 5	Protein
# 4