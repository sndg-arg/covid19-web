# pfam_file_path = "/data/databases/xfam/Pfam-A.hmm"
#
# term = ""
# name = ""
# ident = ""
# ver = 1
# description = ""
# is_obsolete = "F"
# ontology = Ontology.objects.get_or_create(name="Pfam")[0]
#
# with open(pfam_file_path) as pfam_handle:
#     for line in tqdm(pfam_handle):
#         if line.strip().startswith("DESC"):
#             description = line.split("DESC")[1].strip()
#
#         elif line.strip().startswith("ACC"):
#             term = line.split("ACC")[1].strip().lower()
#             ver = int(term.split(".")[1])
#             term = term.split(".")[0]
#         elif line.strip().startswith("NAME"):
#             if term:
#                 dbTerm = Term(name=name, definition=description, identifier=term,
#                               is_obsolete=is_obsolete, ontology=ontology,version=ver).save()
#                 term = ""
#                 ident = ""
#                 description = ""
#
#             name = line.split("NAME")[1].strip()