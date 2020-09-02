import re
from urllib.parse import urlparse, parse_qsl

from django.db.models import Q
from django import template

register = template.Library()

from pdbdb.models import PDB, Residue, ResidueSet, Property, PDBResidueSet

@register.inclusion_tag('tags/structure.html')
def pdb_structure(pdbid):
    pdbobj = PDB.objects.prefetch_related("residues").get(code=pdbid.lower())
    context = {}

    context["chains"] = [{"name": x} for x in set([r.chain for r in pdbobj.residues.all() if r.chain.strip()])]
    context["layers"] = []

    resnames = list(Residue.objects.filter(pdb=pdbobj, type=Residue.HETATOM).values("resname").distinct())
    if "HOH" in resnames or "WAT" in resnames:
        context["layers"].append("water")
    elif len([x for x in resnames if resnames not in ["HOH", "WAT"]]):
        context["layers"].append("hetero")

    from collections import defaultdict
    dna_data = defaultdict(lambda: [])
    for chain_resname in Residue.objects.filter(
            pdb=pdbobj, type="R", resname__in=["DA", "DC", "DG", "DT"]).values("chain", "resname").distinct():
        dna_data[chain_resname["chain"]].append(chain_resname["resname"])
    context["dna"] = []
    for chain, residues in dna_data.items():
        if len(residues) <= 4:
            context["layers"].append("dna")
            context["dna"] += [x for x in context["chains"] if x["name"] == chain]
            context["chains"] = [x for x in context["chains"] if x["name"] != chain]

    ds = Property.objects.get(name="druggability_score")
    rs = ResidueSet.objects.get(name="FPocketPocket")

    # sq = ResidueSetProperty.objects.select_related(pdbresidue_set)\
    #     .filter(property=ds,value__gte=0.2,pdbresidue_set=OuterRef("id"))

    context["pockets"] = PDBResidueSet.objects.prefetch_related("properties__property",
                                                                "residue_set_residue__residue__atoms").filter(
        Q(pdb=pdbobj), Q(residue_set=rs), Q(properties__property=ds) & Q(properties__value__gte=0.2)).all()
    for p in context["pockets"]:
        p.druggability = [x.value for x in p.properties.all() if x.property == ds][0]
        p.atoms = []
        for rsr in p.residue_set_residue.all():
            for a in rsr.residue.atoms.all():
                p.atoms.append(a.serial)
    context["residuesets"] = [{"name": "csa", "residues": range(700, 750)}]
    context["pdbid"] = pdbid.lower()

    return {**context}
