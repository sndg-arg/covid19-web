from django.db.models import Q
from django.http.response import HttpResponse
from django.views.generic import TemplateView

from pdbdb.models import PDB, ResidueSet, PDBResidueSet, Property, Residue
from bioseq.models.Bioentry import Bioentry
from bioseq.models.PDBVariant import PDBVariant
import json
from collections import Counter


class StructureView(TemplateView):
    # http://nglviewer.org/ngl/
    # http://nglviewer.org/ngl/api/class/src/stage/stage.js~Stage.html#instance-method-loadFile
    template_name = "structure_detail.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["pdb"] = self.kwargs["pdbid"].lower()

        context["protein"] = Bioentry.objects.prefetch_related("qualifiers__term").get(
            bioentry_id=self.kwargs["prot_id"])
        pdbobj = PDB.objects.get(code=self.kwargs["pdbid"].lower())

        context["chains"] = [{"name": x} for x in set([r.chain for r in pdbobj.residues.all() if r.chain.strip()])]
        context["layers"] = []

        resnames = list(
            Residue.objects.filter(pdb=pdbobj).exclude(type__in=["R"]).values_list("resname", flat=True).distinct())
        if "HOH" in resnames or "WAT" in resnames:
            context["layers"].append("water")

        if len([x for x in resnames if x not in ["HOH", "WAT", "STP"]]):
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

        context["residuesets"] = []
        for pdb_rs in PDBResidueSet.objects.filter(pdb=pdbobj).exclude(
            residue_set__name=PDBResidueSet.pocket_name).prefetch_related("residue_set_residue__residue"):
            sele = "or".join([f'({rs_res.residue.resid} and :{rs_res.residue.chain})' for rs_res in
                              pdb_rs.residue_set_residue.all()])
            context["residuesets"].append({"name": pdb_rs.name, "description": pdb_rs.description, "sele": sele})

        ds = Property.objects.get(name=Property.druggability)
        rs = ResidueSet.objects.get(name=PDBResidueSet.pocket_name)

        # sq = ResidueSetProperty.objects.select_related(pdbresidue_set)\
        #     .filter(property=ds,value__gte=0.2,pdbresidue_set=OuterRef("id"))

        # .prefetch_related("properties__property",
        #                   "residue_set_residue__residue__atoms")
        raw_pockets = PDBResidueSet.objects.filter(
            Q(pdb=pdbobj), Q(residue_set=rs), Q(properties__property=ds) & Q(properties__value__gte=0.2)).values(
            "properties__value", "residue_set_residue__residue__chain", "residue_set_residue__residue__resid", "name"
        )
        context["pockets"] = []
        curr_pocket = None
        curr_drug = None
        residues = []
        for p in raw_pockets:
            res = f'(:{p["residue_set_residue__residue__chain"]} and {p["residue_set_residue__residue__resid"]})'
            if curr_pocket == p["name"]:
                residues.append(res)
            else:
                if curr_pocket:
                    context["pockets"].append({"druggability": curr_drug, "name": curr_pocket,
                                               "residues": residues})
                curr_pocket = p["name"]
                curr_drug = p["properties__value"]
                residues = [res]
        if curr_pocket:
            context["pockets"].append({"druggability": curr_drug, "name": curr_pocket,
                                       "residues": residues})

        # for rsr in p.residue_set_residue.all():
        # for a in rsr.residue.atoms.all():
        # p.atoms.append(a.serial)

        context["variants"] = list(PDBVariant.objects.prefetch_related(
            "variant", "variant__sample_variants__sample",
            "residue__residue_sets__pdbresidue_set__residue_set",
            "residue__residue_sets__pdbresidue_set__properties__property"
        ).filter(variant__sample_variants__sample__country="Argentina", residue__pdb__code=context["pdb"]))

        if context["variants"]:
            context["layers"].append("variants")
            sele = []
            for v in context["variants"]:
                sele.append(f'({v.residue.resid} and :{v.residue.chain} and .CA)')
            context["variants_sele"] = " or ".join(set(sele))
        return context


class StructureStaticView(TemplateView):
    template_name = "structure_static.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["pdb"] = self.kwargs["pdbid"].lower()

        context["protein"] = Bioentry.objects.prefetch_related("qualifiers__term").get(
            bioentry_id=self.kwargs["prot_id"])
        pdbobj = PDB.objects.get(code=self.kwargs["pdbid"].lower())

        context["chains"] = [{"name": x} for x in set([r.chain for r in pdbobj.residues.all() if r.chain.strip()])]
        context["layers"] = []

        resnames = list(
            Residue.objects.filter(pdb=pdbobj).exclude(type__in=["R"]).values_list("resname", flat=True).distinct())
        if "HOH" in resnames or "WAT" in resnames:
            context["layers"].append("water")

        if len([x for x in resnames if x not in ["HOH", "WAT", "STP"]]):
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

        context["residuesets"] = []
        for pdb_rs in PDBResidueSet.objects.filter(pdb=pdbobj).exclude(
            residue_set__name=PDBResidueSet.pocket_name).prefetch_related("residue_set_residue__residue"):
            sele = "or".join([f'({rs_res.residue.resid} and :{rs_res.residue.chain})' for rs_res in
                              pdb_rs.residue_set_residue.all()])
            context["residuesets"].append({"name": pdb_rs.name, "description": pdb_rs.description, "sele": sele})

        context["variants"] = list(PDBVariant.objects.prefetch_related(
            "variant", "variant__sample_variants__sample",
            "residue__residue_sets__pdbresidue_set__residue_set",
            "residue__residue_sets__pdbresidue_set__properties__property"
        ).filter(variant__sample_variants__sample__country="Argentina", residue__pdb__code=context["pdb"]))
        context["pdbres2pos"] = {x.residue.chain + "_" + str(x.residue.resid): str(x.variant.pos + 1)
                                 for x in context["variants"]}
        if context["variants"]:
            context["layers"].append("variants")
            sele = []
            for v in context["variants"]:
                sele.append(f'({v.residue.resid} and :{v.residue.chain} and .CA)')
            counter = Counter(sele)
            context["variants_sele"] = " or ".join([k for k, v in counter.items() if v > 1])
        return context


def pdb_download(request, pdbid):
    response = HttpResponse(content_type='text/json')
    if "download" in request.GET:
        response['Content-Disposition'] = f'attachment; filename="{pdbid}.json"'

    pdbobj = PDB.objects.prefetch_related("residue_sets__residue_set",
                                          "residue_sets__properties__property"
                                          ).get(code=pdbid.lower())

    alphas = []
    for r in Residue.objects.filter(pdb=pdbobj, resname="STP"):
        alphas += r.lines()

    rss = []
    data = {
        "pdb": pdbid, "residue_sets": rss,
        "alpha_spheres": "\n".join(alphas),
    }

    for rs in pdbobj.residue_sets.all():
        residue_set = {"name": rs.name, "type": rs.residue_set.name, "desc": rs.description,
                       "props": [{"name": prop.property.name, "value": prop.prop_value()} for prop in
                                 rs.properties.all()],
                       "residues": [{"icode": rsr.residue.icode, "resname": rsr.residue.resname,
                                     "chain": rsr.residue.chain, "resid": rsr.residue.resid}
                                    for rsr in rs.residue_set_residue.all()]}
        rss.append(residue_set)

    json.dump(data, response, indent=4, sort_keys=True)
    return response
