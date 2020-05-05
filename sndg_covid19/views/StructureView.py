from django.db.models import Q
from django.views.generic import TemplateView

from pdbdb.models import PDB, ResidueSet, PDBResidueSet, Property, Residue
from bioseq.models.Bioentry import Bioentry


class StructureView(TemplateView):
    # http://nglviewer.org/ngl/api/class/src/stage/stage.js~Stage.html#instance-method-loadFile
    template_name = "structure_detail.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["pdb"] = self.kwargs["pdbid"].lower()

        context["protein"] = Bioentry.objects.get(bioentry_id=self.kwargs["prot_id"])
        pdbobj = PDB.objects.prefetch_related("residues").get(code=self.kwargs["pdbid"].lower())

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
        return context
