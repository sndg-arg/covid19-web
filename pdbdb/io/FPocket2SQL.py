import json
import os
import sys
import warnings

from Bio import BiopythonWarning, BiopythonParserWarning, BiopythonDeprecationWarning, BiopythonExperimentalWarning
from django.db import transaction
from django.db.models import Max

from pdbdb.models import PDB, Property, PDBResidueSet, ResidueSet, ResidueSetProperty, ResidueSetResidue, \
    AtomResidueSet, Residue, Atom

warnings.simplefilter('ignore', RuntimeWarning)
warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', BiopythonParserWarning)
warnings.simplefilter('ignore', BiopythonDeprecationWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)

from tqdm import tqdm


from SNDG.Structure.FPocket import FPocket, fpocket_properties_map

pocket_prop_map = {v: k for k, v in fpocket_properties_map.items()}


class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)


class FPocket2SQL:

    def __init__(self):
        self.pocket_props = None
        self.rsfpocker = None
        self.res_pockets = None
        self.pdb = None


    def create_or_get_pocket_properties(self):
        self.pocket_props = {name: Property.objects.get_or_create(name=name, description=desc)[0]
                             for name, desc in fpocket_properties_map.items()}

        self.rsfpocker = ResidueSet.objects.get_or_create(name="FPocketPocket", description="")[0]

    def load_pdb(self,code):
        self.pdb = PDB.objects.prefetch_related("residues__atoms").get(code=code)
        Residue.objects.filter(pdb=self.pdb, resname="STP").delete()
        PDBResidueSet.objects.filter(pdb=self.pdb,residue_set__name="FPocketPocket").delete()


    def run_fpocket(self,tmp="/tmp/pockets/",pdb_path=None):
        if not os.path.exists(tmp):
            os.makedirs(tmp)
        if not pdb_path:
            pdb_path = tmp + self.pdb.code + ".pdb"

        pockets_path = os.path.abspath(tmp) + "/" + self.pdb.code + ".pockets.json"
        if not os.path.exists(pockets_path):
            if not os.path.exists(pdb_path):
                with open(pdb_path, "w") as h:
                    h.write(self.pdb.text)
            res = FPocket(pdb_path, tmp).hunt_pockets()
            res.save(pockets_path )
        with open(pockets_path) as h:
            self.res_pockets = json.load(h)

    def _process_pocket_alphas(self,pocket,nro_atm):
        res_alpha = {}
        for stp_line in pocket.as_lines:
            nro_atm += 1
            resid = int(stp_line[22:26])
            if resid in res_alpha:
                r = res_alpha[resid]
            else:
                r = Residue(pdb=self.pdb, chain=stp_line[22:23], resid=resid,
                            type="",
                            resname="STP", disordered=1)
                r.save()
                res_alpha[resid] = r
            Atom(residue=r, serial=nro_atm, name=stp_line[12:16],
                 x=float(stp_line[30:38].strip()), y=float(stp_line[38:46].strip()),
                 z=float(stp_line[46:54].strip()),
                 occupancy=float(stp_line[54:60].strip()), bfactor=float(stp_line[60:66].strip()),
                 element="").save()
        return nro_atm

    def load_pockets(self):
        rss = []
        nro_atm = Atom.objects.filter(residue__pdb=self.pdb).aggregate(Max("serial"))["serial__max"]
        with transaction.atomic():
            for pocket in tqdm(self.res_pockets):
                pocket = Struct(**pocket)
                rs = PDBResidueSet(name="%i" % pocket.number, pdb=self.pdb, residue_set=self.rsfpocker)
                rss.append(rs)
                rs.save()
                nro_atm = self._process_pocket_alphas(pocket,nro_atm)


                atoms = Atom.objects.select_related("residue").filter(residue__pdb=self.pdb,
                                                                      serial__in=[int(x) for x in pocket.atoms])
                residues = set([x.residue for x in atoms])
                rs_dict = {}
                for r in residues:
                    rsr = ResidueSetResidue(residue=r, pdbresidue_set=rs)
                    rsr.save()
                    rs_dict[r.id] = rsr

                for atom in atoms:
                    AtomResidueSet(atom=atom, pdb_set=rs_dict[atom.residue.id]).save()

                for k, v in pocket.properties.items():
                    prop = pocket_prop_map[k]
                    prop_model = self.pocket_props[prop]
                    ResidueSetProperty(pdbresidue_set=rs, property=prop_model, value=v).save()