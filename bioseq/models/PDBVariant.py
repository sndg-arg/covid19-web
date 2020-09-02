# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models

from bioseq.models.Variant import Variant
from pdbdb.models import Residue, PDBResidueSet, Property, ResidueSet


class PDBVariant(models.Model):
    variant = models.ForeignKey(Variant, models.CASCADE, "pdb_variants")
    residue = models.ForeignKey(Residue, models.CASCADE, "pdb_variants")

    class Meta:
        managed = True
        db_table = 'pdbvariant'

    def ann(self):
        data = []
        for pdb_rs in self.residue.residue_sets.all():
            # TODO Breaks encapsulation from pocket residue set
            if pdb_rs.pdbresidue_set.residue_set.name == PDBResidueSet.pocket_name:
                rs_ann = f'In pocket {pdb_rs.pdbresidue_set.name} with druggability {pdb_rs.pdbresidue_set.properties_dict()[Property.druggability]}'
            elif ResidueSet.BINDING_SITE in pdb_rs.pdbresidue_set.residue_set.name:
                # bl = pdb_rs.pdbresidue_set.properties_dict().get(Property.bound_ligand, "")
                # bl = f' bound to {bl}' if bl else ""
                # rs_ann = f'{pdb_rs.pdbresidue_set.residue_set.name} {pdb_rs.pdbresidue_set.name} {bl}'
                rs_ann = pdb_rs.pdbresidue_set.description
            else:
                rs_ann = f'{pdb_rs.pdbresidue_set.residue_set.name} {pdb_rs.pdbresidue_set.name}'
            data.append(
                {"name": pdb_rs.pdbresidue_set.name, "type": pdb_rs.pdbresidue_set.residue_set.name, "desc": rs_ann})
        return data

