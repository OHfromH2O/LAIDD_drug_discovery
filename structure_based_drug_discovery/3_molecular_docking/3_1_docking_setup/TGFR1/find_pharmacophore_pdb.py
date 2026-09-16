#!/usr/bin/env python
import numpy as np
from openbabel import pybel

def find_hydrophobic_ring(m):
    smarts_pattern_list = ['a1:a:a:a:1', 'a1:a:a:a:a1', 'a1:a:a:a:a:a:1',
                           'a1:a:a:a:a:a:a:1', 'a1:a:a:a:a:a:a:a:1']
#    smarts_pattern_list = ['*1***1', '*1****1', '*1*****1',
#                           '*1******1', '*1*******1']
    atoms = m.atoms
    feature_coor_list = list()
    for smarts_pattern in smarts_pattern_list:
        smarts = pybel.Smarts(smarts_pattern)
        atom_idx_group_list = smarts.findall(m)
#        print(atom_idx_group_list)
        for atom_idx_group in atom_idx_group_list:
            pseudo_atom = []
            for atom_idx in atom_idx_group:
                coor = np.array(atoms[atom_idx-1].coords)
                pseudo_atom += [coor]
            pseudo_atom = np.array(pseudo_atom)
            pseudo_atom_coor = pseudo_atom.mean(axis=0)
            feature_coor_list += [pseudo_atom_coor]
    return feature_coor_list

def main():

    pdb_chain_id = '3HMMA'
    ligand_id = '855'
    ref_ligand_file = 'pdb/%s_%s.pdb' %(pdb_chain_id, ligand_id)
    ms = pybel.readfile('pdb', ref_ligand_file)
    m_ref = list(ms)[0]
    feature_coor_ref_list = find_hydrophobic_ring(m_ref)
    print(feature_coor_ref_list)


if __name__ == '__main__':
    main()


