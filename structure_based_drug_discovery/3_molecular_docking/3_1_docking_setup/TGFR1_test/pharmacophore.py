#!/usr/bin/env python
import os
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
        for atom_idx_group in atom_idx_group_list:
            pseudo_atom = []
            for atom_idx in atom_idx_group:
                coor = np.array(atoms[atom_idx - 1].coords)
                pseudo_atom += [coor]
            pseudo_atom = np.array(pseudo_atom)
            pseudo_atom_coor = pseudo_atom.mean(axis=0)
            feature_coor_list += [pseudo_atom_coor]
    return feature_coor_list


def find_pharmacophore(ligand_file, cutoff, pharmacophore_coor_ref_list):

    if not os.path.exists(ligand_file):
        return [0]

    ms = pybel.readfile('pdb', ligand_file)
    ms = list(ms)
    prediction_list = list()
    for m in ms:
        feature_coor_list = find_hydrophobic_ring(m)
        prediction = 0.0
        count = 0
        for i in range(3):
            check = False
            for ref_coor in pharmacophore_coor_ref_list[i]:
                dif_min = 100
                for j, feature_coor in enumerate(feature_coor_list):
                    diff = np.linalg.norm(ref_coor - feature_coor)
                    if dif_min > diff:
                        dif_min = diff
                if dif_min < cutoff:
                    check = True
            if check:
                count += 1
        if count >= 3:
            prediction = 1.0
        prediction_list += [prediction]
    return np.array(prediction_list)


def read_score(ligand_file):

    fp = open(ligand_file)
    lines = fp.readlines()
    fp.close()
    score_list = list()
    for line in lines:
        if line.startswith('REMARK VINA RESULT:'):
            score = float(line[19:].strip().split()[0])
            score_list += [score]
    return np.array(score_list)


def main():

    cutoff = 2.5
    pharmacophore_coor_ref_list = [[[12.97766667, 64.211, 4.3995]],
                                   [[17.38733333, 68.97966667,  2.449]],
                                   [[15.542, 66.681,  6.613666], [16.22216667, 67.62283333,  8.757]]]

    list_file = 'ligand_list.smi'
    fp = open(list_file)
    lines = fp.readlines()
    fp.close()

    for line in lines:
        lis = line.strip().split()
        ligand_id = lis[0]
        ligand_file = 'dock/dock_%s.pdb' % (ligand_id)

        score = read_score(ligand_file)
        prediction = find_pharmacophore(ligand_file, cutoff,
                                        pharmacophore_coor_ref_list)
        print(ligand_file, score, prediction)


if __name__ == '__main__':
    main()
