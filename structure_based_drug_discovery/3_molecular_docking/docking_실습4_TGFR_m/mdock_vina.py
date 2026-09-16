#!/usr/bin/env python
import os
import numpy as np
import subprocess
from pbi import ligandtools
from pbi.pdbtools import PDBtools
from openbabel import pybel
from multiprocessing import Manager
from multiprocessing import Process
from multiprocessing import Queue
import argparse


def add_hydrogen(smi0, pH=7.4):
    try:
        m = pybel.readstring("smi", smi0)
        m.OBMol.AddHydrogens(False, True, pH)
        smi_p = m.write("smi").strip()
    except Exception:
        smi_p = smi0
    return smi_p


def docking(vina, config_file, ligand_file, output_file):

    run_line = '%s --config %s --ligand %s --out %s' % (
        vina, config_file, ligand_file, output_file)
    e = None
    try:
        subprocess.check_output(run_line.split(),
                                stderr=subprocess.STDOUT,
                                universal_newlines=True)
    except Exception as e:
        return e
    return e


def read_score(dock_pdb_file, vina_v='vina'):
    fp = open(dock_pdb_file)
    lines = fp.readlines()
    fp.close()
    score_list = list()
    for line in lines:
        if vina_v == 'vina':
            if line.startswith('REMARK VINA RESULT'):
                score = float(line[19:29].strip())
                score_list += [score]
        elif vina_v == 'smina':
            if line.startswith('REMARK minimizedAffinity'):
                score = float(line[24:36].strip())
                score_list += [score]
    score_list = np.array(score_list)
    return score_list


def creator(q, data, num_sub_proc):
    for d in data:
        idx = d[0]
        q.put((idx, d[1]))
    for i in range(0, num_sub_proc):
        q.put('DONE')


def worker(q, param, return_dict):
    vina, config_file, pH, dock_dir, vina_v = param
#    pid = os.getpid()
    while True:
        qqq = q.get()
        if qqq == 'DONE':
            # print('proc =', pid)
            break

        (idx, d) = qqq
        mol_id = d[0]
        smi = d[1]
        mol_id2 = mol_id[0:7]
        dock_dir1 = dock_dir + "/" + mol_id2
        if not os.path.exists(dock_dir1):
            try:
                os.makedirs(dock_dir1)
            except FileExistsError as e:
                print(e, flush=True)

        pdb_file = dock_dir1 + '/' + mol_id + '.pdb'
        pdbqt_file = dock_dir1 + '/' + mol_id + '.pdbqt'
        dock_pdbqt_file = dock_dir1 + '/' + '/dock_' + mol_id + '.pdbqt'
        dock_pdb_file = dock_dir1 + '/' + '/dock_' + mol_id + '.pdb'

        if pH is not None:
            smi_p = add_hydrogen(smi, pH=pH)
        else:
            smi_p = smi
        ligandtools.gen_3d(smi_p, pdb_file, mol_id=mol_id, file_format='pdb',
                           timeout=20)
        PDBtools.ligand_to_pdbqt(pdb_file, pdbqt_file)
        docking(vina, config_file, pdbqt_file, dock_pdbqt_file)
        ligandtools.pdbqt_to_pdb_ref(dock_pdbqt_file, dock_pdb_file, pdb_file)
        dock_score = read_score(dock_pdb_file, vina_v=vina_v)

        return_dict[idx] = dock_score


def main():

    parser = argparse.ArgumentParser(description='docking with multi process')
    parser.add_argument('-v', '--vina', type=str, required=True,
                        default='vina', help='vina run file ')
    parser.add_argument('-c', '--dock_config', type=str, required=True,
                        default='config.txt', help='docking config file ')
    parser.add_argument('-s', '--list_file', type=str, required=True,
                        default='ligand_list.smi', help='molecule list file')
    parser.add_argument('-d', '--out_dir', type=str, required=True,
                        default='dock', help='out_dir')
    parser.add_argument('-o', '--score_file', type=str, required=False,
                        default='docking.txt', help='number of sub proces')
    parser.add_argument('-p', '--ncpu', type=int, required=False,
                        default='1', help='number of sub proces')

    args = parser.parse_args()

    vina = args.vina
    config_file = args.dock_config
    list_file = args.list_file
    out_dir = args.out_dir
    score_file = args.score_file
    num_sub_proc = args.ncpu
    pH = 7.4

    vina_v = 'vina'
    if vina == 'smina':  # vina or smina
        vina_v = 'smina'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    fp = open(list_file)
    lines = fp.readlines()
    fp.close()

    smiles_list = list()
    for line in lines:
        lis = line.strip().split()
        mol_id = lis[0]
        smi = lis[1]
        smiles_list += [[mol_id, smi]]

    data = list(enumerate(smiles_list))
    num_data = len(data)
    num_sub_proc = min(num_sub_proc, num_data)

    param = (vina, config_file, pH, out_dir, vina_v)

    q1 = Queue()
    manager = Manager()
    return_dict = manager.dict()
    proc_master = Process(target=creator,
                          args=(q1, data, num_sub_proc))
    proc_master.start()

    procs = []
    for sub_id in range(0, num_sub_proc):
        proc = Process(target=worker, args=(q1, param, return_dict))
        procs.append(proc)
        proc.start()
    q1.close()
    q1.join_thread()
    proc_master.join()
    for proc in procs:
        proc.join()

    fp = open(score_file, 'w')
    line_out = 'mol_id SMILES docking_score\n'
    fp.write(line_out)
    for dd in data:
        key = dd[0]
        mol_id = dd[1][0]
        smi = dd[1][1]
        line_out = '%s %s' % (mol_id, smi)
        if key in return_dict:
            dscore = return_dict[key]
            line_score = ''
            for score0 in dscore:
                line_score += '%.3f,' % (score0)
            line_out = line_out + ' ' + line_score.strip(',')

        else:
            line_out += ' None'
        line_out += '\n'
        fp.write(line_out)
    fp.close()


if __name__ == "__main__":
    main()
