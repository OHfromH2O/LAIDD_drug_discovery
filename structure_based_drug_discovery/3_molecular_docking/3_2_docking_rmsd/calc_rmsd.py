import numpy as np
from Bio.PDB import PDBParser

def calc_rmsd(dock_file, ref_file):
    parser = PDBParser(QUIET=True)
    s_dock = parser.get_structure('dock', dock_file)
    s_ref  = parser.get_structure('ref',  ref_file)

    coords_dock = [a.get_coord() for a in s_dock.get_atoms()
                   if a.element != 'H']
    coords_ref  = [a.get_coord() for a in s_ref.get_atoms()
                   if a.element != 'H']

    n = min(len(coords_dock), len(coords_ref))
    diff = np.array(coords_dock[:n]) - np.array(coords_ref[:n])
    rmsd = np.sqrt((diff**2).sum() / n)

    print(f"Dock atoms : {len(coords_dock)}")
    print(f"Ref  atoms : {len(coords_ref)}")
    print(f"RMSD       : {rmsd:.3f} Å")
    print("✅ 성공적 재현 (RMSD < 2Å)" if rmsd < 2.0 else "⚠️ 재현 실패 (RMSD ≥ 2Å)")

calc_rmsd('dock_3HMMA_855.pdb', '3HMMA_855.pdb')
