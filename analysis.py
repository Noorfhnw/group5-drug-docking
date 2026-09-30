"""
Shared analysis script. Both students edit THIS file on their own branches.

Student A implements : dock_ligand, pose_rmsd
Student B implements : dock_all, plot_ranking
BOTH implement       : load_box, main, summary_sentence   <- expect a merge conflict here

Run: python analysis.py
"""
import yaml
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CONFIG = yaml.safe_load(open("config.yaml"))

import glob, os
import numpy as np
from vina import Vina

# ---------- BOTH ----------
def load_box(config):
    """Return (center [x,y,z], size [x,y,z]) from data/box.yaml, using CONFIG for overrides."""
    raise NotImplementedError


# ---------- Student A ----------
def dock_ligand(receptor, ligand, center, size, exhaustiveness):
    """Dock one ligand with Vina. Return (best_score, pose_pdbqt_string).
    Hint: v = Vina(sf_name='vina'); v.set_receptor(...); v.set_ligand_from_file(...);
    v.compute_vina_maps(center=..., box_size=...); v.dock(exhaustiveness=..., n_poses=5)"""
    v = Vina(sf_name="vina", verbosity=0)
    v.set_receptor(receptor)
    v.set_ligand_from_file(ligand)
    v.compute_vina_maps(center=list(center), box_size=list(size))
    v.dock(exhaustiveness=exhaustiveness, n_poses=5)

    # energies(): one row per pose, first column = total score, best pose first
    best_score = float(np.atleast_2d(v.energies())[0, 0])
    pose = v.poses(n_poses=1)
    return best_score, pose


def _pdbqt_heavy_atoms(pdbqt):
    """Parse a PDBQT (file path or text) -> (elements list, Nx3 coords) of heavy atoms."""
    if os.path.exists(pdbqt):
        with open(pdbqt) as fh:
            text = fh.read()
    else:
        text = pdbqt

    elements, coords = [], []
    for line in text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        # PDBQT: the autodock type sits in the trailing columns (77-79)
        atype = line[77:79].strip() or line[12:16].strip()
        if atype.upper() in ("CL", "BR"):
            element = atype.capitalize()
        else:
            element = atype[0].upper()
        if element == "H":
            continue  # heavy atoms only
        # AutoDock types: A = aromatic carbon, NA/OA/SA = H-bond acceptor flavours
        element = {"A": "C"}.get(element, element)
        elements.append(element)
        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return elements, np.asarray(coords, dtype=float)


def pose_rmsd(pose_pdbqt, crystal_pdbqt):
    """Heavy-atom RMSD between the docked pose and the crystal ligand (same atom order
    is NOT guaranteed - match by element and nearest neighbour, or use the symmetry-
    corrected RMSD from meeko/rdkit)."""
    pose_el, pose_xyz = _pdbqt_heavy_atoms(pose_pdbqt)
    xtal_el, xtal_xyz = _pdbqt_heavy_atoms(crystal_pdbqt)

    if len(pose_el) == 0 or len(xtal_el) == 0:
        raise ValueError("no heavy atoms parsed from one of the PDBQT inputs")
    if len(pose_el) != len(xtal_el):
        raise ValueError(
            f"atom count mismatch: pose has {len(pose_el)}, crystal has {len(xtal_el)}"
        )

    # Optimal one-to-one assignment, but only between atoms of the same element:
    # forbidden pairs get an infinite cost so the solver never picks them.
    from scipy.optimize import linear_sum_assignment

    d2 = ((pose_xyz[:, None, :] - xtal_xyz[None, :, :]) ** 2).sum(axis=2)
    same = np.equal(np.array(pose_el)[:, None], np.array(xtal_el)[None, :])
    cost = np.where(same, d2, 1e9)
    rows, cols = linear_sum_assignment(cost)
    if not same[rows, cols].all():
        raise ValueError("pose and crystal ligand have different element compositions")

    return float(np.sqrt(d2[rows, cols].mean()))


# ---------- Student B ----------
def dock_all(receptor, ligand_dir, center, size, exhaustiveness):
    """Dock every *.pdbqt in ligand_dir. Return a DataFrame: ligand, score."""
    raise NotImplementedError


def plot_ranking(scores, out="results/ranking.png"):
    """Bar plot of Vina scores sorted best-first, sotorasib highlighted."""
    raise NotImplementedError


# ---------- BOTH ----------
def summary_sentence(rmsd, scores):
    """One sentence: redocking RMSD, and the rank of sotorasib among the six ligands."""
    raise NotImplementedError


def main():
    center, size = load_box(CONFIG)
    # Student A: redock sotorasib + RMSD -> results/redock.csv
    # Student B: dock all + ranking plot
    # After the merge: both, then print(summary_sentence(rmsd, scores))
    raise NotImplementedError


if __name__ == "__main__":
    main()
