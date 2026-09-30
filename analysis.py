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
    raise NotImplementedError


def pose_rmsd(pose_pdbqt, crystal_pdbqt):
    """Heavy-atom RMSD between the docked pose and the crystal ligand (same atom order
    is NOT guaranteed - match by element and nearest neighbour, or use the symmetry-
    corrected RMSD from meeko/rdkit)."""
    raise NotImplementedError


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
